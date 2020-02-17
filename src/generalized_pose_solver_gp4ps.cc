// Copyright (c) 2020, Torsten Sattler
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the name of the copyright holder nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// author: Torsten Sattler, torsten.sattler.de@googlemail.com

#include <iostream>

#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <PoseLib/gp4ps.h>

#include "generalized_pose_solver_gp4ps.h"

namespace multi_camera_pose {

// Non-linear refinement code for Ceres-based optimization.
struct ReprojectionError {
  ReprojectionError(double x, double y, double X, double Y, double Z,
                    double fx, double fy, double aax_x, double aax_y,
                    double aax_z, double trans_x, double trans_y,
                    double trans_z)
  : point2D_x(x),
  point2D_y(y),
  point3D_X(X),
  point3D_Y(Y),
  point3D_Z(Z),
  f_x(fx),
  f_y(fy),
  ax_x(aax_x), ax_y(aax_y), ax_z(aax_z), t_x(trans_x), t_y(trans_y),
  t_z(trans_z){}
  
  template <typename T>
  bool operator()(const T* const camera, T* residuals) const {
    // The entries 3 to 5 are the camera translation.
    T p[3];
    p[0] = static_cast<T>(point3D_X);
    p[1] = static_cast<T>(point3D_Y);
    p[2] = static_cast<T>(point3D_Z);
    
    // The first three entries correspond to the rotation matrix stored in an
    // angle-axis representation.
    T p_rot[3];
    ceres::AngleAxisRotatePoint(camera, p, p_rot);
    
    T p_t[3];
    for (int i = 0; i < 3; ++i) p_t[i] = p_rot[i] + camera[3 + i];
    
    // Transforms to the local camera coordinate system.
    T aax[3];
    aax[0] = static_cast<T>(ax_x);
    aax[1] = static_cast<T>(ax_y);
    aax[2] = static_cast<T>(ax_z);
    ceres::AngleAxisRotatePoint(aax, p_t, p);
    p_t[0] = p[0] + static_cast<T>(t_x) * camera[6];
    p_t[1] = p[1] + static_cast<T>(t_y) * camera[6];
    p_t[2] = p[2] + static_cast<T>(t_z) * camera[6];
    
    T x_proj = static_cast<T>(f_x) * p_t[0] / p_t[2];
    T y_proj = static_cast<T>(f_y) * p_t[1] / p_t[2];
    
    residuals[0] = static_cast<T>(point2D_x) - x_proj;
    residuals[1] = static_cast<T>(point2D_y) - y_proj;
    
    return true;
  }
  
  // Factory function
  static ceres::CostFunction* CreateCost(const double x, const double y,
                                         const double X, const double Y,
                                         const double Z, const double fx,
                                         const double fy, double aax_x,
                                         double aax_y, double aax_z,
                                         double trans_x, double trans_y,
                                         double trans_z) {
    return (new ceres::AutoDiffCostFunction<ReprojectionError, 2, 7>(
                                                                               new ReprojectionError(x, y, X, Y, Z, fx, fy, aax_x, aax_y, aax_z, trans_x, trans_y, trans_z)));
  }
  
  // Assumes that the measurement is centered around the principal point.
  // This camera model does not take any radial distortion into account. If
  // radial distortion is present, one should undistort the measurements first.
  double point2D_x;
  double point2D_y;
  // The 3D point position is fixed as we are only interested in refining the
  // camera parameters.
  double point3D_X;
  double point3D_Y;
  double point3D_Z;
  double f_x;
  double f_y;
  double ax_x, ax_y, ax_z;
  double t_x, t_y, t_z;
};
  
GeneralizedPoseSolverGP4Ps::
    GeneralizedPoseSolverGP4Ps(
        const MultiCameraRig& rig, const double squared_inlier_threshold,
        const Points2D& points2D, const Points3D& points3D,
        const std::vector<int>& camera_indices)
    : rig_(rig),
      squared_inlier_threshold_(squared_inlier_threshold),
      points2D_(points2D),
      points3D_(points3D),
      camera_indices_(camera_indices) {
  num_data_ = static_cast<int>(points2D_.size());
  num_cameras_ = static_cast<int>(rig.cameras.size());
  
  // Initializes the viewing rays.
  rays_.resize(num_data_);
  for (int i = 0; i < num_data_; ++i) {
    const Camera& cam = rig_.cameras[camera_indices[i]];
    Eigen::Vector3d d(points2D[i][0] / cam.focal_x,
                      points2D[i][1] / cam.focal_y, 1.0);
    d.normalize();
    rays_[i] = cam.R.transpose() * d;
  }
}

int GeneralizedPoseSolverGP4Ps::MinimalSolver(
    const std::vector<int>& sample, GenCamPoses* poses) const {
  poses->clear();
  
  std::vector<Vector3d> sample_rays(4);
  std::vector<Vector3d> sample_points3D(4);
  std::vector<Vector3d> sample_positions(4);
  
  for (int i = 0; i < 4; ++i) {
    sample_rays[i] = rays_[sample[i]];
    sample_points3D[i] = points3D_[sample[i]];
    sample_positions[i] = rig_.cameras[camera_indices_[sample[i]]].c;
  }

  pose_lib::CameraPoseVector pl_poses;
  pose_lib::gp4ps(sample_positions, sample_rays, sample_points3D, &pl_poses);

  if (pl_poses.empty()) return 0;
  for (const pose_lib::CameraPose& pose : pl_poses) {
    const double kError = EvaluateModelOnPoint(pose, sample[3]);
    if (kError < squared_inlier_threshold_) {
      poses->push_back(pose);
//      break;
    }
  }

  return static_cast<int>(poses->size());
}

int GeneralizedPoseSolverGP4Ps::NonMinimalSolver(
    const std::vector<int>& sample, GenCamPose* pose) const {
  // Call the minimal solver and refine the pose.
  GenCamPoses poses;
  MinimalSolver(sample, &poses);
  if (poses.empty()) return 0;
  
  *pose = poses[0];
  LeastSquares(sample, pose);
  return 1;
}

double GeneralizedPoseSolverGP4Ps::EvaluateModelOnPoint(
    const GenCamPose& pose, int i) const {
  // Transforms into the coordinate system of the rig.
  Vector3d p_r = (pose.R * points3D_[i] + pose.t);
  
  // Transforms into the coordinate system of the camera.
  const Camera& cam = rig_.cameras[camera_indices_[i]];
  Vector3d p_c = cam.R * p_r + cam.t * pose.alpha;

  // Check whether point projects behind the camera.
  if (p_c[2] < 0.0) return std::numeric_limits<double>::max();

  Eigen::Vector2d p_2d = p_c.head<2>() / p_c[2];
  p_2d[0] *= cam.focal_x;
  p_2d[1] *= cam.focal_y;

  return (p_2d - points2D_[i]).squaredNorm();
}

void GeneralizedPoseSolverGP4Ps::LeastSquares(
    const std::vector<int>& sample, GenCamPose* pose) const {
  Eigen::AngleAxisd aax(pose->R);
  Eigen::Vector3d aax_vec = aax.axis() * aax.angle();
  double camera[7];
  camera[0] = aax_vec[0];
  camera[1] = aax_vec[1];
  camera[2] = aax_vec[2];
  camera[3] = pose->t[0];
  camera[4] = pose->t[1];
  camera[5] = pose->t[2];
  camera[6] = pose->alpha;
  
 
  ceres::Problem refinement_problem;
  const int kSampleSize = static_cast<int>(sample.size());
  for (int i = 0; i < kSampleSize; ++i) {
    const int kIdx = sample[i];
    const Eigen::Vector2d& p_img = points2D_[kIdx];
    const Eigen::Vector3d& p_3D = points3D_[kIdx];
    const Camera& cam = rig_.cameras[camera_indices_[kIdx]];
    Eigen::AngleAxisd aax2(cam.R);
    Eigen::Vector3d aax2_vec = aax2.axis() * aax2.angle();
    
    ceres::CostFunction* cost_function =
    ReprojectionError::CreateCost(p_img[0], p_img[1], p_3D[0], p_3D[1], p_3D[2], cam.focal_x, cam.focal_y, aax2_vec[0], aax2_vec[1], aax2_vec[2], cam.t[0], cam.t[1], cam.t[2]);
    refinement_problem.AddResidualBlock(cost_function, nullptr, camera);
  }
  
  ceres::Solver::Options options;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_progress_to_stdout = false;
  //  options.function_tolerance = 0.000001;
  ceres::Solver::Summary summary;
  ceres::Solve(options, &refinement_problem, &summary);
  
  //  std::cout << summary.BriefReport() << std::endl;
  
  if (summary.IsSolutionUsable()) {
    Eigen::Vector3d axis(camera[0], camera[1], camera[2]);
    double angle = axis.norm();
    axis.normalize();
    aax.axis() = axis;
    aax.angle() = angle;
    
    pose->R = aax.toRotationMatrix();
    pose->t = Eigen::Vector3d(camera[3], camera[4], camera[5]);
    pose->alpha = camera[6];
  }
}

}  // namespace multi_camera_pose
