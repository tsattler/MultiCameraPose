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

#include <PoseLib/gp4ps.h>

#include "generalized_pose_solver_gp4ps.h"

namespace multi_camera_pose {

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
                      points2D[i][0] / cam.focal_y, 1.0);
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
//  std::cout << " Calling gp4s" << std::endl;
  pose_lib::CameraPoseVector pl_poses;
  pose_lib::gp4ps(sample_positions, sample_rays, sample_points3D, &pl_poses);
//  std::cout << " done " << std::endl;
//  std::cout << pl_poses.size() << std::endl;
  if (pl_poses.empty()) return 0;
  for (const pose_lib::CameraPose& pose : pl_poses) {
//    std::cout << pose.alpha << " " << sample.size() << std::endl;
//    std::cout << pose.R << std::endl;
//    std::cout << pose.t.transpose() << std::endl;
    const double kError = EvaluateModelOnPoint(pose, sample[4]);
//    std::cout << kError << std::endl;
    if (kError < squared_inlier_threshold_) {
      poses->push_back(pose);
      break;
    }
  }

  return static_cast<int>(poses->size());
}

int GeneralizedPoseSolverGP4Ps::NonMinimalSolver(
    const std::vector<int>& sample, GenCamPose* pose) const {
  // For now, just call the minimal solver.
  GenCamPoses poses;
  MinimalSolver(sample, &poses);
  if (poses.empty()) return 0;
  
  *pose = poses[0];
  return 1;
}

double GeneralizedPoseSolverGP4Ps::EvaluateModelOnPoint(
    const GenCamPose& pose, int i) const {
  // Transforms into the coordinate system of the rig.
  Vector3d p_r = (pose.R * points3D_[i] + pose.t) / pose.alpha;
  
  // Transforms into the coordinate system of the camera.
  const Camera& cam = rig_.cameras[camera_indices_[i]];
  Vector3d p_c = cam.R * p_r + cam.t;

  // Check whether point projects behind the camera.
  if (p_c[2] < 0.0) return std::numeric_limits<double>::max();

  Eigen::Vector2d p_2d = p_c.head<2>() / p_c[2];
  p_2d[0] *= cam.focal_x;
  p_2d[1] *= cam.focal_y;

  return (p_2d - points2D_[i]).squaredNorm();
}

void GeneralizedPoseSolverGP4Ps::LeastSquares(
    const std::vector<int>& sample, GenCamPose* pose) const {
  // For now, do nothing and simply return the pose, i.e., no refinement done.
}

}  // namespace multi_camera_pose
