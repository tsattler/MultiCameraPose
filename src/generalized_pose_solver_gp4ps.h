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

#ifndef MULTICAMERAPOSE_SRC_GENERALIZED_POSE_SOLVER_GP4PS_H_
#define MULTICAMERAPOSE_SRC_GENERALIZED_POSE_SOLVER_GP4PS_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <random>
#include <vector>

#include <Eigen/Core>
#include <Eigen/StdVector>

#include <PoseLib/types.h>

namespace multi_camera_pose {
using Eigen::aligned_allocator;
using Eigen::Matrix3d;
using Eigen::Vector2d;
using Eigen::Vector3d;
  
// An generalized absolute pose is rotation and translation from world to camera
// coordinates as well as scaling factor alpha that handles the fact that the
// pose is only defined up to scale. The transformation is defined as
//  alpha * p_i + lambda_i * x_i = R * X_i + t
// where R and t are the rotation and translation respectively. p_i is the
// camera center of the i-th observation and x_i is the corresponding
// observation direction as a unit vector. Both are defined in the coordinate
// system of the multi-camera rig. X_i is the world position of the
// corresponding 3D point.
typedef pose_lib::CameraPose GenCamPose;
typedef std::vector<GenCamPose, aligned_allocator<CameraPose>> GenCamPoses;

// Defines a camera with intrinsics and extrinsics.
// We use a simple camera model, where we assume that all 3D points have
// already been centered around the principal point and distortion has been
// removed from the image.
struct Camera {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Eigen::Matrix3d R;
  Eigen::Vector3d t;
  Eigen::Vector3d c;
  double focal_x;
  double focal_y;
};
typedef std::vector<Camera, aligned_allocator<Camera>> Cameras;

// Defines a multi-camera rig as a set of individual cameras. The pose of each
// camera is given by the relative pose of the camera wrt. to a local camera
// coordinate system.
struct MultiCameraRig {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Cameras cameras;
};

typedef std::vector<Vector2d, aligned_allocator<Vector2d>> Points2D;
typedef std::vector<Vector3d, aligned_allocator<Vector3d>> Points3D;
typedef std::vector<Vector3d, aligned_allocator<Vector3d>> ViewingRays;

typedef std::vector<Vector3d, aligned_allocator<Vector3d>> CameraPositions;
typedef std::vector<Matrix3d, aligned_allocator<Matrix3d>> CameraRotations;

// Implements a camera pose solver for generalized calibrated cameras. Uses
// the PoseLib the GP4P+s solver as a minimal solver and Ceres for non-linear
// optimization. The non-minimal solver is simply implemented as first calling
// the minimal solver followed by non-linear optimization.
// To avoid returning multiple models by the minimal solver, a fifth point is
// used to pick at most one model. The output is the generalized absolute pose
// of each pose in the multi-camera rig in world coordinates as defined above.
class GeneralizedPoseSolverGP4Ps {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  // The input to the constructor is a set of cameras that define a multi-camera
  // rig, a set of 2D keypoint positions, the corresponding 3D points, and the
  // camera index for each match. The 2D points are expected to be centered
  // around the principal point (i.e., the principal has already been
  // subtracted) and to be undistorted.
  // In addition, the squared inlier threshold used by *SAC is required as
  // input. It is used to pick at most one of the solutions created by GP4Ps.
  GeneralizedPoseSolverGP4Ps(
      const MultiCameraRig& rig, const double squared_inlier_threshold,
      const Points2D& points2D, const Points3D& points3D,
      const std::vector<int>& camera_indices);

  inline int min_sample_size() const { return 5; }

  inline int non_minimal_sample_size() const { return 5; }

  inline int num_data() const { return num_data_; }

  int MinimalSolver(const std::vector<int>& sample,
                    GenCamPoses* poses) const;

  // Returns 0 if no model could be estimated and 1 otherwise.
  int NonMinimalSolver(const std::vector<int>& sample,
                       GenCamPoses* pose) const;

  // Evaluates the model on the i-th data point.
  double EvaluateModelOnPoint(const GenCamPose& pose, int i) const;

  // Linear least squares solver. Calls NonMinimalSolver.
  void LeastSquares(const std::vector<int>& sample, GenCamPose* pose) const;

 protected:
  MultiCameraRig rig_;
  double squared_inlier_threshold_;
  // Vector holding the 2D point positions.
  Points2D points2D_;
  // Vector holding the corresponding 3D point positions.
  Points3D points3D_;
  // Vector holding the viewing rays, defined in the rig coordinate system.
  ViewingRays rays_;
  // For each match, stores the corresponding camera in the multi-camera rig.
  std::vector<int> camera_indices_;

  int num_data_;
  int num_cameras_;
};

}  // namespace multi_camera_pose

#endif  // MULTICAMERAPOSE_SRC_GENERALIZED_POSE_SOLVER_GP4PS_H_
