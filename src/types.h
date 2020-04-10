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

#ifndef MULTICAMERAPOSE_SRC_TYPES_H_
#define MULTICAMERAPOSE_SRC_TYPES_H_

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
typedef std::vector<GenCamPose, aligned_allocator<GenCamPose>> GenCamPoses;

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
  
struct QueryData {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  std::string name;
  double c_x;
  double c_y;
  
  double focal_x;
  double focal_y;
  
  int width;
  int height;
  
  Eigen::Quaterniond q;
  Eigen::Vector3d c;
  
  std::string camera_type;
  
  std::vector<double> radial;
};

struct PoseWInliers {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Eigen::Matrix3d R;
  Eigen::Vector3d c;
  Eigen::Vector3d t;
  int num_inliers;
  double score;
};

typedef std::vector<QueryData, Eigen::aligned_allocator<QueryData>> Queries;

}  // namespace multi_camera_pose

#endif  // MULTICAMERAPOSE_SRC_TYPES_H_
