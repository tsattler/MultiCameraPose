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

#ifndef MULTICAMERAPOSE_SRC_GENERALIZED_POSE_SOLVER_GP3P_H_
#define MULTICAMERAPOSE_SRC_GENERALIZED_POSE_SOLVER_GP3P_H_

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

#include "types.h"

namespace multi_camera_pose {

// Implements a camera pose solver for generalized calibrated cameras. Uses
// the PoseLib the GP3P solver as a minimal solver and Ceres for non-linear
// optimization. The non-minimal solver is simply implemented as first calling
// the minimal solver followed by non-linear optimization.
// To avoid returning multiple models by the minimal solver, a fourth point is
// used to reject models. The output is the generalized absolute pose of each
// pose in the multi-camera rig in world coordinates as defined in types.h.
class GeneralizedPoseSolverGP3P {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  // The input to the constructor is a set of cameras that define a multi-camera
  // rig, a set of 2D keypoint positions, the corresponding 3D points, and the
  // camera index for each match. The 2D points are expected to be centered
  // around the principal point (i.e., the principal has already been
  // subtracted) and to be undistorted.
  // In addition, the squared inlier threshold used by *SAC is required as
  // input. It is used to pick at most one of the solutions created by GP3P.
  GeneralizedPoseSolverGP3P(
      const MultiCameraRig& rig, const double squared_inlier_threshold,
      const Points2D& points2D, const Points3D& points3D,
      const std::vector<int>& camera_indices);

  inline int min_sample_size() const { return 4; }

  inline int non_minimal_sample_size() const { return 4; }

  inline int num_data() const { return num_data_; }

  int MinimalSolver(const std::vector<int>& sample,
                    GenCamPoses* poses) const;

  // Returns 0 if no model could be estimated and 1 otherwise.
  int NonMinimalSolver(const std::vector<int>& sample,
                       GenCamPose* pose) const;

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

#endif  // MULTICAMERAPOSE_SRC_GENERALIZED_POSE_SOLVER_GP3P_H_
