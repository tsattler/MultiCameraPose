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
//     * Neither the name of Torsten Sattler nor the
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

#ifndef MULTICAMERAPOSE_SRC_COMMON_H_
#define MULTICAMERAPOSE_SRC_COMMON_H_

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/StdVector>

#include "types.h"

namespace multi_camera_pose {

template <typename T>
double ComputeMedian(std::vector<T>* data) {
  T mean = static_cast<T>(0.0);
  for (size_t i = 0; i < data->size(); ++i) {
    mean += (*data)[i];
  }
  mean /= static_cast<T>(data->size());
  std::cout << " mean : " << mean << std::endl;

  std::sort(data->begin(), data->end());
  if (data->size() % 2u == 1u) {
    return static_cast<double>((*data)[data->size() / 2]);
  } else {
    double a = static_cast<double>((*data)[data->size() / 2 - 1]);
    double b = static_cast<double>((*data)[data->size() / 2]);
    return (a + b) * 0.5;
  }
}



// Distortion function, implementing various camera models.
void Distortion(const QueryData& camera, const double u, const double v,
                double* du, double* dv) {
  if (camera.camera_type.compare("SIMPLE_RADIAL") == 0) {
    const double kR2 = u * u + v * v;
    const double kRadial = camera.radial[0] * kR2;
    *du = u * kRadial;
    *dv = v * kRadial;
  } else if (camera.camera_type.compare("RADIAL") == 0) {
    const double kR2 = u * u + v * v;
    const double kRadial = camera.radial[0] * kR2 + camera.radial[1] * kR2 * kR2;
    *du = u * kRadial;
    *dv = v * kRadial;
  } else if (camera.camera_type.compare("BROWN_3_PARAMS") == 0) {
    const double kR2 = u * u + v * v;
    const double kRadial = camera.radial[0] * kR2 + camera.radial[1] * kR2 * kR2 + camera.radial[2] * kR2 * kR2 * kR2;
    *du = u * kRadial;
    *dv = v * kRadial;
  } else if (camera.camera_type.compare("PINHOLE") == 0) {
    *du = 0;
    *dv = 0;    
  } else {
    std::cerr << " ERROR: Distortion function for camera model "
              << camera.camera_type << " not yet implemented" << std::endl;
  }
}

// Undistortion code taken from Colmap.
// TODO(sattler): Replace with own code before release.
// Assumes that principal point has been subtracted and that the coordinates
// have been divided by the focal length.
void IterativeUndistortion(const QueryData& camera, double* u, double* v) {
  // Parameters for Newton iteration using numerical differentiation with
  // central differences, 100 iterations should be enough even for complex
  // camera models with higher order terms.
  const size_t kNumIterations = 100;
  const double kMaxStepNorm = 1e-10;
  const double kRelStepSize = 1e-6;

  Eigen::Matrix2d J;
  const Eigen::Vector2d x0(*u, *v);
  Eigen::Vector2d x(*u, *v);
  Eigen::Vector2d dx;
  Eigen::Vector2d dx_0b;
  Eigen::Vector2d dx_0f;
  Eigen::Vector2d dx_1b;
  Eigen::Vector2d dx_1f;

  for (size_t i = 0; i < kNumIterations; ++i) {
    const double step0 = std::max(std::numeric_limits<double>::epsilon(),
                                  std::abs(kRelStepSize * x(0)));
    const double step1 = std::max(std::numeric_limits<double>::epsilon(),
                                  std::abs(kRelStepSize * x(1)));
    Distortion(camera, x(0), x(1), &dx(0), &dx(1));
    Distortion(camera, x(0) - step0, x(1), &dx_0b(0), &dx_0b(1));
    Distortion(camera, x(0) + step0, x(1), &dx_0f(0), &dx_0f(1));
    Distortion(camera, x(0), x(1) - step1, &dx_1b(0), &dx_1b(1));
    Distortion(camera, x(0), x(1) + step1, &dx_1f(0), &dx_1f(1));
    J(0, 0) = 1 + (dx_0f(0) - dx_0b(0)) / (2 * step0);
    J(0, 1) = (dx_1f(0) - dx_1b(0)) / (2 * step1);
    J(1, 0) = (dx_0f(1) - dx_0b(1)) / (2 * step0);
    J(1, 1) = 1 + (dx_1f(1) - dx_1b(1)) / (2 * step1);
    const Eigen::Vector2d step_x = J.inverse() * (x + dx - x0);
    x -= step_x;
    if (step_x.squaredNorm() < kMaxStepNorm) {
      break;
    }
  }

  *u = x(0);
  *v = x(1);
}

// Loads the list of query images together with their intrinsics and extrinsics.
bool LoadListIntrinsicsAndExtrinsics(const std::string& filename,
                                     Queries* query_images) {
  std::ifstream ifs(filename.c_str(), std::ios::in);
  if (!ifs.is_open()) {
    std::cerr << " ERROR: Cannot read the image list from " << filename
              << std::endl;
    return false;
  }
  std::string line;

  query_images->clear();

  while (std::getline(ifs, line)) {
    std::stringstream s_stream(line);

    std::string camera_type;

    QueryData q;
    s_stream >> q.name >> camera_type >> q.width >> q.height;
    q.camera_type = camera_type;
    if (camera_type.compare("SIMPLE_RADIAL") == 0 ||
        camera_type.compare("VSFM") == 0) {
      q.radial.resize(1);
      s_stream >> q.focal_x >> q.c_x >> q.c_y >> q.radial[0];
      q.focal_y = q.focal_x;
    } else if (camera_type.compare("PINHOLE") == 0) {
      q.radial.clear();
      s_stream >> q.focal_x >> q.focal_y >> q.c_x >> q.c_y;
    } else if (camera_type.compare("OPENCV") == 0) {
      // The OPENCV camera model used in Colmap (see
      // https://github.com/colmap/colmap/blob/master/src/base/camera_models.h
      // for details).
      q.radial.resize(4);
      s_stream >> q.focal_x >> q.focal_y >> q.c_x >> q.c_y >> q.radial[0]
               >> q.radial[1] >> q.radial[2] >> q.radial[3];
    } else if (camera_type.compare("VSFM") == 0) {
      q.radial.resize(1);
      s_stream >> q.focal_x >> q.c_x >> q.c_y >> q.radial[0];
      q.focal_y = q.focal_x;
    } else if (camera_type.compare("BROWN_3_PARAMS") == 0) {
      q.radial.resize(3);
      s_stream >> q.focal_x >> q.focal_y >> q.c_x >> q.c_y >> q.radial[0]
               >> q.radial[1] >> q.radial[2];
    } else if (camera_type.compare("RADIAL") == 0) {
      q.radial.resize(2);
      s_stream >> q.focal_x >> q.c_x >> q.c_y >> q.radial[0] >> q.radial[1];
      q.focal_y = q.focal_x;
    } else {
      std::cerr << " ERROR: Unknown camera model " << camera_type << " for "
                << "image " << q.name << std::endl;
      return false;
    }
    s_stream >> q.q.w() >> q.q.x() >> q.q.y() >> q.q.z() >> q.c[0] >> q.c[1] >>
        q.c[2];
    
    query_images->push_back(q);
  }

  ifs.close();

  return true;
}

// Loads the 2D-3D matches found for that image from a text file.
bool LoadMatches(const std::string& filename, bool invert_Y_Z,
                 Points2D* points2D, Points3D* points3D) {
  points2D->clear();
  points3D->clear();

  std::ifstream ifs(filename.c_str(), std::ios::in);
  if (!ifs.is_open()) {
    std::cerr << " ERROR: Cannot read the matches from " << filename
              << std::endl;
    return false;
  }
  std::string line;

  while (std::getline(ifs, line)) {
    std::stringstream s_stream(line);

    Eigen::Vector2d p2D;
    Eigen::Vector3d p3D;
    s_stream >> p2D[0] >> p2D[1] >> p3D[0] >> p3D[1] >> p3D[2];

    if (invert_Y_Z) {
      p3D[1] *= -1.0;
      p3D[2] *= -1.0;
    }

    points2D->push_back(p2D);
    points3D->push_back(p3D);
  }

  return true;
}

bool AssembleMultiCameraRig(const Queries& queries,
                            const std::vector<int>& indices,
                            MultiCameraRig* rig) {
  const int kNumCams = static_cast<int>(indices.size());
  if (kNumCams == 0) {
    std::cerr << "ERROR: Not enough cameras for a multi-camera system"
              << std::endl;
    return false;
  }
  
  // Uses the first camera as the reference camera.
  rig->cameras.resize(kNumCams);
  rig->cameras[0].R = Eigen::Matrix3d::Identity();
  rig->cameras[0].t = Eigen::Vector3d::Zero();
  rig->cameras[0].c = Eigen::Vector3d::Zero();
  rig->cameras[0].focal_x = queries[indices[0]].focal_x;
  rig->cameras[0].focal_y = queries[indices[0]].focal_y;
  
  Eigen::Matrix3d R0(queries[indices[0]].q);
  Eigen::Vector3d t0 = -R0 * queries[indices[0]].c;
  
  for (int i = 1; i < kNumCams; ++i) {
    Eigen::Matrix3d Ri(queries[indices[i]].q);
    Eigen::Vector3d ti = -Ri * queries[indices[i]].c;
    
    Eigen::Matrix3d R = Ri * R0.transpose();
    Eigen::Vector3d t = ti - R * t0;
    Eigen::Vector3d c = -R.transpose() * t;
    
    rig->cameras[i].R = R;
    rig->cameras[i].t = t;
    rig->cameras[i].c = c;
    
    rig->cameras[i].focal_x = queries[indices[i]].focal_x;
    rig->cameras[i].focal_y = queries[indices[i]].focal_y;
  }
  
  return true;
}

}  // namespace multi_camera_pose

#endif  // MULTICAMERAPOSE_SRC_COMMON_H_
