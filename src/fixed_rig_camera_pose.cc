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

#include <RansacLib/ransac.h>

#include "common.h"
#include "generalized_pose_solver_gp3p.h"
#include "types.h"


int main(int argc, char** argv) {
  using ransac_lib::LocallyOptimizedMSAC;
  using namespace multi_camera_pose;

  std::cout << " usage: " << argv[0] << " images_with_intrinsics outfile "
            << "inlier_threshold num_lo_steps invert_Y_Z points_centered "
            << " undistortion_needed sequence_length [match-file postfix]"
            << std::endl;
  if (argc < 9) return -1;

  bool invert_Y_Z = static_cast<bool>(atoi(argv[5]));
  bool points_centered = static_cast<bool>(atoi(argv[6]));
  bool undistortion_needed = static_cast<bool>(atoi(argv[7]));
  int sequence_length = atoi(argv[8]);

  Queries query_data;
  std::string list(argv[1]);

  if (!LoadListIntrinsicsAndExtrinsics(list, &query_data)) {
    std::cerr << " ERROR: Could not read the data from " << list << std::endl;
    return -1;
  }
  const int kNumQuery = static_cast<int>(query_data.size());
  std::cout << " Found " << kNumQuery << " query images " << std::endl;

  std::ofstream ofs(argv[2], std::ios::out);
  if (!ofs.is_open()) {
    std::cerr << " ERROR: Cannot write to " << argv[2] << std::endl;
    return -1;
  }

  std::string matchfile_postfix = ".individual_datasets.matches.txt";
  if (argc >= 10) {
    matchfile_postfix = std::string(argv[9]);
  }

  std::vector<double> orientation_error(kNumQuery,
                                        std::numeric_limits<double>::max());
  std::vector<double> position_error(kNumQuery,
                                     std::numeric_limits<double>::max());

  std::vector<double> position_thresholds = {0.05, 0.03, 0.02, 0.01};
  std::vector<double> orientation_thresholds = {5.0, 3.0, 2.0, 1.0};
  const int kNumThresholds = 4;
  std::vector<int> num_poses_within_threshold(kNumThresholds, 0);

  double mean_ransac_time = 0.0;

  int num_better_reprojection_error_than_gt = 0;
  int num_reproj_tested = 0;
  
  std::vector<PoseWInliers, Eigen::aligned_allocator<PoseWInliers>> best_poses(kNumQuery);
  for (int i = 0; i < kNumQuery; ++i) {
    best_poses[i].num_inliers = 0;
    best_poses[i].score = std::numeric_limits<double>::max();
  }
  
  for (int i = 0; i < kNumQuery; i += sequence_length) {
    std::cout << std::endl << std::endl;
    
    // Assembles the rig.
    int cams_end = std::min(i + sequence_length, kNumQuery);
    std::vector<int> rig_indices(cams_end - i);
    std::iota(rig_indices.begin(), rig_indices.end(), i);
    MultiCameraRig rig;
    std::cout << " Assembling camera rig" << std::endl;
    if (!AssembleMultiCameraRig(query_data, rig_indices, &rig)) continue;
    
    std::cout << " done" << std::endl;
    
    // Obtains the matching data.
    Points2D points2D_rig;
    Points3D points3D_rig;
    std::vector<int> camera_ids;
    for (int j = i; j < cams_end; ++j) {
      std::string matchfile(query_data[j].name);
      matchfile.append(matchfile_postfix);
      
      Points2D points2D;
      Points3D points3D;
      if (!LoadMatches(matchfile, invert_Y_Z, &points2D, &points3D)) {
        std::cerr << "  ERROR: Could not load matches from " << matchfile
                  << std::endl;
        continue;
      }
      if (points2D.empty()) continue;
      
      std::vector<int> id(points2D.size(), j - i);
      
      points2D_rig.insert(points2D_rig.end(), points2D.begin(), points2D.end());
      points3D_rig.insert(points3D_rig.end(), points3D.begin(), points3D.end());
      camera_ids.insert(camera_ids.end(), id.begin(), id.end());
    }
    
    const int kNumMatches = static_cast<int>(points2D_rig.size());
    std::cout << " Found " << kNumMatches << " matches" << std::endl;
    if (kNumMatches <= 5) {
      std::cout << " Found only " << kNumMatches << " matches for query image "
                << query_data[i].name << " -> skipping image" << std::endl;

      continue;
    }
    

    if (!points_centered) {
      for (int j = 0; j < kNumMatches; ++j) {
        points2D_rig[j][0] -= query_data[camera_ids[j] + i].c_x;
        points2D_rig[j][1] -= query_data[camera_ids[j] + i].c_y;
      }
    }
    
    if (undistortion_needed) {
      for (int j = 0; j < kNumMatches; ++j) {
        double u = points2D_rig[j][0] / query_data[camera_ids[j] + i].focal_x;
        double v = points2D_rig[j][1] / query_data[camera_ids[j] + i].focal_y;
        IterativeUndistortion(query_data[camera_ids[j] + i], &u, &v);
        points2D_rig[j][0] = u * query_data[camera_ids[j] + i].focal_x;
        points2D_rig[j][1] = v * query_data[camera_ids[j] + i].focal_y;
      }
    }

    ransac_lib::LORansacOptions options;
    options.min_num_iterations_ = 200u;
    options.max_num_iterations_ = 10000u;
    options.min_sample_multiplicator_ = 7;
    //    options.num_lsq_iterations_ = 0;
    //    options.num_lo_steps_ = 0;
    options.num_lsq_iterations_ = 4;
    options.num_lo_steps_ = atoi(argv[4]);
    options.lo_starting_iterations_ = 60;
    options.final_least_squares_ = true;
    //    options.threshold_multiplier_ = 2.0;

    std::random_device rand_dev;
    options.random_seed_ = rand_dev();

    const double kInThreshPX = static_cast<double>(atof(argv[3]));
    options.squared_inlier_threshold_ = kInThreshPX * kInThreshPX;

    GeneralizedPoseSolverGP3P solver(rig, kInThreshPX * kInThreshPX,
                                     points2D_rig, points3D_rig, camera_ids);

    LocallyOptimizedMSAC<GenCamPose, GenCamPoses, GeneralizedPoseSolverGP3P> lomsac;
    ransac_lib::RansacStatistics ransac_stats;
    GenCamPose best_model;

    std::cout << "   " << query_data[i].name << " : running LO-MSAC on "
              << kNumMatches << " matches " << std::endl;
    auto ransac_start = std::chrono::system_clock::now();

    int num_ransac_inliers =
        lomsac.EstimateModel(options, solver, &best_model, &ransac_stats);
    auto ransac_end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = ransac_end - ransac_start;
    mean_ransac_time += elapsed_seconds.count();
    std::cout << "   ... LOMSAC found " << num_ransac_inliers << " inliers in "
              << ransac_stats.num_iterations
              << " iterations with an inlier ratio of "
              << ransac_stats.inlier_ratio << std::endl;
    std::cout << "   ... LOMSAC took " << elapsed_seconds.count() << " s"
              << std::endl;
    std::cout << "   ... LOMSAC executed " << ransac_stats.number_lo_iterations
              << " local optimization stages" << std::endl;

    if (num_ransac_inliers < 5) continue;
    
    std::cout << best_model.alpha << std::endl;
    
    // Updates the poses of all cameras in the rig.
    for (int j = i; j < cams_end; ++j) {
      if (ransac_stats.best_model_score < best_poses[j].score) {
        best_poses[j].num_inliers = num_ransac_inliers;
        best_poses[j].score = ransac_stats.best_model_score;
        Eigen::Matrix3d R = rig.cameras[j - i].R * best_model.R;
        Eigen::Vector3d t = rig.cameras[j - i].R * best_model.t + rig.cameras[j - i].t;
        Eigen::Vector3d c = -R.transpose() * t;
        best_poses[j].R = R;
        best_poses[j].t = t;
        best_poses[j].c = c;
      }
    }
  }

  for (int i = 0; i < kNumQuery; ++i) {
//     if (best_poses[i].num_inliers <= 5) continue;
    
    Eigen::Matrix3d R = best_poses[i].R;
    Eigen::Vector3d t = best_poses[i].t;
    Eigen::Vector3d c = best_poses[i].c;
    Eigen::Quaterniond q(R);
    q.normalize();

    // Measures the pose error.
    double c_error = (c - query_data[i].c).norm();
    Eigen::Matrix3d R1 = R.transpose();
    Eigen::Matrix3d R2(query_data[i].q);
    Eigen::AngleAxisd aax(R1 * R2);
    double q_error = aax.angle() * 180.0 / M_PI;
    orientation_error[i] = q_error;
    position_error[i] = c_error;

    for (int k = 0; k < kNumThresholds; ++k) {
      if (c_error <= position_thresholds[k] &&
          q_error <= orientation_thresholds[k]) {
        num_poses_within_threshold[k] += 1;
      }
    }

    ofs << query_data[i].name << " " << q.w() << " " << q.x() << " " << q.y()
        << " " << q.z() << " " << t[0] << " " << t[1] << " " << t[2]
        << std::endl;
  }

  std::sort(orientation_error.begin(), orientation_error.end());
  std::sort(position_error.begin(), position_error.end());

  std::cout << std::endl
            << " Mean RANSAC time: "
            << mean_ransac_time / static_cast<double>(kNumQuery) << " s "
            << std::endl;

  double median_pos = ComputeMedian<double>(&position_error);
  double median_rot = ComputeMedian<double>(&orientation_error);
  std::cout << " Median position error: " << median_pos << "m "
            << " Median orientation error: " << median_rot << " deg"
            << std::endl;
  for (int k = 0; k < kNumThresholds; ++k) {
    std::cout << " % images within " << position_thresholds[k] * 100.0
              << "cm and " << orientation_thresholds[k] << "deg: "
              << static_cast<double>(num_poses_within_threshold[k]) /
                     static_cast<double>(kNumQuery) * 100.0
              << std::endl;
  }
  ofs.close();
  return 0;
}
