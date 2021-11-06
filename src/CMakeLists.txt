cmake_minimum_required (VERSION 3.0)

#if (EXISTS "${CMAKE_SOURCE_DIR}/cmake")
#  set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
#endif ()

SET(CMAKE_CXX_FLAGS "-std=gnu++11")

find_package (Eigen3 REQUIRED)

find_package (Ceres REQUIRED)

add_definitions (-march=native)

include_directories (
  ${CMAKE_SOURCE_DIR}
  ${EIGEN_INCLUDE_DIRS}
  ${CERES_INCLUDE_DIRS}
  ${POSELIB_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/ext/RansacLib
  ${CMAKE_SOURCE_DIR}/ext/PoseLib
  ${CMAKE_SOURCE_DIR}/ext
)

add_executable (multi_camera_pose multi_camera_pose.cc generalized_pose_solver_gp4ps.cc generalized_pose_solver_gp4ps.h types.h common.h ${CMAKE_SOURCE_DIR}/ext/PoseLib/PoseLib/gp4ps.cc ${CMAKE_SOURCE_DIR}/ext/PoseLib/PoseLib/misc/re3q3.cc ${CMAKE_SOURCE_DIR}/ext/PoseLib/PoseLib/misc/univariate.cc)
target_link_libraries (multi_camera_pose ${CERES_LIBRARIES})

add_executable (fixed_rig_camera_pose fixed_rig_camera_pose.cc generalized_pose_solver_gp3p.cc generalized_pose_solver_gp3p.h types.h common.h ${CMAKE_SOURCE_DIR}/ext/PoseLib/PoseLib/gp3p.cc ${CMAKE_SOURCE_DIR}/ext/PoseLib/PoseLib/misc/re3q3.cc)
target_link_libraries (fixed_rig_camera_pose ${CERES_LIBRARIES})

if(MSVC)	
  target_compile_options(multi_camera_pose PRIVATE /bigobj /fp:fast)
  target_compile_options(fixed_rig_camera_pose PRIVATE /bigobj /fp:fast)
else()
  target_compile_options(multi_camera_pose PRIVATE -march=native -Wall -Werror -Wno-sign-compare -Wno-unused-variable -ffast-math)
  target_compile_options(fixed_rig_camera_pose PRIVATE -march=native -Wall -Werror -Wno-sign-compare -Wno-unused-variable -ffast-math)
endif()
