cmake_minimum_required (VERSION 3.0)

if (EXISTS "${CMAKE_SOURCE_DIR}/cmake")
  set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
endif ()

SET(CMAKE_CXX_FLAGS "-std=gnu++11")

find_package (Eigen3 REQUIRED)

find_package (PoseLib REQUIRED)

find_package (Ceres REQUIRED)

add_definitions (-march=native)

include_directories (
  ${CMAKE_SOURCE_DIR}
  ${EIGEN_INCLUDE_DIRS}
  ${CERES_INCLUDE_DIRS}
  ${POSELIB_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/ext/RansacLib
)

#add_executable (multi_camera_pose generalized_pose_solver_gp4ps.cc generalized_pose_solver_gp4ps.h ${POSELIB_INCLUDE_DIR}/PoseLib/types.h ${POSELIB_INCLUDE_DIR}/PoseLib/gp4ps.h ${POSELIB_INCLUDE_DIR}/PoseLib/gp4ps.cc ${POSELIB_INCLUDE_DIR}/PoseLib/re3q3.h ${POSELIB_INCLUDE_DIR}/PoseLib/re3q3.cc)
#target_link_libraries (multi_camera_pose ${POSELIB_LIBRARY}) # ${CERES_LIBRARIES})

add_executable (multi_camera_pose multi_camera_pose.cc generalized_pose_solver_gp4ps.cc generalized_pose_solver_gp4ps.h)
target_link_libraries (multi_camera_pose ${POSELIB_LIBRARY} ${CERES_LIBRARIES})

if(MSVC)	
  target_compile_options(multi_camera_pose PRIVATE /bigobj /fp:fast)
else()
  target_compile_options(multi_camera_pose PRIVATE -march=native -Wall -Werror -Wno-sign-compare -Wno-unused-variable -ffast-math)
endif()