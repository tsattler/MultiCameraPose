# - Try to find the PoseLib libray
# Once done this will define
#  POSELIB_FOUND - system has the PoseLib library
#  POSELIB_INCLUDE_DIR - the PoseLib include directory

if (POSELIB_INCLUDE_DIR)
  # in cache already
  set(POSELIB_FIND_QUIETLY TRUE)
else ()
  find_path(POSELIB_INCLUDE_DIR PoseLib/gp4ps.h
      PATHS
      "/Users/sattlert/code/PoseLib"
    )

  if(POSELIB_INCLUDE_DIR)
    set(POSELIB_FOUND TRUE)
  else(POSELIB_INCLUDE_DIR)
    set(POSELIB_FOUND FALSE)
  endif(POSELIB_INCLUDE_DIR)
endif()

