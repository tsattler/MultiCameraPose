cmake_minimum_required (VERSION 3.0)

if (EXISTS "${CMAKE_SOURCE_DIR}/cmake")
  set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
endif ()

# Same compiler options as PoseLib
if(NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

#if(MSVC)	
#	target_compile_options(poselib PRIVATE /bigobj /fp:fast)
#	target_compile_options(benchmark PRIVATE /bigobj /fp:fast)
#else()
#	target_compile_options(poselib PRIVATE -march=native -Wall -Werror -Wno-#sign-compare -Wno-unused-variable -ffast-math)
#	target_compile_options(benchmark PRIVATE -march=native -Wall -Werror -#Wno-sign-compare -Wno-unused-variable -ffast-math)
#endif()

add_subdirectory (src)