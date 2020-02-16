# Multi Camera Pose Estimation
Uses pre-computed 2D-3D matches together with the known relative poses between multiple images for multi-camera pose estimation.

## Setup
This project has the following dependencies:
* [RansacLib](https://github.com/tsattler/RansacLib)
* [PoseLib](https://github.com/vlarsson/PoseLib)
* Eigen3
* Ceres

RansacLib is included as a submodule. After cloning the repository, run
```
git submodule update --init --recursive
```
