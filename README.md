# Multi Camera Pose Estimation
Uses pre-computed 2D-3D matches together with the known relative poses between multiple images for multi-camera pose estimation.

## Setup
This project depends on two submodules:
* [PoseLib](https://github.com/vlarsson/PoseLib)
* [RansacLib](https://github.com/tsattler/RansacLib)

After cloning the repository, run
```
git submodule update --init --recursive
```
