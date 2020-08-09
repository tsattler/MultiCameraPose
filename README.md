# Multi Camera Pose Estimation
Uses pre-computed 2D-3D matches together with the known relative poses between multiple images for multi-camera pose estimation.

**More details, including licensing terms, will be released soon**

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

## Citing
When using this software for publications, please cite
```
@inproceedings{wald2020,
  title={{Beyond Controlled Environments: 3D Camera Re-Localization in Changing Indoor Scenes}},
  author={Wald, Johanna and Sattler, Torsten and Golodetz, Stuart and Cavallari, Tommaso and Tombari, Federico},
  booktitle = {Proceedings IEEE European Conference on Computer Vision (ECCV)},
  year = {2020}
}
```

When using the pose estimator based on the GP4Ps solver, please cite
```
@InProceedings{Kukelova2016CVPR,
  author = {Kukelova, Zuzana and Heller, Jan and Fitzgibbon, Andrew},
  title = {Efficient Intersection of Three Quadrics and Applications in Computer Vision},
  booktitle = {The IEEE Conference on Computer Vision and Pattern Recognition (CVPR)},
  year = {2016}
}
```
