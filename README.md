# Multi Camera Pose Estimation
Code for estimating the absolute pose of a multi-camera system from a set of 2D-3D matches.

## Setup
This project has the following dependencies:
* [RansacLib](https://github.com/tsattler/RansacLib)
* [PoseLib](https://github.com/vlarsson/PoseLib)
* Eigen3
* Ceres

RansacLib and PoseLib are included as a submodule. After cloning the repository, run
```
git submodule update --init --recursive
```

## Running the executables
There are two executables" `fixed_rig_camera_pose` and `multi_camera_pose`. 
`fixed_rig_camera_pose` assumes that the absolute scale of the transformation between the images is known. 
`multi_camera_pose` does not require the scale to be known, e.g., when the poses are estimated by SLAM, but rather 
estimates the scale. Both executables define multi-camera rigs define from sequences of images and expect a list of 2D-3D 
matches to be given for each image in the multi-camera system. Both executables share the same command line parameters:
* `images_with_intrinsics` is the file name of a text file that contains image names, camera intrinsics, and camera extrinsics. Each line consists of the following information:
  * `image_name`: the name of the image.
  * `intrinsics`: the intrinsics of the image, consisting of a camera type and the parameters. We use Colmaps camera definition. Please see [Colmap's camera definitions](https://colmap.github.io/cameras.html). An example for this part is `SIMPLE_RADIAL 1024 1024 640.0 512 512 0.2`. 
  * The camera extrinsics in the form `qw qx qy qz cx cy cz`, where `qw qx qy qz` is a unit quaternion defining a rotation from world to camera coordinates for this image and `cx cy cz` is the position of the image in world coordinates. I.e., a point `Xw` in world coordinates is transformed into the local camera coordinate system of the image as `Xc = R * (Xw - c)`, where `R` is the rotation matrix defined by the quaternion and `c` is the position of the image in the world coordinate system. These poses can be defined in an arbitrary coordinate frame. The executables will automatically extract relative poses between the images.
* `outfile` is the file name of a text file into which the estimated poses will be written. For each image in `images_with_extrinsics`, a pose will be written (if a corresponding pose can be estimated) in a line of the output file. The format of that line is `image_name qw qx qy qz tx ty tz`. Here `image_name` is the name of the image (as specified in `images_with_extrsincis`, `qw qx qy qz` again is a unit quaternion defining the rotation from the coordinate system of the 3D points (see below) to the local camera coordinate system and `tx ty tz` is the corresponding translation. If `R` is the rotation matrix corresponding to the quaternion, then a point `Xw` in world coordinates is transformed into the local camera coordinate system as `Xc = R * Xw + t`, where `t` is the translation vector given by `tx ty tz`.
* `inlier_threshold`: the inlier threshold to be used in RANSAC, given in pixels.
* `num_lo_steps`: the number of local optimization steps performed in RANSAC whenever a new best minimal pose is found.
* `invert_Y_Z`: the 2D-3D matches for each image are read from text files, where each line has the format `x y X Y Z`. Here, `x y` defines a 2D keypoint. Set this variable to `1` if `x y` is given in a coordinate system where the y-axis is pointing upwards.
* `points_centered`: set to `1` if the 2D keypoint coordinates `x y` are already centered around the principal point. If set to `0`, the executables will center the keypoints before using them.
* `undistortion_needed`: set to `1` if the 2D keypoint coordinates need to be undistorted and to `0` if the keypoints that are read from the text files are already undistorted.
* `sequence_length`: both executables assume that the images specified in `images_with_intrinsics` are given in sequential order. If set to `k`, e.g., `3`, `fixed_rig_camera_pose` will use the first `k` images to define a multi-camera system and attempt to localize them jointly. It will then use the next `k` images to define the next multi-camera rig, etc. `multi_camera_pose` will use the first `k` images to define the first multi-camera system. It will then use images `2, ..., k+1` to define the next multi-camera system, then`3, ..., k + 2`, etc.
* `[match-file postfix]` (optional): the postfix of the match files, set to `.individual_datasets.matches.txt` per default. For a given image name `a.jpg` in `images_with_intrinsics`, both executables will attempt to load 2D-3D matches from a text file called `a.jpg.individual_datasets.matches.txt` (for the default value). The text file stores each 2D-3D match in a single line, with the format `x y X Y Z`. Here, `x y` is the 2D coordinate of the matching 2D point and `X Y Z` is the corresponding 3D point.

Simply calling one of the programs without parameters will give you a list of command line arguments.

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

When using the pose estimators based on the GP3P or GP4Ps solvers, please cite
```
@InProceedings{Kukelova2016CVPR,
  author = {Kukelova, Zuzana and Heller, Jan and Fitzgibbon, Andrew},
  title = {Efficient Intersection of Three Quadrics and Applications in Computer Vision},
  booktitle = {The IEEE Conference on Computer Vision and Pattern Recognition (CVPR)},
  year = {2016}
}
```
