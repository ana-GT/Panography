# CMake for Assignment 3 
cmake_minimum_required( VERSION 2.8 )

project( Assignment3 )

include_directories( . /usr/local/include
  /usr/local/include/opencv /usr/include/eigen3 )
link_directories( /usr/local/lib )
set( SRC_FILES getCameraMatrix.cpp )
set( OpenCV_LIBS opencv_core
  opencv_ml
  opencv_video
  opencv_objdetect
  opencv_photo
  opencv_legacy
  opencv_stitching
  opencv_videostab
  opencv_highgui
  opencv_imgproc
  opencv_calib3d
  opencv_contrib
  opencv_nonfree
  opencv_features2d
  opencv_flann )

#add_executable( getCameraMatrix ${SRC_FILES} )
#target_link_libraries( getCameraMatrix ${OpenCV_LIBS} )  

# OpenCV code
add_executable( calibration calibration.cpp )
target_link_libraries( calibration ${OpenCV_LIBS} )

add_executable( imagelist_creator imagelist_creator.cpp )
target_link_libraries( imagelist_creator ${OpenCV_LIBS} )

add_executable( camera_calibration camera_calibration.cpp )
target_link_libraries( camera_calibration ${OpenCV_LIBS} )

add_executable( Assignment3_Ransac Assignment3_Ransac.cpp Ransac.cpp )
target_link_libraries( Assignment3_Ransac ${OpenCV_LIBS} )
