/**
 * @file Assignment3_Ransac.cpp
 */

#include "Ransac.h"
#include <stdio.h>

char rgbWindow[50] = "RGB Window";
char matchWindow[50] = "Match Window";
char keypointWindow[50] = "Keypoint Window";
char warpedWindow[50] = "Warped Window";
/**
 * @function main
 */
int main( int argc, char* argv[] )  {
  
  Ransac ransac( argv[1] );
  ransac.matchAllFrames();
  ransac.getAllTransforms();
  ransac.getTransformedImages();
  ransac.saveTransformedImages();

  cv::namedWindow( rgbWindow, CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO );
  cv::namedWindow( matchWindow, CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO );
  cv::namedWindow( warpedWindow, CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO );
  cv::imshow( rgbWindow, ransac.getRgb( 0 ) ); 
    
  int key;
  int ind = 0;
  while( true ) {
    key = cv::waitKey(30);

    // If a key is pressed
    if( key != -1 ) {
      
      if( key == 27 ) {
	break;
      }
      
      if( key == 'a' ) {
	ind++;
	ind = ind % ransac.getNumFrames();
	printf("Showing image %d out of %d \n", ind, ransac.getNumFrames() );
	cv::imshow( rgbWindow, ransac.getRgb( ind ) );
	cv::imshow( warpedWindow, ransac.getWarped( ind ) );
	cv::imshow( matchWindow, ransac.getMatchesDraw( ind, 0 ) );
	//cv::imshow( matchWindow, ransac.Ransac_Homography2D( ind, 0 ) );
      }
   
    } // end if( key != -1 )
  }

  return 0;
}
