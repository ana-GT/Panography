/**
 * @file Ransac.h
 * @author A. Huaman :D
 * @date  2012/10/07
 */

#ifndef __RANSAC_H__
#define __RANSAC_H__

//-- OpenCV headers
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
// For SURF
#include "opencv2/nonfree/features2d.hpp"

#include <vector>

/**
 * @class Ransac
 */
class Ransac {

 public:

  // Constructor and destructor
  Ransac( char *_filenames );
  ~Ransac();

  //-- RANSAC
  cv::Mat Ransac_Homography2D( int _ind1, int _ind2 );
  cv::Mat getModel_Homography2D( int _ind1,
				 int _ind2, 
				 std::vector<int> _indices );  
  float getErrorModel( int _ind1, int _ind2, int _ind,
		       const cv::Mat &_params );
  void getPointsFromMatch( int _ind1, int _ind2, 
			   const cv::DMatch &_match,
			   float &_x1, float &_y1, 
			   float &_x2, float &_y2 );
  std::vector<int> getRandomIndices( int _numSamples, 
				     int _totalSamples );
  bool isInSet( int _index, 
		std::vector<int> _currentIndices );
  int getRandom( int _max );
  cv::Mat applyHomography2DToImage( const cv::Mat &_homography2D,
				    const cv::Mat &_image );

  void getAllTransforms();
  void getTransformedImages();
  void saveTransformedImages();

  //-- Matching
  void matchAllFrames();
  cv::Mat getMatchesDraw( int _ind1, int _ind2  );
  cv::Mat getSomeMatchesDraw( int _ind1, int _ind2,
			      const std::vector<int> &matchesIndices );
  cv::Mat getKeypointsDraw( int _ind );

  //-- Data acquisition  
  bool readImagesData( char* _filenames, 
		       std::vector<cv::Mat> &_images );
  bool getGrayImages( const std::vector<cv::Mat> &_rgbImages,
		      std::vector<cv::Mat> &_grayImages );

  //-- Get functions
  inline cv::Mat getRgb( int _ind );
  inline cv::Mat getWarped( int _ind );
  inline int getNumFrames();

  //-- Public variables
  std::vector<cv::Mat> mLocalTransforms;
  std::vector<cv::Mat> mGlobalTransforms;

 private:
  // Image  data variables
  std::vector<cv::Mat> mRgbImages;
  std::vector<cv::Mat> mGrayImages;
  std::vector<cv::Mat> mWarpedImages;

  int mNumFrames;
  int mHeight;
  int mWidth;
  
  // SURF Keypoint detection and descriptor
  int mMinHessian;
  float mRadiusFactor;

  std::vector< std::vector<cv::KeyPoint> > mKeypoints; 
  std::vector< cv::Mat > mDescriptors;
  std::vector< std::vector< std::vector< cv::DMatch > > > mMatches;
  

  // RANSAC variables
  int mN; /**< Num Trials */
  int mS; /**< Model instantiation size */
  float mThresh;

  // SVD 
  cv::SVD svd;
  cv::Mat mA;
  std::vector<cv::Mat> mParams;
};
  
////////////////// INLINE FUNCTIONS /////////////////////////////////////// 

/**
 * @function getRgb
 */
inline cv::Mat Ransac::getRgb( int _ind ) {
  
  if( _ind < mNumFrames ) {
      return mRgbImages[_ind];
  }
}

/**
 * @function getWarped
 */
inline cv::Mat Ransac::getWarped( int _ind ) {
  
  if( _ind < mNumFrames ) {
      return mWarpedImages[_ind];
  }
}

/**
 * @function getNumFrames
 */
inline int Ransac::getNumFrames() {
  return mNumFrames;
}


#endif /** __RANSAC_H__ */
  
