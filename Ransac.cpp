/**
 * @file Ransac.cpp
 * @author A. Huaman
 */

#include "Ransac.h"
#include <time.h>
#include <Eigen/Core>
#include <stdio.h>
#include <iostream>
#include <fstream>

/**
 * @function Ransac
 * @brief Constructor
 */
Ransac::Ransac( char *_filenames ) {
  
  // Initialize a few things
  mMinHessian = 400;
  mRadiusFactor = 3.0; // 3.0
  mThresh = 5; // pixels.
  mFactor = 2;

  // For random generation
  srand( time(NULL) );

  // Read images
  if( !readImagesData( _filenames, mRgbImages ) ) {
    printf( "[X] No Rgb Data read! \n" );
  }

  mNumFrames = mRgbImages.size();
  printf( "Loaded %d images \n", mNumFrames );

  // Get gray
  getGrayImages( mRgbImages, mGrayImages );
  printf( "Loaded %d Gray images \n", mGrayImages.size() );  

  mHeight = mRgbImages[0].rows;
  mWidth = mRgbImages[0].cols;

}

/**
 * @function ~Ransac
 * @brief Destructor
 */
Ransac::~Ransac() {

}

////////////////////////////// RANSAC ////////////////////////////////////////////

/**
 * @function Ransac_Homography2D2
 */
cv::Mat Ransac::Ransac_Homography2D2( int _ind1, int _ind2 ) {
  
  int totalSamples = mMatches[_ind1][_ind2].size();
  std::vector<cv::Point2f> points1, points2;
  cv::Point2f p1, p2;

  // Get corresponding points
  for( int i = 0; i < totalSamples; ++i ) {
    getPointsFromMatch( _ind1, _ind2, 
			mMatches[_ind1][_ind2][i],
			p1, p2 );
    points1.push_back( p1 );
    points2.push_back( p2 );
  }

  //-- Find homography
  cv::Mat H = cv::findHomography( cv::Mat(points1), cv::Mat( points2), 
				  CV_RANSAC );
  return H;
}

/**
 * @function Ransac_Homography2D
 */
cv::Mat Ransac::Ransac_Homography2D( int _ind1, int _ind2 ) {

  mN = 500;
  mS = 4;

  int totalSamples;
  std::vector<int> randomSampleIndices;
  cv::Mat params; cv::Mat bestParams;
  std::vector<int> bestInitSet;
  int count; int bestCount; float error;

  totalSamples = mMatches[_ind1][_ind2].size();
  bestCount = 0;
  printf("[Ransac] Model: %d total samples: %d \n", mS, totalSamples );

  // Loop
  for( int i = 0; i < mN; ++i ) {

    randomSampleIndices = getRandomIndices( mS, totalSamples );
    params = getModel_Homography2D( _ind1, _ind2, randomSampleIndices  );

    // Calculate inliers / outliers
    count = 0;
    for( int j = 0; j < totalSamples; ++j ) {
      error = getErrorModel( _ind1, _ind2, j, params );
      if( error < mThresh ) {
	count++;
      }
    }
    // ..compare with best
    if( count > bestCount ) {
      bestCount = count;
      bestParams = params;
      bestInitSet = randomSampleIndices;
    }
  }
  printf( "** Best params  with count %d / %d (%d - %d )  \n", bestCount, totalSamples, _ind1, _ind2 ); 
  std::cout<< bestParams << std::endl;

  //return getSomeMatchesDraw( _ind1, _ind2, bestInitSet );
  return bestParams;
}


/**
 * @function getModel_Homography2D 
 */
cv::Mat Ransac::getModel_Homography2D( int _ind1,
				       int _ind2, 
				       std::vector<int> _indices ) {
 
  cv::Mat param(9, 1, CV_32FC1 );
  cv::Mat A;
  float x1, y1; 
  float x2, y2;
  
  // Check that there are four params
  if( _indices.size() != 4 ) {
    printf("[X] [geModel_Homography2D] What the heck? Give me 4 freaking points! Exiting \n");
    return param;
  }

  int ind;
  A = cv::Mat::zeros( 8, 9, CV_32FC1 );
  cv::Mat D, U, Vt, V;
  
  for( int i = 0; i < _indices.size(); ++i ) {
    getPointsFromMatch( _ind1, _ind2, mMatches[_ind1][_ind2][ _indices[i] ],
			  x1, y1, x2, y2 );
    // Fill A and B
    ind = i*2;
    A.at<float>(ind,3) = x1; 
    A.at<float>(ind,4) = y1; 
    A.at<float>(ind,5) = 1; 
    A.at<float>(ind,6) = -y2*x1; 
    A.at<float>(ind,7) = -y2*y1; 
    A.at<float>(ind,8) = -y2*1;

    ind = i*2 + 1;
    A.at<float>(ind,0) = x1; 
    A.at<float>(ind,1) = y1; 
    A.at<float>(ind,2) = 1; 
    A.at<float>(ind,6) = -x2*x1; 
    A.at<float>(ind,7) = -x2*y1; 
    A.at<float>(ind,8) = -x2*1;
  }

  // Compute the wonderful SVD :)
  cv::SVD::compute( A, D, U, Vt );
  V = Vt.t();
  printf("V cols: %d rows: %d \n", V.cols, V.rows);
  // H is the last column of V (nullspace of A)
  cv::Mat transf = cv::Mat::zeros( 3, 3, CV_32FC1 );
  
  float h33 = V.at<float>(8,7);

  transf.at<float>(0,0) = V.at<float>(0,7) / h33;
  transf.at<float>(0,1) = V.at<float>(1,7) / h33;
  transf.at<float>(0,2) = V.at<float>(2,7) / h33;

  transf.at<float>(1,0) = V.at<float>(3,7) / h33;
  transf.at<float>(1,1) = V.at<float>(4,7) / h33;
  transf.at<float>(1,2) = V.at<float>(5,7) / h33;

  transf.at<float>(2,0) = V.at<float>(6,7) / h33;
  transf.at<float>(2,1) = V.at<float>(7,7) / h33;
  transf.at<float>(2,2) = V.at<float>(8,7) / h33;

  return transf;
}

/**
 * @function getErrorModel
 * @brief Outputs error of Homography2D model 
 */
float Ransac::getErrorModel( int _ind1, int _ind2, 
			     int _ind,
			     const cv::Mat &_params ) {

  float x1, y1;
  float x2, y2;
  float px2, py2, pw2;

  getPointsFromMatch( _ind1, _ind2, 
		      mMatches[_ind1][_ind2][_ind],
		      x1, y1, 
		      x2, y2 );

  // Apply model
  px2 = _params.at<float>(0,0)*x1 +  _params.at<float>(0,1)*y1 +
    _params.at<float>(0,2)*1.0;
 
  py2 = _params.at<float>(1,0)*x1 +  _params.at<float>(1,1)*y1 +
    _params.at<float>(1,2)*1.0;

  pw2 = _params.at<float>(2,0)*x1 +  _params.at<float>(2,1)*y1 +
    _params.at<float>(2,2)*1.0;

  px2 = px2 / pw2;
  py2 = py2 / pw2;

  return( sqrt( (x2 - px2)*(x2-px2) +  (y2 - py2)*(y2-py2)  ) );
}


/**
 * @function getPointsFromMatch
 */
void Ransac::getPointsFromMatch( int _ind1, int _ind2, 
				const cv::DMatch &_match,
				float &_x1, float &_y1, 
				float &_x2, float &_y2 ) {

  int idx1; int idx2;
  cv::Point2f p1; cv::Point2f p2;
  int h1, w1; int h2, w2;
  
  // Get both 2D Points from keypoint info
  idx1 = _match.queryIdx;
  idx2 = _match.trainIdx;
  
  p1 = mKeypoints[_ind1][idx1].pt; 
  _y1 = p1.y; _x1 = p1.x;
  
  p2 = mKeypoints[_ind2][idx2].pt; 
  _y2 = p2.y; _x2 = p2.x;
  
}

/**
 * @function getPointsFromMatch
 */
void Ransac::getPointsFromMatch( int _ind1, int _ind2, 
				 const cv::DMatch &_match,
				 cv::Point2f &_p1, 
				 cv::Point2f &_p2 ) {

  int idx1; int idx2;
  
  // Get both 2D Points from keypoint info
  idx1 = _match.queryIdx;
  idx2 = _match.trainIdx;
  
  _p1 = mKeypoints[_ind1][idx1].pt; 
  _p2 = mKeypoints[_ind2][idx2].pt; 
  
}


/**
 * @function getRandomIndices
 */
std::vector<int> Ransac::getRandomIndices( int _numSamples, 
					   int _totalSamples ) {

  std::vector<int> randomIndices;
  int currentIndice;

  if( _numSamples > _totalSamples ) {
    printf( "[X] [getRandomIndices] Model needs more than available keypoints. Exiting! \n");
  }

  else {
    for( int i = 0; i < _numSamples; ++i ) {  
      do {
	currentIndice = getRandom( _totalSamples );
      } while( isInSet( currentIndice, randomIndices ) == true );
      randomIndices.push_back( currentIndice );
    }
  }

  return randomIndices;
}

/**
 * @function isInSet
 */
bool Ransac::isInSet( int _index, 
		      std::vector<int> _currentIndices ) {

  for( int i = 0; i < _currentIndices.size(); ++i ) {
    if( _index == _currentIndices[i] ) {
      return true;
    }
  }
  return false;
}

/**
 * @function getRandom
 */
int Ransac::getRandom( int _max ) {
  return( rand() % _max );
}


/**
 * @function applyHomography2DToImage
 */ 
cv::Mat Ransac::applyHomography2DToImage( const cv::Mat &_homography2D, 
					  const cv::Mat &_image ) {

  cv::Mat transformed;
  return transformed;
}

/**
 * @function getAllTransforms
 */
void Ransac::getAllTransforms() {

  printf( "... Getting all homographies \n" );
  
  mGlobalTransforms.resize( mNumFrames );
  
  cv::Mat tf_global( 3, 3, CV_32FC1 );
  
  // Take global frame as 0
  tf_global = cv::Mat::eye( 3, 3, CV_32FC1 );
  mGlobalTransforms[0] = tf_global.clone();

  // Get the transforms to 0 frame only
  for( int i = 1; i < mNumFrames; ++i ) {
    tf_global = Ransac_Homography2D2( i, 0 );
    mGlobalTransforms[i] = tf_global.clone();
  }

  // Print these bad guys 
  for( int i = 0; i < mNumFrames; ++i ) {
    std::cout << "Transform "<<i<<" to 0: \n"<<mGlobalTransforms[i]<<std::endl;
  }
}

/**
 * @function getTransformedImages
 */
void Ransac::getTransformedImages() {

  // Apply the transformations to the images
  mWarpedImages.resize( mNumFrames );
  int warp_h = mHeight*2;
  int warp_w = mWidth*2;

  //  mWarpedImages[0] = mRgbImages[0].clone();

  cv::Mat temp_warped = cv::Mat( warp_h, warp_w, CV_8UC3 );
  for( int i = 0; i < mNumFrames; ++i ) {
    
    cv::warpPerspective( mRgbImages[i], 
			 temp_warped, 
			 mGlobalTransforms[i], 
			 cv::Size( warp_w, warp_h ),
			 cv::INTER_LINEAR,
			 cv::BORDER_CONSTANT,
			 cv::Scalar() );
    mWarpedImages[i] = temp_warped.clone();
  }
  
}

/**
 * @function saveTransformedPointClouds
 */
void Ransac::saveTransformedImages() {

  // Reset some values
  mOrigX.resize(mNumFrames);
  mOrigY.resize(mNumFrames);

  // Create ROIs to use
  std::vector<cv::Mat> ROIs;
  ROIs.resize( mNumFrames );

  // Create a big image (Panography)
  cv::Mat pan = cv::Mat( 4*mHeight, 8*mWidth, CV_8UC3 );

  // Locate center image at the middle
  int pan_w = pan.cols;
  int pan_h = pan.rows;
  int pan_w2 = pan_w / 2;
  int pan_h2 = pan_h / 2;
  int warp_w; int warp_h;
  int warp_w2; int warp_h2;

  //-- [0] Create ROI in the middle
  warp_w = mWarpedImages[0].cols; 
  warp_h = mWarpedImages[0].rows;

  mOrigX[0] = pan_w2 - mWidth/2;
  mOrigY[0] = pan_h2 - mHeight/2;
 
  ROIs[0] = cv::Mat( pan, cv::Rect( mOrigX[0], mOrigY[0], warp_w, warp_h) );
	
  //-- [1-4] Create ROI for image i
  /*
  for( int i = 1; i < mNumFrames; ++i ) {

    warp_w = mWarpedImages[i].cols;
    warp_h = mWarpedImages[i].rows;

    mOrigX[i] = mOrigX[0] + 0*(int) ( mGlobalTransforms[i].at<double>(0, 2) );
    mOrigY[i] = mOrigY[0] + 0*(int) ( mGlobalTransforms[i].at<double>(1, 2) );
    printf("Warp[%d]: Ox: %d Oy: %d w: %d h: %d \n", i, mOrigX[i], mOrigY[i], warp_w, warp_h);
    
    ROIs[i] = cv::Mat( pan, cv::Rect( mOrigX[i], mOrigY[i], warp_w, warp_h) );
  }
  */
  //-- Put them all together
  /*
  for( int i = mNumFrames - 1; i >= 0; --i ) {
    mWarpedImages[i].copyTo( ROIs[i] );
  }
  */

  ROIs[4] = cv::Mat( pan, cv::Rect( mOrigX[0], mOrigY[0], 
				    mWarpedImages[4].cols,
				    mWarpedImages[4].rows ));
  mWarpedImages[4].copyTo( ROIs[4] );
  cv::Mat( pan, cv::Rect(0,0,mWidth, mHeight ));

  ROIs[3] = cv::Mat( pan, cv::Rect( mOrigX[0], mOrigY[0], 
				    mWarpedImages[3].cols,
				    mWarpedImages[3].rows ));
  mWarpedImages[3].copyTo( ROIs[3] );
  cv::Mat( pan, cv::Rect(0,0,mWidth, mHeight ));


  mWarpedImages[0].copyTo(ROIs[0]);
		  
  cv::imwrite("PanographyYeah.png", pan );


  //-- Write them individually
  for( int i = 0; i < mNumFrames; ++i ) {
    char imgName[50];
    sprintf( imgName, "warped%03d.png", i );
    cv::imwrite( imgName, mWarpedImages[i] );
  }
  
}

///////////////////////////// MATCHING //////////////////////////////////////////////////

/**
 * @function matchAllFrames
 */
void Ransac::matchAllFrames() {
 
  printf("Match all frames \n");
  cv::SURF surf( mMinHessian );
  
  // Detect keypoints and generate descriptors
  for( int i = 0; i < mNumFrames; ++i ) {

    std::vector<cv::KeyPoint> keypoint;
    cv::Mat descriptor;
    surf( mGrayImages[i], cv::Mat(), keypoint, descriptor, false );
    mKeypoints.push_back( keypoint );
    mDescriptors.push_back( descriptor );
    printf("[%d] Num keypoints: %d \n", i, keypoint.size() );
  }

  // Create a FlannBasedMatcher
  cv::FlannBasedMatcher matcher;
  std::vector< std::vector< std::vector< cv::DMatch > > > tempMatches;

  // Match each descriptor successively
  printf( " Start matching process \n" );
  for( int i = 0; i < mNumFrames; ++i ) {

    std::vector< std::vector< cv::DMatch > > matches_iToAll;
    for( int j = 0; j < mNumFrames; ++j ) {
      std::vector< cv::DMatch> matches_iToj;
      matcher.match( mDescriptors[i], mDescriptors[j], matches_iToj );
      matches_iToAll.push_back( matches_iToj );
    }
    tempMatches.push_back( matches_iToAll );
  }

  // Resize
  mMatches.resize( mNumFrames );
  for( int i = 0; i < mNumFrames; ++i ) {
    mMatches[i].resize( mNumFrames );
  }

  // Quick calculation of max and min distances
  double max_dist; 
  double min_dist;
  double dist;

  for( int i = 0; i < mNumFrames; ++i ) {
    for( int j = 0; j < mNumFrames; ++j ) {

      max_dist = 0; min_dist = 100;
      for( int k = 0; k < tempMatches[i][j].size(); ++k ) {
	dist = tempMatches[i][j][k].distance;
	if( dist < min_dist ) { min_dist = dist; }
	if( dist > max_dist ) { max_dist = dist; }
      }

      // Get only good matches
      for( int k = 0; k < tempMatches[i][j].size(); ++k ) {
	if( tempMatches[i][j][k].distance <= mRadiusFactor*min_dist ) {
	  mMatches[i][j].push_back( tempMatches[i][j][k] );
	}
      }

    }
  }

  printf("End matching process with these results: \n");
  for( int i = 0; i < mNumFrames; ++ i ) {
    printf("<%d> ", i );
    for( int j = 0; j < mNumFrames; ++j ) {
      printf(":[%d] %d  ", j, mMatches[i][j].size() );
    }
    printf("\n");
  }

}

/**
 * @function getMatchesDraw
 */
cv::Mat Ransac::getMatchesDraw( int _ind1, int _ind2  ) {

  //-- Draw matches
  cv::Mat matchesImage;

  cv::drawMatches( mRgbImages[_ind1], 
		   mKeypoints[_ind1], 
		   mRgbImages[_ind2], 
		   mKeypoints[_ind2], 
		   mMatches[_ind1][_ind2], 
		   matchesImage,
		   cv::Scalar::all(-1),
		   cv::Scalar::all(-1),
		   std::vector<char>(),
		   2 ); // 2: NOT_DRAW_SINGLE_POINTS  
  
  return matchesImage;
}

/**
 * @function getSomeMatchesDraw
 */
cv::Mat Ransac::getSomeMatchesDraw( int _ind1, int _ind2,
				    const std::vector<int> &matchesIndices ) {

  cv::Mat matchesImage;
  std::vector<cv::DMatch> matches;
  float x1, y1, x2, y2;
  for( int i = 0; i < matchesIndices.size(); ++i ) {
    matches.push_back( mMatches[_ind1][_ind2][ matchesIndices[i] ] );

    getPointsFromMatch( _ind1, _ind2, matches[i],
			x1, y1, x2, y2 ); 
    printf("Drawing match between: (%.3f, %.3f) and  (%.3f, %.3f) \n", x1, y1, x2, y2 );
  }

  cv::drawMatches( mRgbImages[_ind1],
		   mKeypoints[_ind1],
		   mRgbImages[_ind2],
		   mKeypoints[_ind2],
		   matches,
		   matchesImage,
		   cv::Scalar::all(-1),
		   cv::Scalar::all(-1),
		   std::vector<char>(),
		   2 ); // 2: NOT_DRAW_SINGLE_POINTS 
  return matchesImage;
}

/**
 * @function getKeypointsDraw
 */
cv::Mat Ransac::getKeypointsDraw( int _ind ) {

  cv::Mat keypointsDraw;
  cv::drawKeypoints( mRgbImages[_ind],  mKeypoints[_ind], 
		     keypointsDraw, cv::Scalar::all(-1), 
		     cv::DrawMatchesFlags::DEFAULT );
  return keypointsDraw;
}
    

//////////////////////////// DATA ACQUISITION ///////////////////////////////////////////

/**
 * @function readImagesData
 */
bool Ransac::readImagesData( char* _filenames, 
			     std::vector<cv::Mat> &_images ) {

  std::vector<std::string> imageFiles;
  
  // Open Stream
  _images.resize(0);
  std::ifstream ifs( _filenames );
  
  // Read
  std::string temp;
  while( getline( ifs, temp ) ) {
    imageFiles.push_back( temp );
  }
  
  // Load them 
  for( int i = 0; i < imageFiles.size(); ++i ) {

    cv::Mat image;

    // Read 3-channel image
    image = cv::imread( imageFiles[i], 1 );
    if( !image.data ) {
      std::cout<<"[X] Could not read file "<< imageFiles[i]<< std::endl;
     _images.resize(0);
      return false;
    }    

    // Save a bigger image
    cv::Mat bigImage = cv::Mat( image.rows*mFactor, image.cols*mFactor, CV_8UC3 );
    cv::Mat roi = cv::Mat( bigImage, cv::Rect(image.cols*mFactor/2 - image.cols/2,
					      image.rows*mFactor/2 - image.rows/2,
					      image.cols, image.rows) );
    image.copyTo( roi );
    
    _images.push_back( bigImage );
  }

  return true;
 
}

/**
 * @function getGrayImages
 * @brief 
 */
bool Ransac::getGrayImages( const std::vector<cv::Mat> &_rgbImages,
			    std::vector<cv::Mat> &_grayImages ) {

  _grayImages.resize(0);

  for( int i = 0; i < _rgbImages.size(); ++i ) {

    cv::Mat grayImage;
    _rgbImages[i].convertTo( grayImage, CV_8UC1 );
    _grayImages.push_back( grayImage );
  }

  return true;
}

