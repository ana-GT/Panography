/**
 * @fie imagelist_creator
 * @brief Creates a .yaml or .xml list of files from the command line args
 */


#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <string>
#include <iostream>

/**
 * @function help
 */
static void help( char** av )
{
  std::cout << "\nThis creates a yaml or xml list of files from the command line args" << std::endl;
  std::cout<< "Usage:\n./" << av[0] << " imagelist.yaml *.png" << std::endl;
  std::cout<< "Try using different extensions.(e.g. yaml yml xml xml.gz etc...)"<< std::endl;
  std::cout<< "This will serialize this list of images or whatever with opencv's FileStorage framework" << std::endl;
}

/**
 * @function main
 */
int main(int ac, char** av)
{
  if (ac < 3)
    {
      help(av);
      return 1;
    }
  
  std::string outputname = av[1];
  
  // Check if the output is an image - prevent overwrites!
  cv::Mat m = cv::imread(outputname); 
  if( !m.empty() ){
    std::cerr << "[X] Fail! Please specify an output file, don't want to overwrite you images!" << std::endl;
    help(av);
    return 1;
  }

  // Do the storage
  cv::FileStorage fs( outputname, cv::FileStorage::WRITE );
  fs << "images" << "[";
  for(int i = 2; i < ac; i++){
    fs << std::string(av[i]);
  }
  fs << "]";
  return 0;
}
