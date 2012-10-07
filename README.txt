
Hi, here some brief instructions of how to run my code:

1. Run CMake and you will obtain three executables:

* showPCD : To show an existing PCD (./showPCD <textfile> <PCD-index>
* grabPCD: To grab PCDs from Kinect (press Space to capture each frame)
* Reconstruct PCD: This is the project itself (RANSAC + Matching  +...)

The first two files were used to grab data and check how it looked before
processing it. The third one is the main one.

2. To actually run the executable: 

./Reconstruct3D <textfile-with-PCD-files-location>

The argument must be a text file that has the location of the PCDs files
grabbed with the Kinect. For instance you might run

./Reconstruct3D pcdFilenames.txt

in pcdFilenames.txt there are a couple of dummy PCDfile names, just put your files in the folder and run this thing. The output will be a visualization
of the reconstructed 3D area + each PCD frame converted to the global frame (chosen to be frame 0 ) and a big PCD with all the individual global PCDs bundled together.

If you want to see the OpenCV matching stuff, check Reconstruct3D.cpp
and uncomment following the instructions (it is commented now so the PCD
viewer run its loop). Additionally the program outputs in the terminal
the local transformation obtained with RANSAC, if you want to verify :)

Thanks for reading!

Ana
