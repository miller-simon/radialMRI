README.txt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ABOUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The Oakland Reconstruction Toolbox is designed to offer multiple methods of
reconstruction of images from radially sampled NMR data collected on Bruker AVII or AVIIIHD systems in 
MATLAB. It could easily be adapted for other scanner data. Reconstruction methods include Standard Gridding, Gridding with Oversampling,
Gridding with Oversampling and Deapodization, and Filtered Back Projection.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INSTRUCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Before using the toolbox, make sure to add the toolbox (the folder 
OaklandReconstructionToolbox) and its subfolders to the MATLAB path. The current folder 
should be 'OaklandReconstructionToolbox' when the code is run.

When inserting scan data into the 'Expts' file, be sure to insert it into the folder 
associated with the correct system. Also, the folder containing the scan data should
simply be labeled with the correct scan number. It should contain a method file that 
includes Kx and Ky coordinates.

The function oaklandReconstruct.m is intended to provide access to all reconstruction
techniques. When running this file, it will ask for which system the scan was generated
on and a scan number. 

If you would like to use your own trajectory coordinates, you can call oaklandReconstructForceK.m and
pass trajectory coordinates as arguments.

If the scan folder is organized correctly and the file 'method'
is undisturbed, the Kx and Ky coordinates will be read successfully. If an error is
thrown, the method file was likely disturbed by a user. Try to add a space to every line
where Kx and Ky are identified in the method document.

Next, it will ask which reconstruction algorithm is desired. It will then call the
algorithm and output the data into the folder 'Images.' 

The techniques are basically characterized as a gridding approach  or a filtered back 
projection approach. Gridding approaches include Standard Gridding, Gridding with 
Oversampling, and Gridding with Oversampling and Deapodization. Filtered back projection 
approaches include Standard Filtered Back Projection and Filtered Back Projection with Spatial 
Filtering. Through experimentation, we have determined Gridding with Oversampling and
Deapodization to be the most accurate method.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stdGridding,m, griddingOS.m, and griddingDeapodize.m important variables
    complexMatrix - matrix that contains 3 tuples of coordinates and values for every sample
    fineVals - values interpolated on to the ultrafine grid
    imageMat - the image resulting from the reconstruction thus far
    imaginaryData - the imaginary components of the FID
    Kx, newKx - K space x coordinates from a radial sample
    Ky, newKy - K space y coordinates from a radial sample
    neighbors - the neighborhood around a standrad Cartesian grid point used to asssign its values
    radVec - vector of approximate radii at which samples are taken
    realData - the real components of the FID 
    stdVals - values interpolated on to the standard grid
    thetaVec - angular values at which radial samples are taken
    transformedData - complex data resulting from 2DFFT
    x - meshgrid of approximate x locations at which samples were taken
    xScale - how much x coordinates are scaled by to achieve 256x256 K-space
    y - meshgrid of approximate x locations at which samples were taken
    yScale - how much y coordinates are scaled by to achieve 256x256 K-space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stdFBP.m and fbpSpatialFilt.m important variables
    data - nonuniform complex K space data
    fil - the spatial filter that corrects for the parabolic structure in reconstructed images
    imageMat - the image resulting from the reconstruction thus far
    imaginaryData - the imaginary components of the FID
    interpData - K space data interpolated to evenly spaced points along diameter
    p - the polynomial coeffecients fit for the single tube scan
    projections - matrix of projection data
    radVec - vector of approximate radii at which samples are taken
    realData - the real components of the FID
    thetaVec - angular values at which radial samples are taken
    transformedData - complex data resulting from FFT


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OTHER GENERAL REMARKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The functions testFBP.m and testOS.m have been included to allow for experimentation in reconstruction 
techniques. The variables are the same as those from stdFBP.m and stdGridding.m, respectively.

There was a discrepancy between the original code and the naming conventions used in the corresponding documents. 
Thus, variable names were changed. If, for some reason, a reconstruction fails, go to the 
ReconstructionMethods folder, then the OriginalMethods folder. You can replace the nonworking code 
with this, in hopes that the issue was related to some change in variables used.

This software is for research purposes only. It was developed in the Xia Lab  at Oakland University
with assistance from Farid Badar, Simon Miller, and Dr. Yang Xia.

For any questions about how to operate the toolbox, email simonmiller@oakland.edu.










