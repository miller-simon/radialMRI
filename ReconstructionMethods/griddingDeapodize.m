function griddingDeapodize(Kx,Ky,main,class,scan)
%This is the gridding with oversampling and deapodization function for the
%Oakland Toolbox.

%This sets the file path for a given scan's FID.
data = strcat(main,'/',class,'/',scan,'/fid');

%If the length of Kx is 131, we run the following to construct a 256x256
%image.
if length(Kx) == 131
    %This opens the data.
    fileID = fopen(data, 'r','l');
    fidData = fread(fileID,[512,804],'int32')';
    fclose(fileID);
   
    %Anything beyond the 262 column is 0. We do not need this information.
    fidData = fidData(:,1:262);
   
    %Separate the FID into real and imaginary, then combine into complex data.
    realData = fidData(:,1:2:262-1);
    imaginaryData = fidData(:,2:2:262);
    complexData = realData + 1i*imaginaryData;
    
    %This shifts the Kx and Ky vectors such that they begin at 0.
    newKx = (-1*Kx(1)) + Kx;
    newKy = (-1*Ky(1)) + Ky;
   
    %This creates a vector of radii of a given readout.
    radVec = sqrt((newKx.^2) + (newKy.^2));
   
    %This creates vectors with evenly spaced angles from [0,360). Each one
    %corresponds to a readout.
    thetaVec = linspace(0,360-(180/402),804);
    thetaVec = deg2rad(thetaVec);

    %We must convert polar coordinates to Cartesian coordinates.
    [r,t] = meshgrid(radVec,thetaVec);
    x = r.*cos(t);
    y = r.*sin(t);

    %We put all of the information for a sample in one place for ease of
    %use.
    complexMatrix = zeros(804,131,3);
    complexMatrix(:,:,1) = x;
    complexMatrix(:,:,2) = y;
    complexMatrix(:,:,3) = complexData;
    
    %Our data is initally organized as radial readouts. We would like each
    %column to span the diameter of K-space.
    completeMatrix = zeros(262,402,3);
    for l = 1:3
        A = complexMatrix(1:402,:,l);
        B = complexMatrix(403:804,:,l);
        B = flip(B,2);
        completeMatrix(:,:,l) = [B A]';
    end

    %We elect to scale the K-space to a 256x256 area with the bottom left
    %as (0,0) to easily correspond to matrix indices.
    xScale = 256/(2*max(max(completeMatrix(:,:,1))));
    yScale = 256/(2*max(max(completeMatrix(:,:,2))));
    completeMatrix(:,:,1) = (xScale*completeMatrix(:,:,1)) + 128;
    completeMatrix(:,:,2) = (yScale*completeMatrix(:,:,2)) + 128;
    
    %These will act as inputs for the interpolation below.
    sampX = completeMatrix(:,:,1);
    sampY = completeMatrix(:,:,2);
    sampV = completeMatrix(:,:,3);
    
    %Cubic Interpolation
    %The vector a will contain information about the grid fineness. We can
    %use multiple values in a to obtain different reconstructions.
    a = (16);
    for step = 1:length(a)
        %We define an utrafine grid to interpolate on to.
        [fineX, fineY] = meshgrid(0:1/a(step):256,0:1/a(step):256);

        %We interpolate on to the grid.
        fineVals = griddata(sampX,sampY,sampV,fineX,fineY,'cubic');

        %If any values are NaN, we set them to 0.
        fineVals(isnan(fineVals)) = 0;
        
        %Next, we interpolate on to the twice fine Cartesian grid with .5
        %spacing. This is called oversampling.
        stdVals = zeros(512,512);
        for i = 1:512
            for j = 1:512
                %We use neighboring points to assign values to the
                %standard grid. We currently use a simple average.
                neighbors = fineVals((a(step)*i*.5)-(a(step)/4):(a(step)*i*.5),(a(step)*j*.5)-(a(step)/4):(a(step)*j*.5));  
                stdVals(i,j)=sum(sum(neighbors))/(((a(step)/4)+1)^2);
            end
        end
        
        %We take the modulus of the 2DFFT to obtain image data
        transformedData = fftshift(fft2(stdVals));
        imageMat = sqrt((real(transformedData).^2) + (imag(transformedData).^2))';
        imageMat = rot90(imageMat, 2);


        %This is the deapodization process.
        %We generate the 2D FFT of the Gaussian kernel to divide against
        %the image.
        X = linspace(max(radVec)*-1,max(radVec),512);
        sig = .25;
        y = FFTgaussian(X,0,sig);
        ker = y'*y;
        imageMat = imageMat./ker;
        
        %We select only the center portion of the image because we
        %oversampled in the frequency domain.
        imageMat = imageMat(129:384,129:384);
        
        %This makes the unknown parts of the image black.
        black = makeBowl(256);
        imageMat = imageMat.*black;

        %This saves the image.
        binTitle = sprintf('oaklandDeapodize%sScan%s.bin',class,scan);
        binTitle = fullfile('/Users/yangxia/Desktop/SimonMiller/OaklandReconstructionToolbox/Images',binTitle);
        fid=fopen(binTitle,'w');
        fwrite(fid,imageMat,'double');
        fclose(fid);
    end
else %-------------------------------------------------------------------------------------------------
%If the length of Kx is not 131, we assume that we are constructing a
%128x128 image.
%This opens the data.
    fileID = fopen(data, 'r','l');
    fidData = fread(fileID,[256,402],'int32')';
    fclose(fileID);
   
    %Anything beyond the 134 column is 0. We do not need this information.
    fidData = fidData(:,1:134);
    
    %Separate FID into real and imaginary, then combine into complex data.
    realData = fidData(:,1:2:134-1);
    imaginaryData = fidData(:,2:2:134);
    complexData = realData + 1i*imaginaryData;
    
    %This shifts the Kx and Ky vectors such that they begin at 0.
    newKx = (-1*Kx(1)) + Kx;
    newKy = (-1*Ky(1)) + Ky;
    
    %This creates a vector of radii of a given readout.
    radVec = sqrt((newKx.^2)+(newKy.^2));

    %This creates vectors with evenly spaced angles from [0,360). Each one
    %corresponds to a readout.
    thetaVec = linspace(0,360-(180/201),402);
    thetaVec = deg2rad(thetaVec);

    %We must convert polar coordinates to Cartesian coordinates.
    [r,t] = meshgrid(radVec,thetaVec);
    x = r.*cos(t);
    y = r.*sin(t);
    
    %We put all of the information for a sample in one place for ease of
    %use.
    complexMatrix = zeros(402,67,3);
    complexMatrix(:,:,1) = x;
    complexMatrix(:,:,2) = y;
    complexMatrix(:,:,3) = complexData;
    
    %Our data is initally organized as radial readouts. We would like each
    %column to span the diameter of K-space.
    completeMatrix =zeros(134,201,3);
    for l =1:3
        A = complexMatrix(1:201,:,l);
        B = complexMatrix(202:402,:,l);
        B = flip(B,2);
        completeMatrix(:,:,l) = [B A]';
    end


    %We elect to scale the K-space to a 128x128 area with the bottom left
    %as (0,0) to easily correspond to matrix indices.
    xScale = 128/(2*max(max(completeMatrix(:,:,1))));
    yScale = 128/(2*max(max(completeMatrix(:,:,2))));
    completeMatrix(:,:,1) = (xScale*completeMatrix(:,:,1)) + 64;
    completeMatrix(:,:,2) = (yScale*completeMatrix(:,:,2)) + 64;

    %These will act as inputs for the interpolation below.
    sampX = completeMatrix(:,:,1);
    sampY = completeMatrix(:,:,2);
    sampV = completeMatrix(:,:,3);
    
    %128 Cubic Interpolation
    %The vector a will contain information about the grid fineness. We can
    %use multiple values in a to obtain different reconstructions.
    a = (16);
    for step = 1:length(a)
        %We define an utrafine grid to interpolate on to.
        [fineX, fineY] = meshgrid(0:1/a(step):128,0:1/a(step):128);

        %We interpolate on to the grid.
        fineVals = griddata(sampX,sampY,sampV,fineX,fineY,'cubic');
        
        %If any values are NaN, we set them to 0.
        fineVals(isnan(fineVals)) = 0;
        
        %Next, we interpolate on to the twice fine Cartesian grid with .5
        %spacing. This is called oversampling.
        stdVals = zeros(256,256);
        for i = 1:256
            for j = 1:256
                %We use neighboring points to assign values to the
                %standard grid. We currently use a simple average.
                neighbors = fineVals((a(step)*i*.5)-(a(step)/4):(a(step)*i*.5),(a(step)*j*.5)-(a(step)/4):(a(step)*j*.5));  
                stdVals(i,j) = sum(sum(neighbors))/(((a(step)/4)+1)^2);
            end
        end
        
        %We take the modulus of the 2DFFT to obtain image data
        transformedData = fftshift(fft2(stdV));
        imageMat = sqrt((real(transformedData).^2) + (imag(transformedData).^2))';
        imageMat = rot90(imageMat, 2);
        
        %This is the deapodization process.
        %We generate the 2D FFT of the Gaussian kernel to divide against
        %the image.
        X = linspace(max(radVec)*-1,max(radVec),256);
        sig = .25;
        y = FFTgaussian(X,0,sig);
        ker = y'*y;
        imageMat = imageMat./ker;
        
        %We select only the center portion of the image because we
        %oversampled in the frequency domain.
        imageMat = imageMat(65:192,65:192);

        %This makes the unknown parts of the image black.
        black = makeBowl(128);
        imageMat = imageMat.*black;

        %This saves the image.
        binTitle = sprintf('oaklandDeapodize%sScan%s.bin',class,scan);
        binTitle = fullfile('/Users/yangxia/Desktop/SimonMiller/OaklandReconstructionToolbox/Images',binTitle);
        fid=fopen(binTitle,'w');
        fwrite(fid,imageMat,'double');
        fclose(fid);
    end
end