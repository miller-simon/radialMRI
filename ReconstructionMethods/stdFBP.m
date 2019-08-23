function imageMat = stdFBP(Kx,Ky,main,class,scan)
%This is the filtered back projection function for the Oakland Toolbox.

%This shifts the Kx and Ky vectors such that they begin at approximately 0.
Kx = 1*Kx(1) + Kx + 1E-10;
Ky = -1*Ky(1) + Ky + 1E-10;

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

    %Separate the FID into real and imaginary.
    realData = fidData(:,1:2:262-1);
    imaginaryData = fidData(:,2:2:262);

    %Then combine them to get the complex data for each sample.
    complexData = realData + 1i*imaginaryData;
    
    %Our data is initally organized as radial readouts. We would like each
    %column to span the diameter of K-space.
    A = complexData(1:402,:);
    B = complexData(403:804,:);
    B = flip(B,2);
    data = [B A]';
     
    %This creates vectors with evenly spaced angles from [0,360). Each one
    %corresponds to a readout.
    thetaVec = linspace(0,360-(180/402),804);

    %This creates a vector of radii of a given readout.
    radVec = sqrt(Kx.^2 + Ky.^2);
    
    %We then transform it to act as a number line centered at 0.
    x = [-1*(flip(radVec)) radVec];
    
    %This defines evenly spaced query points along a diameter of K-space.
    xQuery = linspace(-1*max(rad),max(rad),256);
    
    %We use the sample data to interpolate to the query points for each
    %K-space diameter.
    interpData = zeros(256,402);
    for i = 1:402
        y = data(:,i);
        interpData(:,i) = interp1(x,y,xQuery,'linear');
    end
    
    %We take the modulus of the 1DFFT to obtain projection information.
    transformedData = fftshift(fft(interpData));
    projections = sqrt(real(transformedData).^2 + imag(transformedData).^2);
    
    %We transform the sinogram such that it is continuous and smooth.
    projections = [projections(:,202:402) projections(:,1:201)];
    projections(:,1) = projections(:,2);
    
    %This creates 360 degrees of projection data.
    projections = [projections flipud((projections))];
    rightProjections = [projections(256,403:804); projections(1:255,403:804)];
    projections = [projections(:,1:402) rightProjections];
    
    %This performs the filtered back projection.
    imageMat = iradon(abs(projections),thetaVec,'linear','Ram-Lak',1,256);

    %This makes the unknown parts of the image black.
    black = makeBowl(256);
    imageMat = imageMat.*black; 
    
    %This sets values less than 0, produced by iradon, equal to 0.
    imageMat(imageMat < 0)=0;
    
    %This saves the image.
    binTitle = sprintf('oaklandFBP%sScan%s.bin',class,scan);
    binTitle = fullfile('/Users/yangxia/Desktop/SimonMiller/OaklandReconstructionToolbox/Images',binTitle);
    fid=fopen(binTitle,'w');
    fwrite(fid,imageMat,'double');
    fclose(fid);
    
else
%If the length of Kx is not 131, we assume that we are constructing a
%128x128 image.
    
    %This opens the data.
    fileID = fopen(data, 'r','l');
    fidData = fread(fileID,[256,402],'int32')';
    fclose(fileID);

    %Anything beyond the 134 column is 0. We do not need this information.
    fidData = fidData(:,1:134);

    %Separate the FID into real and imaginary.
    realData = fidData(:,1:2:134-1);
    imaginaryData = fidData(:,2:2:134);
    
    %Then combine them to get the complex data for each sample.
    complexData = realData + 1i*imaginaryData;   
   
    %Our data is initally organized as radial readouts. We would like each
    %column to span the diameter of K-space.
    A = complexData(1:201,:);
    B = complexData(202:402,:);
    B = flip(B,2);
    data = [B A]';

    %This creates vectors with evenly spaced angles from [0,360). Each one
    %corresponds to a readout.
    thetaVec = linspace(0,360-(180/201),402);
    
    %This creates a vector of radii of a given readout.
    radVec = sqrt(Kx.^2+Ky.^2);
    
    %We then transform it to act as a number line centered at 0.
    x = [-1*(flip(radVec)) radVec];
    
    %This defines evenly spaced query points along a diameter of K-space.
    xQuery = linspace(-1*max(rad),max(rad),128);
    
    %We use the sample data to interpolate to the query points for each
    %K-space diameter.
    interpData = zeros(128,201);
    for i = 1:201
        y = data(:,i);
        interpData(:,i) = interp1(x,y,xQuery,'linear');
    end
    
    %We take the modulus of the 1DFFT to obtain projection information.
    projections = fftshift(fft(interpData));
    projections = sqrt(real(projections).^2+imag(projections).^2);
    
    %We transform the sinogram such that it is continuous and smooth.
    projections = [projections(:,101:201) projections(:,1:100)];
    projections(:,1) = projections(:,2);
    
    %This creates 360 degrees of projection data.
    projections = [projections flipud((projections))];
    rightProjections = [projections(128,202:402); projections(1:127,202:402)];
    projections = [projections(:,1:201) rightProjections];

    %This performs the filtered back projection.
    imageMat = iradon(abs(projections),thetaVec,'linear','Ram-Lak',1,128);

    %This makes the unknown parts of the image black.
    black = makeBowl(128);
    imageMat = imageMat.*black;
    
    %This sets values less than 0, produced by iradon, equal to 0.
    imageMat(imageMat < 0)=0;
    
    %This saves the image.
    binTitle = sprintf('oaklandFBP%sScan%s.bin',class,scan);
    binTitle = fullfile('/Users/yangxia/Desktop/SimonMiller/OaklandReconstructionToolbox/Images',binTitle);
    fid = fopen(binTitle,'w');
    fwrite(fid,imageMat,'double');
    fclose(fid);
end

    

