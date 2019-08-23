function griddingOSOld(Kx,Ky,main,class,scan)
%This is the gridding with oversampling function for the Oakland Toolbox.

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
    Real = fidData(:,1:2:262-1);
    Imaginary = fidData(:,2:2:262);
   
    rrr = Real;
    iii = Imaginary;
    
    %This shifts the Kx and Ky vectors such that they begin at 0.
    newKx=(-1*Kx(1))+Kx;
    newKy=(-1*Ky(1))+Ky;
   
    %This creates a vector of radii of a given readout.
    radVec=sqrt((newKx.^2)+(newKy.^2));
   
    %This creates vectors with evenly spaced angles from [0,360). Each one
    %corresponds to a readout.
    thetaVec = linspace(0,360-(180/402),804);
    thetaVec=deg2rad(thetaVec);

    %We must convert polar coordinates to Cartesian coordinates.
    [r,t]=meshgrid(radVec,thetaVec);
    x = r.*cos(t);
    y = r.*sin(t);

    %We put all of the information for a sample in one place for ease of
    %use.
    betterMatrix = zeros(804,131,4);
    betterMatrix(:,:,1)=x;
    betterMatrix(:,:,2)=y;
    betterMatrix(:,:,3)=rrr;
    betterMatrix(:,:,4)=iii;

    %Our data is initally organized as radial readouts. We would like each
    %column to span the diameter of K-space.
    betterMatrix2=zeros(262,402,4);
    for l =1:4
        bMA = betterMatrix(1:402,:,l);
        bMB = betterMatrix(403:804,:,l);
        bMB = flip(bMB,2);
        betterMatrix2(:,:,l) = [bMB bMA]';
    end

    betterMatComplex(:,:,1) = betterMatrix2(:,:,1);
    betterMatComplex(:,:,2) = betterMatrix2(:,:,2);
    
    %This adds the real and imaginary data to give is complex data for each
    %sample.
    betterMatComplex(:,:,3) = betterMatrix2(:,:,3) + 1i*betterMatrix2(:,:,4);


    %We elect to scale the K-space to a 256x256 area with the bottom left
    %as (0,0) to easily correspond to matrix indices.
    xScale = 256/(2*max(max(betterMatComplex(:,:,1))));
    yScale = 256/(2*max(max(betterMatrix2(:,:,2))));
    scaledMat(:,:,1)=xScale*betterMatComplex(:,:,1);
    scaledMat(:,:,2)=yScale*betterMatComplex(:,:,2);

    scaledMat(:,:,1)=scaledMat(:,:,1)+(256/2);
    scaledMat(:,:,2)=scaledMat(:,:,2)+(256/2);

    scaledMat(:,:,3)=betterMatComplex(:,:,3);

    %These will act as inputs for the interpolation below.
    sampX = scaledMat(:,:,1);
    sampY = scaledMat(:,:,2);
    sampV = scaledMat(:,:,3);
    
    
    assignin('base','sampX',sampX);
    assignin('base','sampV',sampV);
    %Cubic Interpolation
    %The vector a will contain information about the grid fineness. We can
    %use multiple values in a to obtain different reconstructions.
    a=(16);
    for step=1:length(a)
        %We define an utrafine grid to interpolate on to.
        [fineX, fineY]=meshgrid(0:1/a(step):256,0:1/a(step):256);

        %We interpolate on to the grid.
        fineV=griddata(sampX,sampY,sampV,fineX,fineY,'cubic');

        %If any values are NaN, we set them to 0.
        fineV(isnan(fineV))=0;
        
        %Next, we interpolate on to the twice fine Cartesian grid with .5
        %spacing. This is called oversampling.
        stdV=zeros(512,512);
        for i=1:512
            for j=1:512
                %We use neighboring points to assign values to the
                %standard grid. We currently use a simple average.
                z=fineV((a(step)*i*.5)-(a(step)/4):(a(step)*i*.5),(a(step)*j*.5)-(a(step)/4):(a(step)*j*.5));  
                stdV(i,j)=sum(sum(z))/(((a(step)/4)+1)^2);
            end
        end
        
        
        assignin('base','stdV',stdV);
        
        %We take the modulus of the 2DFFT to obtain image data
        Co = fftshift(fft2(stdV));
        R = abs(real(Co));
        I = abs(imag(Co));
        ImageMat = sqrt((R).^2 + (I).^2)';
        ImageMat = rot90(ImageMat, 2);
        
        %We select only the center portion of the image because we
        %oversampled in the frequency domain.
        ImageMat=ImageMat(129:384,129:384);
        
        %This makes the unknown parts of the image black.
        black = makeBowl(256);
        ImageMat=ImageMat.*black;
        assignin('base','ImageMat',ImageMat);
        %This saves the image.
        binTitle = sprintf('oaklandOversampling%sScan%s.bin',class,scan);
        binTitle = fullfile('/Users/yangxia/Desktop/SimonMiller/OaklandReconstructionToolbox/Images',binTitle);
        fid=fopen(binTitle,'w');
        fwrite(fid,ImageMat,'double');
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
    
    %Separate fid into real and imaginary
    Real = fidData(:,1:2:134-1);
    Imaginary = fidData(:,2:2:134);
   
    rrr = Real;
    iii = Imaginary;
    
    %This shifts the Kx and Ky vectors such that they begin at 0.
    newKx=(-1*Kx(1))+Kx;
    newKy=(-1*Ky(1))+Ky;
    
    %This creates a vector of radii of a given readout.
    radVec=sqrt((newKx.^2)+(newKy.^2));

    %This creates vectors with evenly spaced angles from [0,360). Each one
    %corresponds to a readout.
    thetaVec = linspace(0,360-(180/201),402);
    thetaVec=deg2rad(thetaVec);

    %We must convert polar coordinates to Cartesian coordinates.
    [r,t]=meshgrid(radVec,thetaVec);
    x = r.*cos(t);
    y = r.*sin(t);
    
    %We put all of the information for a sample in one place for ease of
    %use.
    betterMatrix = zeros(402,67,4);
    betterMatrix(:,:,1)=x;
    betterMatrix(:,:,2)=y;
    betterMatrix(:,:,3)=rrr;
    betterMatrix(:,:,4)=iii;

    %Our data is initally organized as radial readouts. We would like each
    %column to span the diameter of K-space.
    betterMatrix2=zeros(134,201,4);
    for l =1:4
        bMA = betterMatrix(1:201,:,l);
        bMB = betterMatrix(202:402,:,l);
        bMB = flip(bMB,2);
        betterMatrix2(:,:,l) = [bMB bMA]';
    end

    betterMatComplex(:,:,1) = betterMatrix2(:,:,1);
    betterMatComplex(:,:,2) = betterMatrix2(:,:,2);
    %This adds the real and imaginary data to give is complex data for each
    %sample.
    betterMatComplex(:,:,3) = betterMatrix2(:,:,3) + 1i*betterMatrix2(:,:,4);

    %We elect to scale the K-space to a 128x128 area with the bottom left
    %as (0,0) to easily correspond to matrix indices.
    xScale = 128/(2*max(max(betterMatComplex(:,:,1))));
    yScale = 128/(2*max(max(betterMatrix2(:,:,2))));
    scaledMat(:,:,1)=xScale*betterMatComplex(:,:,1);
    scaledMat(:,:,2)=yScale*betterMatComplex(:,:,2);

    scaledMat(:,:,1)=scaledMat(:,:,1)+(128/2);
    scaledMat(:,:,2)=scaledMat(:,:,2)+(128/2);
    scaledMat(:,:,3)=betterMatComplex(:,:,3);

    %These will act as inputs for the interpolation below.
    sampX = scaledMat(:,:,1);
    sampY = scaledMat(:,:,2);
    sampV = scaledMat(:,:,3);

    
    
    
    sampV(67,:)=sampV(67,:)+sampV(66,:);
    sampV(68,:)=sampV(68,:)+sampV(69,:);
    
    
    
    
    assignin('base','sampX',sampX);
    assignin('base','sampV',sampV);
    
    %128 Cubic Interpolation
    %The vector a will contain information about the grid fineness. We can
    %use multiple values in a to obtain different reconstructions.
    a=(16);
    for step=1:length(a)
        %We define an utrafine grid to interpolate on to.
        [fineX, fineY]=meshgrid(0:1/a(step):128,0:1/a(step):128);

        %We interpolate on to the grid.
        fineV=griddata(sampX,sampY,sampV,fineX,fineY,'cubic');
        
        %If any values are NaN, we set them to 0.
        fineV(isnan(fineV))=0;
        assignin('base','fineV',fineV);
        %Next, we interpolate on to the twice fine Cartesian grid with .5
        %spacing. This is called oversampling.
        stdV=zeros(256,256);
        for i=1:256
            for j=1:256
                %We use neighboring points to assign values to the
                %standard grid. We currently use a simple average.
                z=fineV((a(step)*i*.5)-(a(step)/4):(a(step)*i*.5),(a(step)*j*.5)-(a(step)/4):(a(step)*j*.5));  
                stdV(i,j)=sum(sum(z))/(((a(step)/4)+1)^2);
            end
        end
              
        assignin('base','stdV',stdV);
        
        %We take the modulus of the 2DFFT to obtain image data
        Co = fftshift(fft2(stdV));
        ImageMat = sqrt((real(Co).^2) + (imag(Co).^2))';
        ImageMat = rot90(ImageMat, 2);
        
        %We select only the center portion of the image because we
        %oversampled in the frequency domain.
        ImageMat=ImageMat(65:192,65:192);

        %This makes the unknown parts of the image black.
        black = makeBowl(128);
        ImageMat=ImageMat.*black;

        %This saves the image.
        assignin('base','im',ImageMat);
        binTitle = sprintf('oaklandOversampling%sScan%s.bin',class,scan);
        binTitle = fullfile('/Users/yangxia/Desktop/SimonMiller/OaklandReconstructionToolbox/Images',binTitle);
        fid=fopen(binTitle,'w');
        fwrite(fid,ImageMat,'double');
        fclose(fid);
    end
end