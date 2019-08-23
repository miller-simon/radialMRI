function im = fbpSpatialFiltOld(KX,KY,main,class,scan)
%This is the filtered back projection function with spatial filtering for the Oakland Toolbox.
%It is recommended for the AVIIIHD scans over the standard filtered back
%projection function.

%This shifts the KX and KY vectors such that they begin at approximately 0.
KX=-1*KX(1)+KX+1E-10;
KY=-1*KY(1)+KY+1E-10;

%This sets the file path for a given scan's FID.
data = strcat(main,'/',class,'/',scan,'/fid');

%If the length of Kx is 131, we run the following to construct a 256x256
%image.
if length(KX) == 131
    %This opens the data.
    fileID = fopen(data, 'r','l');
    fidData = fread(fileID,[512,804],'int32')';
    fclose(fileID);

    %Anything beyond the 262 column is 0. We do not need this information.
    fidData = fidData(:,1:262);

    %Separate the FID into real and imaginary.
    Real = fidData(:,1:2:262-1);
    Imaginary = fidData(:,2:2:262);

    %Then combine them to get the complex data for each sample.
    DATA=Real+1i*Imaginary;
    
    %Our data is initally organized as radial readouts. We would like each
    %column to span the diameter of K-space.
    bMA = DATA(1:402,:);
    bMB = DATA(403:804,:);
    bMB = flip(bMB,2);
    DATA= [bMB bMA]';
    
    
    %This creates vectors with evenly spaced angles from [0,360). Each one
    %corresponds to a readout.
    THETA=linspace(0,360-(180/402),804);

    %This creates a vector of radii of a given readout.
    rad = sqrt(KX.^2+KY.^2);
    
    %We then transform it to act as a number line centered at 0.
    x=[-1*(flip(rad)) rad];
    
    %This defines evenly spaced query points along a diameter of K-space.
    xq=linspace(-1*max(rad),max(rad),256);
    
    %We use the sample data to interpolate to the query points for each
    %K-space diameter.
    interpData=zeros(256,402);
    for i=1:402
        v=DATA(:,i);
        interpData(:,i)=interp1(x,v,xq,'linear');
    end

    DATA=interpData;

    %We take the modulus of the 1DFFT to obtain projection information.
    PR=fftshift(fft(DATA));
    PR=sqrt(real(PR).^2+imag(PR).^2);
    
    %We transform the sinogram such that it is continuous and smooth.
    PR=[PR(:,202:402) PR(:,1:201)];
    PR(:,1)=PR(:,2);
    
    %This creates 360 degrees of projection data.
    PR=[PR flipud((PR))];
    PRright=[PR(256,403:804); PR(1:255,403:804)];
    PR=[PR(:,1:402) PRright];

    %This performs the filtered back projection.
    im=iradon(abs(PR),THETA,'linear','Ram-Lak',1,256);

    
    %This polynomial was fit based on a single tube scan.
    p =1.0e+03 *[ -0.0012    0.2985   -2.0157];
       
    %The polynomail is used to construct a spatial domain filter.
    y=polyval(p,1:256);
    fil=abs(y'*y);
    fil=fil+min(min(fil));
    fil(fil < fil(50,50))=fil(50,50);
    im=im./fil;
    im(im==inf)=0;
    
    %This makes the unknown parts of the image black.
    blk=makeBowl(256);
    im=im.*blk;
   
    %This sets values less than 0, produced by iradon, equal to 0.
    im(im<0)=0;
    
    %This saves the image.
    binTitle = sprintf('oaklandFBPSpatialFilt%sScan%s.bin',class,scan);
    binTitle = fullfile('/Users/yangxia/Desktop/SimonMiller/OaklandReconstructionToolbox/Images',binTitle);
    fid=fopen(binTitle,'w');
    fwrite(fid,im,'double');
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
    Real = fidData(:,1:2:134-1);
    Imaginary = fidData(:,2:2:134);
    
    %Then combine them to get the complex data for each sample.
    DATA=Real+1i*Imaginary;
    
    %Our data is initally organized as radial readouts. We would like each
    %column to span the diameter of K-space.
    bMA = DATA(1:201,:);
    bMB = DATA(202:402,:);
    bMB = flip(bMB,2);
    DATA= [bMB bMA]';

    %This creates vectors with evenly spaced angles from [0,360). Each one
    %corresponds to a readout.
    THETA=linspace(0,360-(180/201),402);
    
    %This creates a vector of radii of a given readout.
    rad = sqrt(KX.^2+KY.^2);
    
    %We then transform it to act as a number line centered at 0.
    x=[-1*(flip(rad)) rad];
    
    %This defines evenly spaced query points along a diameter of K-space.
    xq=linspace(-1*max(rad),max(rad),128);
    
    %We use the sample data to interpolate to the query points for each
    %K-space diameter.
    interpData=zeros(128,201);
    for i=1:201
        v=DATA(:,i);
        interpData(:,i)=interp1(x,v,xq,'linear');
    end

    DATA=interpData;
    
    %We take the modulus of the 1DFFT to obtain projection information.
    PR=fftshift(fft(DATA));
    PR=sqrt(real(PR).^2+imag(PR).^2);
    
    %We transform the sinogram such that it is continuous and smooth.
    PR=[PR(:,101:201) PR(:,1:100)];
    PR(:,1)=PR(:,2);
    %PR(:,101:201)=flip(PR(:,101:201));
    
    %This creates 360 degrees of projection data.
    PR=[PR flipud((PR))];
    PRright=[PR(128,202:402); PR(1:127,202:402)];
    PR=[PR(:,1:201) PRright];

    assignin('base','PR',PR);
    %This performs the filtered back projection.
    im=iradon(abs(PR),THETA,'linear','Ram-Lak',1,128);
    
    %This polynomial was fit based on a single tube scan.
    p =1.0e+03 *[-0.0047    0.5782   -1.5414];
    
    %The polynomail is used to construct a spatial domain filter.
    y=polyval(p,1:128);
    fil=abs(y'*y);
    fil=fil+min(min(fil));
    fil(fil < fil(25,25))=fil(25,25);
    im=im./fil;
    im(im==inf)=0;
    
    %This makes the unknown parts of the image black.
    blk=makeBowl(128);
    im=im.*blk;
    
     %This sets values less than 0, produced by iradon, equal to 0.
    im(im<0)=0;
    
    %This saves the image.
    binTitle = sprintf('oaklandFBPSpatialFilt%sScan%s.bin',class,scan);
    binTitle = fullfile('/Users/yangxia/Desktop/SimonMiller/OaklandReconstructionToolbox/Images',binTitle);
    fid=fopen(binTitle,'w');
    fwrite(fid,im,'double');
    fclose(fid);
end