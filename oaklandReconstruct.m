function oaklandReconstruct()
%This function provides access to each reconstruction method in the
%toolbox.

%This sets the file path for the folder where each scan's information will
%be stored.
main = '/Users/yangxia/Desktop/SimonMiller/OaklandReconstructionToolbox/Expts';
%The user will enter the input for system used and associated scan number.
disp(' ');
class = input('Enter 2 for AVII or 3 for AVIIIHD: ','s');
disp(' ');
scan = input('Enter Scan number: ', 's');

if class == '2'
    class='AVII';
else
    class='AVIIIHD';
end

disp(' ');

%This sets the path which contains the Kx,Ky information for the selected
%scan.
data = strcat(main,'/',class,'/', scan,'/method');

%We create vectors with the Kx and Ky coordinates for  a readout from the scan.
fid = fopen(data, 'rt');
TextAsCell=textscan(fid, '%s','Delimiter','%n');
fclose(fid);
TextAsCells = TextAsCell{1};

indexKx = find(contains(TextAsCells,'Kx'));
indexKx=indexKx(1)+1;

indexKy = find(contains(TextAsCells,'Ky'));
indexKy=indexKy(1)+1;

Kx=[TextAsCells{indexKx:indexKy-4,1}];

lengthK=indexKy-4-indexKx;

Kx=str2num(Kx);   %#ok<*ST2NM>

Ky=[TextAsCells{indexKy:indexKy+lengthK,1}];
Ky=str2num(Ky);

assignin('base','KX',Kx);
assignin('base','KY',Ky);

%This prompts the user to select a reconstruction method.
disp(' ');
disp('Select a method of reconstruction. Enter G for Standard Gridding,');
disp('O for Gridding with Oversampling, D for Gridding with Oversampling and Deapodization,')
type=input('F for Filtered Back Projection (FBP), or S for FBP with Spatial Filtering: ','s');

%This calls the selected method.
if type == 'G' || type == 'g'
    stdGriddingFIX(Kx,Ky,main,class,scan);
elseif type == 'O' || type == 'o'
    griddingOS(Kx,Ky,main,class,scan);
elseif type == 'D' || type == 'd'
    griddingDeapodize(Kx,Ky,main,class,scan);
elseif type == 'F' || type == 'f'
    stdFBP(Kx,Ky,main,class,scan);
elseif type == 'S' || type == 's'
    fbpSpatialFilt(Kx,Ky,main,class,scan);
elseif type =="to" || type=="TO"
    fov = input('Enter FOV size: ', 's');
    testOS(Kx,Ky,main,class,scan, fov);
elseif type =="tf" || type=="TF"
    fov = input('FOV size', 's');
    testFBP(Kx,Ky,main,class,scan, fov);
end