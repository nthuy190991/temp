% Load datasets from ISPRS Benchmark datasets
% system=1: Macbook
% system=2: Windows PC
function [points,imgRGB,ri,rx,ry]=load_datasets_isprs(iReg,system,verbose)

dir=pwd;

if system==1
    dir_code='/Users/thanhhuynguyen/Dropbox/Thesis/dev/';
    dir_data='/Users/thanhhuynguyen/Storage/DATA/ISPRS_Benchmark/';
elseif system==2
    dir_code='C:\Users\thnguyen\Dropbox\Thesis\dev\';
    dir_data='D:\workingData\ISPRS_Vaihingen\';
end

%% Load datasets
% Load the optical image
A=imread(strcat(dir_data,'Ortho/TOP_Mosaic_09cm.tif'));
 
R = maprefcells();
R.CellExtentInWorldX=0.09;
R.CellExtentInWorldY=0.09;

X0=496274.985000;
Y0=5420850.075000;
RasterSize=[size(A,1) size(A,2)];

R.XWorldLimits=[X0 X0+RasterSize(2)*R.CellExtentInWorldX];
R.YWorldLimits=[Y0-RasterSize(1)*R.CellExtentInWorldY Y0];
R.RasterSize=[size(A,1) size(A,2)];
R.ColumnsStartFrom='north';

% figure
% mapshow(A,R)

% [x,y]=ginput(2)

%% Values from .tfw 
% Line 1: pixel size in the x-direction in map units (GSD).
% Line 2: rotation about y-axis. 
% Line 3: rotation about x-axis.
% Line 4: pixel size in the y-direction in map in map units (GSD).
% Line 5: x-coordinate of the upper left corner of the image.
% Line 6: y-coordinate of the upper left corner of the image.

% 0.090000
% 0.000000
% 0.000000
% -0.090000
% 496274.985000
% 5420850.075000


%% Load LiDAR point cloud
if iReg<4 || iReg==5
    fileID = fopen(strcat(dir_data,'ALS/area', num2str(iReg),'_txyzic.txt'),'r');

    formatSpec = '%f,%f,%f,%f,%d,%d';
    size_pts=[6 Inf];

    tmp=fscanf(fileID, formatSpec, size_pts);
    fclose(fileID);

    % ptCloudDBMat = tmp(2:5,:)'; %[XYZI]
    % ptCloudDBClass = tmp(6,:)'; %[C]

    points=tmp(2:6,:)'; %[XYZIC]
elseif iReg==6
    fileID = fopen(strcat(dir_data,'ALS/building49.txt'),'r');

    formatSpec = '%f,%f,%f,%f,%d,%d';
    size_pts=[6 Inf];

    tmp=fscanf(fileID, formatSpec, size_pts);
    fclose(fileID);
    
    points=tmp(2:6,:)'; %[XYZIC]
else
    points=[];
end
%% Get datasets on the considered region

row=[9800 11950; 14300 16550; 15820 18320; 9200 19000; 13800 14300; 8195 8706];
col=[9200 10570; 6700 8640; 8800 10450; 5000 12000; 12150 12550; 13678 14189];
     
% iReg=2;
imgRGB = A(row(iReg,1):row(iReg,2), col(iReg,1):col(iReg,2), :);
rx=R.XWorldLimits(1)+(col(iReg,1):col(iReg,2))*R.CellExtentInWorldX;
ry=R.YWorldLimits(2)-(row(iReg,1):row(iReg,2))*R.CellExtentInWorldY;

%small raster
ri=R;
ri.XWorldLimits = [rx(1), rx(end)];
ri.YWorldLimits = [ry(end), ry(1)];
ri.RasterSize = size(imgRGB);
if verbose
figure
    mapshow(rgb2gray(A),R)
    hold on
    mapshow(imgRGB,ri)
end

% cd(dir_code)
% [pts, ~]=findPoints(points,[rx(1), rx(end)],[ry(1), ry(end)],0);
% if verbose
% figure
%     mapshow(imgRGB,ri)
%     hold on
%     pcshow(pts(:,1:3),'MarkerSize',20)
%     view(2)

if verbose
figure
    mapshow(imgRGB,ri)
    hold on
    pcshow(points(:,1:3),'MarkerSize',20)
    view(2)
end


end