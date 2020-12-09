% Load datasets associated with a region 
% system=1: Macbook
% system=2: Windows PC
function [pts,imgRGB,ri,rx,ry]=load_datasets(iReg,system,verbose)

dir=pwd;

% if ~exist('iReg')
%     iReg=2;
% end
% 
% if ~exist('system')
%     system=1;
% end


if system==1
    dir_code='/Users/thanhhuynguyen/Dropbox/Thesis/dev/';
    dir_data='/Users/thanhhuynguyen/Storage/DATA/';
elseif system==2
    dir_code='C:\Users\thnguyen\Dropbox\Thesis\dev\';
    dir_data='D:\workingData\Quebec_data\';
end

%% Load datasets
% Load the optical image
[A, R] = geotiffread(strcat(dir_data,'Orthos_2016/2016_CMQ_243500-5181500.tif'));

% Load LiDAR point cloud
if iReg~=3
    fileID = fopen(strcat(dir_data,'LiDAR/nineth2p_txyzic.txt'),'r');
else
    fileID = fopen(strcat(dir_data,'LiDAR/nineth1_txyzic.txt'),'r');
end
formatSpec = '%f,%f,%f,%f,%d,%d';
size_pts=[6 Inf];

tmp=fscanf(fileID, formatSpec, size_pts);
fclose(fileID);

% ptCloudDBMat = tmp(2:5,:)'; %[XYZI]
% ptCloudDBClass = tmp(6,:)'; %[C]

points=tmp(2:6,:)'; %[XYZIC]

%% Get datasets on the considered region, i.e. reg2

% 1st row: residential, 2nd: avg-sized buildings, 3rd: big and complex buildings
row=[1 2000; 1 1550; 3300 5700];
col=[1 2000; 2350 4000; 1500 3500];
% iReg=2;

imgRGB = A(row(iReg,1):row(iReg,2), col(iReg,1):col(iReg,2), :);
rx=R.XWorldLimits(1)+(col(iReg,1):col(iReg,2))*R.CellExtentInWorldX;
ry=R.YWorldLimits(1)+(size(A,1)-fliplr(row(iReg,1):row(iReg,2)))*R.CellExtentInWorldY; %note: directions of row and y are opposite

%small raster
ri=R;
ri.XWorldLimits = [rx(1), rx(end)];
ri.YWorldLimits = [ry(1), ry(end)];
ri.RasterSize = size(imgRGB);
if verbose
figure
    mapshow(rgb2gray(A),R)
    hold on
    mapshow(imgRGB,ri)
end

cd(dir_code)
[pts, ~]=findPoints(points,[rx(1), rx(end)],[ry(1), ry(end)],0);
if verbose
figure
    mapshow(imgRGB,ri)
    hold on
    pcshow(pts(:,1:3),'MarkerSize',20)
    view(2)
end

cd(dir)

end