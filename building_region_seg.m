% disk_location = 'Z';

% Read LiDAR data
% txtfile = strcat(disk_location,':\workingData\Quebec data\LiDAR\nineth1_txyzic.txt');
txtfile = '/Users/thanhhuynguyen/Storage/DATA/LiDAR/nineth1_txyzic.txt';
        
fileID = fopen(txtfile,'r');
formatSpec = '%f,%f,%f,%f,%d,%d';
sizeA=[6 Inf];

A=fscanf(fileID, formatSpec, sizeA);
fclose(fileID);

ptCloudDBMat = A(2:5,:)'; %[XYZI]
ptCloudDBClass = A(6,:)'; %[C]

figure
pcshow(ptCloudDBMat(:,1:3),ptCloudDBClass,'MarkerSize',20)
view(2)
% colormap('jet')
title('1: Unclassified, 2: Ground, 3: Low vege., 4: Med. vege., 5: High vege.')


figure
pcshow(ptCloudDBMat(ptCloudDBClass==3|ptCloudDBClass==4|ptCloudDBClass==5,1:3),...
    'MarkerSize',20)
view(2)

points=ptCloudDBMat;
classes=ptCloudDBClass;
xmin=min(points(:,1));
xmax=max(points(:,1));
ymin=min(points(:,2));
ymax=max(points(:,2));

% Thresholding on z values
win = 20; % window of 20 meters
points_non_gnd = [];
points_bd = [];
points_non_bd = [];
nb_bin=200;
for x=xmin:win:xmax-win
    for y=ymin:win:ymax-win
        
        %% Thresholding using a value from the height of gnd points
        [pts,idx] = findPoints(points,[x x+win],[y y+win],0);
        c=classes(idx,:);
        
        if ~isempty(pts(c==2,:))
            zmean_gnd_pts = mean(pts(c==2,3));
            zstd_gnd_pts = std(pts(c==2,3));
            thres=zmean_gnd_pts+max(3,3*zstd_gnd_pts); 
            % cf. [M. Awrangjeb et al., 2013, Rule-based Segmentation...]
        else
            thres=0;
        end
        
%         display([num2str(x) ', ' num2str(y) ', ' num2str(thres)])
        idx2=find(pts(:,3)>thres);
        pts_non_ground = pts(idx2,:);
        c_ng = c(idx2,:);
        
        points_non_gnd = [points_non_gnd; pts_non_ground];
        
        %% remove points that are close to a tree
        pts_bd=[]; % to stock building points
        pts_non_bd=[]; % to stock remaining non building points
        if (~isempty(pts_non_ground))
            w=0.3; % 0.3 meters of tree allowed
            for i=1:size(pts_non_ground,1)
                if c_ng(i)==1
                    xyz=pts_non_ground(i,1:3);

                    [~,idx] = findPoints(pts_non_ground((c_ng==3)|(c_ng==4)|(c_ng==5),1:3),...
                        [xyz(1)-w, xyz(1)+w],[xyz(2)-w, xyz(2)+w],0);
                    if (isempty(idx))
                        pts_bd = [pts_bd; pts_non_ground(i,:)];
                    else
                        pts_non_bd = [pts_non_bd; pts_non_ground(i,:)];
                    end
                else
                    pts_non_bd = [pts_non_bd; pts_non_ground(i,:)];
                end
            end
        end
        points_bd = [points_bd; pts_bd];
        points_non_bd = [points_non_bd; pts_non_bd];
    end
end


figure
    pcshow(points_non_bd(:,1:3),'g')
    hold on
    plot3(points_bd(:,1),points_bd(:,2),points_bd(:,3),'ro')
    view(2)
    legend('Non-ground tree', 'Building')
    
% Use the built-in pointcloud denoising function
ptCloudIn = pointCloud(points_bd(:,1:3));
[~, inl] = pcdenoise(ptCloudIn);

pts_non_gnd_denoised=points_bd(inl,:);

figure
    pcshow(pts_non_gnd_denoised(:,1:3))
    view(2)
    
    
%% Create building mask (binary image)

% without camera pose parameters (i.e. using vertical projection)

xmax=max(pts_non_gnd_denoised(:,1));
xmin=min(pts_non_gnd_denoised(:,1));
ymax=max(pts_non_gnd_denoised(:,2));
ymin=min(pts_non_gnd_denoised(:,2));

resolution = [1 1];% [1 1];
Ncol=round((xmax-xmin)/resolution(1));
Nrow=round((ymax-ymin)/resolution(2));
step=[(xmax-xmin)/Ncol, (ymax-ymin)/Nrow];

mask = zeros(Nrow,Ncol);
for k=1:size(pts_non_gnd_denoised,1)
    P=pts_non_gnd_denoised(k,:);
    
    col=ceil((P(1)-xmin)/step(1))+1;
    row=ceil((P(2)-ymin)/step(2))+1;
    mask(row,col)=1;    
end

figure
    subplot(121)
    imagesc(mask)
    colormap(gray)
    axis equal
    
xWorldLimits = [xmin xmax];
yWorldLimits = [ymin ymax];
mask_ref=maprefcells(xWorldLimits, yWorldLimits, [Nrow,Ncol]);

%% applying a opening morphological operation
% se = [0 1 0; 1 1 1; 0 1 0];
se = strel('diamond',2);
mask2 = imopen(mask,se);

    subplot(122)
    imagesc(mask2)
    colormap(gray)
    axis equal
    
se = strel('diamond',1);
mask3=imclose(mask2,se);

%% Labeling the building binary mask
[L,n]=bwlabel(mask3,8);
% [L,n]=bwlabel(mask2,8);

figure
imagesc(L)

L2=zeros(size(L));
n2=zeros(size(n));
k=1;

for i=1:n
    [r, c] = find(L==i);
    if (numel(r)>20) %only take regions that have more than 20 pixels
        for j=1:numel(r)
            L2(r(j),c(j))=k;
        end
        k=k+1;
    end
end
n2=k-1;

figure
imagesc(L2)
axis equal
axis image
axis xy

%% saving the results

building_mask = mask2; %mask3;
building_label = L2;
nb_labels = n2;
save('building_mask_nineth2p.mat','building_mask','	','pts_non_gnd_denoised',...
    'building_label','nb_labels')

