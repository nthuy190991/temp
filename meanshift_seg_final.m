%% Orthoimage
filename = '/Users/thanhhuynguyen/Storage/DATA/Orthos_2016/2016_CMQ_243500-5181500.tif';

[A, R] = geotiffread(filename);

% 1st row: residential, 2nd: avg-sized buildings, 3rd: big and complex buildings
row=[1 2000; 1 1550; 3300 5700; 3950 6100; 3940 5620];
col=[1 2000; 2350 4000; 1500 3500; 1750 3600; 1750 3490];
% bandwidth param
tb_bw=[2.5, 3.5, 4, 4, 4];

% flag of using 2 or 3-dimensional image pixel values
tb_flag_3d = [0 1];

% region removal thresholds
S_min=300;
Sr_thres=7;
Se_thres=0.85;
Sa_thres=[700 Inf];

close all
for iColorSpace = 1:1%2%length(tb_flag_3d)
    flag_3d=tb_flag_3d(iColorSpace)
    
    for iReg=5:5%3:3%size(row,1)
        iReg
        I = A(row(iReg,1):row(iReg,2), col(iReg,1):col(iReg,2), :);

        bandwidth=tb_bw(iReg)

        % L*a*b color space
        lab_he = rgb2lab(I);
        nrows = size(lab_he,1);
        ncols = size(lab_he,2);
        if (flag_3d)
            ab = lab_he(:,:,1:3);
            dat = reshape(ab,nrows*ncols,3);
        else
            ab = lab_he(:,:,2:3);
            dat = reshape(ab,nrows*ncols,2);
        end
        
        % meanshift clustering
        [clustCent,data2cluster,cluster2dataCell] = MeanShiftCluster(dat',bandwidth);
        pixel_labels = reshape(data2cluster,nrows,ncols);

        % segment and segment boundary extraction       
        bound_all=[];
        mask_all=[];
        seg_area_all=[];
        for l=1:numel(unique(pixel_labels))
            [bound, mask, seg_area]=segmentFiltering(pixel_labels,l,Sa_thres);
            
            bound_all=[bound_all; bound];
            mask_all=[mask_all; mask];
            seg_area_all=[seg_area_all; seg_area];
        end
        
        
        mat_filename=strcat('reg',num2str(iReg),'_bw',num2str(bandwidth),'_',num2str(flag_3d),'.mat');
        save(mat_filename,'pixel_labels','bound_all','mask_all','seg_area_all')

        figure%(iReg)
        imagesc(pixel_labels), axis equal
        
        figure%(10+iReg)
        imshow(rgb2gray(I))
        hold on
        for i=1:size(bound_all,1)
            b=bound_all{i,:};
            if ~isempty(b)
            p2=plot(b(:,2),b(:,1),'-','Color',rand(3,1),'LineWidth',2);
            end
        end
        title('Boundary extraction result')
        
        figure%(20+iReg)
        imshow(rgb2gray(I))
        hold on
        for i=1:size(mask_all,1)
            p=mask_all{i,:};
            plot(p(:,2),p(:,1),'.','Color',rand(3,1),'LineWidth',2);
        end
        title('Segmentation result')
        
    end
end
