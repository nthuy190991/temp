dir='/Users/thanhhuynguyen/Dropbox/Thesis/dev/';
dir2='/Users/thanhhuynguyen/Storage/DATA/';

for iReg=1:3
    
    cd ..
    [~,imgRGB,ri,rx,ry]=load_datasets_isprs(iReg,1,0);
    cd Snakes
    %% Metric result
    filename=strcat(dir2,'ISPRS_Benchmark/Reference_3d_reconstruction/Reference_Buildings/');
    [S,A]=shaperead(strcat(filename,'building_outline_area_',num2str(iReg),'.shp'));
    
    I=rgb2gray(imgRGB);
    
    % Ground truth building mask
    pt=zeros(numel(I),2);
    for i=1:size(I,1)
        for j=1:size(I,2)
            pt((i-1)*size(I,2)+j,:)=[i,j];
        end
    end
    
    ref_mask=zeros(size(I,1),size(I,2));
    for i=1:length(S)
        if (mean(S(i).X(1:end-1))<=ri.XWorldLimits(2) && mean(S(i).X(1:end-1))>=ri.XWorldLimits(1) &&...
                mean(S(i).Y(1:end-1))<=ri.YWorldLimits(2) && mean(S(i).Y(1:end-1))>=ri.YWorldLimits(1))
            
            [x,y]=convertGeoCoorMatCoor(ri.XWorldLimits,ri.YWorldLimits,...
                S(i).X(1:end-1),S(i).Y(1:end-1),ri.RasterSize);
            inBox = reshape(inpolygon(pt(:,1),pt(:,2),y,x),size(I,2),size(I,1));
            ref_mask = ref_mask | inBox';
        end
    end
    
    
    load(strcat('ISPRS_area',num2str(iReg),'_result.mat'),'snake_mask','snake_results','acc_snake')
    
    ref_mask_color=zeros(size(I,1),size(I,2));
    ref_mask_color(ref_mask==1&snake_mask==1)=1;
    ref_mask_color(ref_mask==0&snake_mask==1)=2;
    ref_mask_color(ref_mask==1&snake_mask==0)=3;
%     figure
%     imagesc(ref_mask_color)
    
    
    figure
    imshow(snake_mask)
    hold on
    for i=1:length(S)
        [x,y]=convertGeoCoorMatCoor(ri.XWorldLimits,ri.YWorldLimits,...
                S(i).X(1:end-1),S(i).Y(1:end-1),ri.RasterSize);

%         if (polyarea(S(i).X(1:end-1),S(i).Y(1:end-1))>50)
            plot(x,y,'LineWidth',3)
            text(mean(x,'omitnan'),mean(y,'omitnan'),num2str(i),'Color','g')
%         end
    end
    



    figure
    imagesc(~snake_mask)
    colormap(gray)
    hold on
    
    
    
    [y,x]=find(ref_mask_color==2);
    plot(x,y,'r.')
    
    [y,x]=find(ref_mask_color==3);
    plot(x,y,'b.')
    
    [y,x]=find(ref_mask_color==1);
    plot(x,y,'y.')
    
    axis equal
    axis off
    
end


