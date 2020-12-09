% building boundary extraction from meanshift segmentation

clear
close all
clc


%% Load image

filename = '/Users/thanhhuynguyen/Storage/DATA/Orthos_2016/2016_CMQ_243500-5181500.tif';
[A, R] = geotiffread(filename);

% 1st row: residential, 2nd: avg-sized buildings, 3rd: big and complex buildings
row=[1 2000; 1 1550; 3300 5700];
col=[1 2000; 2350 4000; 1500 3500];

%% Load meanshift result
matfile = 'reg1_bw3.5_1.mat';
load(matfile)

% seg_npts_all = seg_area_all;
% clear seg_area_all
nb_segs = size(bound_all,1);

seg_npts_all=zeros(1,nb_segs);
for i=1:nb_segs
    seg_npts_all(i)=size(mask_all{i},1);
end

s=strsplit(matfile,'_');
if strcmp(s{1},'reg1')
    iReg=1;
elseif strcmp(s{1},'reg2')
    iReg=2;
else
    iReg=3;
end

I = A(row(iReg,1):row(iReg,2), col(iReg,1):col(iReg,2), :);

figure
hist(seg_npts_all,100)

k=0;
c='rbymgc'
figure
imshow(rgb2gray(I))
hold on
for i=[216,217,247,248,255,276]%[1182,1315,3421,3873]%1:nb_segs%
    k=k+1;
%     if seg_npts_all(i)>300% && seg_npts_all(i)<2e4
        b=bound_all{i,:};
        p2=plot(b(:,2),b(:,1),'-','Color',c(k),'LineWidth',2.5);
%         p=mask_all{i,:};
%         plot(p(:,2),p(:,1),'.','Color',rand(3,1),'LineWidth',2);
%         center=[mean(b(:,2)), mean(b(:,1))];
%         text(center(1),center(2),num2str(i),'Color','r')
%     end
end
xlim([221 775])
ylim([995 1585])
%% Minimum Bounding Rectangle/Box (MBR) and Largest Empty Rectangle (LER)
segment_thres=[1000 7e4];

k=1;
tb_segs_inrange=[];
mbr={};
mbr_area=[];
seg_area=[]; 
seg_area_from_mask=[];
mbr_fill_percent=[];
mbr_fill_percent2=[];

% ler={};
ler_area=[];
ler_fill_percent=[];

ler_to_mbr=[];

for i=1:nb_segs
    if seg_npts_all(i)>segment_thres(1) && seg_npts_all(i)<segment_thres(2)
        
        tb_segs_inrange=[tb_segs_inrange i];
        
        %% Min bounding box
        % get min bounding box
        p=mask_all{i};
        bb = minBoundingBox([p(:,2)'; p(:,1)']);
        bb2=[bb bb(:,1)];
        mbr(k,:)={bb2};
%         bb2=mbr{k,:};

        % calculate bounding box area
        mbr_area = [mbr_area, polyarea(bb2(1,:), bb2(2,:))];
        
        % calculate segment area
        b=bound_all{i};
        seg_area=[seg_area, polyarea(b(:,2),b(:,1))];
        seg_area_from_mask=[seg_area_from_mask size(mask_all{i},1)];
%         
%         % calculate the % of filling of each segment inside its bounding box
        mbr_fill_percent=[mbr_fill_percent, seg_area(k)/mbr_area(k)];
        mbr_fill_percent2=[mbr_fill_percent2, seg_area_from_mask(k)/mbr_area(k)];
        
        %% Largest empty box
        % get largest empty box
        p2=mask_all{i};
        bb3 = maxEmptyBox([p2(:,2)'; p2(:,1)'],0);
        bb4=[bb3 bb3(:,1)];
        ler(k,:)={bb4};
%         bb4=ler{k,:};
        
         % calculate bounding box area
        ler_area = [ler_area, polyarea(bb4(1,:), bb4(2,:))];
        
        % calculate the % of filling of each segment inside its bounding box
        ler_fill_percent=[ler_fill_percent, ler_area(k)/seg_area(k)];
        
        
        ler_to_mbr=[ler_to_mbr ler_area(k)/mbr_area(k)];

        k=k+1;
    end
end

nb_segs_inrange=size(tb_segs_inrange,2);

save('reg2_Lab_res.mat','tb_segs_inrange','mbr','ler','ler_to_mbr')

%% Display results

figure
imshow(rgb2gray(I))
hold on
for i=1:nb_segs_inrange
    p=mask_all{tb_segs_inrange(i),:};
%     b=boundary(p,1);
    p0=plot(p(:,2),p(:,1),'.','Color',rand(3,1));
    
    bb=mbr{i};
    p1=plot(bb(1,:), bb(2,:), 'y-','LineWidth',1.5);
    
    bb2=ler{i,1};
    p2=plot(bb2(1,:), bb2(2,:), 'g-','LineWidth',1.25);
    
%     center=[mean(bb(1,1:4)), mean(bb(2,1:4))];
%     text(center(1),center(2),num2str(fill_percent(i),'%.2f'),'Color','r')
end
title('Segments with MBR and LER')
legend([p0 p1 p2],'segments','MBR','LER')

%% Display N highest filling percentage segments
n=10;
[~,idx]=maxk(mbr_fill_percent,n);
figure(90)
subplot(131)
imshow(rgb2gray(I))
hold on
for i=idx
    p=mask_all{tb_segs_inrange(i),:};
    plot(p(:,2),p(:,1),'.','Color',rand(3,1));

    bb=mbr{i};
    plot(bb(1,:), bb(2,:), 'r-','LineWidth',1)
end
title([num2str(n) ' segments that have highest MBR filling percentage'])

[~,idx]=maxk(ler_fill_percent,n);
figure(90)
subplot(132)
imshow(rgb2gray(I))
hold on
for i=idx
    p=mask_all{tb_segs_inrange(i),:};
    plot(p(:,2),p(:,1),'.','Color',rand(3,1));

    bb=ler{i};
    plot(bb(1,:), bb(2,:), 'r-','LineWidth',1)  
end
title([num2str(n) ' segments that have highest LER filling percentage'])

[~,idx]=maxk(ler_to_mbr,n);
figure(90)
subplot(133)
imshow(rgb2gray(I))
hold on
for i=idx
    p=mask_all{tb_segs_inrange(i),:};
    plot(p(:,2),p(:,1),'.','Color',rand(3,1));

    bb=mbr{i};
    plot(bb(1,:), bb(2,:), 'r-','LineWidth',1)
end
title([num2str(n) ' segments that have highest LER/MBR percentage'])

%%
% Display segments that have more than 70% of MBR filling percentage
figure(100)
subplot(131)
imshow(rgb2gray(I))
hold on
for i=1:nb_segs_inrange
    if mbr_fill_percent(i)>0.7
        p=mask_all{tb_segs_inrange(i),:};
        plot(p(:,2),p(:,1),'.','Color',rand(3,1));

        bb=mbr{i};
        plot(bb(1,:), bb(2,:), 'r-','LineWidth',1)
    end
end
title('Segments with S/MBR>70%')



% Display segments that have more than 70% of LER filling percentage
figure(100)
subplot(132)
imshow(rgb2gray(I))
hold on
for i=1:nb_segs_inrange
    if ler_fill_percent(i)>0.7
        p=mask_all{tb_segs_inrange(i),:};
        plot(p(:,2),p(:,1),'.','Color',rand(3,1));

        bb=ler{i};
        plot(bb(1,:), bb(2,:), 'r-','LineWidth',1)
    end
end
title('Segments with LER/S>70%')

% Display segments that have more than 50% of LER/MBR
figure(100)
subplot(133)
imshow(rgb2gray(I))
hold on
for i=1:nb_segs_inrange
    if ler_to_mbr(i)>0.5
        p=mask_all{tb_segs_inrange(i),:};
        plot(p(:,2),p(:,1),'.','Color',rand(3,1));

%         bb=ler{i};
%         plot(bb(1,:), bb(2,:), 'r-','LineWidth',1)
    end
end
title('Segments with LER/MBR>50%')

