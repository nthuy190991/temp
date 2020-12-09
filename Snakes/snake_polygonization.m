%% Polygonization
load building_idx202_snakes_res
% load building_idx195_snakes_res
% load building_idx207_snakes_res


figure, imshow(I);
hold on; 
plot([P(:,2);P(1,2)],[P(:,1);P(1,1)],'cv-');
plot([O_rve_df(:,2);O_rve_df(1,2)],[O_rve_df(:,1);O_rve_df(1,1)],'r*-','LineWidth',2);

% mbr_P = minBoundingBox([P(:,2)'; P(:,1)']);
% mbr_O = minBoundingBox([O_rve_df(:,2)'; O_rve_df(:,1)']);
% P=fliplr(P);
O=fliplr(O_rve_df);


mbr_P = minBoundingBox([P(:,1)'; P(:,2)']);
mbr_O = minBoundingBox([O(:,1)'; O(:,2)']);

plot([mbr_P(1,1:4) mbr_P(1,1)],[mbr_P(2,1:4) mbr_P(2,1)],'c-');
plot([mbr_O(1,1:4) mbr_O(1,1)],[mbr_O(2,1:4) mbr_O(2,1)],'r-','LineWidth',2);

%% Get the dominant angle from the longer edge of MBR
if norm([mbr(2,2)-mbr(2,1),mbr(1,2)-mbr(1,1)])>norm([mbr(2,3)-mbr(2,2),mbr(1,3)-mbr(1,2)])
%     ang1=atan((mbr(2,2)-mbr(2,1))./(mbr(1,2)-mbr(1,1))); % dominant angle
%     ang2=atan((mbr(2,3)-mbr(2,2))./(mbr(1,3)-mbr(1,2)));
    line=[mbr(1,1:2); mbr(2,1:2)];
else
%     ang1=atan((mbr(2,3)-mbr(2,2))./(mbr(1,3)-mbr(1,2))); % dominant angle
%     ang2=atan((mbr(2,2)-mbr(2,1))./(mbr(1,2)-mbr(1,1)));
    line=[mbr(1,2:3); mbr(2,2:3)];
end

%% Look for the edges that are parallel to the dominant angle edge
[p,q]=get_major_line_mbr(O,mbr_O,0.5,50,1);


%% Try DP simplification (to be replaced by reducem built-in function)
Ps = dpsimplify(P,30);

figure, imshow(I);
hold on; 
plot([P(:,2);P(1,2)],[P(:,1);P(1,1)],'c-');
plot([Ps(:,2);Ps(1,2)],[Ps(:,1);Ps(1,1)],'r-');

Os = dpsimplify(O,10);

figure, imshow(I);
hold on; 
plot([O_rve_df(:,2);O_rve_df(1,2)],[O_rve_df(:,1);O_rve_df(1,1)],'c-');
plot([Os(:,2);Os(1,2)],[Os(:,1);Os(1,1)],'r*-','LineWidth',2);

%% Dutter's method on Snake result
cd Snakes
load building_idx207_snakes_res

figure
title('Snake models result')
imshow(I);
hold on; 
plot([P(:,2);P(1,2)],[P(:,1);P(1,1)],'c.-','LineWidth',2,'MarkerSize',10);
plot([O_basic(:,2);O_basic(1,2)],[O_basic(:,1);O_basic(1,1)],'y.-','LineWidth',2,'MarkerSize',10);
plot([O_gvf(:,2);O_gvf(1,2)],[O_gvf(:,1);O_gvf(1,1)],'b.-','LineWidth',2,'MarkerSize',10);
plot([O_rve(:,2);O_rve(1,2)],[O_rve(:,1);O_rve(1,1)],'g.-','LineWidth',2,'MarkerSize',10);
plot([O_df(:,2);O_df(1,2)],[O_df(:,1);O_df(1,1)],'r.-','LineWidth',4,'MarkerSize',10);
plot([O_rve_df(:,2);O_rve_df(1,2)],[O_rve_df(:,1);O_rve_df(1,1)],'m.-','LineWidth',2,'MarkerSize',10);
legend('Initial seed','Basic','GVF','GVF+RVE','GVF+DF','GVF+RVE+DF')

O=(O_df);
Os = dpsimplify(O,5); % Douglas-Peucker simplification

cd ..
Vspec=20;
Mspec=100;
[vertices,~,S,Psub]= polygonization_dutty(O,O,Vspec,Mspec,0,1);
[vertices_new,~,S,Psub]= polygonization_dutty(O,P,Vspec,Mspec,0,1);

figure, imshow(I);
hold on; 
plot([O(:,2);O(1,2)],[O(:,1);O(1,1)],'y+','MarkerSize',10);
plot([Os(:,2);Os(1,2)],[Os(:,1);Os(1,1)],'gv-','LineWidth',2);
plot([vertices(2,:) vertices(2,1)],[vertices(1,:) vertices(1,1)],'bo-','LineWidth',2)
plot([vertices_new(2,:) vertices_new(2,1)],[vertices_new(1,:) vertices_new(1,1)],'rs-','LineWidth',2)
legend('Snake points','Douglas-Peucker','Dutter''s algorithm','Improved algorithm')

control=ginput(2)
angle_gt=atand((control(2,2)-control(1,2))/(control(2,1)-control(1,1)))

angle_bad=atand((vertices(2,1)-vertices(2,2))/(vertices(1,1)-vertices(1,2)))
angle_good=atand((vertices_new(2,1)-vertices_new(2,2))/(vertices_new(1,1)-vertices_new(1,2)))
abs(angle_good)-angle_gt
abs(angle_bad)-angle_gt
% Try increase snake point density
figure
hold on; 
plot([O(:,2);O(1,2)],[O(:,1);O(1,1)],'rx');
grid on, axis equal
hold on
for i=2:size(O,1)
    text(O(i,2),O(i,1),num2str(i))
end

nb_pts_add=2;
O_inc=increase_snake_density(O,nb_pts_add);

figure
hold on; 
plot([O(:,2);O(1,2)],[O(:,1);O(1,1)],'rx');
grid on, axis equal
plot(O_inc(:,2),O_inc(:,1),'b+');
for i=2:size(O_inc,1)
    text(O_inc(i,2),O_inc(i,1),num2str(i))
end

[vertices_inc,~,~,~]= polygonization_dutty(O_inc,Vspec,Mspec,0,1);

figure, imshow(I);
hold on; 
plot([O(:,2);O(1,2)],[O(:,1);O(1,1)],'c-');
plot([vertices(2,:) vertices(2,1)],[vertices(1,:) vertices(1,1)],'rs-','LineWidth',2)
plot([vertices_inc(2,:) vertices_inc(2,1)],[vertices_inc(1,:) vertices_inc(1,1)],'b+-','LineWidth',2)

legend('Snake result','Dutter''s method on snake','Dutter''s method on density-increased snake')


%% Use Dutter's method on LiDAR filtered building points
cd Snakes
load building_idx207_snakes_res

cd ..
Vspec=10;
Mspec=50;
[vertices,points,~,~]=polygonization_dutty(P,Vspec,Mspec,10,0);

figure, imshow(I);
hold on; 
plot([P(:,2);P(1,2)],[P(:,1);P(1,1)],'c-');
plot([vertices(2,:) vertices(2,1)],[vertices(1,:) vertices(1,1)],'bs-','LineWidth',2)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Use Dutter's method on LiDAR filtered building points and snake results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
cd Snakes
load building_idx202_snakes_res

figure
title('Snake models result')
imshow(I);
hold on; 
plot([P(:,2);P(1,2)],[P(:,1);P(1,1)],'c.-','LineWidth',2,'MarkerSize',10);
plot([O_basic(:,2);O_basic(1,2)],[O_basic(:,1);O_basic(1,1)],'y.-','LineWidth',2,'MarkerSize',10);
plot([O_gvf(:,2);O_gvf(1,2)],[O_gvf(:,1);O_gvf(1,1)],'b.-','LineWidth',2,'MarkerSize',10);
plot([O_rve(:,2);O_rve(1,2)],[O_rve(:,1);O_rve(1,1)],'g.-','LineWidth',2,'MarkerSize',10);
plot([O_df(:,2);O_df(1,2)],[O_df(:,1);O_df(1,1)],'r.-','LineWidth',2,'MarkerSize',10);
% plot([O_rve_df(:,2);O_rve_df(1,2)],[O_rve_df(:,1);O_rve_df(1,1)],'k.-','LineWidth',2,'MarkerSize',10);
legend('Initial seed','Basic','GVF','GVF+RVE','GVF+DF')%,'GVF+RVE+DF')

cd ..
%get the best snake result
O=(O_df);

Vspec=10;
Mspec=50;
[vertices_li,points_li,~,~]=polygonization_dutty(P,Vspec,Mspec,10,0);
% [vertices_op,points_op,~,~]=polygonization_dutty(O,Vspec,Mspec,10,0);

O_inc=increase_snake_density(O,2); % increase snake point density to improve dutter polygonization result
[vertices_op,points_op,~,~]= polygonization_dutty(O_inc,Vspec,Mspec,10,0);



figure
imshow(I)
hold on
plot([vertices_li(2,:) vertices_li(2,1)],[vertices_li(1,:) vertices_li(1,1)],'bs-','LineWidth',2)
plot([P(:,2);P(1,2)],[P(:,1);P(1,1)],'b.','MarkerSize',10);

plot([vertices_op(2,:) vertices_op(2,1)],[vertices_op(1,:) vertices_op(1,1)],'rv-','LineWidth',2)
plot([O(:,2);O(1,2)],[O(:,1);O(1,1)],'r.','MarkerSize',10);

legend('LiDAR polygones','LiDAR projected points','Optical image polygones','Snake points')

%% Vertex matching (by nearest rule)

% R, tol, col, row are from snake_building_extraction.m code
[x_li,y_li]=convertMatCoorGeoCoor(R.XWorldLimits,R.YWorldLimits,...
    vertices_li(2,:)+min(col)-tol,vertices_li(1,:)+min(row)-tol,R.RasterSize);

[x_op,y_op]=convertMatCoorGeoCoor(R.XWorldLimits,R.YWorldLimits,...
    vertices_op(2,:)+min(col)-tol,vertices_op(1,:)+min(row)-tol,R.RasterSize);

vertices_li_ref=[x_li; y_li];
vertices_op_ref=[x_op; y_op];

figure
mapshow(A,R)
hold on
plot([vertices_li_ref(1,:) vertices_li_ref(1,1)],[vertices_li_ref(2,:) vertices_li_ref(2,1)],'bs-','LineWidth',2)
plot([vertices_op_ref(1,:) vertices_op_ref(1,1)],[vertices_op_ref(2,:) vertices_op_ref(2,1)],'rv-','LineWidth',2)

z_avg=mean(b_proj(:,3)); % ??????
% Find the nearest vertex
if size(vertices_li_ref,2)==size(vertices_op_ref,2)
    N=size(vertices_li_ref,2);
    points2d=vertices_op_ref;
    
%     points2d=ginput(6);
%     points2d=points2d';
    points3d=zeros(3,N);
    for i=1:N
        d=zeros(1,N);
        for j=1:N
            d(j)=norm(points2d(:,i)-vertices_li_ref(:,j));
        end
        [~,idx]=min(d);
        x=vertices_li_ref(1,idx);
        y=vertices_li_ref(2,idx);
        
        points3d(:,i)=[x,y,z_avg];
    end
end

if size(points2d,2)<6
    new_corr=
    points2d=
end

figure
mapshow(A,R)
hold on
plot(points2d(1,:),points2d(2,:),'b*')
plot(points3d(1,:),points3d(2,:),'r*')
for i=1:N
    plot([points2d(1,i),points3d(1,i)],[points2d(2,i), points3d(2,i)],'g-')
end

%% Estimate the local transformation model
[P_local,b_local_proj]=TME_GS_algo(b,points2d,points3d);
[P_local,p_local_proj]=TME_GS_algo(p,points2d,points3d);

figure
mapshow(A,R)
hold on
plot(p_local_proj(:,1),p_local_proj(:,2),'r*')
%plot(b_proj(:,1),b_proj(:,2),'c*')

%% get global transformation model
[A, R] = geotiffread('/Users/thanhhuynguyen/Storage/DATA/Orthos_2016/2016_CMQ_243500-5181500.tif');
row=[1 2000; 1 1550; 3300 5700];
col=[1 2000; 2350 4000; 1500 3500];

iReg=2;
imgRGB = A(row(iReg,1):row(iReg,2), col(iReg,1):col(iReg,2), :);
rx=R.XWorldLimits(1)+(col(iReg,1):col(iReg,2))*R.CellExtentInWorldX;
ry=R.YWorldLimits(1)+(size(A,1)-fliplr(row(iReg,1):row(iReg,2)))*R.CellExtentInWorldY; %note: directions of row and y are opposite

%small raster
ri=R;
ri.XWorldLimits = [rx(1), rx(end)];
ri.YWorldLimits = [ry(1), ry(end)];
ri.RasterSize = size(imgRGB);

load corr4TME_reg2_byGTM
[~,b_global_proj]=TME_GS_algo_IGARSS(imgRGB,ri,b,point2d',point3d');
[~,p_global_proj]=TME_GS_algo_IGARSS(imgRGB,ri,p,point2d',point3d');



%% Compare global vs local transformation model
figure
mapshow(A,R)
hold on
plot(p_local_proj(:,1),p_local_proj(:,2),'r.')
plot(p_global_proj(:,1),p_global_proj(:,2),'c.')
legend('Local','Global')

