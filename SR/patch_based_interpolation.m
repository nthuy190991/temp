%% Load datasets
iReg=2;
[pts,imgRGB,ri,rx,ry]=load_datasets(iReg,1,0);

opt_img=rgb2gray(imgRGB);

szPatch=[500 550];
nbPatchy=floor(size(imgRGB,1)/szPatch(1));
nbPatchx=floor(size(imgRGB,2)/szPatch(2));
nbPatch=nbPatchy*nbPatchx;

nb_row=size(imgRGB,1);

patch={};
patch_ref={};
patch_pts={};
patch_row_col={};

load corr4TME_reg2_byGTM.mat
[~,pts_proj]=TME_GS_algo_wo_centering(pts,point2d',point3d',0);

for iPatch=1:nbPatch
    iPatch
    row=(floor((iPatch-1)/nbPatchx)*szPatch(1)+1):(floor((iPatch-1)/nbPatchx)+1)*szPatch(1);
    if mod(iPatch,nbPatchx)==0
        col=(nbPatchx-1)*szPatch(2)+1:nbPatchx*szPatch(2);
    else
        col=(max(1,(mod(iPatch,nbPatchx)-1)*szPatch(2)+1):mod(iPatch,nbPatchx)*szPatch(2));
    end
    patch(iPatch)={imgRGB(row,col,:)};
    
    patch_row_col(iPatch)={[min(row) min(col)]};

    ri_small=ri;
    ri_small.RasterSize=[length(row) length(col)];
    ri_small.XWorldLimits=[rx(col(1)) rx(col(end))];
    ri_small.YWorldLimits=[ry(nb_row-row(end)) ry(nb_row-row(1))];
    patch_ref(iPatch)={ri_small};
    
    [~,index]=findPoints(pts_proj,[rx(col(1)) rx(col(end))],[ry(nb_row-row(end)) ry(nb_row-row(1))], 0);
    pts_small=pts(index,:);
    patch_pts(iPatch)={pts_small};
end

%% Check if patches are well divided
% figure
% k=1;
% for i=1:nbPatchx
%     for j=1:nbPatchy
%         subplot(nbPatchy,nbPatchx,k)
%         mapshow(patch{k},patch_ref{k})
%         title(['patch ' num2str(k)])
%         k=k+1;
%     end
% end
% 
% figure
% k=1;
% for i=1:nbPatchx
%     for j=1:nbPatchy
%         subplot(nbPatchy,nbPatchx,k)
%         pts_small=patch_pts{k};
%         scatter(pts_small(:,1),pts_small(:,2),10,pts_small(:,3),'filled')
%         title(['patch ' num2str(k)])
%         k=k+1;
%     end
% end
% 
% k=1; Ix=[];
% for i=1:nbPatchy
%     Iy=[];
%     for j=1:nbPatchx
%         display(size(patch{k}))
%         Iy=[Iy rgb2gray(patch{k})];
%         k=k+1;
%     end
%     Ix=[Ix; Iy];
% end
% 
%     
% 
% figure
% k=1;
% for i=1:nbPatchx
%     for j=1:nbPatchy
%         subplot(nbPatchy,nbPatchx,k)
%         mapshow(patch{k},patch_ref{k})
%         hold on
%         pts_small=patch_pts{k};
%         pcshow(pts_small(:,1:3))
%         view(2)
%         axis off
%         title(['patch ' num2str(k)])
%         k=k+1;
%     end
% end


%%


%% Get the transformation model

[c,r]=convertGeoCoorMatCoor(ri.XWorldLimits,ri.YWorldLimits,point2d(:,1),point2d(:,2),ri.RasterSize);

figure
imshow(rgb2gray(imgRGB))
hold on
for i=1:length(c)
    plot(c(i),r(i),'rx')
end

xy=[c r]';
XYZ=point3d';
[P, K, R, t, C, error] = runGoldStandard(xy, XYZ);

theta_x = atan2d(R(3,2), R(3,3));
theta_y = atan2d(-R(3,1), sqrt(R(3,2)^2 + R(3,3)^2));
theta_z = atan2d(R(2,1), R(1,1));
R-rotz(theta_z)*roty(theta_y)*rotx(theta_x)

[~,N] = size(xy);

xy = [xy;ones(1,N)];
XYZ = [XYZ;ones(1,N)];

xy_est = K*[R t]*XYZ;

for i = 1:N
    xy_est(:,i) = xy_est(:,i)/xy_est(end,i);
end
error=norm(xy-xy_est,'fro');

figure
imshow(rgb2gray(imgRGB))
hold on
for i=1:length(c)
    plot(c(i),r(i),'rx')
    plot(xy_est(1,i),xy_est(2,i),'gx')
end



[x,y]=convertMatCoorGeoCoor(ri.XWorldLimits,ri.YWorldLimits,xy_est(1,:),xy_est(2,:),ri.RasterSize);
norm([x',y']-point2d,'fro')

cd SR


%% coarse registration parameters
figure
imshow(imgRGB)
hold on
% for iPatch=1:nbPatch
%     iPatch
%     pts_small=patch_pts{iPatch};
%     opt_img=rgb2gray(patch{iPatch});
%     row_col=patch_row_col{iPatch};
%     load('temp_6params_9patches.mat','K')
%     x0=patch_result{iPatch,1};
    x0=[theta_x theta_y theta_z; C(1:3)'];
    R=rotz(x0(1,3))*roty(x0(1,2))*rotx(x0(1,1)); %first line of x0: yaw, pitch, roll
    C=x0(2,:)'; %second line of x0: X, Y, Z from translation
    P=K*R*[eye(3) -C(1:3)];
    
    XYZ=pts(:,1:3)';
    n=size(XYZ,2);
    xy_est=P*[XYZ; ones(1,n)];
    for i = 1:n
        xy_est(:,i) = xy_est(:,i)/xy_est(end,i);
    end
    pts_proj=zeros(n,3);
    pts_proj(:,1:2)=xy_est(1:2,:)';
    pts_proj(:,3)=pts(:,3);
    pcshow(pts_proj)
% end
view(2)

% figure
% imshow(imgRGB)
% hold on
% pcshow()

%% bounded fminsearch

load fminsearch_theta_6params_9patches_MI
figure
imshow(imgRGB)
hold on
for iPatch=1:nbPatch
    iPatch
    pts_small=patch_pts{iPatch};
    opt_img=rgb2gray(patch{iPatch});
    row_col=patch_row_col{iPatch};
    
%     load('temp_6params_9patches.mat','K')
    
    x0=patch_result{iPatch,1};
    % x0=[theta_x theta_y theta_z; C(1:3)'];
    R=rotz(x0(1,3))*roty(x0(1,2))*rotx(x0(1,1)); %first line of x0: yaw, pitch, roll
    C=x0(2,:)'; %second line of x0: X, Y, Z from translation
    P=K*R*[eye(3) -C(1:3)];
    % disp(P)
    XYZ=pts_small(:,1:3)';
    n=size(XYZ,2);
    xy_est=P*[XYZ; ones(1,n)];
    for i = 1:n
        xy_est(:,i) = xy_est(:,i)/xy_est(end,i);
    end
    pts_proj=zeros(n,3);
    pts_proj(:,1:2)=xy_est(1:2,:)';
%     pts_proj(:,1)=pts_proj(:,1)-300;
%     pts_proj(:,2)=-pts_proj(:,2)-250;
    pts_proj(:,3)=pts_small(:,3);
    pcshow(pts_proj)
end
view(2)


%% weighted

x=zeros(1,nbPatch);
y=zeros(1,nbPatch);
for iPatch=1:nbPatch
    ri_small=patch_ref{iPatch};
    x(iPatch)=mean(ri_small.XWorldLimits);
    y(iPatch)=mean(ri_small.YWorldLimits);
end

nb_refs=9; % number of patches considered as references

pts_proj=[];

for iPatch=1:nbPatch
    iPatch
    pts_small=patch_pts{iPatch};
    opt_img=rgb2gray(patch{iPatch});
    row_col=patch_row_col{iPatch};
%     load('temp_6params_9patches.mat','K')
        
    for i=1:length(pts_small)
        dist=zeros(1,nbPatch);
        for j=1:nbPatch
            dist(j)=norm(pts_small(i,1:2)-[x(j) y(j)]);
        end
        [dk,idx]=mink(dist,nb_refs);
        
        weight=1./dk.^2;
        
        tb_x0_theta=zeros(nb_refs,3);
        tb_x0_C=zeros(nb_refs,3);
        for j=1:nb_refs
            temp=patch_result{idx(j),1};
            tb_x0_theta(j,:)=temp(1,:);
            tb_x0_C(j,:)=temp(2,:);
        end
        
        
        x0_theta=zeros(1,3);
        x0_C=zeros(1,3);
        if any(dk==0)
            x0_theta=tb_x0_theta(idx(1),:);
            x0_C=tb_x0_C(idx(1),:);
        else
            for j=1:3
                x0_theta(j)=sum(weight'.*tb_x0_theta(:,j))/sum(weight);
                x0_C(j)=sum(weight'.*tb_x0_C(:,j))/sum(weight);
            end
        end
        
        x0=[x0_theta; x0_C];
        
        R=rotz(x0(1,3))*roty(x0(1,2))*rotx(x0(1,1)); %first line of x0: yaw, pitch, roll
        C=x0(2,:)'; %second line of x0: X, Y, Z from translation
        P=K*R*[eye(3) -C(1:3)];
        
        % disp(P)
        XYZ=pts_small(i,1:3)';
        n=size(XYZ,2);
        xy_est=P*[XYZ; ones(1,n)];
        for j = 1:n
            xy_est(:,j) = xy_est(:,j)/xy_est(end,j);
        end
        p=zeros(n,3);
        p(:,1:2)=xy_est(1:2,:)';
        p(:,3)=pts_small(i,3);
        
        pts_proj=[pts_proj; p];
    end
end

pts_IDW_proj=pts_proj;

figure
imshow(imgRGB)
hold on
pcshow(pts_IDW_proj)
view(2)


%% some results

figure
hold on
for iPatch=1:nbPatch
    row_col=patch_row_col{iPatch};
    plot([row_col(2),row_col(2)+szPatch(2)],[row_col(1),row_col(1)],'r--','LineWidth',3)
    plot([row_col(2),row_col(2)],[row_col(1),row_col(1)+szPatch(1)],'r--','LineWidth',3)
    plot([row_col(2)+szPatch(2),row_col(2)+szPatch(2)],[row_col(1),row_col(1)+szPatch(1)],'r--','LineWidth',3)
    plot([row_col(2),row_col(2)+szPatch(2)],[row_col(1)+szPatch(1),row_col(1)+szPatch(1)],'r--','LineWidth',3)
end

load fminsearch_theta_6params_9patches_MI

for iPatch=1:nbPatch
    iPatch
    pts_small=patch_pts{iPatch};
    opt_img=rgb2gray(patch{iPatch});
    row_col=patch_row_col{iPatch};
    
%     load('temp_6params_9patches.mat','K')
    
    x0=patch_result{iPatch,1};
    % x0=[theta_x theta_y theta_z; C(1:3)'];
    R=rotz(x0(1,3))*roty(x0(1,2))*rotx(x0(1,1)); %first line of x0: yaw, pitch, roll
    C=x0(2,:)'; %second line of x0: X, Y, Z from translation
    P=K*R*[eye(3) -C(1:3)];
    % disp(P)
    XYZ=pts_small(:,1:3)';
    n=size(XYZ,2);
    xy_est=P*[XYZ; ones(1,n)];
    for i = 1:n
        xy_est(:,i) = xy_est(:,i)/xy_est(end,i);
    end
    pts_proj=zeros(n,3);
    pts_proj(:,1:2)=xy_est(1:2,:)';
%     pts_proj(:,1)=pts_proj(:,1)-300;
%     pts_proj(:,2)=-pts_proj(:,2)-250;
    pts_proj(:,3)=pts_small(:,3);
    pcshow(pts_proj,'MarkerSize',20)
end
axis ij
axis off
view(2)

xlim([243 516])
ylim([430 676])


%%

figure
hold on
for iPatch=1:nbPatch
    row_col=patch_row_col{iPatch};
    plot([row_col(2),row_col(2)+szPatch(2)],[row_col(1),row_col(1)],'r--','LineWidth',1)
    plot([row_col(2),row_col(2)],[row_col(1),row_col(1)+szPatch(1)],'r--','LineWidth',1)
    plot([row_col(2)+szPatch(2),row_col(2)+szPatch(2)],[row_col(1),row_col(1)+szPatch(1)],'r--','LineWidth',1)
    plot([row_col(2),row_col(2)+szPatch(2)],[row_col(1)+szPatch(1),row_col(1)+szPatch(1)],'r--','LineWidth',1)
end

pcshow(pts_IDW_proj)
view(2)