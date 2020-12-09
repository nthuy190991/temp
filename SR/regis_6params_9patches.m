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
    
    pts_small=findPoints(pts,[rx(col(1)) rx(col(end))],[ry(nb_row-row(end)) ry(nb_row-row(1))], 0);
    patch_pts(iPatch)={pts_small};
end

%% Check if patches are well divided
figure
k=1;
for i=1:nbPatchx
    for j=1:nbPatchy
        subplot(nbPatchy,nbPatchx,k)
        mapshow(patch{k},patch_ref{k})
        title(['patch ' num2str(k)])
        k=k+1;
    end
end

figure
k=1;
for i=1:nbPatchx
    for j=1:nbPatchy
        subplot(nbPatchy,nbPatchx,k)
        pts_small=patch_pts{k};
        scatter(pts_small(:,1),pts_small(:,2),10,pts_small(:,3),'filled')
        title(['patch ' num2str(k)])
        k=k+1;
    end
end

k=1; Ix=[];
for i=1:nbPatchy
    Iy=[];
    for j=1:nbPatchx
        display(size(patch{k}))
        Iy=[Iy rgb2gray(patch{k})];
        k=k+1;
    end
    Ix=[Ix; Iy];
end


% nb_row=size(imgRGB,1);
% opt_img=rgb2gray(imgRGB(row,col,:));
% pts_small=findPoints(pts,[rx(col(1)) rx(col(end))],[ry(nb_row-row(end)) ry(nb_row-row(1))], 0);
%
% ri_small=ri;
% ri_small.RasterSize=[length(row) length(col)];
% ri_small.XWorldLimits=[rx(col(1)) rx(col(end))];
% ri_small.YWorldLimits=[ry(nb_row-row(end)) ry(nb_row-row(1))];
%
% figure
% mapshow(opt_img,ri_small)
% hold on
% scatter(pts_small(:,1),pts_small(:,2),10,pts_small(:,3),'filled')

%% Get the transformation model
load corr4TME_reg2_byGTM
% [P0,pts_proj,avg]=TME_GS_algo(pts_small,point2d',point3d',0);

[c,r]=convertGeoCoorMatCoor(ri.XWorldLimits,ri.YWorldLimits,point2d(:,1),point2d(:,2),ri.RasterSize);

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

[x,y]=convertMatCoorGeoCoor(ri.XWorldLimits,ri.YWorldLimits,xy_est(1,:),xy_est(2,:),ri.RasterSize);
norm([x',y']-point2d,'fro')

cd SR

%% fminsearch on R and C (i.e. exterior orientation params)
patch_result={};
for iPatch=1:nbPatch
    iPatch
    pts_small=patch_pts{iPatch};
    opt_img=rgb2gray(patch{iPatch});
    row_col=patch_row_col{iPatch};
    
    save('temp_6params_9patches.mat','opt_img','pts_small','K','row_col')
    
    x0=[theta_x theta_y theta_z; C(1:3)']
    
    fun = @calc_mi_with_theta_6params;
    tic
    [x,fval,exitflag,output] = fminsearch(fun,x0)
    toc
    patch_result(iPatch,1)={x};
    patch_result(iPatch,2)={fval};
%     save(strcat('fminsearch_theta_6params_9patches_MI_p', num2str(iPatch),'.mat'),'patch','patch_ref','patch_pts','patch_result')
end

%% try bounded fminsearch
% patch_result={};

% load fminsearch_theta_6params_9patches_MI_all_bnd
% patch_result
patch_result={};
rangexyz=[2.5 2.5 2.5];
rangetheta=[0.25 0.25 1];
% rangexyz=[1 1 1];
% rangetheta=[0.1 0.1 0.5];
for iPatch=[1 2]%nbPatch“
    iPatch
    pts_small=patch_pts{iPatch};
    opt_img=rgb2gray(patch{iPatch});
    row_col=patch_row_col{iPatch};
    save('temp_6params_9patches.mat','opt_img','pts_small','K','row_col')
    x0=[theta_x theta_y theta_z; C(1:3)']
    x1=[[theta_x theta_y theta_z]-rangetheta; C(1:3)'-rangexyz];
    x2=[[theta_x theta_y theta_z]+rangetheta; C(1:3)'+rangexyz];
    fun = @calc_mi_with_theta_6params;
    tic
    [x,fval,exitflag,output] = fminsearchbnd(fun,x0,x1,x2)
    toc
    patch_result(iPatch,1)={x};
    patch_result(iPatch,2)={fval};
%     save(strcat('fminsearch_theta_6params_9patches_MI_p', num2str(iPatch),'.mat'),'patch','patch_ref','patch_pts','patch_result')
end

% save(strcat('fminsearch_theta_6params_9patches_MI_all_bnd.mat'),'patch_result')

%% NCMI 8-parameters
patch_result={};
for iPatch=1:nbPatch
    iPatch
    pts_small=patch_pts{iPatch};
    opt_img=rgb2gray(patch{iPatch});
    ri_small=patch_ref{iPatch};
    save('temp3.mat','opt_img','pts_small','avg','P0','ri_small')
    
    x0 = P(1:2,:);
    
    fun = @calc_ncmi_with_P;
    tic
    [x,fval,exitflag,output] = fminsearch(fun,x0)
    toc
    patch_result(iPatch,1)={x};
    patch_result(iPatch,2)={fval};
    %     save(strcat('fminsearchP_8params_all_patches_NCMI_p', num2str(iPatch),'.mat'),'patch','patch_ref','patch_pts','patch_result')
end
% save('fminsearchP_8params_all_patches_NCMI.mat','patch','patch_ref','patch_pts','patch_result')


%%
figure
imshow(rgb2gray(imgRGB(1:nbPatchx*szPatch(1)+1,1:nbPatchy*szPatch(2)+1,:)))
hold on
plot(xy(1,:), xy(2,:), 'gx', 'MarkerSize', 20, 'LineWidth', 5)

set(gcf, 'Position',  [100, 100, 700, 700])

alpha(.25)

hold on
for iPatch=1:nbPatch
    row_col=patch_row_col{iPatch};
    plot([row_col(2),row_col(2)+szPatch(2)],[row_col(1),row_col(1)],'r--','LineWidth',3)
    plot([row_col(2),row_col(2)],[row_col(1),row_col(1)+szPatch(1)],'r--','LineWidth',3)
    plot([row_col(2)+szPatch(2),row_col(2)+szPatch(2)],[row_col(1),row_col(1)+szPatch(1)],'r--','LineWidth',3)
    plot([row_col(2),row_col(2)+szPatch(2)],[row_col(1)+szPatch(1),row_col(1)+szPatch(1)],'r--','LineWidth',3)
end

spacex=180;
spacey=65;
fontsz=16;
for iPatch=1:nbPatch
    x=patch_result{iPatch,1};
    fval=patch_result{iPatch,2};
    
    %     if iPatch<5
    %         x=patch_result{iPatch,1};
    %         fval=patch_result{iPatch,2};
    %     end
    
    if isempty(x)
        x=x0;
        fval=0;
    end
    
    row_col=patch_row_col{iPatch};
    txt=strcat('$\textnormal{max MI}$=', sprintf('%3.3f',-fval));
    txt2=strcat('$\Delta X=', sprintf('%3.2f',x(2,1)-C(1)), '$ m');
    txt3=strcat('$\Delta Y=', sprintf('%3.2f',x(2,2)-C(2)), '$ m');
    
    t4=sprintf('%3.2f',x(2,3)-C(3));
    if strcmp(t4,'-0.00')
        t4='0.00';
    end
    txt4=strcat('$\Delta Z=', t4, '$ m');
    
    txt5=strcat('$\Delta \omega=', sprintf('%3.2f',x(1,1)-theta_x), '^\circ$');
    
    t6=sprintf('%3.2f',x(1,2)-theta_y);
    if strcmp(t6,'-0.00')
        t6='0.00';
    end
    txt6=strcat('$\Delta \phi=', t6, '^\circ$');
    
    t7=sprintf('%3.2f',x(1,3)-theta_z);
    if strcmp(t7,'-0.00')
        t7='0.00';
    end
    txt7=strcat('$\Delta \kappa=', t7, '^\circ$');
    
    text(row_col(2)+szPatch(1)/2-spacex,row_col(1)+szPatch(2)/2-3*spacey-20,txt,'Color','b','FontSize',fontsz,'Interpreter','latex')
    text(row_col(2)+szPatch(1)/2-spacex,row_col(1)+szPatch(2)/2-2*spacey-20,txt2,'Color','b','FontSize',fontsz,'Interpreter','latex')
    text(row_col(2)+szPatch(1)/2-spacex,row_col(1)+szPatch(2)/2-spacey-20,txt3,'Color','b','FontSize',fontsz,'Interpreter','latex')
    text(row_col(2)+szPatch(1)/2-spacex,row_col(1)+szPatch(2)/2-20,txt4,'Color','b','FontSize',fontsz,'Interpreter','latex')
    text(row_col(2)+szPatch(1)/2-spacex,row_col(1)+szPatch(2)/2+spacey-20,txt5,'Color','b','FontSize',fontsz,'Interpreter','latex')
    text(row_col(2)+szPatch(1)/2-spacex,row_col(1)+szPatch(2)/2+2*spacey-20,txt6,'Color','b','FontSize',fontsz,'Interpreter','latex')
    text(row_col(2)+szPatch(1)/2-spacex,row_col(1)+szPatch(2)/2+3*spacey-20,txt7,'Color','b','FontSize',fontsz,'Interpreter','latex')
end

%% avg MI score with P0
for iPatch=1:nbPatch
    iPatch
    pts_small=patch_pts{iPatch};
    opt_img=rgb2gray(patch{iPatch});
    row_col=patch_row_col{iPatch};
    
    x0=[theta_x theta_y theta_z; C(1:3)'];
    
    save('temp_6params_9patches.mat','opt_img','pts_small','K','row_col')
    res=calc_mi_with_theta_6params(x0);
    mi_score(iPatch)=-res
end

mean(mi_score) 

pts_small=pts;
opt_img=rgb2gray(imgRGB);
row_col=[1 1];

x0=[theta_x theta_y theta_z; C(1:3)'];

save('temp_6params_9patches.mat','opt_img','pts_small','K','row_col')
res=calc_mi_with_theta_6params(x0);
mi=-res


%%
figure
imshow(rgb2gray(imgRGB(1:nbPatchx*szPatch(1)+1,1:nbPatchy*szPatch(2)+1,:)))
hold on

set(gcf, 'Position',  [100, 100, 500, 500])

alpha(.25)

hold on
for iPatch=1:nbPatch
    row_col=patch_row_col{iPatch};
    plot([row_col(2),row_col(2)+szPatch(2)],[row_col(1),row_col(1)],'r--','LineWidth',3)
    plot([row_col(2),row_col(2)],[row_col(1),row_col(1)+szPatch(1)],'r--','LineWidth',3)
    plot([row_col(2)+szPatch(2),row_col(2)+szPatch(2)],[row_col(1),row_col(1)+szPatch(1)],'r--','LineWidth',3)
    plot([row_col(2),row_col(2)+szPatch(2)],[row_col(1)+szPatch(1),row_col(1)+szPatch(1)],'r--','LineWidth',3)
end

spacex=140;
spacey=30;
fontsz=15;
for iPatch=1:nbPatch
    x=patch_result{iPatch,1};
    fval=patch_result{iPatch,2};
    fval2=mi_score(iPatch);
    
    %     if iPatch<5
    %         x=patch_result{iPatch,1};
    %         fval=patch_result{iPatch,2};
    %     end
    
    if isempty(x)
        x=x0;
        fval=0;
    end
    
    row_col=patch_row_col{iPatch};
    txt=strcat('$\textnormal{MI local }$ $\theta$=', sprintf('%3.3f',-fval));
    txt2=strcat('$\textnormal{MI global }$ $\theta$=', sprintf('%3.3f',fval2));
    
    text(row_col(2)+szPatch(1)/2-spacex,row_col(1)+szPatch(2)/2-1*spacey-20,txt,'Color','b','FontSize',fontsz,'Interpreter','latex')
    text(row_col(2)+szPatch(1)/2-spacex,row_col(1)+szPatch(2)/2+1*spacey-20,txt2,'Color','b','FontSize',fontsz,'Interpreter','latex')
end

% for iPatch=1:nbPatch
%     row_col=patch_row_col{iPatch};
% end
