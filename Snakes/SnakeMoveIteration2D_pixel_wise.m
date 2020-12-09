function P=SnakeMoveIteration2D_pixel_wise(I,A,B,P,Fext,gamma,kappa,delta_mat_flag,delta,delta_mat,rve_flag,epsilon,ss_flag,ss_refs,nu)
% This function will calculate one iteration of contour Snake movement
%
% P=SnakeMoveIteration2D(S,P,Fext,gamma,kappa)
%
% inputs,
%   B : Internal force (smoothness) matrix
%   P : The contour points N x 2;
%   Fext : External vector field (from image)

%   gamma : Time step
%   kappa : External (image) field weight
%   delta : Balloon Force weight
%
% outputs,
%   P : The (moved) contour points N x 2;
%
% Function is written by D.Kroon University of Twente (July 2010)

%% Clamp contour to boundary
P(:,1)=min(max(P(:,1),1),size(Fext,1));
P(:,2)=min(max(P(:,2),1),size(Fext,2));

% Get image force on the contour points
Fext1(:,1)=kappa*interp2(Fext(:,:,1),P(:,2),P(:,1));
Fext1(:,2)=kappa*interp2(Fext(:,:,2),P(:,2),P(:,1));

% display('Using user-input E_img')
% dir='/Users/thanhhuynguyen/Storage/Dev/DSAC-master/';
% textfile=strcat(dir,'mapE_building49_80epochs.txt');
% fileID = fopen(textfile,'r');
% A = fscanf(fileID,'%f',[256 256]);
% Fext = A';
% fclose(fileID);
% Fext1=kappa*interp2(Fext,P(:,2),P(:,1));


% Interp2, can give nan's if contour close to border
Fext1(isnan(Fext1))=0;

% Calculate the baloon force on the contour points
if (~delta_mat_flag)
%     disp('Scalar delta')
    N=GetContourNormals2D(P,4);
    Fext2=delta*N;
% display([num2str(mean(Fext2(:,1))), ', ', num2str(mean(Fext2(:,2)))])
else
%     disp('Pixel-wise delta')
    
%     Fext2=delta/2*interp2(delta_mat,P(:,2),P(:,1));
    N=GetContourNormals2D(P,4);
    Fext2=delta*(N.*interp2(delta_mat,P(:,2),P(:,1)));
end
%% New external term, based on Regional Variance Energy (proposed by Kabolizade, 2010)
rve=0;
if (rve_flag)
    val_snake_pts = interp2(I,P(:,2),P(:,1));
    mi=mean(val_snake_pts);
    rve=epsilon*mean((val_snake_pts-mi).^2);
end

%% New external term, based on shape similarity
shape_sim=0;
if (ss_flag)    
    [shape_sim,~]=calc_shape_sim(P,ss_refs,nu);
%     display(['SS=' num2str(Fext4), ', hd=' num2str(hd)])
end

%% Update contour positions
% if (~delta_mat_flag)
    ssx = gamma*P(:,1) + Fext1(:,1) + Fext2(:,1) + rve + shape_sim;
    ssy = gamma*P(:,2) + Fext1(:,2) + Fext2(:,2) + rve + shape_sim;
    P(:,1) = B * ssx;
    P(:,2) = B * ssy;
% else
%     ssx = gamma*P(:,1) + Fext1(:,1) + Fext2 + rve + shape_sim;
%     ssy = gamma*P(:,2) + Fext1(:,2) + Fext2 + rve + shape_sim;
%     P(:,1) = B * ssx;
%     P(:,2) = B * ssy;
% end

%% New 2019/11/08: Update contour position
% ssx = A*P(:,1) - (Fext1(:,1) + Fext2(:,1) + rve + shape_sim);
% ssy = A*P(:,2) - (Fext1(:,2) + Fext2(:,2) + rve + shape_sim);
% P(:,1) = P(:,1) - gamma * ssx;
% P(:,2) = P(:,2) - gamma * ssy;


% figure(100)
% plot(Fext2(:,1),Fext2(:,2))

% % Clamp contour to boundary
% P(:,1)=min(max(P(:,1),1),size(Fext,1));
% P(:,2)=min(max(P(:,2),1),size(Fext,2));

end


%% Shape similarity energy term, calculated by Hausdorff distance
function [ss,hd]= calc_shape_sim(P,refs,nu)
    [hd,~]=HausdorffDist(P,refs);
%     display(hd)
    ss=nu(1)*(1-exp(-hd.^2/nu(2)));
end