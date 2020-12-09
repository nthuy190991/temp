function P=SnakeMoveIteration2D(B,P,Fext,I,gamma,kappa,delta,rve_flag,epsilon,df_flag,df_seeds,nu)
% This function will calculate one iteration of contour Snake movement
%
% P=SnakeMoveIteration2D(S,P,Fext,gamma,kappa)
%
% inputs,
%   B : Internal force (smoothness) matrix
%   P : The contour points N x 2;
%   Fext : External vector field (from image)

% Huy's modif:
%   I : the grayscale image
%   epsilon : RVE coefficient


%   gamma : Time step
%   kappa : External (image) field weight
%   delta : Balloon Force weight
%
% outputs,
%   P : The (moved) contour points N x 2;
%
% Function is written by D.Kroon University of Twente (July 2010)

% % Clamp contour to boundary
% P(:,1)=min(max(P(:,1),1),size(Fext,1));
% P(:,2)=min(max(P(:,2),1),size(Fext,2));

% Get image force on the contour points
Fext1(:,1)=kappa*interp2(Fext(:,:,1),P(:,2),P(:,1));
Fext1(:,2)=kappa*interp2(Fext(:,:,2),P(:,2),P(:,1));
% Interp2, can give nan's if contour close to border
Fext1(isnan(Fext1))=0;

% Calculate the baloonforce on the contour points
N=GetContourNormals2D(P,4);
Fext2=delta*N;


% New external term, based on Regional Variance Energy (proposed by Kabolizade, 2010)
Fext3=0;
if (rve_flag)
%     display(size(I))
%     display([round(P)])
        
    idx1=find(round(P(:,1))<=0);
    idx2=find(round(P(:,2))<=0);
    idx3=find(round(P(:,1))>size(I,1));
    idx4=find(round(P(:,2))>size(I,2));

    Q=P;
    Q([idx1;idx2;idx3;idx4],:)=[];
%     round(Q)
%     idx = sub2ind(size(I),floor(P(:,1)),floor(P(:,2)));

    idx = sub2ind(size(I),round(Q(:,1)),round(Q(:,2)));
    mi=sum(I(idx))/(length(Q)-1);
    Fext3=epsilon*sum((I(idx)-mi).^2);
end
% display(['RVE=' num2str(Fext3)])

% New external term, based on shape similarity
Fext4=0;

if (df_flag)    
    [Fext4,hd]=calc_shape_sim(P,df_seeds,nu);
%     display(['SS=' num2str(Fext4), ', hd=' num2str(hd)])
end

% New external term, based on corner of seed point set
Fext5=0;
% if (df_flag)    
%     [hd,~]=HausdorffDist(P,df_seeds);
% %     display(hd)
%     Fext5=nu(1)*(1-exp(-hd.^2/nu(2)));
% end

% Update contour positions
ssx = gamma*P(:,1) + Fext1(:,1) + Fext2(:,1) + Fext3 + Fext4 + Fext5;
ssy = gamma*P(:,2) + Fext1(:,2) + Fext2(:,2) + Fext3 + Fext4 + Fext5;
P(:,1) = B * ssx;
P(:,2) = B * ssy;

% Fext4=calc_shape_sim(P,df_seeds,nu);
%     display(['SS=' num2str(Fext4)])
%     display(['Fext1_x=' num2str(mean(Fext1(:,1))) ', Fext1_y=' num2str(mean(Fext1(:,2)))])
%     display(['Fext2_x=' num2str(mean(Fext2(:,1))) ', Fext2_y=' num2str(mean(Fext2(:,2)))])

% % Clamp contour to boundary
% P(:,1)=min(max(P(:,1),1),size(Fext,1));
% P(:,2)=min(max(P(:,2),1),size(Fext,2));

end


function [ss,hd]= calc_shape_sim(P,seeds,nu)
    [hd,~]=HausdorffDist(P,seeds);
%     display(hd)
    ss=nu(1)*(1-exp(-hd.^2/nu(2)));
end