function [P, K, R, t, C, error] = runGoldStandard(xy, XYZ)

[~,N] = size(xy);
%normalize data points
[xy_normalized,XYZ_normalized,T,U] = normalization(xy, XYZ);

%compute DLT
[Pn] = dlt(xy_normalized, XYZ_normalized);

%minimize geometric error
pn = [Pn(1,:) Pn(2,:) Pn(3,:)];
for i=1:20
    [pn] = fminsearch(@fminGoldStandard, pn, [], xy_normalized, ...
        XYZ_normalized, i/5);
end
% [pn] = fminsearch(@fminGoldStandard, pn, [], xy_normalized, ...
%         XYZ_normalized, 0);
    
Pn = [pn(1:4);pn(5:8);pn(9:12)];

%denormalize camera matrix
P=inv(T)*Pn*U;

%factorize camera matrix in to K, R and t
M=inv(P(:,1:3));
p4=P(:,4);
[QQ,RR]=qr(M);
R=inv(QQ);
K=inv(RR);

%specific for my implementation
K=-K;
R=-R;
K(1,1)=-K(1,1);
R(1,:)=-R(1,:);

t=inv(K)*p4;

%compute reprojection error
xy = [xy;ones(1,N)];
XYZ = [XYZ;ones(1,N)];

xy_est = K*[R t]*XYZ;

for i = 1:N
    xy_est(:,i) = xy_est(:,i)/xy_est(end,i);
end

%computation of the center
[~,~,VP] = svd(P);
C=VP(:,end)/VP(end,end);

error=norm(xy-xy_est,'fro');
end

%%
function f = fminGoldStandard(p, xy, XYZ, w)
[~,N] = size(xy);

%reassemble P
P = [p(1:4);p(5:8);p(9:12)];


%factorize camera matrix in to K, R and t
%this is necessary to get the additional
%penality
M=inv(P(:,1:3));
p4=P(:,4);
[QQ,RR]=qr(M);
R=inv(QQ);
K=inv(RR);

%specific for my implementation with 24 points
K(1,1)=-K(1,1);
R(1,:)=-R(1,:);

t=inv(K)*p4;

%get value for s
s=K(1,2);
alphax=K(1,1);
alphay=K(2,2);

%get the error
xy_est = K*[R t]*XYZ;
for i = 1:N
    xy_est(:,i) = xy_est(:,i)/xy_est(end,i);
end

%compute cost
f=norm(xy-xy_est,'fro')+w*s^2+w*(alphax-alphay)^2;
end


function [P] = dlt(xy, XYZ)
%computes DLT, xy and XYZ should be normalized before calling this function
[~,N] = size(xy);

%Staking of matrices together
DLT=[];
for i = 1:N
    DLT=[DLT;SubDLT(xy(:,i), XYZ(:,i))];
end
[UU, SS, VV] = svd(DLT(1:end,:));
P=reshape(VV(:,end),4,3)';

end

%This function computes DLT matrix for a single couple of points
function [Pi] = SubDLT(xy, XYZ)
Pi=[xy(3)*XYZ', zeros(size(XYZ')), -xy(1)*XYZ';
    zeros(size(XYZ')), -xy(3)*XYZ', xy(2)*XYZ'];
end
