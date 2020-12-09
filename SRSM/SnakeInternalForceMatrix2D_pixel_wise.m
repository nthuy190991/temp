function [A,B]=SnakeInternalForceMatrix2D_pixel_wise(P,nPoints,alpha,tab_beta,gamma)
%
% B=SnakeInternalForceMatrix2D(nPoints,alpha,beta,gamma)
%
% inputs,
%   nPoints : The number of snake contour points
%   alpha : membrame energy  (first order)
%   beta : thin plate energy (second order)
%   gamma : Step Size (Time)
%
% outputs,
%   B : The Snake Smoothness regulation matrix
%
% Function is written by D.Kroon University of Twente (July 2010)

beta=zeros(1,nPoints+2);
for i=1:nPoints
    beta(i+1)=interp2(tab_beta,P(i,2),P(i,1));
end
beta(1)=mean([beta(2),beta(nPoints+1)]);
beta(nPoints+2)=beta(1);

if any(isnan(beta))
    beta=zeros(1,nPoints+2);
    for i=1:nPoints
        beta(i+1)=interp2(tab_beta,P(i,2),P(i,1),'spline');
    end
    beta(1)=mean([beta(2),beta(nPoints+1)]);
    beta(nPoints+2)=beta(1);
end

% Penta diagonal matrix, one row:
s=2:nPoints+1;
b1=beta(s-1);
b2=-(alpha + 2*(beta(s-1)+beta(s)));
b3=2*alpha + beta(s-1)+4*beta(s)+beta(s+1);
b4=-(alpha + 2*(beta(s+1)+beta(s)));
b5=beta(s+1);
    
% Make the penta matrix (for every contour point)
A=diag(b1)*circshift(eye(nPoints),2);
A=A+diag(b2)*circshift(eye(nPoints),1);
A=A+diag(b3)*circshift(eye(nPoints),0);
A=A+diag(b4)*circshift(eye(nPoints),-1);
A=A+diag(b5)*circshift(eye(nPoints),-2);

% Calculate the inverse
B=inv(A + gamma.* eye(nPoints));



