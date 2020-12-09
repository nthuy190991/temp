function [NCMI, CMI]=ncmi(A,B,C,varargin) 
% created on 07 May 2019

if nargin>=4
    L=varargin{1};
else
    L=256;
end

A=double(A); 
B=double(B); 
C=double(C); 

% disp(L)

% na = hist(A(:),L);
% na = na/sum(na);
% 
% nb = hist(B(:),L); 
% nb = nb/sum(nb);

nc = hist(C(:),L); 
nc = nc/sum(nc);

n2ab = hist2(A(:),B(:),L); 
n2ab = n2ab/sum(n2ab(:));

% n2ac = hist2(A,C,L); 
% n2ac = n2ac/sum(n2ac(:));
% 
% n2bc = hist2(B,C,L); 
% n2bc = n2bc/sum(n2bc(:));

% Iac=sum(minf(n2ac,na'*nc)); 
% Ibc=sum(minf(n2bc,nb'*nc)); 

[n3,~,~,~] = histcn([A(:),B(:),C(:)],L-1,L-1,L-1);
n3 = n3/sum(n3(:));

% disp(entro(na)+entro(nb)-entro(n2ab))
% disp(sum(minf(n2ab,na'*nb)))
% 
% disp(entro(na)+entro(nb)-entro(n2ab)-sum(minf(n2ab,na'*nb)))

CMI=entro(n2ab)+entro(nc)-entro(n3);

NCMI=(entro(n2ab)+entro(nc))/entro(n3);

% -----------------------

% Mutual info
function y=minf(pab,papb)

I=find(papb(:)>1e-12 & pab(:)>1e-12); % function support 
y=pab(I).*log2(pab(I)./papb(I));
% y=pab(:).*log2(pab(:)./papb(:));

% Entropy
function h=entro(px)

I=find(px(:)>1e-12); % function support 
h=-sum(px(I).*log2(px(I)));
% h=-sum(px(:).*log2(px(:)));


