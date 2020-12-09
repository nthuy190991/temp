function [delta_f,gx,gy]=gradient_Dxy_fast(x,Ny)

%     x=reshape(x,1,[]);
%     Dx=zeros(length(x),2);
%     for n=1:length(x)
%         Dx(n,1)=x(n+Ny)-x(n);
%         Dx(n,2)=x(n+1)-x(n);
%     end

% gx=grad_Dx(x,Ny);
% gy=grad_Dy(x);
% [gx, gy]=grad_Dx_Dy(x,Ny);

gx=zeros(length(x),1);
gx(1:Ny)=0.5*(x(1:Ny)-x((1:Ny)+Ny));
gx(end-Ny+1:end)=(x(end-Ny+1:end)-x((end-Ny+1:end)-Ny));
gx(Ny+1:end-Ny)=2*x(Ny+1:end-Ny)-x((Ny+1:end-Ny)-Ny) -x((Ny+1:end-Ny)+Ny);

gy=zeros(length(x),1);
gy(1)=(x(2)-x(1));
gy(end)=x(end);
gy(2:end-1)=2*x(2:end-1)-x((2:end-1)-1)-x((2:end-1)+1);

delta_f=gx+gy;
end

