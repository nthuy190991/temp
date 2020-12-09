function [delta_f,gx,gy]=gradient_Dxy(x,Ny)

%     x=reshape(x,1,[]);
%     Dx=zeros(length(x),2);
%     for n=1:length(x)
%         Dx(n,1)=x(n+Ny)-x(n);
%         Dx(n,2)=x(n+1)-x(n);
%     end

gx=grad_Dx(x,Ny);
gy=grad_Dy(x);
delta_f=gx+gy;

end

function res=grad_Dx(x,Ny)
    res=zeros(length(x),1);
    for i=1:Ny
        res(i)=2*(x(i)-x(i+Ny));
    end
%     res(length(x)-Ny:end)=2*x(length(x)-Ny:end);
    for i=length(x)-Ny+1:length(x)
        res(i)=2*(x(i)-x(i-Ny));
    end
    for i=Ny+1:length(x)-Ny
        res(i)=2*(x(i)-x(i-Ny)) + 2*(x(i)-x(i+Ny));
    end
end

function res=grad_Dy(x)
    res=zeros(length(x),1);
    res(1)=2*(x(2)-x(1));
    res(end)=2*x(end);
    for i=2:length(x)-1
        res(i)=2*(x(i)-x(i-1)) + 2*(x(i)-x(i+1));
    end
end