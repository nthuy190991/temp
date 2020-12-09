function calc_f=cost_func_Dxy(x,Ny)

Dx=zeros(length(x),2);
for n=1:length(x)-Ny
    Dx(n,1)=x(n+Ny)-x(n);
    Dx(n,2)=x(n+1)-x(n);
end
for n=length(x)-Ny+1:length(x)-1
    Dx(n,1)=x(n);
    Dx(n,2)=x(n+1)-x(n);
end
Dx(length(x),2)=x(n);
calc_f=sum(Dx(:,1).^2+Dx(:,2).^2);
end