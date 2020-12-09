function res=calc_ncmi_with_P(x0)
load('temp3.mat','opt_img','pts','avg','P0','ri')

P=zeros(3,4);
P(1:2,:)=x0;
P(3,4)=1;

rx2=linspace(ri.XWorldLimits(1),ri.XWorldLimits(2),ri.RasterSize(2));
ry2=linspace(ri.YWorldLimits(1),ri.YWorldLimits(2),ri.RasterSize(1));

[int_img,z_img]=SR_for_NCMI(P,pts,avg,ri.RasterSize(1),ri.RasterSize(2),rx2,ry2,500);
[NCMI,~]=ncmi(z_img,int_img,opt_img);
res=-NCMI;
end