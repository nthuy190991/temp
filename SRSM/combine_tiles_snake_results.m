filename='res_BE_LiDAR2017_4784_Feb17.mat'
load(filename)

str=strsplit(filename,'_LiDAR2017_');
tile_no=str{2}(1:4);

strx=tile_no(1:2);
stry=tile_no(3:4);

xlims=[2e5+str2double(strx)*1e3 2e5+str2double(strx)*1e3+1e3];
ylims=[51e5+str2double(stry)*1e3 51e5+str2double(stry)*1e3+1e3];
reso_grid=[0.15 0.15];

tol=0.5;
near_xlims=xlims+[-tol tol]';
near_ylims=ylims+[-tol tol]';

figure, hold on, grid on
for i=1:length(snake_results)
    s=snake_results{i};
    x=(s(:,2)+patch_corner(i,1))*reso_grid(1)+xlims(1);
    y=ylims(2)-(s(:,1)+patch_corner(i,2))*reso_grid(2);
    plot(x,y,'b','LineWidth',1)
    
    idx1=find(x>near_xlims(1,1)&x<near_xlims(2,1));
    idx2=find(x>near_xlims(1,2)&x<near_xlims(2,2));
    idx3=find(y>near_ylims(1,1)&y<near_ylims(2,1));
    idx4=find(y>near_ylims(1,2)&y<near_ylims(2,2));
    idx=[idx1; idx2; idx3; idx4];
    plot(x(idx),y(idx),'r','LineWidth',2)
end

%%

list=dir('res_BE_LiDAR2017_4*83_Feb17.mat');
% list(1)=[]; % not account the tile 4381

minx=Inf;
maxx=0;
miny=Inf;
maxy=0;
for idx=1:length(list)
    load(list(idx).name,'snake_mask')
    
    str=strsplit(list(idx).name,'_LiDAR2017_');
    tile_no=str{2}(1:4);

    strx=tile_no(1:2);
    stry=tile_no(3:4);

    xlims=[2e5+str2double(strx)*1e3 2e5+str2double(strx)*1e3+1e3];
    ylims=[51e5+str2double(stry)*1e3 51e5+str2double(stry)*1e3+1e3];
    
    minx=min(minx,xlims(1));
    maxx=max(maxx,xlims(2));
    miny=min(miny,ylims(1));
    maxy=max(maxy,ylims(2));
end

n_cols=round((maxx-minx)/reso_grid(1));
n_rows=round((maxy-miny)/reso_grid(2));
A=zeros(n_rows,n_cols);

for idx=1:length(list)
    load(list(idx).name,'snake_mask')
    
    str=strsplit(list(idx).name,'_LiDAR2017_');
    tile_no=str{2}(1:4);

    strx=tile_no(1:2);
    stry=tile_no(3:4);

    xlims=[2e5+str2double(strx)*1e3 2e5+str2double(strx)*1e3+1e3];
    ylims=[51e5+str2double(stry)*1e3 51e5+str2double(stry)*1e3+1e3];
    
    r0=max(1,round((ylims(1)-miny)/reso_grid(2)));
    c0=max(1,round((xlims(1)-minx)/reso_grid(1)));

    A(r0:r0+size(snake_mask,1)-1,c0:c0+size(snake_mask,2)-1)=(snake_mask);
end


figure
imshow(A)