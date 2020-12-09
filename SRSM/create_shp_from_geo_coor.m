function [shape,dbfspec]=create_shp_from_geo_coor(filename,xlims,ylims,resolution,write_flag,verbose)

load(strcat(filename,'.mat'),'snake_mask','snake_results','patch_corner')

nbBuildings=size(snake_results,1);

tb_empty=[];
for i=1:nbBuildings
    s=snake_results{i};
    if isempty(s)
        tb_empty=[tb_empty, i];
    end
end

snake_results2=snake_results;
snake_results2(tb_empty)=[];

patch_corner2=patch_corner;
patch_corner2(tb_empty,:)=[];

A={};
for i=1:nbBuildings-length(tb_empty)
    A(i).FID_=0;
    A(i).Entity='3DPolyline';
    A(i).Layer='3D_Gebaeude_Dissolve_Haeuser';
    A(i).Color=7;
    A(i).Linetype='Continuous';
    A(i).Elevation=0;
    A(i).LineWt=25;
    A(i).RefName='';
end


S={};
for i=1:nbBuildings-length(tb_empty)
    s=snake_results2{i};
    %     row=snake_results{i}(:,1)+patch_corner(i,2);
    %     col=snake_results{i}(:,2)+patch_corner(i,1);
    %
    %     [x,y]=convertMatCoorGeoCoor(rx,ry,col,row,ri.RasterSize);


    x=(s(:,2)+patch_corner2(i,1))*resolution(1)+xlims(1);
    y=ylims(2)-(s(:,1)+patch_corner2(i,2))*resolution(2);

    S(i).Geometry='Polygon';
    S(i).BoundingBox=[min(x) min(y); max(x) max(y)];
    S(i).X=x';
    S(i).Y=y';
end

% figure
% hold on
% for i=1:nbBuildings
%     x=snake_results{i}(:,1)+patch_corner(i,2);
%     y=snake_results{i}(:,2)+patch_corner(i,1);
%     plot(x,y)
% end


% shape = mapshape({S.X}, {S.Y}, A);
shape = mapshape({S.X}, {S.Y}, A);
shape.Geometry = S(1).Geometry;
dbfspec = makedbfspec(shape);


if write_flag
    shapewrite(shape,strcat(filename,'.shp'),'DbfSpec',dbfspec)
end

if verbose
    shapeinfo(strcat(filename,'.shp'))
    
    [S3,A3]=shaperead(strcat(filename,'.shp'));

    figure
    hold on
    for i=1:length(S3)
        X=S3(i).X;
        Y=S3(i).Y;
        plot(X,Y)

        b=S3(i).BoundingBox;
        plot([b(1,1) b(1,1) b(2,1) b(2,1) b(1,1)],...
             [b(1,2) b(2,2) b(2,2) b(1,2) b(1,2)])
    end

end
end