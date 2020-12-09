% Run building *3D* boundary extraction on a large scale
function [bounds,points] = buildingBoundaryExtraction(filename,saveFlag)
    if strcmp(filename,'')
        filename='building_mask.mat';
    end

    load(filename,'building_label','mask_ref','nb_labels','pts_non_gnd_denoised')
    figure
    mapshow(building_label,mask_ref)

    tol=1; 
    tolz=1.5;
    bounds={};
    points={};
    for i=1:nb_labels
        i
        [row,col]=find(building_label==i);
        r = [mask_ref.CellExtentInWorldX*(col-1) + mask_ref.XWorldLimits(1),...
             mask_ref.CellExtentInWorldY*(row-1) + mask_ref.YWorldLimits(1)];
        % rx = [min(r(:,1)) max(r(:,1))];
        % ry = [min(r(:,2)) max(r(:,2))];

        % get boundary from the mask
        b = boundary(r(:,1),r(:,2),1);

        % find points inside mask-based boundary
        [pts]=findPointsInBoundaries(pts_non_gnd_denoised(:,1:3),r(b,:),tol,tolz);

        % get boundary from 3D points
        b2 = boundary(pts(:,1),pts(:,2),1);
        bounds(i,:) = {pts(b2,:)};
        points(i,:) = {pts};
    end

    clear b b2
    figure
    pcshow(pts_non_gnd_denoised(:,1:3))
    hold on
    for i=1:nb_labels
        b=bounds{i,:};
        plot3(b(:,1),b(:,2),b(:,3),'r','LineWidth',3)
    %     p=points{i,:};
    %     pcshow(p,'g')
    end
    grid on
    view(2)

    if saveFlag
        s=strsplit(filename,'_');

        filename2=strcat('building_boundary_',s(end))
        save(filename2{1},'bounds','points')
    end
end

%% Find points within a given preliminary set of boundaries
function pts=findPointsInBoundaries(points,bound,tol,tolz)
    pts=findPtsInBound(points,bound,tol,tolz);
    k=0;
    while size(pts,1) ~= size(findPtsInBound(points,pts,tol,tolz),1)
        k=k+1;
        bound=pts;
        pts=findPtsInBound(points,bound,tol,tolz);
    end
end

%% Sub-function finding points inside boundaries
function pts=findPtsInBound(points,bound,tol,tolz)
bound=unique(bound,'row');
pts=[];
for i=1:size(bound,1)

    p=bound(i,:);
    if size(p,2)<3
        pts = [pts; findPoints(points,p(1)+[-tol tol],p(2)+[-tol tol],0)];
    else
        pts = [pts; find3DPoints(points,p(1)+[-tol tol],p(2)+[-tol tol],...
                                        p(3)+[-tolz tolz],0)];
    end
end
end

%%
function [pts, inliers]=find3DPoints(points,rx,ry,rz,epsilon)

x=points(:,1);
y=points(:,2);
z=points(:,3);

inRangeX = x >= rx(1)-epsilon & x <= rx(2)+epsilon;
inRangeY = y >= ry(1)-epsilon & y <= ry(2)+epsilon;
inRangeZ = z >= rz(1)-epsilon & z <= rz(2)+epsilon;

inliers=[];
    for i=1:size(points,1)
        if (inRangeX(i)==1)&&(inRangeY(i)==1)&&(inRangeZ(i)==1)
            %i
            inliers=[inliers i];
        end
    end
    pts=points(inliers,:);
end