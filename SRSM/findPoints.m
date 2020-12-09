function [pts, inliers]=findPoints(points,rx,ry,epsilon)

x=points(:,1);
y=points(:,2);

inRangeX = x >= rx(1)-epsilon & x <= rx(2)+epsilon;
inRangeY = y >= ry(1)-epsilon & y <= ry(2)+epsilon;

% display([num2str(x) ',' num2str(y) ',' num2str(rx) ',' num2str(ry)])

inliers=[];
for i=1:size(points,1)
    if (inRangeX(i)==1)&&(inRangeY(i)==1)
        %i
        inliers=[inliers i];
    end
end
pts=points(inliers,:);
end