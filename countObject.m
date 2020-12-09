% A count function to see if the object is well detected (det) based on its
% reference boundary (ref). 
% Tolerance of area percentage is 50%, i.e. if an object is detected 50% of
% its area, it is counted as a true detection, else false detection.
function [flag,com]=countObject(det, ref, tol, verbose)
    
    xmin=min([det(:,1); ref(:,1)]);
    ymin=min([det(:,2); ref(:,2)]);
    sizex=max([det(:,1); ref(:,1)])-xmin;
    sizey=max([det(:,2); ref(:,2)])-ymin;
    I=zeros(round(sizex),round(sizey));
    
    % Ground truth building mask
    pt=zeros(numel(I),2);
    for i=1:size(I,1)
        for j=1:size(I,2)
            pt((i-1)*size(I,2)+j,:)=[i,j];
        end
    end

    det_mask = reshape(inpolygon(pt(:,1),pt(:,2),det(:,1)-xmin,det(:,2)-ymin),size(I,2),size(I,1));
    ref_mask = reshape(inpolygon(pt(:,1),pt(:,2),ref(:,1)-xmin,ref(:,2)-ymin),size(I,2),size(I,1));

    if verbose
        figure
        subplot(121), imshow(det_mask)
        subplot(122), imshow(ref_mask)
    end
    
    com=sum(sum(det_mask&ref_mask))/sum(sum(ref_mask));

    if com>tol
        flag=1;
    else
        flag=0;
    end
    
end