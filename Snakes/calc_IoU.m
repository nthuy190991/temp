
%% Intersection Over Union

% Src: https://zhenye-na.github.io/2018/05/22/intersection-over-union-for-object-detection.html
% def bb_intersection_over_union(boxA, boxB):
%     # determine the (x, y)-coordinates of the intersection rectangle
%     xA = max(boxA[0], boxB[0])
%     yA = max(boxA[1], boxB[1])
%     xB = min(boxA[2], boxB[2])
%     yB = min(boxA[3], boxB[3])
% 
%     # compute the area of intersection rectangle
%     interArea = (xB - xA) * (yB - yA)
% 
%     # compute the area of both the prediction and ground-truth
%     # rectangles
%     boxAArea = (boxA[2] - boxA[0]) * (boxA[3] - boxA[1])
%     boxBArea = (boxB[2] - boxB[0]) * (boxB[3] - boxB[1])
% 
%     # compute the intersection over union by taking the intersection
%     # area and dividing it by the sum of prediction + ground-truth
%     # areas - the interesection area
%     iou = interArea / float(boxAArea + boxBArea - interArea)
% 
%     # return the intersection over union value
%     return iou

% function iou=calc_IoU(boxA, boxB, shape)
%     if strcmp(shape,'L')
%         iou=0;
%     elseif strcmp(shape,'R') % rectangular
%         
%         % determine the (x, y)-coordinates of the intersection rectangle
%         xA=max(boxA(1), boxB(1));
%         yA=max(boxA(2), boxB(2));
%         xB=min(boxA(3), boxB(3));
%         yB=min(boxA(4), boxB(4));
%         
%         % compute the area of intersection rectangle
%         interArea = (xB - xA) * (yB - yA);
%         
%         % compute the area of both the prediction and ground-truth rectangles
%         boxAArea = (boxA(3) - boxA(1)) * (boxA(4) - boxA(2));
%         boxBArea = (boxB(3) - boxB(1)) * (boxB(4) - boxB(2));
%         
%         % compute the intersection over union by taking the intersection
%         % area and dividing it by the sum of prediction + ground-truth 
%         % areas - the interesection area
%         iou = interArea / (boxAArea + boxBArea - interArea);
%     end
% end

% boxA: estimated box
% boxB: ground truth box
function [iou,com,corr]=calc_IoU(I, boxA, boxB, verbose)
    
    pt=zeros(numel(I),2);
    for i=1:size(I,1)
        for j=1:size(I,2)
            pt((i-1)*size(I,2)+j,:)=[i,j];
        end
    end
    
    inBoxA = inpolygon(pt(:,1),pt(:,2),boxA(1,:),boxA(2,:));
    inBoxB = inpolygon(pt(:,1),pt(:,2),boxB(1,:),boxB(2,:));
    
    iou=sum(inBoxA&inBoxB)/sum(inBoxA|inBoxB);
    
    com=sum(inBoxA&inBoxB)/sum(inBoxB);
    corr=sum(inBoxA&inBoxB)/sum(inBoxA);
    
    if verbose
        ptInter=pt(inBoxA&inBoxB,:);
        ptUnion=pt(inBoxA|inBoxB,:);
        figure 
        subplot(121)
        imshow(I)
        hold on
        plot([boxA(2,:),boxA(2,1)], [boxA(1,:),boxA(1,1)], 'b')
        plot([boxB(2,:),boxB(2,1)], [boxB(1,:),boxB(1,1)], 'r')
        plot(ptInter(:,2),ptInter(:,1),'y.')
        
        subplot(122)
        imshow(I)
        hold on
        plot([boxA(2,:),boxA(2,1)], [boxA(1,:),boxA(1,1)], 'b')
        plot([boxB(2,:),boxB(2,1)], [boxB(1,:),boxB(1,1)], 'r')
        plot(ptUnion(:,2),ptUnion(:,1),'y.')
    end
    
end
