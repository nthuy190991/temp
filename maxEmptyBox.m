function [bb]=maxEmptyBox(X,verbose)

xmin=min(X(1,:));
xmax=max(X(1,:));
ymin=min(X(2,:));
ymax=max(X(2,:));

xsize=xmax-xmin;
ysize=ymax-ymin;

offset=[20, 20];

%% create a binary image representing the points X
J=zeros(ysize+offset(2)*2,xsize+offset(1)*2);
for i=1:size(X,2)
    x=X(1,i)-xmin+offset(1);
    y=X(2,i)-ymin+offset(2);
    J(y,x)=1;
end

%% hole filling for binary image J
% for y=1:ysize
%     idx=find(J(y,:)==1);
%     J(y,min(idx):max(idx))=1;
% end
% 
% for x=1:xsize
%     idx=find(J(:,x)==1);
%     J(min(idx):max(idx),x)=1;
% end
% 
% idx=find(J==1);
% px=zeros(2,length(idx));
% k=1;
% for y=1:ysize
%     for x=1:xsize
%         if J(y,x)==1
%             px(1,k)=x+xmin-1;
%             px(2,k)=y+ymin-1;
%             k=k+1;
%         end
%     end
% end

SE=strel('disk',7);
J=imclose(J,SE);
% [r,c]=find(J==1);
% px=[c+xmin-offset(1),r+ymin-offset(2)];

%% for-loop to look for the optimal rectangle area
ang=1:0.5:90;
area=zeros(1,length(ang));
for i=1:length(ang)
    
    rotated = imrotate(J,ang(i),'nearest','crop');
    [~,~,~, MM] = FindLargestRectangles(rotated);%, [0 0 1]);
    area(i)=sum(MM(:));
end
[~,iMax]=max(area);
im=imrotate(J,ang(iMax),'nearest','crop');

[C, H, W, ~] = FindLargestRectangles(im);
[~, pos] = max(C(:));
[r, c] = ind2sub(size(im), pos);

if (verbose)
    figure
    imshow(J)
    figure
    imshow(im)
    hold on
    rectangle('Position',[c,r,W(r,c),H(r,c)], 'EdgeColor','r', 'LineWidth',3);
end

%% determine rectangle corners of the original image
theta=ang(iMax);
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
corners = [c r; c+W(r,c) r; c+W(r,c) r+H(r,c); c r+H(r,c)]';
% center=mean(corners,2);
center=[fliplr(size(im))/2]';
rotpoints = R*(corners-center) + center;

if (verbose)
    figure
    imagesc(J)
    hold on
    plot(rotpoints(1,:), rotpoints(2,:), 'r-','LineWidth',1)
end

%% offset with X coordinates
bb=zeros(size(rotpoints));
bb(1,:)=rotpoints(1,:)+xmin-offset(1);
bb(2,:)=rotpoints(2,:)+ymin-offset(2);

end