function [bound, mask, seg_area]=segmentFiltering(all_labels,label,Sa_thres)

%removal of segments based on their size

bw=zeros(size(all_labels));
bw(all_labels==label)=1;

[L,n] = bwlabel(bw);

% dont keep components that are too small
l=1;
bound={};
mask={};
seg_area=[];
for i=1:n
    idx=find(L==i);
    Sa=length(idx);
    if Sa>Sa_thres(1) && Sa<Sa_thres(2)
        [r,c]=ind2sub(size(L),idx);
        b=boundary(r,c,1);
        bound(l,:)={[r(b), c(b)]};
        mask(l,:)={[r, c]};
        seg_area=[seg_area; length(r)];
        l=l+1;
    end
end
end