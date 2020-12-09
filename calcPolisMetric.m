%% Calculate PoLiS metric
% A Metric for Polygon Comparison and Building Extraction Evaluation
% Janja Avbelj, Rupert Müller, and Richard Bamler
% IEEE GRSL, vol. 12, no. 1, 2015
% https://ieeexplore.ieee.org/abstract/document/6849454
%
% Example test code: 
% polyA=polyshape([0 0 3 3],[1 0 0 3]);
% polyB=polyshape([0 1 4 2],[2 3 1 1]);
% figure
% plot(polyA)
% hold on, grid on
% plot(polyB)
% d=calcPolisMetric(polyA.Vertices,polyB.Vertices)
%  
% Written by Thanh Huy Nguyen, 2020
% based on the pseudo-code provided in the paper, and following
% https://github.com/most-ieeta/preprocessing_geometry/blob/master/src/polygon.cpp
function d=calcPolisMetric(polyA,polyB)
    d12=0; % directed distance between 1 and 2
    d21=0; % directed distance between 2 and 1
    q=length(polyA);
    r=length(polyB);
    
    for j=1:q
        a_j=polyA(j,:);
        d12=d12+minDistPt2Poly(a_j,polyB);
    end
    for k=1:r
        b_k=polyB(k,:);
        d21=d21+minDistPt2Poly(b_k,polyA);
    end
    d=(d12/q+d21/r)/2;
end

%% Minimum distance between a point to a polygon
function d=minDistPt2Poly(pt,poly)
    d=0;
    min_dist=Inf;
    for i=1:length(poly)-1
        closest_pt=getClosestPt(pt,poly(i,:),poly(i+1,:)); 
        dp1=norm([pt-poly(i,:)]);
        dp2=norm([pt-poly(i+1,:)]);
        if ((pt(1)<poly(i,1) && pt(1)<poly(i+1,1)) ||...  
           (pt(1)>poly(i,1) && pt(1)>poly(i+1,1)))
            dline=Inf;
        else
            dline=norm([pt-closest_pt]);
        end
        if (dp1<min_dist); min_dist=dp1; end
        if (dp2<min_dist); min_dist=dp2; end
        if (dline<min_dist); min_dist=dline; end
    end
    d=d+min_dist;
end

%% Get the point closest to p on the line defined by p1 and p2
% 2D point p=[x,y];
function pc=getClosestPt(p,p1,p2) 
    if p1(1)~=p2(1)
        a=(p1(2)-p2(2))/(p1(1)-p2(1));
        b=-1;
        c=p1(2)-a*p1(1);
    else
        a=1;
        b=0;
        c=-p1(1);
    end
    pc=zeros(size(p));
    pc(1)=(b*(b*p(1)-a*p(2))-a*c)/(a^2+b^2);
    pc(2)=(-a*(b*p(1)-a*p(2))-b*c)/(a^2+b^2);
end