function [MKNNG, val] = MedianKNearestNeighborGraph(dist,K,med)

% KNearestNeighborGraph - Computes the K Nearest Neighbor graph 
%
% Computes the K Nearest Neighbor Graph (KNNG) of a set of points 
% KNNG is a DIRECTED graph
%
% CALL:
% [KNNG]=KNearestNeighborGraph(dist,K)
% 
% INPUT:
% dist: NxN Euclidean distance matrix between each pair of points
%               dist(i,k)=norm(data(i,:)-data(k,:), (dist(k,k)=0)
%              
% OUPUT:
% KNNG: N cells KNNG{i}=[a b...f] set of index in data rows of the KNNG
% neighbors of i
% 
% Author   : Michael Aupetit
%            Qatar Computing Research Institute (QCRI)
%            Hamad Bin Khalifa University (HBKU)
%            Doha, Qatar
%            maupetit@qf.org.qa/michael.aupetit@gmail.com
% Last Rev : Feb 08 2016 by Michael Aupetit

N=size(dist,1);
dist(eye(N)==1)=inf;

for i = 1:N
    for j = 1:N
        if (dist(i,j)>med)
            dist(i,j) = Inf;
        end
    end
end

[val, indSort] = sort(dist);
MKNNG = cell(N,1);

for i = 1:N
%     MKNNG{i}=indSort(1:K,i)';
    for j = 1:min(K,size(val,1))
        if val(j,i)<Inf
            MKNNG{i}(j)=indSort(j,i);
        end
    end
end

%% Normal KNN Graph

% [N,N]=size(dist);
% dist(eye(N)==1)=inf;
% [val, indSort]=sort(dist);
% KNNG=cell(N,1);
% 
% for i=1:N
%     KNNG{i}=indSort(1:K,i)';
% end
% 
% for i=1:N
%     if (MKNNG{i}~=KNNG{i})
%         fprintf('different\n')
%     else
%         fprintf('same\n')
%     end
% end