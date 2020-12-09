function [xyn, XYZn, T, U] = normalization(xy, XYZ)

%get number of points
[~,N] = size(xy);
%data normalization
%first compute centroid
xy_centroid = mean(xy,2);
XYZ_centroid = mean(XYZ,2);

%then, compute scale
xym = repmat(xy_centroid,1,N);
XYZm = repmat(XYZ_centroid,1,N);
sigma2 = mean(sqrt(sum((xy-xym).^2,1)))/sqrt(2);
sigma3 = mean(sqrt(sum((XYZ-XYZm).^2,1)))/sqrt(3);
%create T and U transformation matrices
T = [1/sigma2 0 -xy_centroid(1)/sigma2;
0 1/sigma2 -xy_centroid(2)/sigma2;
0 0 1];
U = [1/sigma3 0 0 -XYZ_centroid(1)/sigma3;
0 1/sigma3 0 -XYZ_centroid(2)/sigma3;
0 0 1/sigma3 -XYZ_centroid(3)/sigma3;
0 0 0 1];

%and normalize the points according to the transformations
xyn = T*[xy;ones(1,N)];
XYZn = U*[XYZ;ones(1,N)];

end
