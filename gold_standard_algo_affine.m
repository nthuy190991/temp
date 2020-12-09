function [P] = gold_standard_algo_affine(x,X)
% algorithm page 185 (Multiple View Geometry in Computer Vision, 2004)

disp('Gold Standard algorithm, for estimating P of an affine camera')
n=4; 
disp(['Taking ' num2str(n) ' point correspondences'])

% Normalization (translation + scaling)
[x_tilde, T] = normalize_pts(x(:,1:n), 0);
[X_tilde, U] = normalize_pts(X(:,1:n), 0);

% Add array of 1 to be homogeneous coordinates
X_tilde=[X_tilde; ones(1, size(X_tilde,2))];

b=[];
for i=1:n
    b=[b; x_tilde(1,i); x_tilde(2,i)];
end

A=[];
for i=1:n
    A=[A; X_tilde(:,i)' 0 0 0 0; 0 0 0 0 X_tilde(:,i)'];
end

% Pseudo-inverse of A
Ap = (A'*A)\A';

% Estimate the matrix for x_tilde <-> X_tilde
p=Ap*b;
P_tilde=[p(1:4)'; p(5:8)'; 0 0 0 1];

figure
    plot(x_tilde(1,:), x_tilde(2,:), 's')
    hold on
    x2=P_tilde*X_tilde;
    plot(x2(1,:), x2(2,:), 'r*')
    title('If the squares and stars overlap, then it is correct')
    grid on

% denormalization (matrix for x <-> X)
P = inv(T)*P_tilde*U;

figure
    plot(x(1,1:n), x(2,1:n), 's')
    hold on
    x2=P*[X(:,1:n); ones(1,n)];
    plot(x2(1,:), x2(2,:), 'r*')
    title('If the squares and stars overlap, then it is correct')
    grid on
end