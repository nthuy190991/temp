function [P,MSE] = gold_standard_algo(x,X,verbose)
% algorithm page 181 (Multiple View Geometry in Computer Vision, 2004)
n=size(x,2);

if (verbose)
    disp('Gold Standard algorithm, for estimating P knowing point correspondences')
    disp(['Taking ' num2str(n) ' point correspondences'])
end

% n=6;
% x=10*randn(2,N)+[20 10]'*ones(1,N);
% X=10*randn(3,N)+[10 60 40]'*ones(1,N);

%% Linear solution (compute an initial estimate of P using a linear method)
% Normalization (translation + scaling)
[x_tilde, T] = normalize_pts(x(:,1:n), 0);
[X_tilde, U] = normalize_pts(X(:,1:n), 0);

% [x_tilde, X_tilde, T, U] = normalization(x(:,1:n), X(:,1:n));

% Add array of 1 to be homogeneous coordinates
x_tilde = [x_tilde; ones(1, size(x_tilde,2))];
X_tilde = [X_tilde; ones(1, size(X_tilde,2))];

A=[];
for i=1:n
    A=[A; 0 0 0 0,                      -x_tilde(3,i)*X_tilde(:,i)',    x_tilde(2,i)*X_tilde(:,i)'; 
          x_tilde(3,i)*X_tilde(:,i)',   0 0 0 0,                        -x_tilde(1,i)*X_tilde(:,i)'];
end

% SVD of matrix A
[~,D,V] = svd(A);
p=V(:,end);

% Estimate the matrix for x_tilde <-> X_tilde
P_tilde=[p(1:4)'; p(5:8)'; p(9:12)'];

if (verbose)
figure
    plot(x_tilde(1,:), x_tilde(2,:), 's')
    hold on
    x2=P_tilde*X_tilde;
    plot(x2(1,:), x2(2,:), 'r*')
    title('Initial estimate of P~')
    grid on
end

%% Minimize geometric error
x0=p; % initial estimation of P
xdata=X_tilde;
ydata=x_tilde;
fun = @(x,xdata)[x(1:4)'; x(5:8)'; x(9:12)']*xdata;

options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt',...
    'OptimalityTolerance',1e-20,'FunctionTolerance',1e-20,'StepTolerance',1e-20);

p_min = lsqcurvefit(fun,x0,xdata,ydata,[],[],options);
P_tilde_min = [p_min(1:4)'; p_min(5:8)'; p_min(9:12)'];

if (verbose)
figure
    plot(x_tilde(1,:), x_tilde(2,:), 's')
    hold on
    x2=P_tilde*X_tilde;
    plot(x2(1,:), x2(2,:), 'r*')
    x3=P_tilde_min*X_tilde;
    plot(x3(1,:), x3(2,:), 'g*')
    title('After minization')
    grid on
    legend('x~','P~_{init}X~','P~_{final}X~')
end

%% denormalization (matrix for x <-> X)
P = inv(T)*P_tilde_min*U;

if (verbose)
figure
    plot(x(1,1:n), x(2,1:n), 's')
    hold on
    x2=P*[X(:,1:n); ones(1,n)];
    plot(x2(1,:), x2(2,:), 'r*')
    title('If the squares and stars overlap, then it is correct')
    grid on
    legend('x','PX')
end

error=zeros(1,n);
x2=P*[X(:,1:n); ones(1,n)];
for i=1:n
    error(i) = norm([x(:,i);1] - x2(:,i));
end
MSE=mean(error);
end