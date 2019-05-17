function [x,status] = l1tf_cvx1(y,lambda)
% [x,status] = l1tf_cvx(y,lambda)
%
% finds the solution of the l1 trend estimation problem
%
%  minimize    (1/2)||y-x||^2+lambda*||Dx||_1,
%  plus the codition that xt is less or equal than xt+1
%
% with variable x, and problem data y and lambda, with lambda >0.
% D is the second difference matrix, with rows [0... -1 2 -1 ...0]
%
% uses CVX, a package for generic convex optimization
% www.stanford.edu/~boyd/cvx/
%
% CVX is not optimized for the l1 trend filtering problem, and
% can be slow or give low accuracy results for large problems
%
% Input arguments:
%
% - y:          n-vector; original signal
% - lambda:     scalar; positive regularization parameter
%
% Output arguments:
%
% - x:          n-vector; primal optimal point
% - status:     string;
%               'solved', 'maxiter exceeded'
%
% for more details,
% see "l1 Trend Filtering", S. Kim, K. Koh, ,S. Boyd and D. Gorinevsky
% www.stanford.edu/~boyd/l1_trend_filtering.html
%
 
%----------------------------------------------------------------------
%               INITIALIZATION
%----------------------------------------------------------------------

% DIMENSIONS
n   = length(y);    % length of signal x
m   = n-2;          % length of Dx

% OPERATOR MATRICES
I2  = speye(n-2,n-2);
O2  = zeros(n-2,1);
D   = [I2 O2 O2]+[O2 -2*I2 O2]+[O2 O2 I2];
% changes here
temp= speye(n,n);
temp(n,n)=0;
D2 =  temp + [zeros(n, 1), [-speye(n - 1); zeros(1, (n - 1))]];
F = D2(1:(n-1),:);


%disp('Solving l1 trend filtering problem using CVX.');
%disp('Can be slow, or give inaccurate results for large problems.');

cvx_begin quiet
    cvx_precision best
    variable x(n)
    minimize( 0.5*sum_square(y-x)+lambda*norm(D*x,1) );
    subject to 
        F*x <= zeros((n-1),1); 
cvx_end
status = cvx_status;