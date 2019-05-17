function[x,status] = l1Mono(y,lambda)

% last version 
% check the dimension and change the differential term

% [x,status] = l1tf(y,lambda)
%
% finds the solution of the l1 trend and monotonically estimation problem
%
%  minimize    (1/2)||y-x||^2+lambda*||Dx||_1, 
%  subject to Fx is smaller or equal to 0
%
% with variable x, and problem data y and lambda, with lambda >0.
% D is the second difference matrix, with rows [0... -1 2 -1 ...0]
% F is the monotonicity constraint matrix, with rows [1 -1 0 ... 0 ]
%
% So we transfer the optimization problem to
% (1/2)||y-x||^2+lambda*||z||_1+w*Fx
% subject to w is larger or equal to 0 and Dx =  z

% Then
% (1/2)||y-x||^2+lambda*||z||_1+w*Fx+v*(Dx-z)

% and the dual problem
%
%  minimize (1/2) v'DD'v + (1/2)v'DF'w + (1/2)w'FD'v + (1/2)w'FF'v
%  -y'D'v - y'F'w
%  subject to  norm(v,inf) <= lambda and -w <= 0
%
% with variable v and w.
%
% Input arguments:
%
% - y:          n-vector; original signal
% - lambda:     scalar; positive regularization parameter
%
% Output arguments:
%
% - x:          n-vector; primal optimal point
% - z:          2(n-1)-vector; dual optimal point
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

% PARAMETERS
ALPHA     = 0.45;   % backtracking linesearch parameter (0,0.5]
BETA      = 0.95;    % backtracking linesearch parameter (0,1)
MU        = 2;      % IPM parameter: t update
%MAXITER   = 80;     % IPM parameter: max iteration of IPM
MAXITER   = 80;
%MAXLSITER = 20;     % IPM parameter: max iteration of line search
MAXLSITER = 20;
TOL       = 1e-4;   % IPM parameter: tolerance

% DIMENSIONS
n   = length(y);    % length of signal x
m1  = n-2;          % length of Dx
m2  = n-1;          % length of Fx

%warning off all;

% OPERATOR MATRICES
% dim for D:(n-2)*n    
% dim for F: (n-1)*n
I2  = speye(n-2,n-2);
O2  = zeros(n-2,1);
D   = [I2 O2 O2]+[O2 -2*I2 O2]+[O2 O2 I2];
temp= speye(n,n);
temp(n,n)=0;

D2 =  temp + [zeros(n, 1), [-speye(n - 1); zeros(1, (n - 1))]];
F = D2(1:(n-1),:);

DDT = D*D';
DFT = sparse(D*F');
FDT = sparse(F*D');
FFT = F*F';
Dy  = D*y;
Fy  = F*y;

% VARIABLES
v   = zeros(m1,1);   % dual variable
w   = zeros(m2,1)+0.00001;
mu1 = ones(m1,1);    % dual of dual variable
mu2 = ones(m1,1);    % dual of dual variable
mu3 = ones(m2,1);    % dual of dual variable

t    = 1e-10; 
pobj =  Inf;
dobj =  0;
step =  Inf;
f1   =  v-lambda;
f2   = -v-lambda;
f3   = -w;

disp('--------------------------------------------');
disp('l1 trend filtering method plus constraint via primal-dual algorithm');

disp('--------------------------------------------');
%disp(sprintf('\n%s %13s %12s %8s %9s %17s \n',...
    %'Iteration','Primal obj.','Dual obj.','Gap','t','Step size'));
disp(sprintf('\n%s %13s %12s %8s\n',...
    'Iteration','Primal obj.','Dual obj.','Gap'));

%----------------------------------------------------------------------
%               MAIN LOOP
%----------------------------------------------------------------------
for iters = 0:MAXITER
    
    DTv  = (v'*D)';
    DDTv = D*DTv;
    FTw  = (w'*F)'; 
    DFTw = D*FTw;
    FDTv = F*DTv;
    FFTw = F*FTw;
    
    w1    = Dy-(mu1-mu2);
    w2    = Fy+mu3;

    
    % two ways to evaluate primal objective:
    % 1) using dual variable of dual problem
    % 2) using optimality condition
    
    
    %%%%%% question part
   % pobj1 = 0.5*w'*(DDT\w)+lambda*sum(mu1+mu2);
   % pobj1 = 0.5*w1'*(DDT\w1)+0.5*w2'*(FFT\w2)+lambda*sum(mu1+mu2);
    pobj1 = -0.5*v'*w1-0.5*w'*w2+w1'*v+w2'*w+lambda*sum(mu1+mu2);
    %pobj2 = 0.5*DTz'*DTz+lambda*sum(abs(Dy-DDTz));
    pobj2 = 0.5*(FTw+DTv)'*(FTw+DTv)+lambda*sum(abs(Dy-DFTw-DDTv))+(FTw)'*y-(FTw)'*FTw-(FTw)'*DTv;
    pobj = min(pobj1,pobj2);
    %pobj=pobj1;
    dobj = -0.5*(DTv)'*DTv-0.5*(DTv)'*FTw-0.5*(FTw)'*DTv-0.5*(FTw)'*FTw+Dy'*v+Fy'*w;
%     pobj = min(min(pobj1,pobj2), pobj);
%     dobj = max(-0.5*DTz'*DTz+Dy'*z, dobj);
    gap  =  pobj - dobj;
    %%%%%% question part
    
    
    

    %disp(sprintf('%6d %16.4e %13.5e %10.2e %11.2e %13.2e',...
        %iters, pobj, dobj, gap, t, step));
    disp(sprintf('%6d %15.4e %13.5e %10.2e',...
        iters, pobj, dobj, gap));

    % STOPPING CRITERION
    if (gap <= TOL)
        status = 'solved';
        disp(status);
        x = y-F'*w-D'*v;
        return;
    end

    if (step >= 0.2) 
        t =max(2*m2*MU/gap, 1.2*t);
    end

    % CALCULATE NEWTON STEP
    
    rv     =  DDTv + DFTw - w1;
    rw     =  FDTv + FFTw - w2;
    S1     =  DDT-sparse(1:m1,1:m1,mu1./f1+mu2./f2);
    S2     =  FFT-sparse(1:m2,1:m2,mu3./f3);
    r1     = -DDTv - DFTw + Dy + (1/t)./f1 - (1/t)./f2;
    r2     = -FDTv - FFTw + Fy - (1/t)./f3;
    
   % DFTFDTinv = (DFT*FDT)\eye(m1);
   % FDTinv = DFTFDTinv*DFT;
   % DFTinv = FDT*DFTFDTinv;
   
   %%%% use the block inverse formula
   invS1 = S1\eye(size(S1));
   invS2FDTinvS1DFT = (S2-FDT*invS1*DFT)\eye(size(S2));
   Ginv = [invS1+invS1 * DFT*invS2FDTinvS1DFT*FDT*invS1,-invS1*DFT*invS2FDTinvS1DFT;-invS2FDTinvS1DFT*FDT*invS1, invS2FDTinvS1DFT];
   %Ginv = [inv(S1)+S1\DFT*(S2-FDT*S1\DFT)\FDT*inv(S1),-inv(S1)*DFT*inv(S2-FDT*S1\DFT);-inv(S2-FDT*S1\DFT)*FDT*inv(S1), inv(S2-FDT*S1\DFT)];
     
    result = Ginv * [r1;r2];
    
    dv = result(1);
    dw = result(2);
   
    %dv  =  S1\r1;
    %dw  =  S2\r2;
    
    
    %rightv = S2*DFTinv*r1-r2;
    %leftv  = S2*DFTinv*S1-FDT;
    %dv     = leftv\rightv;
    %rightw = S1*FDTinv*r2-r1;
    %leftw  = S1*FDTinv*S2-DFT;
    %dw     = leftw\rightw;
    dmu1    = -(mu1+((1/t)+dv.*mu1)./f1);
    dmu2    = -(mu2+((1/t)-dv.*mu2)./f2);
    dmu3    = -(mu3+((1/t)-dw.*mu3)./f3);

    
    resDual = [rv;rw];
    resCent = [-mu1.*f1-1/t; -mu2.*f2-1/t; -mu3.*f3-1/t];
    residual= [resDual; resCent];

    % BACKTRACKING LINESEARCH
    negIdx1 = (dmu1 < 0); 
    negIdx2 = (dmu2 < 0);
    negIdx3 = (dmu3 < 0);
    step = 1;
    if (any(negIdx1))
        step = min( step, 0.99*min(-mu1(negIdx1)./dmu1(negIdx1)) );
    end
    if (any(negIdx2))
        step = min( step, 0.99*min(-mu2(negIdx2)./dmu2(negIdx2)) );
    end
    if (any(negIdx3))
        step = min( step, 0.99*min(-mu3(negIdx3)./dmu3(negIdx3)) );
    end

    for liter = 1:MAXLSITER
        newv    =  v  + step*dv;
        neww    =  w  + step*dw;
        newmu1  =  mu1 + step*dmu1;
        newmu2  =  mu2 + step*dmu2;
        newmu3  =  mu3 + step*dmu3;
        newf1   =  newv - lambda;
        newf2   = -newv - lambda;
        newf3   = -neww;

        % UPDATE RESIDUAL
        
        newResDual  = [DDT*newv + DFT*neww - Dy + newmu1 - newmu2; FDT*newv + FFT*neww - Fy - newmu3];
        %newResDual  = [DDT*newv  - Dy + newmu1 - newmu2; FFT*neww - Fy - newmu3];
        newResCent  = [-newmu1.*newf1-1/t; -newmu2.*newf2-1/t;-newmu3.*newf3-1/t];
        newResidual = [newResDual; newResCent];
        
        if ( max([max(newf1),max(newf2),max(newf3)]) < 0 && ...
            norm(newResidual) <= (1-ALPHA*step)*norm(residual) )
            break;
        end
        step = BETA*step;
    end
    % UPDATE PRIMAL AND DUAL VARIABLES
    v  = newv; w = neww;  
    mu1 = newmu1; mu2 = newmu2; mu3 = newmu3; 
    f1 = newf1; f2 = newf2;  f3 = newf3;
end

% The solution may be close at this point, but does not meet the stopping
% criterion (in terms of duality gap).
x = y-D'*v-F'*w;
if (iters >= MAXITER)
    status = 'maxiter exceeded';
    disp(status);
    return;
end
