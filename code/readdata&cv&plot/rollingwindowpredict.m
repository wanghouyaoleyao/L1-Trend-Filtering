% for one asset, predict the next-day volatility for every rolling window
% question: rolling window should be day 1-5 ,2-6,3-7 or 1-5,6-10,11-15
% !!! need to revise(because of holidays or missing data)
%

% set maximum regularization parameter lambda_max
% question: here, should lambda_max change in every rolling window? or it
% should be the one for whole data
% here we choose one for every rolling window
%lambda_max = l1tf_lambdamax(S);
function [predvollist,vollist,lamfraclist] = rollingwindowpredict(S)

%parameter1=[];
%parameter2=[];
%parameter3=[];
%parameter4=[];
%parameter5=[];
%parameter6=[];
%for i =1:9
%   parameter1 = [parameter1,i*1e-6];
%   parameter2 = [parameter2,i*1e-5];
%   parameter3 = [parameter3,i*1e-4];
%   parameter4 = [parameter4,i*1e-3];
%   parameter5 = [parameter5,i*1e-2];
%   parameter6 = [parameter6,i*1e-1];
%end
%parameter = [parameter1,parameter2,parameter3,parameter4,parameter5,parameter6,1];

parameter = exp(-10:0.5:log(1));

%parameter = parameter1(1:10);
cvx_precision best
rw = 1;
predvollist = [];
lamfraclist = [];
vollist = [];
while (rw+4)*390<length(S)
    disp([' rolling window = ',num2str(rw),' out of ',num2str(floor(length(S)/390 -4))])
    data = S((rw-1)*390+1:(rw+4)*390);
    lambda_max = l1tf_lambdamax(data);
    lambdaset = parameter*lambda_max;
    % choose the best lambda for this rolling window based on 5-fold cv
    rwlambda = chooselambda(data,lambdaset);
    % estimate the volatility for whole rolling window based on cv result
    [z,status] = l1tf_cvx1(data, rwlambda);
    vol = diff(z);
    vollist = [vollist,vol];
    predvol = vol(length(vol));
    predvollist = [predvollist,predvol];
    rw = rw+1;
    lam_frac = rwlambda/lambda_max;
    lamfraclist = [lamfraclist,lam_frac];
    

%----------------------------------------------------------------------
% 	l1 trend filtering
%----------------------------------------------------------------------
%[z1,status] = l1tf(y, 0.001*lambda_max);
%[z2,status] = l1tf(y, 0.01*lambda_max);
%[z3,status] = l1tf_revised(y, 0.01*lambda_max);
end
