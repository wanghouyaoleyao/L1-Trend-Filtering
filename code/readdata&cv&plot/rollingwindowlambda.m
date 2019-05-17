% for one asset, find the optimal lambda for every rolling window
% rolling window is day 1-5 ,2-6,3-7 
% rolling window exercise
% we have a total of 781997 observations
% we take a rolling window of size wS = 390*5 (one week of data)
% train the model on this data
% predict the next day volatility as equal to the last vol form the model
% move by one day (2-6)
% 
% here we choose 1-5,2-6,3-7
% day 1 to 5 (a week)
% !!! need to revise(because of holidays or missing data)
%

% set maximum regularization parameter lambda_max
% question: here, should lambda_max change in every rolling window? or it
% should be the one for whole data
% here we choose one for every rolling window
%lambda_max = l1tf_lambdamax(S);
function [lambdalist] = rollingwindowlambda(S)

parameter1=[];
parameter2=[];
parameter3=[];
parameter4=[];
parameter5=[];
parameter6=[];
for i =1:9
   parameter1 = [parameter1,i*1e-6];
   parameter2 = [parameter2,i*1e-5];
   parameter3 = [parameter3,i*1e-4];
   parameter4 = [parameter4,i*1e-3];
   parameter5 = [parameter5,i*1e-2];
   parameter6 = [parameter6,i*1e-1];
end
parameter = [parameter1,parameter2,parameter3,parameter4,parameter5,parameter6,1];


rw = 1;
lambdalist = [];
while (rw+4)*390<length(S)
    data = S((rw-1)*390+1:(rw+4)*390);
    lambda_max = l1tf_lambdamax(data);
    lambdaset = parameter*lambda_max;
    rwlambda = chooselambda(data,lambdaset);
    lambdalist = [lambdalist,rwlambda];
    rw = rw+1;
    
    
%lambda_max = l1tf_lambdamax(data);

%----------------------------------------------------------------------
% 	l1 trend filtering
%----------------------------------------------------------------------
%[z1,status] = l1tf(y, 0.001*lambda_max);
%[z2,status] = l1tf(y, 0.01*lambda_max);
%[z3,status] = l1tf_revised(y, 0.01*lambda_max);
end
lambdalist = lambdalist;
