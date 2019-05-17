
% for every rolling window, choose one optimal lambda
function [finallambda] = chooselambda(y,lambdaset)
mselist = [];
for i = 1:length(lambdaset)
    disp([ 'lambda i =',num2str(i),' out of ',num2str(length(lambdaset))])
    mse = cv(y,lambdaset(i));
    mselist = [mselist,mse];
end
[c,i] = min(mselist);
finallambda = lambdaset(i);

    