%
%   Simple usage example of l1tf
%
% rand('state',0); randn('state',0); % make data reproducible

%----------------------------------------------------------------------
% 	data generation
%----------------------------------------------------------------------
%y = 1000*S(1:300000);
y = 100*S(1:2400);
% eps =randn(length(S),1);

% set maximum regularization parameter lambda_max
lambda_max = l1tf_lambdamax(y);

%----------------------------------------------------------------------
% 	l1 trend filtering
%----------------------------------------------------------------------
%[z1,status] = l1tf(y, 0.001*lambda_max);
%[z2,status] = l1tf(y, 0.01*lambda_max);
%[z3,status] = l1tf_revised(y, 0.01*lambda_max);
lambda_used = 1e-6*lambda_max;
[z4,status] = l1tf_cvx1(y, lambda_used);
disp(sprintf('  lambda_max : %e', lambda_max)) ;
disp(sprintf('  lambda used : %e', lambda_used)) ;
subplot(2,1,1)
plot(diff(y),'.-')
title('ret^2')
subplot(2,1,2)
plot(1:size(z4,1)-1,diff(z4),'.-r');
title('estvol')


% uncomment line below to solve l1 trend filtering problem using CVX
%[z1,status] = l1tf_cvx(y, 0.001*lambda_max);
%[z2,status] = l1tf_cvx(y, 0.005*lambda_max);
%[z3,status] = l1tf_cvx(y, 0.025*lambda_max);
%plot(1:size(y,1),y, 'r',1:size(z3,1),z3,'b'); %ylim([minx maxx]); 

%hold on
%plot(1:size(z3,1)-1,y(2:end),'r');
%figure(1);
%plot(1:size(z2,1)-1,diff(z2),'b');
%hold on
%plot(1:size(z3,1)-1,diff(z3),'r');

%plot(1:size(z3,1)-1,diff(z3),'b'); %ylim([minx maxx]); 
%hold off

%----------------------------------------------------------------------
% 	plot results
%----------------------------------------------------------------------
%xyzs = [y z1 z2 z3];
%maxx = max(max(xyzs));
%minx = min(min(xyzs));

%plot(1:n,z1,'b',1:n,y,'r'); ylim([minx maxx]); title('lambda = 0.001');
%plot(1:n,z2,'b',1:n,y,'r'); ylim([minx maxx]); title('lambda = 0.005');
%plot(1:n,z3,'b',1:n,y,'r'); ylim([minx maxx]); title('lambda = 0.025');
