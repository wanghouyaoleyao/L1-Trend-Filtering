clear all;
clc;
% read data from txt file
fid = fopen("E:\HF-PCA\data\DJ30\Copy of AAPL.txt");
data = textscan(fid,"%f%f:%f%f%f%f%f%f","delimiter",",");
fclose(fid);
date = data{1};
hour = data{2};
minute = data{3};
close1 = data{7};
%find out all the missing value
%tryi = [0;0;0];
%for i = 1:781486
%    if (minute(i+1) ~= minute(i)+1 && minute(i+1) ~= minute(i)-59 && minute(i+1) ~= minute(i)-29 )
%        tryi = [tryi(1,:),date(i),date(i+1);tryi(2,:),hour(i),hour(i+1);tryi(3,:),minute(i),minute(i+1)];
%    end
         
%end
%tryi = tryi';
%fill the NA
closefull = close1;
logprice = log(closefull);
% log return
logreturn =  diff(logprice);
lrs = logreturn.^2;
S = cumsum(lrs);
S = [0;S];






    

