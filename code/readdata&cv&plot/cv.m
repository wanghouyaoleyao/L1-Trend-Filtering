function[MSE] = cv(y,lambda)
MSE = [];
for i =1:5    
   traindata = [y(1:(i-1)*390);y(i*390+1:1950)];
   validatedata = y((i-1)*390+1:i*390);
   
   [z,status] = l1Mono(traindata, lambda);
   [z,status] = l1Mono(traindata, lambda);
   vol = diff(z);
   
   % question: what should be the volatility for validate set
   % last volatility of the previous folder 
   % for the first one, first volatility of the next folder
   if i ==1
       predvol = vol(1);
   else
       predvol = vol((i-1)*390-1);
   end
   
   actualvol = diff(validatedata); 
   
   MSEeach = sum((actualvol-predvol).^2)/length(actualvol);
   
   % calculate MSE in validatedata
   MSE = [MSE,MSEeach];
   
    
end
MSE = mean(MSE);

