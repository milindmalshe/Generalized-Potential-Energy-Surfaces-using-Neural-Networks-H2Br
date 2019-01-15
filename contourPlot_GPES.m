
load('net_fc_2.mat');

load('net_fr_HH_2.mat');
load('net_fr_HBr_2.mat');

load('net_ftheta_HHBr_2.mat');


r_start= 0.5;
r_end= 4.0;

stepSize= 0.1;

numPoints= (r_end-r_start)/stepSize + 1;

for i= 1:numPoints
%     r_H1H2(1,i)= r_start + (i-1)*stepSize;
    
   for j= 1:numPoints
       for k= 1:numPoints
           
           r_HHBr(i,j,k)= r_start + (i-1)*stepSize + (j-1)*stepSize + (k-1)*stepSize;
           
       end
   end
end