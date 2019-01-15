

% function [PE_tot,force] = zNNG98_H2BR(in,minin,maxin,minout,maxout)

function [Vhat_Milind,DPEDxyz] = zNNG98_H2BR(in,coord)


% load('MILIND_net_fromOhm_lm_df1_f7.mat');
% load('MILIND_inOut_fromOhm_minmax.mat');
% 
% net_object=newff([minin maxin],[150 1],{'tansig' 'purelin'});
% net_object.IW{1,1}=net.IW{1,1};
% net_object.LW{2,1}=net.LW{2,1};
% net_object.b{1}=net.b{1};
% net_object.b{2}=net.b{2};
% 
% 
% 
% inn = tramnmx(in,minin,maxin); % inn is the scaled input between 0-1
% 
% % [PE_tot,force] = postmnmx(sim(net,inn),minout,maxout);
% [PE_tot] = postmnmx(sim(net_object,inn),minout,maxout);
% 
% 
% n1 = net_object.IW{1,1}*inn + net_object.b{1};
% a1 = tansig(n1);
% 
% df1 = (ones(length(net_object.IW{1,1}),1) - a1.*a1);
% 
% deriv = net_object.IW{1,1}'*(df1.*[net_object.LW{2,1}]');
% 
% DPEDin = deriv.*((maxout-minout)./(maxin-minin));
% 
% force= DPEDin; 


% % % % % 
x=coord(:,1); y=coord(:,2); z=coord(:,3);
% z(2)= z(2) - 0.003;

r1 = sqrt( (x(1)-x(2))^2 + (y(1)-y(2))^2 + (z(1)-z(2))^2 );
r2 = sqrt( (x(1)-x(3))^2 + (y(1)-y(3))^2 + (z(1)-z(3))^2 );
r3 = sqrt( (x(2)-x(3))^2 + (y(2)-y(3))^2 + (z(2)-z(3))^2 );
cosTh = (r1^2 + r2^2 - r3^2)/(2*r1*r2);
theta = acos(cosTh);

% in = [r1;r2;theta];
in= [r1;r2;r3];

% in= [3.6196; 1.5143; 2.5011];
% in(1)= in(1) - 0.001;



load('net_fc_2.mat');

load('net_fr_HH_2.mat');
load('net_fr_HBr_2.mat');

load('net_ftheta_HHBr_2.mat');


r_H1H2= in(1,1);
r_H1Br= in(2,1);
r_H2Br= in(3,1);


fr_H1H2= 0; fr_H1Br= 0; fr_H1H2= 0; ftheta_H1H2Br= 0;

fc=1;%sim(net_fc,rij);

% for iQ=1:1:Q
    fr_H1H2= sim(net_fr_HH,r_H1H2);
    fr_H1Br= sim(net_fr_HBr,r_H1Br);
    fr_H2Br= sim(net_fr_HBr,r_H2Br);
    
    ftheta_H1H2Br= sim(net_ftheta_HHBr,[r_H1H2; r_H1Br; r_H2Br]);
    
% end
Vhat_Milind= (fr_H1H2+fr_H1Br+fr_H2Br+ftheta_H1H2Br);
% % % % % % ----------

W1= net_fr_HH.IW{1,1};  W2= net_fr_HH.LW{2,1};  b1= net_fr_HH.b{1};  b2= net_fr_HH.b{2};

n1 = W1*r_H1H2 + b1;
a1 = logsig(n1);

df1 = a1.*(ones(length(W1),1) - a1);

deriv_fr_H1H2 = W1'*(df1.*[W2]');
% % % % % % ----------

W1= net_fr_HBr.IW{1,1};  W2= net_fr_HBr.LW{2,1};  b1= net_fr_HBr.b{1};  b2= net_fr_HBr.b{2};

n1 = W1*r_H1Br + b1;
a1 = logsig(n1);

df1 = a1.*(ones(length(W1),1) - a1);

deriv_fr_H1Br = W1'*(df1.*[W2]');
% % % % % % ----------

W1= net_fr_HBr.IW{1,1};  W2= net_fr_HBr.LW{2,1};  b1= net_fr_HBr.b{1};  b2= net_fr_HBr.b{2};

n1 = W1*r_H2Br + b1;
a1 = logsig(n1);

df1 = a1.*(ones(length(W1),1) - a1);

deriv_fr_H2Br = W1'*(df1.*[W2]');
% % % % % % ----------

W1= net_ftheta_HHBr.IW{1,1};  W2= net_ftheta_HHBr.LW{2,1};  b1= net_ftheta_HHBr.b{1};  b2= net_ftheta_HHBr.b{2};

n1 = W1*[r_H1H2; r_H1Br; r_H2Br] + b1;
a1 = logsig(n1);

df1 = a1.*(ones(length(W1),1) - a1);

deriv_ftheta_H1H2Br = W1'*(df1.*[W2]');
% % % % % % ----------


DPEDxyz(1,1)= deriv_fr_H1H2*((x(1)-x(2))/r_H1H2) + deriv_fr_H1Br*((x(1)-x(3))/r_H1Br) + ...
                deriv_ftheta_H1H2Br(1,1)*((x(1)-x(2))/r_H1H2) + deriv_ftheta_H1H2Br(2,1)*((x(1)-x(3))/r_H1Br);
            
DPEDxyz(1,2)= deriv_fr_H1H2*((y(1)-y(2))/r_H1H2) + deriv_fr_H1Br*((y(1)-y(3))/r_H1Br) + ...
                deriv_ftheta_H1H2Br(1,1)*((y(1)-y(2))/r_H1H2) + deriv_ftheta_H1H2Br(2,1)*((y(1)-y(3))/r_H1Br);

DPEDxyz(1,3)= deriv_fr_H1H2*((z(1)-z(2))/r_H1H2) + deriv_fr_H1Br*((z(1)-z(3))/r_H1Br) + ...
                deriv_ftheta_H1H2Br(1,1)*((z(1)-z(2))/r_H1H2) + deriv_ftheta_H1H2Br(2,1)*((z(1)-z(3))/r_H1Br);

            

            
DPEDxyz(2,1)= deriv_fr_H1H2*((x(2)-x(1))/r_H1H2) + deriv_fr_H2Br*((x(2)-x(3))/r_H2Br) + ...
                deriv_ftheta_H1H2Br(1,1)*((x(2)-x(1))/r_H1H2) + deriv_ftheta_H1H2Br(3,1)*((x(2)-x(3))/r_H2Br);
            
DPEDxyz(2,2)= deriv_fr_H1H2*((y(2)-y(1))/r_H1H2) + deriv_fr_H2Br*((y(2)-y(3))/r_H2Br) + ...
                deriv_ftheta_H1H2Br(1,1)*((y(2)-y(1))/r_H1H2) + deriv_ftheta_H1H2Br(3,1)*((y(2)-y(3))/r_H2Br);

DPEDxyz(2,3)= deriv_fr_H1H2*((z(2)-z(1))/r_H1H2) + deriv_fr_H2Br*((z(2)-z(3))/r_H2Br) + ...
                deriv_ftheta_H1H2Br(1,1)*((z(2)-z(1))/r_H1H2) + deriv_ftheta_H1H2Br(3,1)*((z(2)-z(3))/r_H2Br);
            
 

            
DPEDxyz(3,1)= deriv_fr_H1Br*((x(3)-x(1))/r_H1Br) + deriv_fr_H2Br*((x(3)-x(2))/r_H2Br) + ...
                deriv_ftheta_H1H2Br(2,1)*((x(3)-x(1))/r_H1Br) + deriv_ftheta_H1H2Br(3,1)*((x(3)-x(2))/r_H2Br);

DPEDxyz(3,2)= deriv_fr_H1Br*((y(3)-y(1))/r_H1Br) + deriv_fr_H2Br*((y(3)-y(2))/r_H2Br) + ...
                deriv_ftheta_H1H2Br(2,1)*((y(3)-y(1))/r_H1Br) + deriv_ftheta_H1H2Br(3,1)*((y(3)-y(2))/r_H2Br);

DPEDxyz(3,3)= deriv_fr_H1Br*((z(3)-z(1))/r_H1Br) + deriv_fr_H2Br*((z(3)-z(2))/r_H2Br) + ...
                deriv_ftheta_H1H2Br(2,1)*((z(3)-z(1))/r_H1Br) + deriv_ftheta_H1H2Br(3,1)*((z(3)-z(2))/r_H2Br);
            

DPEDxyz;

% DPEDin = deriv.*((maxout-minout)./(maxin-minin));


% % % % % 