function [Force]=calc_force_GPES(net_fc, net_fr_HH, net_fr_HBr,net_ftheta_HHBr,rs, clust_size, type)


% %%% ONLY FOR TESTING PURPOSE %%%%%%%%%%%%%%%%%%
% data = xlsread('H2Br_1Point.xlsx');
% 
% %no of data pts
% % Q=length(data);
% Q=1;
% 
% %in future clust_size and type will be read from the data file
% total_columns=length(data(1,:));
% 
% % column 1 in the data file -> cluster size
% max_clust_size=max(data(1:Q,1)); %5
% 
% 
% % For now, all the atoms are of same type. In future it will be for different types
% type(1:max_clust_size)=[1,1,35];
% 
% 
% for iQ=1:1:Q
%         
%     clust_size(iQ)=data(iQ,1);     %#ok<AGROW>
%     colcount=2;
%     for it=1:1:clust_size(iQ)
%         for jt=it:1:clust_size(iQ)
%             if jt>it
%                 
%                 rs(it,jt,iQ)=data(iQ,colcount); %#ok<AGROW>
%                 rs(jt,it,iQ)=rs(it,jt,iQ); %#ok<AGROW>
%                 colcount=colcount+1;
%             end
%             
%         end
%     end
%     V(iQ)=data(iQ,total_columns); %#ok<AGROW>
% end
% load('net_fc_2'); load('net_fr_HH_2'); load('net_fr_HBr_2');load('net_ftheta_HHBr_2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Force(1:3)=0.0;
HHIWeights=net_fr_HH.iw{1,1};
HHLWeights=net_fr_HH.lw{2,1};
HHIBiases=net_fr_HH.b{1,1};
HHLBiases=net_fr_HH.b{2,1};

HBrIWeights=net_fr_HBr.iw{1,1};
HBrLWeights=net_fr_HBr.lw{2,1};
HBrIBiases=net_fr_HBr.b{1,1};
HBrLBiases=net_fr_HBr.b{2,1};

HHBrIWeights=net_ftheta_HHBr.iw{1,1};
HHBrLWeights=net_ftheta_HHBr.lw{2,1};
HHBrIBiases=net_ftheta_HHBr.b{1,1};
HHBrLBiases=net_ftheta_HHBr.b{2,1};


iQ=1;
%converting rs to my format
for i=1:1:clust_size(iQ)
    for j=1:1:clust_size(iQ)
        if i~=j
            r(i,j)=rs(i,j,iQ); %#ok<AGROW>
        end
    end
end
for i=1:1:clust_size(iQ)
    for j=1:1:clust_size(iQ)
        if j>i
            
            identify2=sprintf('%d %d', type(i),type(j));
            switch identify2
                case {'1 1'}

                    Product1=HHIWeights.*logsig('dn',HHIWeights*r(i,j)+HHIBiases);
                    Product2= HHLWeights*Product1;
                    Product3=Product2*purelin('dn',HHLWeights*logsig(HHIWeights*r(i,j)+HHIBiases)+HHLBiases);
                    Force(i)=Force(i)-Product3;
                    Force(j)=Force(j)+Product3;

                case {'1 35', '35 1'}
                    Product1=HBrIWeights.*logsig('dn',HBrIWeights*r(i,j)+HBrIBiases);
                    Product2= HBrLWeights*Product1;
                    Product3=Product2*purelin('dn',HBrLWeights*logsig(HBrIWeights.*r(i,j)+HBrIBiases)+HBrLBiases);
                    Force(i)=Force(i)-Product3;
                    Force(j)=Force(j)+Product3;
                otherwise
                    disp(identify2)

            end

        end
        for k=1:1:clust_size(iQ)
            if (j>i && k>i && k>j)
              
                identify3=sprintf('%d %d %d', type(i),type(j),type(k));
                switch identify3
                    case {'1 1 35','1 35 1','35 1 1'}
                        Product12ij_1=(HHBrIWeights*[1;0;0]).*logsig('dn',HHBrIWeights*[r(i,j);r(i,k);r(j,k)]+HHBrIBiases);
                        Product12ik_1=(HHBrIWeights*[0;1;0]).*logsig('dn',HHBrIWeights*[r(i,j);r(i,k);r(j,k)]+HHBrIBiases);
                        Product12jk_1=(HHBrIWeights*[0;0;1]).*logsig('dn',HHBrIWeights*[r(i,j);r(i,k);r(j,k)]+HHBrIBiases);
                        
%                         Product12ij_2=(HHBrIWeights*[0;0;1]).*logsig('dn',HHBrIWeights*[r(i,k);r(j,k);r(i,j)]+HHBrIBiases);
%                         Product12ik_2=(HHBrIWeights*[1;0;0]).*logsig('dn',HHBrIWeights*[r(i,k);r(j,k);r(i,j)]+HHBrIBiases);
%                         Product12jk_2=(HHBrIWeights*[0;1;0]).*logsig('dn',HHBrIWeights*[r(i,k);r(j,k);r(i,j)]+HHBrIBiases);
%                          
%                         Product12ij_3=(HHBrIWeights*[0;1;0]).*logsig('dn',HHBrIWeights*[r(j,k);r(i,j);r(i,k)]+HHBrIBiases);
%                         Product12ik_3=(HHBrIWeights*[0;0;1]).*logsig('dn',HHBrIWeights*[r(j,k);r(i,j);r(i,k)]+HHBrIBiases);
%                         Product12jk_3=(HHBrIWeights*[1;0;0]).*logsig('dn',HHBrIWeights*[r(j,k);r(i,j);r(i,k)]+HHBrIBiases);
                        
                        
                        Product22ij_1=HHBrLWeights*Product12ij_1;
                        Product22ik_1=HHBrLWeights*Product12ik_1;
                        Product22jk_1=HHBrLWeights*Product12jk_1;
                        
%                         Product22ij_2=HHBrLWeights*Product12ij_2;
%                         Product22ik_2=HHBrLWeights*Product12ik_2;
%                         Product22jk_2=HHBrLWeights*Product12jk_2;
%                          
%                         Product22ij_3=HHBrLWeights*Product12ij_3;
%                         Product22ik_3=HHBrLWeights*Product12ik_3;
%                         Product22jk_3=HHBrLWeights*Product12jk_3;
                        
                        
                        Product33ij_1=Product22ij_1*purelin('dn',HHBrLWeights*logsig(HHBrIWeights*[r(i,j);r(i,k);r(j,k)]+HHBrIBiases)+HHBrLBiases);
                        Product33ik_1=Product22ik_1*purelin('dn',HHBrLWeights*logsig(HHBrIWeights*[r(i,j);r(i,k);r(j,k)]+HHBrIBiases)+HHBrLBiases);
                        Product33jk_1=Product22jk_1*purelin('dn',HHBrLWeights*logsig(HHBrIWeights*[r(i,j);r(i,k);r(j,k)]+HHBrIBiases)+HHBrLBiases);
                        
                        
%                         Product33ij_2=Product22ij_1*purelin('dn',HHBrLWeights*logsig(HHBrIWeights*[r(i,k);r(j,k);r(i,j)]+HHBrIBiases)+HHBrLBiases);
%                         Product33ik_2=Product22ik_1*purelin('dn',HHBrLWeights*logsig(HHBrIWeights*[r(i,k);r(j,k);r(i,j)]+HHBrIBiases)+HHBrLBiases);
%                         Product33jk_2=Product22jk_1*purelin('dn',HHBrLWeights*logsig(HHBrIWeights*[r(i,k);r(j,k);r(i,j)]+HHBrIBiases)+HHBrLBiases);
%                         
%                         
%                         Product33ij_3=Product22ij_1*purelin('dn',HHBrLWeights*logsig(HHBrIWeights*[r(j,k);r(i,j);r(i,k)]+HHBrIBiases)+HHBrLBiases);
%                         Product33ik_3=Product22ik_1*purelin('dn',HHBrLWeights*logsig(HHBrIWeights*[r(j,k);r(i,j);r(i,k)]+HHBrIBiases)+HHBrLBiases);
%                         Product33jk_3=Product22jk_1*purelin('dn',HHBrLWeights*logsig(HHBrIWeights*[r(j,k);r(i,j);r(i,k)]+HHBrIBiases)+HHBrLBiases);
%                         
                        
                        
                        Force(i)=Force(i)-Product33ij_1;
                        Force(j)=Force(j)+Product33ij_1;
                        
                        Force(i)=Force(i)-Product33ik_1;
                        Force(k)=Force(k)+Product33ik_1;
%                         
                        Force(j)=Force(j)-Product33jk_1;
                        Force(k)=Force(k)+Product33jk_1;
                                               
                        
%                         i;
                    otherwise
                        disp(identify3)
                end

            end
            
        end

    end
end

i=1; %#ok<NASGU> % dummy line




