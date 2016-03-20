

function [q]= calcweight1(xp,y)
M_N=1.0; 
for r = 1 : 500
                                          %% weight calculation %%
      N=(xp(:,r)^2)/20;
      q(1,r)=exp(-0.5*(M_N^(-1))*((y(r,1)-N)^2));
 end
 q(1,:)=q(1,:)./sum(q(1,:));            % weight normalization 
return;