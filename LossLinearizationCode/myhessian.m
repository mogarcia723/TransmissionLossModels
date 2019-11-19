function [Lxx]=myhessian(states,lam)
%this function outputs the hessian of the lagrangian for a given
%approximation.  

mu=lam.eqnonlin;

load('MPCtemp.mat','rs','xs','spA','n','m','D','spAdot','tiesidx','genbusidx','gencost','g');

dF=spdiags(2*rs.*(1./xs).^2 .*(abs(spA)*mu),0,size(tiesidx,1),size(tiesidx,1));

Lxx=[spdiags(2*gencost(:,1),0,g,g),sparse(g,n-1);sparse(n-1,g),.5*100*spAdot'*dF*spAdot];