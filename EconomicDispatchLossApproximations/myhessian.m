function [Lxx]=myhessian(states,lam,eta,rlx,approxtype)
%this function outputs the hessian of the lagrangian for a given
%approximation.  If eta is non-empty, then eta represents the loss
%distribution factor vector.  

if rlx==1
    mu=lam.ineqnonlin;
else
    mu=lam.eqnonlin;
end

load('MPCtemp.mat','mpc','Vopt','ys','spA','TR','theta','n','m','D','dgmat','spAdot','tiesidx','genbusidx','gencost');
gs=real(ys);
bs=imag(ys);
rs=real(1./ys);
g=length(genbusidx);
thetadot=states((g+1):(g+n-1));
Pg=states(1:g);


if isempty(eta)
    if approxtype==0
        dF=spdiags(2*gs .* (Vopt(tiesidx(:,1))./TR) .* Vopt(tiesidx(:,2)) .*cos(spAdot*thetadot-theta*pi/180).*(abs(spA)*mu),0,size(tiesidx,1),size(tiesidx,1));
        dH=spdiags(bs.*Vopt(tiesidx(:,1))./TR.*Vopt(tiesidx(:,2)).*sin(spAdot*thetadot-theta*pi/180).*(spA*mu),0,size(tiesidx,1),size(tiesidx,1));
    elseif approxtype==1
        dF=spdiags(2*gs .* (Vopt(tiesidx(:,1))./TR) .* Vopt(tiesidx(:,2)) .*(abs(spA)*mu),0,size(tiesidx,1),size(tiesidx,1));
        dH=sparse(size(tiesidx,1),size(tiesidx,1));
    elseif approxtype==2
        dF=spdiags(2*gs .*(abs(spA)*mu),0,size(tiesidx,1),size(tiesidx,1));
        dH=sparse(size(tiesidx,1),size(tiesidx,1));
    elseif approxtype==3
        dF=spdiags(2*rs.*bs.^2 .*(abs(spA)*mu),0,size(tiesidx,1),size(tiesidx,1));
        dH=sparse(size(tiesidx,1),size(tiesidx,1));
    end
    
    Lxx=[spdiags(2*gencost,0,g,g),sparse(g,n-1);sparse(n-1,g),100*spAdot'*(.5*dF+dH)*spAdot];%+.005*speye(n-1+g);
else
    if approxtype==0
        dF=spdiags(2*gs .* (Vopt(tiesidx(:,1))./TR) .* Vopt(tiesidx(:,2)) .*cos(spAdot*thetadot-theta*pi/180).*(2*ones(m,1)*eta'*mu),0,size(tiesidx,1),size(tiesidx,1));
        dH=spdiags(bs.*Vopt(tiesidx(:,1))./TR.*Vopt(tiesidx(:,2)).*sin(spAdot*thetadot-theta*pi/180).*(spA*mu),0,size(tiesidx,1),size(tiesidx,1));
    elseif approxtype==1
        
    elseif approxtype==2
        
    elseif approxtype==3
        %dF=spdiags(2*rs.*bs.^2 .*(2*ones(m,1)*eta'*mu),0,size(tiesidx,1),size(tiesidx,1));
        dF=spdiags(4*rs.*bs.^2 *(eta'*mu),0,size(tiesidx,1),size(tiesidx,1));
        dH=sparse(size(tiesidx,1),size(tiesidx,1));
    end
    
    Lxx=[spdiags(2*gencost,0,g,g),sparse(g,n-1);sparse(n-1,g),100*spAdot'*(.5*dF+dH)*spAdot];
    %Lxx=[spdiags(2*gencost,0,g,g),sparse(g,n-1);sparse(n-1,g),50*spAdot'*dF*spAdot];
end
    




