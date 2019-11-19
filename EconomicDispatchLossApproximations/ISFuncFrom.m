function [ISFrom,dISFrom]=ISFuncFrom(dt,idcs,approxtype)

load('MPCtemp.mat','mpc','Vopt','ys','ysh','tiesidx','genbus','TR','theta','n','D','dgmat','m','ties');
V=Vopt;
alpha=abs(ysh./ys+1);
phi=angle(ysh./ys+1);

if approxtype==0
    ISFrom=abs(ys(idcs)).^2./TR(idcs).^2.*((alpha(idcs).*V(tiesidx(idcs,1))./TR(idcs)).^2+V(tiesidx(idcs,2)).^2-2*alpha(idcs).*V(tiesidx(idcs,1))./TR(idcs).*V(tiesidx(idcs,2)).*cos(phi(idcs)+dt-theta(idcs)*pi/180));  
    if nargout==2
        dISFrom=spdiags(2*alpha(idcs).*abs(ys(idcs)).^2./TR(idcs).^2.*V(tiesidx(idcs,1))./TR(idcs).*V(tiesidx(idcs,2)).*sin(phi(idcs)+dt-theta(idcs)*pi/180),0,length(idcs),length(idcs));
    end
elseif approxtype==1
    ISFrom=abs(ys(idcs)).^2./TR(idcs).^2.*((alpha(idcs).*V(tiesidx(idcs,1))./TR(idcs)).^2+V(tiesidx(idcs,2)).^2-2*alpha(idcs).*V(tiesidx(idcs,1))./TR(idcs).*V(tiesidx(idcs,2)).*(1-.5*(phi(idcs)+dt-theta(idcs)*pi/180).^2));  
    if nargout==2
        dISFrom=spdiags(2*alpha(idcs).*abs(ys(idcs)).^2./TR(idcs).^2.*V(tiesidx(idcs,1))./TR(idcs).*V(tiesidx(idcs,2)).*(phi(idcs)+dt-theta(idcs)*pi/180),0,length(idcs),length(idcs));
    end
% elseif approxtype==2
%     ISFrom=abs(ys(idcs)).^2.*((alpha(idcs)).^2+1-2*alpha(idcs).*(1-.5*(phi(idcs)+dt).^2));  
%     if nargout==2
%         dISFrom=spdiags(2*alpha(idcs).*abs(ys(idcs)).^2.*(phi(idcs)+dt),0,length(idcs),length(idcs));
%     end
elseif approxtype==2
    ISFrom=abs(ys(idcs)).^2.*(dt).^2;
    if nargout==2
        dISFrom=spdiags(2*abs(ys(idcs)).^2.*dt,0,length(idcs),length(idcs));
    end
elseif approxtype==3
    zs=1./ys;
    xs=imag(zs);
    ISFrom=(1./xs(idcs).*dt).^2; 
    if nargout==2
        dISFrom=spdiags(2./xs(idcs).^2.*dt,0,length(idcs),length(idcs));
    end
end