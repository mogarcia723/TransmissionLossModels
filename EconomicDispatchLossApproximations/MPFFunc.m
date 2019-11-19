function [MPF,dMPF]=MPFFunc(dt,V,approxtype,ys,ties,TR,theta)
%this function represents the mid-line power flow function and its derivative

if nargin==6
   theta=zeros(size(TR)); 
end

if approxtype==0
    bs=imag(ys);
    gs=real(ys);
    MPF=.5*gs.*((V(ties(:,1))./TR).^2-V(ties(:,2)).^2)-bs.*V(ties(:,1))./TR.*V(ties(:,2)).*sin(dt-theta);
    if nargout==2
        dMPF=spdiags(-bs.*V(ties(:,1))./TR.*V(ties(:,2)).*cos(dt-theta),0,size(ties,1),size(ties,1));
    end
elseif approxtype==1
    bs=imag(ys);
    gs=real(ys);
    MPF=.5*gs.*((V(ties(:,1))./TR).^2-V(ties(:,2)).^2)-bs.*V(ties(:,1))./TR.*V(ties(:,2)).*(dt-theta);
    dMPF=spdiags(-bs.*V(ties(:,1))./TR.*V(ties(:,2)),0,size(ties,1),size(ties,1));
elseif approxtype==2
    bs=imag(ys);
    MPF=-bs.*(dt);
    dMPF=spdiags(-bs,0,size(ties,1),size(ties,1));
elseif approxtype==3
    zs=1./ys;
    xs=imag(zs);
    MPF=1./xs.*(dt);
    dMPF=spdiags(1./xs,0,size(ties,1),size(ties,1));
end