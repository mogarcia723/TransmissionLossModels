function [losses,dlosses]=LossFunc(dt,V,approxtype,ys,ties,TR,theta)
%this function represents the loss function and its derivative

if approxtype==0
    gs=real(ys);
    losses=gs.*((V(ties(:,1))./TR).^2+V(ties(:,2)).^2)-2*gs.*V(ties(:,1))./TR.*V(ties(:,2)).*cos(dt-theta);
    if nargout==2
        dlosses=spdiags(2*gs .* (V(ties(:,1))./TR) .* V(ties(:,2)) .*sin(dt-theta),0,size(ties,1),size(ties,1));
    end
elseif approxtype==1
    gs=real(ys);
    losses=gs.*((V(ties(:,1))./TR).^2+V(ties(:,2)).^2)-2*gs.*V(ties(:,1))./TR.*V(ties(:,2)).*(1-.5*(dt-theta).^2);
    if nargout==2
        dlosses=spdiags(2*gs .* (V(ties(:,1))./TR) .* V(ties(:,2)) .*(dt-theta),0,size(ties,1),size(ties,1));
    end
elseif approxtype==2
    gs=real(ys);
    losses=gs.*(dt).^2;
    if nargout==2
        dlosses=spdiags(2*gs.*(dt),0,size(ties,1),size(ties,1));
    end
elseif approxtype==3
    %bs=imag(ys);
    rs=real(1./ys);
    zs=1./ys;
    xs=imag(zs);
    losses=rs.*(1./xs).^2.*(dt).^2;
    if nargout==2
        dlosses=spdiags(2*rs.*(1./xs).^2.*(dt),0,size(ties,1),size(ties,1));
    end
end