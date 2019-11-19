function [losses,dlosses]=LossFunc(dt,rs,xs,ties)
%this function represents the loss function and its derivative

losses=rs.*(1./xs).^2.*(dt).^2;
if nargout==2
    dlosses=spdiags(2*rs.*(1./xs).^2.*(dt),0,size(ties,1),size(ties,1));
end
