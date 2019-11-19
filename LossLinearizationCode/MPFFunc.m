function [MPF,dMPF]=MPFFunc(dt,xs,ties)
%this function represents the mid-line power flow function and its derivative

MPF=(1./xs).*(dt);
dMPF=spdiags(1./xs,0,size(ties,1),size(ties,1));
