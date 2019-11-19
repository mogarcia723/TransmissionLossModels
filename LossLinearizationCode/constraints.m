function [fout,fouteq,dfout,dfouteq]=constraints(states)


load('MPCtemp.mat','spAdot','rs','xs','spA','n','D','M','tiesidx','genbusidx','m');
g=length(genbusidx);
thetadot=states((g+1):(g+n-1));
inj=-sparse(D);
for i=1:g
    inj(genbusidx(i))=inj(genbusidx(i))+states(i);
end

if nargout>2%enter this section if the gradient output is required
    [LossVal,dLossVal]=LossFunc(spAdot*thetadot,rs,xs,tiesidx);
    [MPFVal,dMPFVal]=MPFFunc(spAdot*thetadot,xs,tiesidx);
    
    
    fouteq=.5*abs(spA)'*LossVal*100+spA'*MPFVal*100-inj;
    dfouteq=-M';
    dfouteq=[dfouteq;100*spAdot'*(.5*dLossVal*abs(spA)+dMPFVal*spA)];
    
    fout=[];
    dfout=[];
    
else%enter this section of code if the gradient output is not required
    [LossVal]=LossFunc(spAdot*thetadot,rs,xs,tiesidx);
    [MPFVal]=MPFFunc(spAdot*thetadot,xs,tiesidx);
    
    fouteq=.5*abs(spA)'*LossVal*100+spA'*MPFVal*100-inj;
    fout=[];
end