function [fout,fouteq,dfout,dfouteq]=constraints(states,eta,rlx,approxtype)


%length of states is g+n-1


load('MPCtemp.mat','spAdot','Vopt','ys','spA','TR','theta','n','D','dgmat','tiesidx','genbusidx','m');
g=length(genbusidx);
thetadot=states((g+1):(g+n-1));
inj=-sparse(D);
for i=1:g
    inj(genbusidx(i))=inj(genbusidx(i))+states(i);
end

if nargout>2%enter this section if the gradient output is required
    [LossVal,dLossVal]=LossFunc(spAdot*thetadot,Vopt,approxtype,ys,tiesidx,TR,theta*pi/180);
    [MPFVal,dMPFVal]=MPFFunc(spAdot*thetadot,Vopt,approxtype,ys,tiesidx,TR,theta*pi/180);
    
    if isempty(eta)
        fout=.5*abs(spA)'*LossVal*100+spA'*MPFVal*100-inj;
        dfout=dgmat;
        dfout=[dfout;100*(.5*spAdot'*dLossVal*abs(spA)+spAdot'*dMPFVal*spA)];
    else
        fout=eta*sum(LossVal)*100+spA'*MPFVal*100-inj;
        dfout=dgmat;
        dfout=[dfout;100*(spAdot'*sum(dLossVal,2)*eta'+spAdot'*dMPFVal*spA)];
    end
    
    fouteq=[];
    dfouteq=[];
    
    if rlx==0
        fouteq=fout;
        dfouteq=dfout;
        fout=[];
        dfout=[];
    end
else%enter this section of code if the gradient output is not required
    [LossVal]=LossFunc(spAdot*thetadot,Vopt,approxtype,ys,tiesidx,TR,theta*pi/180);
    [MPFVal]=MPFFunc(spAdot*thetadot,Vopt,approxtype,ys,tiesidx,TR,theta*pi/180);
    if isempty(eta)
        fout=.5*abs(spA)'*LossVal*100+spA'*MPFVal*100-inj;
    else
        fout=eta*sum(LossVal)*100+spA'*MPFVal*100-inj;
    end
    
    fouteq=[];
    
    if rlx==0
        fouteq=fout;
        fout=[];
    end
end