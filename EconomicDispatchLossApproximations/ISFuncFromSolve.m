function [fout,dfout]=ISFuncFromSolve(dt,idcs,approxtype,solsign)

load('MPCtemp.mat','mpc','Vopt','ys','ysh','tiesidx','genbus','TR','theta','n','D','dgmat','m','ties');


if approxtype==0 | approxtype==1
    psi=theta(idcs)*pi/180;
    phi=angle(ysh(idcs)./ys(idcs)+1);
else
    psi=0;
    phi=0;
end


sgns=sign(dt-psi+phi);
sgns(find(sgns==0))=1;

mtemp=length(idcs);


if solsign==1
    if nargout==1
        [ISFrom]=ISFuncFrom(abs(dt-psi+phi)+psi-phi,idcs,approxtype);
        fout=100^2*ISFrom-mpc.branch(idcs,6).^2;   
    else
        [ISFrom,dISFrom]=ISFuncFrom(abs(dt-psi+phi)+psi-phi,idcs,approxtype);
        fout=100^2*ISFrom-mpc.branch(idcs,6).^2;  
        dfout=100^2*repmat(sgns,1,mtemp).*dISFrom;
    end
else
    if nargout==1
        [ISFrom]=ISFuncFrom(-abs(dt-psi+phi)+psi-phi,idcs,approxtype);
        fout=100^2*ISFrom-mpc.branch(idcs,6).^2;   
    else
        [ISFrom,dISFrom]=ISFuncFrom(-abs(dt-psi+phi)+psi-phi,idcs,approxtype);
        fout=100^2*ISFrom-mpc.branch(idcs,6).^2;  
        dfout=-100^2*repmat(sgns,1,mtemp).*dISFrom;
    end
end
