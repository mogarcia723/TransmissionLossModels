function [cout,dcout]=cost(states)

%length of states is g+n-1

load('MPCtemp.mat','genbus','n','g','gencost');

if nargout>1
    cout=sum(gencost(:,1).*states(1:g).^2+gencost(:,2).*states(1:g)+gencost(:,3));
    dcout=[2*gencost(:,1).*states(1:g)+gencost(:,2);sparse(n-1,1)];
else 
    cout=sum(gencost(:,1).*states(1:g).^2+gencost(:,2).*states(1:g)+gencost(:,3));
    dcout=[];
end
    