function [cout,dcout]=cost(states)

%length of states is g+n-1

load('MPCtemp.mat','mpc','genbus');
n=size(mpc.bus,1);
g=length(genbus);

cout=sum(mpc.gencost(:,5).*states(1:g).^2+mpc.gencost(:,6).*states(1:g)+mpc.gencost(:,7));
dcout=[2*mpc.gencost(:,5).*states(1:g)+mpc.gencost(:,6);sparse(n-1,1)];