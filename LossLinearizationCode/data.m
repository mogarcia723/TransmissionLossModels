function [TR,Y,injbus,n,m,nullbus,ties,vlim,ys,ysh,yshb,theta,tiesidx,genbus,genbusidx,gencost,A]= data(mpc)

j=sqrt(-1);
rs=mpc.branch(:,3);
xs=mpc.branch(:,4);
bc=mpc.branch(:,5)./2;
t=mpc.branch(:,9);
t(find(t==0))=1;
theta=mpc.branch(:,10);
ties=mpc.branch(:,[1,2]); % mx2 matrix.  Each row is a line. Shows shich busses each line is connected to 
ys=1./(rs+j*xs);


n=size(mpc.bus,1);%number of busses
m=size(mpc.branch,1);%number of lines

tiesidx=zeros(size(ties));
for ell=1:size(ties,1)
    tiesidx(ell,1)=find(ties(ell,1)==mpc.bus(:,1));
    tiesidx(ell,2)=find(ties(ell,2)==mpc.bus(:,1));
end
genbus=mpc.gen(:,1);
genbusidx=zeros(size(genbus));
for g=1:length(mpc.gen(:,1))
    genbusidx(g)=find(genbus(g)==mpc.bus(:,1));
end

gencost=mpc.gencost(:,5);

%% Formulating the admittance matrix Y=G+j*B.  
% reference page 17 of the matpower manual
Y=sparse(n,n);

busnum=mpc.bus(:,1);% list of bus numbers (listed 1-9 in 9 bus case)

for i=1:m
    Y(find(busnum==ties(i,1)),find(busnum==ties(i,1)))=Y(find(busnum==ties(i,1)),find(busnum==ties(i,1)))+(ys(i)+j*bc(i)/2)/t(i)^2;
    Y(find(busnum==ties(i,1)),find(busnum==ties(i,2)))=Y(find(busnum==ties(i,1)),find(busnum==ties(i,2)))-ys(i)/(t(i)*exp(-j*theta(i)));
    Y(find(busnum==ties(i,2)),find(busnum==ties(i,1)))=Y(find(busnum==ties(i,2)),find(busnum==ties(i,1)))-ys(i)/(t(i)*exp(j*theta(i)));
    Y(find(busnum==ties(i,2)),find(busnum==ties(i,2)))=Y(find(busnum==ties(i,2)),find(busnum==ties(i,2)))+(ys(i)+j*bc(i)/2);
end
%% Finding the series admittance in each branch ys=gs+j*bs
ys=1./(rs+j*xs);
%% finding shunt admittance ysh=gsh+j*bsh
ysh=j*bc;
%% finding bus shunt admittance yshb=gshb+j*bshb
yshb=mpc.bus(:,5)/100+j*mpc.bus(:,6)/100;
%% finding the voltage limits for each bus
vlim=mpc.bus(:,[13,12]);

TR=t*exp(j*theta(i));
%% lnode
% all nodes with a real or reactive load and the slack bus
ind=zeros(size(mpc.bus,1),1);
for i=1:size(mpc.bus,1)
    ind(i)=mpc.bus(i,3)||mpc.bus(i,4)||(mpc.bus(i,2)==3);
end
lnode=busnum(find(ind));


%% gnode
gnode=mpc.gen(:,1);

%% injbus
%finding the injection buses 
for i=1:(n)
    if isempty(find(busnum(i)==lnode))    &&    isempty(find(busnum(i)==gnode))
        isinj(i)=0;
    else
        isinj(i)=1;
    end
end
injbus=busnum(find(isinj));
nullbus=busnum(find(~isinj));


%create the branch bus incidence matrix
A=sparse(m,n);
for ell=1:size(ties,1)
    A(ell,tiesidx(ell,1))=1;
    A(ell,tiesidx(ell,2))=-1;
end

















