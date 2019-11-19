function [dthetalb,dthetaub]=deltathetalimits(approxtype)
%this function computes the upper and lower limits on the voltage angle
%differences for each transmission line.  The limits represent squared
%current magnitude limits as stated in the paper.

load('MPCtemp.mat','mpc','ys','ysh','theta','n','m');

unlimited=find(mpc.branch(:,6)==0);%indices of lines with no specified limit
limited=find(mpc.branch(:,6)~=0);%indices of lines with specified limits
if approxtype==0 | approxtype==1
    psi=theta*pi/180;%psi from the paper in radians
    phi=angle(ysh./ys+1);%phi from the paper in radians
else
    psi=zeros(m,1);%psi from the paper in radians
    phi=zeros(m,1);%phi from the paper in radians
end

%% Finding DeltaTheta limits on the from side of the line
%First find upper bound on DeltaTheta
ISFrom=ISFuncFrom(psi(limited)-phi(limited)+pi/2,limited,approxtype);%Maximum possible squared current magnitude
alwayssat=find((100^2*ISFrom)<(mpc.branch(limited,6).^2));%check to see which limits are always satisfied
ISFrom=ISFuncFrom(psi(limited)-phi(limited),limited,approxtype);%Minimum possible squared current magnitude
neversat=find((100^2*ISFrom)>mpc.branch(limited,6).^2);%check to see which limits are never satisfied

%check to see that neversat is empty, otherwise the problem will be infeasible
alwaysidcs=[unlimited;limited(alwayssat)];%indices of lines that are always feasible
solveidcs=1:m;
solveidcs(alwaysidcs)=[];%indices of lines that we need to find bounds for

%solve for bounds by finding solving equation gfunc(DeltaTheta)=0
gfunc=@(dt) ISFuncFromSolve(dt,solveidcs,approxtype,1);%define gfunc. 1 refers to the positive direction of flow
options = optimoptions(@fsolve,'SpecifyObjectiveGradient',true);%set options for fsolve.  Note that gfunc outputs two things including derivative.
%options = optimoptions(@fsolve,'Display','iter','SpecifyObjectiveGradient',true);%set options for fsolve.  Note that gfunc outputs two things including derivative.
[x,Fub,exitflagub,outputub,JAC]=fsolve(gfunc,pi/2*ones(length(solveidcs),1),options);%solve the equations
%Save upper bounds
dthetafromub=psi-phi+ones(m,1)*pi/2;%conservative limits for unbounded lines
dthetafromub(solveidcs)=psi(solveidcs)-phi(solveidcs)+abs(x-psi(solveidcs)+phi(solveidcs));%limits found by solving equations

%Second find lower bound on DeltaTheta
ISFrom=ISFuncFrom(psi(limited)-phi(limited)+pi/2,limited,approxtype);%Maximum possible squared current magnitude
alwayssat=find((100^2*ISFrom)<(mpc.branch(limited,6).^2));%Check to see which limits are always satisfied
ISFrom=ISFuncFrom(psi(limited)-phi(limited),limited,approxtype);%Minimum possible squared current magnitude
neversat=find((100^2*ISFrom)>mpc.branch(limited,6).^2);%Check to see which limits are never satisfied

%check to see that neversat is empty, otherwise the problem is infeasible
alwaysidcs=[unlimited;limited(alwayssat)];%indices of lines that are always feasible
solveidcs=1:m;
solveidcs(alwaysidcs)=[];%Indices of lines that require us to find lower bounds

gfunc=@(dt) ISFuncFromSolve(dt,solveidcs,approxtype,-1);%define gfunc. -1 refers to the negative direction of flow
options = optimoptions(@fsolve,'SpecifyObjectiveGradient',true);
%options = optimoptions(@fsolve,'Display','iter','SpecifyObjectiveGradient',true);
[x,Flb,exitflaglb,output,JAC]=fsolve(gfunc,-pi/2*ones(length(solveidcs),1),options);
%Save lower bounds
dthetafromlb=psi-phi-ones(m,1)*pi/2;%conservative limits for unbounded lines
dthetafromlb(solveidcs)=psi(solveidcs)-phi(solveidcs)-abs(x-psi(solveidcs)+phi(solveidcs));%limits found by solving equations

%output data to see if everything is working properly
%norm(Fub,Inf) and norm(Flb,Inf)  should be small
%exitflags should be 1
[norm(Fub,Inf),exitflagub,max(dthetafromub),norm(Flb,Inf),exitflaglb,min(dthetafromlb)];

%% Finding DeltaTheta limits on the TO side of the line
%this section is nearly identical to the previous section
%except we use the function "ISFuncToSolve" instead of "ISFuncFromSolve"

%First find upper bound on DeltaTheta
ISTo=ISFuncTo(psi(limited)-phi(limited)+pi/2,limited,approxtype);
alwayssat=find((100^2*ISTo)<(mpc.branch(limited,6).^2));
ISTo=ISFuncTo(psi(limited)-phi(limited),limited,approxtype);
neversat=find((100^2*ISTo)>mpc.branch(limited,6).^2);
%check to see that neversat is empty

alwaysidcs=[unlimited;limited(alwayssat)];
solveidcs=1:m;
solveidcs(alwaysidcs)=[];

gfunc=@(dt) ISFuncToSolve(dt,solveidcs,approxtype,1);
options = optimoptions(@fsolve,'SpecifyObjectiveGradient',true);
%options = optimoptions(@fsolve,'Display','iter','SpecifyObjectiveGradient',true);
[x,Fub,exitflagub,output,JAC]=fsolve(gfunc,pi/2*ones(length(solveidcs),1),options);
dthetatoub=psi-phi+ones(m,1)*pi/2;
dthetatoub(solveidcs)=psi(solveidcs)-phi(solveidcs)+abs(x-psi(solveidcs)+phi(solveidcs));

%Second find lower bound on DeltaTheta
ISTo=ISFuncTo(psi(limited)-phi(limited)+pi/2,limited,approxtype);
alwayssat=find((100^2*ISTo)<(mpc.branch(limited,6).^2));
ISTo=ISFuncTo(psi(limited)-phi(limited),limited,approxtype);
neversat=find((100^2*ISTo)>mpc.branch(limited,6).^2);

%check to see that neversat is empty
alwaysidcs=[unlimited;limited(alwayssat)];
solveidcs=1:m;
solveidcs(alwaysidcs)=[];

gfunc=@(dt) ISFuncToSolve(dt,solveidcs,approxtype,-1);
options = optimoptions(@fsolve,'SpecifyObjectiveGradient',true);%,'StepTolerance',1e-18);
%options = optimoptions(@fsolve,'Display','iter','SpecifyObjectiveGradient',true);
[x,Flb,exitflaglb,output,JAC]=fsolve(gfunc,pi/2*ones(length(solveidcs),1),options);
dthetatolb=psi-phi-ones(m,1)*pi/2;
dthetatolb(solveidcs)=psi(solveidcs)-phi(solveidcs)-abs(x-psi(solveidcs)+phi(solveidcs));

%output info to see if everything is working
[norm(Fub,Inf),exitflagub,max(dthetafromub),norm(Flb,Inf),exitflaglb,min(dthetafromlb)];
%% Most restricting limits
%output the tighter of the two bounds found
dthetalb=max(dthetatolb,dthetafromlb);
dthetaub=min(dthetatoub,dthetafromub);