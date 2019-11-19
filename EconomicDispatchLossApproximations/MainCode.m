%% Description
%This code produces the results used in Table I of the ACC paper entitled
%"A General Economic Dispatch Problem with Marginal Losses"
%
%The table compares the "Exact" Economic Dispatch Problem to Approximations
%1-4 as well as the AC OPF problem. Approximations 1-4 are explicitly
%outlined in column 3 of Table I in the paper.
%
%Transmission line limits are enforced as current magnitude limits as
%explained in Section III-C2 in the paper
%
%the Matpower toolbox is used to solve the AC OPF problem and the test
%cases are taken from the NESTA archive of test cases.
%
%The standard MATLAB function fmincon() is used to solve the economic
%dispatch problems.  To improve convergence speed analytical hessians of
%the Lagrangian function are provided.

%% Set MATLAB path
%The MATLAB directory should be set to the main folder "EconomicDispatchLossApproximations"
clc
clear
restoredefaultpath
%add matpower path
addpath(genpath('./matpower6.0'))
%add nesta path
addpath(genpath('./nesta/opf'))
%% Load Test Case Data and Run Power Flow and Run AC OPF
%Choose the test case you want to use, e.g. case300.
%Simple Test Cases: 300 118 9 30 39 57 118   
%Large Complicated test cases: 3375wp 6515rte 2869pegase 9241pegase 2746wop 3120sp 2737sop case2868rte
mpc = loadcase(case2869pegase);%set original mpc case

%Do you want to use the load oversatisfaction relaxation? (only works if all prices are positive)
rlx=0;%load oversatisfaction relaxation 1=use relaxation 0= dont use relaxation


busnums=mpc.bus(:,1);
pfr=runpf(mpc);%find power flow results
opt=mpoption('OPF_VIOLATION',5*10^(-6));%set options for AC OPF
opt=mpoption(opt, 'opf.ac.solver','FMINCON');%use fmincon (Note: other solvers may be faster, e.g. Knitro)
opt=mpoption(opt, 'opf.flow_lim','I','opf.init_from_mpc',1);%set options for AC OPF including current limit interpretation
ACOPFsol=runopf(pfr,opt);%Run AC OPF problem

%% Creating temporary data files
%the function ``data'' outputs parameters associated with the test case
%from mpc.  ie. the admittance matrix Y and the transformer ratios TR
[TR,Y,injbus,n,m,nullbus,ties,vlim,ys,ysh,yshb,theta,tiesidx,genbus,genbusidx,gencost,A]= data(mpc,pfr);

Adot=A(:,2:end);%branch bus incidence matrix with reference bus row removed
spAdot=sparse(Adot);%create sparse version of Adot
spA=sparse(A);%create sparse version of A
D=sparse(mpc.bus(:,3));%demand vector

%dgmat is a matrix used to compute the gradient in the file constraints.m
dgmat=[];
g=size(mpc.gen,1);
for i = 1:g
    dgmat=[dgmat;sparse([sparse(1,genbusidx(i)-1),-1,sparse(1,n-genbusidx(i))])];
end

% Creating Cost Functions:  I assume mpc.gencost(:,1) is equal to 2 so cost functions are polynomial
order=mpc.gencost(:,4);%order of polynomial
g=size(mpc.gen,1);%number of generators
genlims=[mpc.gen(:,10),mpc.gen(:,9)];%generator limits

%voltage magnitudes are fixed for all economic dispatch problems considered (the exact problem and approximations 1-4)
Vopt=ACOPFsol.bus(:,8);%fix voltages to AC OPF solution
%Vopt=pfr.bus(:,8);%Another option: fix voltages to power flow solution
%Vopt=ones(size(ACOPFsol.bus(:,8))); %Another option: fix voltages to nominal values

%MPCtemp.mat contains data related to the test case that will be accessed
%in other functions (e.g. "LossFunc.m" and "myhessian.m")
save('MPCtemp.mat','mpc','Adot','Vopt','ys','ysh','ties','genbus','A','TR','theta','n','D','dgmat','spAdot','spA','g','tiesidx','genbusidx','gencost','m');
%% Results pertaining to the AC OPF problem
resAC.LMPs=ACOPFsol.bus(:,14);%Save LMPs for the AC OPF problem 
resAC.dispatch=100*ACOPFsol.var.val.Pg;%generation dispatch
resAC.time=ACOPFsol.et;%cpu time
resAC.f=ACOPFsol.f;%optimal value
resAC.pnorm0=length(find(resAC.dispatch>.0001));%number of dispatched generators
resAC.psum=sum(resAC.dispatch);%total generation
resAC.meanlmp=mean(resAC.LMPs);%average LMP
resAC.lmprange=[min(resAC.LMPs),max(resAC.LMPs)];%range of LMPs
resAC.dtinfnorm=norm(spAdot*(ACOPFsol.var.val.Va(2:end)-ACOPFsol.var.val.Va(1)),Inf);%max angle difference

%% Setting Economic Dispatch Problem Parameters that are the same for each approximation
x0=[ACOPFsol.gen(:,2);(ACOPFsol.bus(2:end,9)-ACOPFsol.bus(1,9))*pi/180];%Setting initial guess for interior point algorithm
%x0=[pfr.gen(:,2);(pfr.bus(2:end,9)-pfr.bus(1,9))*pi/180];%another option in setting initial guess
lb=[genlims(:,1);-Inf*ones(n-1,1)];%lower bound on decision variable x
ub=[genlims(:,2);Inf*ones(n-1,1)];%upper bound on decision variable x

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Exact Economic Dispatch Problem %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
approxtype=0;% The exact economic dispatch problem is designated approximation 0
%% Line Limits
%Current magnitude limits are used as explained in section III-C2 of the paper
[dthetalb,dthetaub]=deltathetalimits(approxtype);%finding line limits using approximation type 0

%% Solving Optimization Problem
eta=[];
myhessiantemp=@(state,lam) myhessian(state,lam,eta,rlx,approxtype);
constraintstemp=@(state) constraints(state,eta,rlx,approxtype);
options = optimoptions('fmincon','Algorithm','interior-point',...
    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
    'HessianFcn',myhessiantemp,'Display','iter','ConstraintTolerance', 1e-06);%,'CheckGradients',true,'FiniteDifferenceType','central');
tic
[xEXACT,fvalEXACT,exitflagEXACT,outputEXACT,lambdaEXACT,gradEXACT,hessianEXACT]=fmincon(@cost,x0,[sparse(m,g),spAdot;sparse(m,g),-spAdot],sparse([dthetaub;-dthetalb]),[],[],sparse(lb),sparse(ub),constraintstemp,options);
% Collecting results pertaining to the Exact Economic Dispatch problem
resEXACT.time=toc;%cpu time
resEXACT.f=fvalEXACT;
resEXACT.dispatch=xEXACT(1:g);
resEXACT.pnorm0=length(find(resEXACT.dispatch>.0001));
resEXACT.psum=sum(resEXACT.dispatch);
if rlx==1
    resEXACT.LMPs=lambdaEXACT.ineqnonlin;
else
    resEXACT.LMPs=lambdaEXACT.eqnonlin;
end
resEXACT.meanlmp=mean(resEXACT.LMPs);
resEXACT.lmprange=[min(resEXACT.LMPs),max(resEXACT.LMPs)];

resEXACT.dtinfnorm=norm(spAdot*xEXACT((g+1):end),Inf);

%% More results pertaining to the AC OPF problem
resAC.pnorm1=norm(resAC.dispatch-resEXACT.dispatch,1);
resAC.pnorminf=norm(resAC.dispatch-resEXACT.dispatch,Inf);
resAC.fdiff=resEXACT.f-resAC.f;
resAC.lmpnormdiff=max(abs(resAC.LMPs-resEXACT.LMPs)./resEXACT.LMPs);
resAC.lmpmeannormdiff=mean(abs(resAC.LMPs-resEXACT.LMPs)./resEXACT.LMPs);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Approx 1 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
approxtype=1;
%% Line Limits
[dthetalb1,dthetaub1]=deltathetalimits(approxtype);
%% Approx 1
eta=[];
myhessiantemp=@(state,lam) myhessian(state,lam,eta,rlx,approxtype);
constraintstemp=@(state) constraints(state,eta,rlx,approxtype);
options = optimoptions('fmincon','Algorithm','interior-point',...
    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
    'HessianFcn',myhessiantemp,'Display','iter','ConstraintTolerance', 1e-06);
tic;
[x1,fval1,exitflag1,output1,lambda1,grad1,hessian1]=fmincon(@cost,x0,[sparse(m,g),spAdot;sparse(m,g),-spAdot],[dthetaub1;-dthetalb1],[],[],lb,ub,constraintstemp,options);

% Collecting results pertaining to the Approx 1
res1.time=toc;%cpu time
res1.fdiff=resEXACT.f-fval1;
res1.pnorm1=norm(resEXACT.dispatch-x1(1:g),1);
res1.pnorminf=norm(resEXACT.dispatch-x1(1:g),Inf);
res1.pnorm0=length(find(x1(1:g)>.0001));
res1.psum=sum(x1(1:g));
if rlx==1
    res1.meanlmp=mean(lambda1.ineqnonlin);
    res1.lmprange=[min(lambda1.ineqnonlin),max(lambda1.ineqnonlin)];
    res1.lmpnormdiff=max(abs(resEXACT.LMPs-lambda1.ineqnonlin)./resEXACT.LMPs);
    res1.lmpmeannormdiff=mean(abs(resEXACT.LMPs-lambda1.ineqnonlin)./resEXACT.LMPs);
else
    res1.meanlmp=mean(lambda1.eqnonlin);
    res1.lmprange=[min(lambda1.eqnonlin),max(lambda1.eqnonlin)];
    res1.lmpnormdiff=max(abs(resEXACT.LMPs-lambda1.eqnonlin)./resEXACT.LMPs);
    res1.lmpmeannormdiff=mean(abs(resEXACT.LMPs-lambda1.eqnonlin)./resEXACT.LMPs);
end
res1.dtinfnorm=norm(spAdot*x1((g+1):end),Inf);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Approx 2 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
approxtype=2;
%% Line Limits
[dthetalb2,dthetaub2]=deltathetalimits(2);
%% Approx 2
eta=[];
myhessiantemp=@(state,lam) myhessian(state,lam,eta,rlx,approxtype);
constraintstemp=@(state) constraints(state,eta,rlx,approxtype);
options = optimoptions('fmincon','Algorithm','interior-point',...
    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
    'HessianFcn',myhessiantemp,'Display','iter','ConstraintTolerance', 1e-06,'MaxIterations',100);%,'CheckGradients',true,'FiniteDifferenceType','central');
tic
[x2,fval2,exitflag2,output2,lambda2,grad2,hessian2]=fmincon(@cost,x0,[sparse(m,g),spAdot;sparse(m,g),-spAdot],[dthetaub2;-dthetalb2],[],[],lb,ub,constraintstemp,options);

% Collecting results pertaining to the Approx 2
res2.time=toc;%cpu time
res2.fdiff=resEXACT.f-fval2;
res2.pnorm1=norm(resEXACT.dispatch-x2(1:g),1);
res2.pnorminf=norm(resEXACT.dispatch-x2(1:g),Inf);
res2.pnorm0=length(find(x2(1:g)>.0001));
res2.psum=sum(x2(1:g));
if rlx==1
    res2.meanlmp=mean(lambda2.ineqnonlin);
    res2.lmprange=[min(lambda2.ineqnonlin),max(lambda2.ineqnonlin)];
    res2.lmpnormdiff=max(abs(resEXACT.LMPs-lambda2.ineqnonlin)./resEXACT.LMPs);
    res2.lmpmeannormdiff=mean(abs(resEXACT.LMPs-lambda2.ineqnonlin)./resEXACT.LMPs);
else
    res2.meanlmp=mean(lambda2.eqnonlin);
    res2.lmprange=[min(lambda2.eqnonlin),max(lambda2.eqnonlin)];
    res2.lmpnormdiff=max(abs(resEXACT.LMPs-lambda2.eqnonlin)./resEXACT.LMPs);
    res2.lmpmeannormdiff=mean(abs(resEXACT.LMPs-lambda2.eqnonlin)./resEXACT.LMPs);
end
res2.dtinfnorm=norm(spAdot*x2((g+1):end),Inf);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Approx 3 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
approxtype=3;
clc
%% Line Limits
[dthetalb3,dthetaub3]=deltathetalimits(approxtype);
%% Approx 3
eta=[];
myhessiantemp=@(state,lam) myhessian(state,lam,eta,rlx,approxtype);
constraintstemp=@(state) constraints(state,eta,rlx,approxtype);
options = optimoptions('fmincon','Algorithm','interior-point',...
    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
    'HessianFcn',myhessiantemp,'Display','iter','ConstraintTolerance', 1e-06);%,'CheckGradients',true,'FiniteDifferenceType','central');
tic
[x3,fval3,exitflag3,output3,lambda3,grad3,hessian3]=fmincon(@cost,x0,[sparse(m,g),spAdot;sparse(m,g),-spAdot],[dthetaub3;-dthetalb3],[],[],lb,ub,constraintstemp,options);
% Collecting results pertaining to the Approx 3
res3.time=toc;%cpu time
res3.fdiff=resEXACT.f-fval3;
res3.pnorm1=norm(resEXACT.dispatch-x3(1:g),1);
res3.pnorminf=norm(resEXACT.dispatch-x3(1:g),Inf);
res3.pnorm0=length(find(x3(1:g)>.0001));
res3.psum=sum(x3(1:g));
if rlx==1
    res3.meanlmp=mean(lambda3.ineqnonlin);
    res3.lmprange=[min(lambda3.ineqnonlin),max(lambda3.ineqnonlin)];
    res3.lmpnormdiff=max(abs(resEXACT.LMPs-lambda3.ineqnonlin)./resEXACT.LMPs);
    res3.lmpmeannormdiff=mean(abs(resEXACT.LMPs-lambda3.ineqnonlin)./resEXACT.LMPs);
else
    res3.meanlmp=mean(lambda3.eqnonlin);
    res3.lmprange=[min(lambda3.eqnonlin),max(lambda3.eqnonlin)];
    res3.lmpnormdiff=max(abs(resEXACT.LMPs-lambda3.eqnonlin)./resEXACT.LMPs);
    res3.lmpmeannormdiff=mean(abs(resEXACT.LMPs-lambda3.eqnonlin)./resEXACT.LMPs);
end
res3.dtinfnorm=norm(spAdot*x3((g+1):end),Inf);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Approx 4 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Note: Approximation 4 in the paper introduces LDFs (denoted eta) to Approximation 3
approxtype=3;
clc
%%
eta=sparse(n,1);%loss distribution factors
eta(find(mpc.bus(:,2)==3))=1;%assign load allocation to slack bus
myhessiantemp=@(state,lam) myhessian(state,lam,eta,rlx,approxtype);
constraintstemp=@(state) constraints(state,eta,rlx,approxtype);
options = optimoptions('fmincon','Algorithm','interior-point',...
    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
    'HessianFcn',myhessiantemp,'Display','iter','ConstraintTolerance', 1e-06);%,'CheckGradients',true,'FiniteDifferenceType','central');
tic
[x3LDF,fval3LDF,exitflag3LDF,output3LDF,lambda3LDF,grad3LDF,hessian3LDF]=fmincon(@cost,x0,[sparse(m,g),spAdot;sparse(m,g),-spAdot],[dthetaub3;-dthetalb3],[],[],lb,ub,constraintstemp,options);
% Collecting results pertaining to the Approx 4 (aka Approx 3LDF)
res3LDF.time=toc;%cpu time
res3LDF.fdiff=resEXACT.f-fval3LDF;
res3LDF.pnorm1=norm(resEXACT.dispatch-x3LDF(1:g),1);
res3LDF.pnorminf=norm(resEXACT.dispatch-x3LDF(1:g),Inf);
res3LDF.pnorm0=length(find(x3LDF(1:g)>.0001));
res3LDF.psum=sum(x3LDF(1:g));
if rlx==1
    res3LDF.meanlmp=mean(lambda3LDF.ineqnonlin);
    res3LDF.lmprange=[min(lambda3LDF.ineqnonlin),max(lambda3LDF.ineqnonlin)];
    res3LDF.lmpnormdiff=max(abs(resEXACT.LMPs-lambda3LDF.ineqnonlin)./resEXACT.LMPs);
    res3LDF.lmpmeannormdiff=mean(abs(resEXACT.LMPs-lambda3LDF.ineqnonlin)./resEXACT.LMPs);
else
    res3LDF.meanlmp=mean(lambda3LDF.eqnonlin);
    res3LDF.lmprange=[min(lambda3LDF.eqnonlin),max(lambda3LDF.eqnonlin)];
    res3LDF.lmpnormdiff=max(abs(resEXACT.LMPs-lambda3LDF.eqnonlin)./resEXACT.LMPs);
    res3LDF.lmpmeannormdiff=mean(abs(resEXACT.LMPs-lambda3LDF.eqnonlin)./resEXACT.LMPs);
end
res3LDF.dtinfnorm=norm(spAdot*x3LDF((g+1):end),Inf);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Printing Results  %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
%% Results that match Table I in the paper
Approx={'Exact';'Approx 1';'Approx 2';'Approx 3';'Approx 4';'AC OPF'};
OptValDiff=[0;res1.fdiff;res2.fdiff;res3.fdiff;res3LDF.fdiff;resAC.fdiff];
TotGen=[resEXACT.psum;res1.psum;res2.psum;res3.psum;res3LDF.psum;resAC.psum];
Norm0=[resEXACT.pnorm0;res1.pnorm0;res2.pnorm0;res3.pnorm0;res3LDF.pnorm0;resAC.pnorm0];
Norm1=[0;res1.pnorm1;res2.pnorm1;res3.pnorm1;res3LDF.pnorm1;resAC.pnorm1];
NormInf=[0;res1.pnorminf;res2.pnorminf;res3.pnorminf;res3LDF.pnorminf;resAC.pnorminf];
MeanLMP=[resEXACT.meanlmp;res1.meanlmp;res2.meanlmp;res3.meanlmp;res3LDF.meanlmp;resAC.meanlmp];
MinLMP=[resEXACT.lmprange(1);res1.lmprange(1);res2.lmprange(1);res3.lmprange(1);res3LDF.lmprange(1);resAC.lmprange(1)];
MaxLMP=[resEXACT.lmprange(2);res1.lmprange(2);res2.lmprange(2);res3.lmprange(2);res3LDF.lmprange(2);resAC.lmprange(2)];
DiffLMP=[0;res1.lmpnormdiff;res2.lmpnormdiff;res3.lmpnormdiff;res3LDF.lmpnormdiff;resAC.lmpnormdiff];
Dtheta=[resEXACT.dtinfnorm;res1.dtinfnorm;res2.dtinfnorm;res3.dtinfnorm;res3LDF.dtinfnorm;resAC.dtinfnorm];
time=[resEXACT.time;res1.time;res2.time;res3.time;res3LDF.time;resAC.time];

table(Approx,OptValDiff,TotGen,Norm0,Norm1,NormInf,MeanLMP,MinLMP,MaxLMP,DiffLMP,Dtheta,time)

%% Latex code to insert in table (remove semicolon to print)

latexwriting=['& Exact &(\ref{LossFuncDef}) / (\ref{midlinefloweq}) / (\ref{ExactCurrent}) / (\ref{PIdef})& 0 & ',sprintf('%.0f',resEXACT.psum),' & ',sprintf('%.0f',resEXACT.pnorm0),' & 0 & 0 & ',sprintf('%.2f',resEXACT.meanlmp),' & ',sprintf('%.2f',resEXACT.lmprange(1)),' / ',sprintf('%.2f',resEXACT.lmprange(2)),' & 0 & ',sprintf('%.4f',resEXACT.dtinfnorm),' & ',sprintf('%.0f',resEXACT.time),'\\ \cline{2-13}',...
   sprintf('\n'),'& 1 &(\ref{LossFuncFirstApprox}) / (\ref{midlineflowfirstapprox}) / (\ref{ISFirstApprox}) / (\ref{PIdef})&',sprintf('%.2f',res1.fdiff),' & ',sprintf('%.0f',res1.psum),' & ',sprintf('%.0f',res1.pnorm0),' & ',sprintf('%.2f',res1.pnorm1),' & ',sprintf('%.2f',res1.pnorminf),' & ',sprintf('%.2f',res1.meanlmp),' & ',sprintf('%.2f',res1.lmprange(1)),' / ',sprintf('%.2f',res1.lmprange(2)),' & ',sprintf('%.4f',res1.lmpnormdiff),' & ',sprintf('%.4f',res1.dtinfnorm),' & ',sprintf('%.0f',res1.time),'\\ \cline{2-13}',...
   sprintf('\n'),'& 2 &(\ref{approx1}) / (\ref{midlinefloweqtaylor}) / (\ref{IStaylor}) / (\ref{PIdef})&',sprintf('%.2f',res2.fdiff),' & ',sprintf('%.0f',res2.psum),' & ',sprintf('%.0f',res2.pnorm0),' & ',sprintf('%.2f',res2.pnorm1),' & ',sprintf('%.2f',res2.pnorminf),' & ',sprintf('%.2f',res2.meanlmp),' & ',sprintf('%.2f',res2.lmprange(1)),' / ',sprintf('%.2f',res2.lmprange(2)),' & ',sprintf('%.4f',res2.lmpnormdiff),' & ',sprintf('%.4f',res2.dtinfnorm),' & ',sprintf('%.0f',res2.time),'\\ \cline{2-13}',...
   sprintf('\n'),'& 3 &(\ref{AnotherApprox}) / (\ref{midlinefloweqfinal}) / (\ref{ISfinal}) / (\ref{PIdef})&',sprintf('%.2f',res3.fdiff),' & ',sprintf('%.0f',res3.psum),' & ',sprintf('%.0f',res3.pnorm0),' & ',sprintf('%.2f',res3.pnorm1),' & ',sprintf('%.2f',res3.pnorminf),' & ',sprintf('%.2f',res3.meanlmp),' & ',sprintf('%.2f',res3.lmprange(1)),' / ',sprintf('%.2f',res3.lmprange(2)),' & ',sprintf('%.4f',res3.lmpnormdiff),' & ',sprintf('%.4f',res3.dtinfnorm),' & ',sprintf('%.0f',res3.time),'\\ \cline{2-13}',...
   sprintf('\n'),'& 4 &(\ref{AnotherApprox}) / (\ref{midlinefloweqfinal}) / (\ref{ISfinal}) / (\ref{LDFs})&',sprintf('%.2f',res3LDF.fdiff),' & ',sprintf('%.0f',res3LDF.psum),' & ',sprintf('%.0f',res3LDF.pnorm0),' & ',sprintf('%.2f',res3LDF.pnorm1),' & ',sprintf('%.2f',res3LDF.pnorminf),' & ',sprintf('%.2f',res3LDF.meanlmp),' & ',sprintf('%.2f',res3LDF.lmprange(1)),' / ',sprintf('%.2f',res3LDF.lmprange(2)),' & ',sprintf('%.4f',res3LDF.lmpnormdiff),' & ',sprintf('%.4f',res3LDF.dtinfnorm),' & ',sprintf('%.0f',res3LDF.time),'\\\cline{2-13}',...
   sprintf('\n'),'& AC OPF & &',sprintf('%.2f',resAC.fdiff),' & ',sprintf('%.0f',resAC.psum),' & ',sprintf('%.0f',resAC.pnorm0),' & ',sprintf('%.2f',resAC.pnorm1),' & ',sprintf('%.2f',resAC.pnorminf),' & ',sprintf('%.2f',resAC.meanlmp),' & ',sprintf('%.2f',resAC.lmprange(1)),' / ',sprintf('%.2f',resAC.lmprange(2)),' & ',sprintf('%.4f',resAC.lmpnormdiff),' & ',sprintf('%.4f',resAC.dtinfnorm),' & ',sprintf('%.0f',resAC.time),'\\'];
