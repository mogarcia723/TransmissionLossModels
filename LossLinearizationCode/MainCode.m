%% Description
%This code produces the results used in Table II of the IEEE TPWRS paper entitled
%"Approximating Economic Dispatch by Linearizing Transmission Losses"
%
%Table II compares the performance of the Common LCED problem using
%different choices of base-case state.  It provides the dispatch
%approximation error when using a base-case state defined by 
%1)a small perturbation of the state, 2)the solution of the DC OPF problem, 
%and 3)the typical operating point provided by the test case.
%
%Note: The base-case state is chosen in the block of code beginning on line 159.
%      The code is initially set up to use the DC OPF problem when finding the base-case state
%
%The original code used mosek to solve most of the optimization problems.
%To simplify the setup procedures, this code uses MATLAB Optimization
%Toolbox.  For this reason, the results might be slightly different than in
%the paper.
%
%The results in Table I are not provided in this code. I was using the
%Knitro interior point solver to produce the results in Table I.  The
%Knitro software requires you purchase a license, so I removed that code.
%
%the Matpower toolbox is used to solve the AC and DC OPF problems.
%
%The standard MATLAB function fmincon() is used to solve the economic
%dispatch problems.  To improve convergence speed analytical hessians of
%the Lagrangian function are provided.

%% Quick note on units
%-LossFunc takes an argument in radians and outputs an arguement in p.u.
%-The vector P can be written as Pring=Aring'*B*Adot*thetadot where thetadot is in radians and Pring is in p.u.
%-The vector thetadot is written as thetadot=inv(Aring'*B*Adot)*Pring where thetadot is in radians and Pring is in p.u.
%-The base power is 100MW
%-The decision variables are in MW
%% Set MATLAB Path
%The MATLAB directory should be set to the main folder "LossLinearizationCode"
clc
clear
restoredefaultpath
%add you matpower path
addpath(genpath('./matpower6.0'))
%add nesta path
addpath(genpath('./nesta/opf'))

%% Load Test Case Data and Run Power Flow and Run AC OPF
%Choose the test case you want to use.
%Simple Test Cases: 118 9 30 39 57 118   
%Large Complicated test cases: 2383wp 3375wp 9241pegase 6515rte 3375wp 2746wop 3120sp 2737sop 2869pegase 300 case2868rte 2869pegase
%The paper uses 2383wp
mpc = loadcase(case2383wp);%set original mpc case
pfr=runpf(mpc);

opt=mpoption('opf.dc.solver','OT');%use OptimizationToolbox (Note: other solvers may be faster, e.g. MOSEK)
dcopf=rundcopf(mpc,opt);

if mpc.gencost(1,4)==2
    mpc.gencost(:,4)=3;
    mpc.gencost=[mpc.gencost(:,1:4),zeros(size(mpc.gencost,1),1),mpc.gencost(:,5:end)];
end
busnums=mpc.bus(:,1);

%% Creating temporary data files
%the function ``data'' outputs parameters associated with the test case
%from mpc.  ie. the admittance matrix Y and the transformer ratios TR
[TR,Y,injbus,n,m,nullbus,ties,vlim,ys,ysh,yshb,theta,tiesidx,genbus,genbusidx,temp,A]= data(mpc);


spA=sparse(A);%create sparse version of A
Aring=spA(:,2:end);
Adot=Aring;
spAdot=sparse(Adot);%create sparse version of Adot


D=sparse(mpc.bus(:,3));%demand vector
Dring=D(2:end,:);

% generator parameters
g=size(mpc.gen,1);%number of generators
genlims=[mpc.gen(:,10),mpc.gen(:,9)];%generator limits
gencost=mpc.gencost(:,5:7);

%Voltage magnitude fixed to Nominal Values
Vopt=ones(n,1); 

%series impedance, resistance, and reactance
xs=imag(1./ys);
rs=real(1./ys);
B=spdiags(1./xs,0,m,m);

InvMat=inv(Aring'*B*Adot);

M=zeros(n,g);% M is used as a cofficient matrix in Equality Constraints
for i=1:n
    for j=1:g
        if busnums(i)==mpc.gen(j,1)
            M(i,j)=1;
        end
    end
end

Mring=M(2:end,:);


%MPCtemp.mat contains data related to the test case that will be accessed
%in other functions (e.g. "LossFunc.m" and "myhessian.m")
save('MPCtemp.mat','M','Mring','InvMat','A','Adot','spA','spAdot','rs','xs','ties','tiesidx','D','gencost','genbus','genbusidx','n','m','g');

%% Setting the line limits
%line limits
unlimited=find(mpc.branch(:,6)==0);%lines with no limits
limited=find(mpc.branch(:,6)~=0);%lines with limits
%upper line limits
Fbar=zeros(m,1);
Fbar(limited)=mpc.branch(limited,6);
%lower line limits
Fubar=zeros(m,1);
Fubar(limited)=-mpc.branch(limited,6);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% TCED Problem with Angles %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is Opt. A in the paper
clc

%initial guess where x=[G;thetadot]
x0=zeros(g+n-1,1);

%upper and lower bounds on x
lb=[genlims(:,1);-Inf*ones(n-1,1)];
ub=[genlims(:,2);Inf*ones(n-1,1)];

%Linear Inequality Constraints (Ax<=b)
Aieq=[sparse(length(limited),g),100*B(limited,:)*spAdot;sparse(length(limited),g),-100*B(limited,:)*spAdot];
bieq=[Fbar(limited);-Fubar(limited)];


%Non-Linear Equality Constraints
constraintstemp=@(state) constraints(state);

%Hessian of the Lagrangian function
myhessiantemp=@(state,lam) myhessian(state,lam);

%fmincon options
options = optimoptions('fmincon','Algorithm','interior-point',...
    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
    'HessianFcn',myhessiantemp,'Display','iter','ConstraintTolerance', 1e-8,'OptimalityTolerance',1e-8,'StepTolerance',1e-8);%,'CheckGradients',true,'FiniteDifferenceType','central');



%running the TCED problem with angles
tic
[xTCED1,fvalTCED1,exitflagTCED1,outputTCED1,lambdaTCED1,gradTCED1,hessianTCED1]=fmincon(@cost,x0,Aieq,bieq,[],[],lb,ub,constraintstemp,options);
TCED1time=toc;


%% Choosing the Base-Case State
%Here are three options for choice of base-case state

%%Typical operating point as stated in the test case
%thetadotzero=(pfr.bus(2:end,9)-pfr.bus(1,9))*pi/180;%current operating point
%Gzero=pfr.gen(:,2);

%%Solution of the DC OPF problem
thetadotzero=(dcopf.bus(2:end,9)-dcopf.bus(1,9))*pi/180;%DCOPF
Gzero=dcopf.gen(:,2);

%%solution of the TCED problem with angles (Opt. A)
%thetadotzero=xTCED1((g+1):end);%solution to TCED
%Gzero=xTCED1(1:g);


%% Base-Case Values
%determining the base-case values P^0 and N^0 from the base-case state
[losses,dlosses]=LossFunc(Adot*thetadotzero,rs,xs,ties);
dlzero=dlosses*100;% derivative of lossfunc in basecase
lzero=losses*100; % losses in basecase
Nzero=0.5*abs(A)'*lzero;% base case nodal loss allocation 
Pzero=zeros(n,1);% Pzero is generation in each bus in basecase
for i=1:g
    Pzero(genbusidx(i))=Pzero(genbusidx(i))+xTCED1(i);
end
Pringzero=Pzero(2:end,:);  



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common LCED problem (same as LDF LCED problem with N variables reduced to single variable ell)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is Opt. F in the paper
% Decision variable is [G;ell]

%% Compute the sensitivity matrix as in Theorem 1 in the paper 
tic;
dNtilderingzero=(0.5/100)*InvMat'*(Adot'*dlzero*abs(Aring))*inv(eye(n-1)+(0.5/100)*InvMat'*(Adot'*dlzero*abs(Aring)));
dNtilde1zero=0.5*(eye(n-1)-dNtilderingzero)*InvMat'*(Adot'*dlzero*abs(A(:,1))/100);
dNtildezero=[dNtilde1zero,dNtilderingzero];% loss sensitive matrix Ntilde in basecase
senstoc=toc;

%% Formulate the Common LCED problem
eta=Nzero/(ones(1,n)*Nzero);
etaring=eta(2:end);

% Equality Constraints
Aeq1Common=[-ones(1,n)*(dNtildezero')*Mring,1];
beq1Common=ones(1,n)*Nzero-ones(1,n)*dNtildezero'*Pringzero;
Aeq2Common=[ones(1,n)*M,-1];
beq2Common=ones(1,n)*D;

AeqCommon=[Aeq1Common;Aeq2Common];
beqCommon=[beq1Common;beq2Common];

% Inequality Constraints
Aieq1Common=[B(limited,:)*Adot*InvMat*Mring,-B(limited,:)*Adot*InvMat*etaring];
Aieq2Common=-Aieq1Common;
Aieq3Common=[eye(g,g),zeros(g,1)];
Aieq4Common=-Aieq3Common;
AieqCommon=[Aieq1Common;Aieq2Common;Aieq3Common;Aieq4Common];

b1=Fbar(limited)+B(limited,:)*Adot*InvMat*Dring;
b2=-Fubar(limited)-B(limited,:)*Adot*InvMat*Dring;
Pmax=mpc.gen(:,9);
Pmin=mpc.gen(:,10);

bieqCommon=[b1;b2;Pmax;-Pmin];

% Objective Function
C1=mpc.gencost(:,5);
C2=mpc.gencost(:,6);
C3=mpc.gencost(:,7);
H11=2*diag(C1);
HCommon=sparse([H11,zeros(g,1);zeros(1,g),zeros(1,1)]);
fCommon=[C2;zeros(1,1)];

%solving the problem
if isempty(find(HCommon))%if it is a linear problem use linprog
    disp('solving linear LDF problem')
    tic
    [xCommon,fvalCommon,exitCommon,outputCommon,lambdaCommon] = linprog(sparse(fCommon),sparse(AieqCommon),sparse(bieqCommon),sparse(AeqCommon),sparse(beqCommon),[],[]);
    Commontime=toc;
else%if it is a quadratic problem use quadprog
    disp('solving quadratic LDF problem')
    tic
    [xCommon,fvalCommon,exitflagCommon,outputCommon,lambdaCommon] = quadprog(sparse(HCommon),sparse(fCommon),sparse(AieqCommon),sparse(bieqCommon),sparse(AeqCommon),sparse(beqCommon),[],[],[xTCED1(1:g);Nzero]);
    Commontime=toc;
end

%% Results
clc
disp(['time to solve TCED problem with angles (Opt. A): ',num2str(TCED1time)])
disp(['time to solve Common LCED (Opt. F): ',num2str(Commontime)])
disp(['time to compute sensitivity matrix: ',num2str(senstoc)])
disp(['number of binding line constraints: ',num2str(length(find(lambdaTCED1.ineqlin>.0001)))])

format bank
rowoftableII=[norm(thetadotzero-InvMat*(Mring*xCommon(1:g)-Dring-etaring*xCommon(g+1))/100,1),norm(M*xCommon(1:g)-M*xTCED1(1:g),1),norm(M*xCommon(1:g)-M*xTCED1(1:g),inf),abs(fvalTCED1-fvalCommon-sum(gencost(:,3))),Commontime]%variable cost results
