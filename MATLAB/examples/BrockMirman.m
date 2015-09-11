% Borck & Mirman model
clear

%set model parameters
alf = .35;
bet = .98;
sig = .02;
rho = .95;
% set up parameter vector to pass to DSGE function file
param = [alf bet sig rho];

%set numerical parameters
nx = 1;
ny = 0;
nz = 1;
options = optimset('Display','iter');
nobs = 250;
logX = 0;
DO_QZ = 0;
do3 = 1;


Zbar = 0;
% find SS numerically
XYbar = LinApp_FindSS(@BrockMirman_dyn,param,.1,Zbar,nx,ny);
Xbar = XYbar(1:nx);
Ybar = XYbar(nx+1:nx+ny);
theta0 = [Xbar; Xbar; Xbar; Zbar; Zbar];

NN = rho;
%find derivatives and coefficients numerically

[AA, BB, CC, DD, FF, GG, HH, JJ, KK, LL, MM, WW, TT] = ...
    LinApp_Deriv(@BrockMirman_dyn,param,theta0,nx,ny,nz,logX);

[PP, QQ, UU, RR, SS, VV] = ...
    LinApp_Solve(AA,BB,CC,DD,FF,GG,HH,JJ,KK,LL,MM,WW,TT,NN,Zbar);


%generate a history of Z's
Z = zeros(nobs,nz);
% uncomment for simulation
eps = sig*randn(nobs,nz);
% uncomment for IRF
% eps = zeros(nobs,nz);
% eps(3,1) = sig;
for t=1:nobs-1
    Z(t+1,:) = Z(t,:)*NN + eps(t+1,:);
end

% set starting values and simulate
XYbar = Xbar;
X0 = Xbar;

%  steady state linarization
tic;
[XSSL, ~] = LinApp_SSL(X0,Z,XYbar,logX,PP,QQ,UU);
toc

%  current state linarization
tic;
[XCSL, ~] = LinApp_CSL(@BrockMirman_dyn,param,X0,Z,NN,logX);
toc

%  current state linarization with steady state PP & QQ
if do3
    tic;
    [XCSL2, ~] = LinApp_CSL(@BrockMirman_dyn,param,X0,Z,NN,logX,[],1,PP,QQ);
    toc
end

%  exact solution
Xexact = zeros(nobs,nx);
Xexact(1,:) = X0;
for t=1:nobs-1
    Xexact(t+1,:) = alf*bet*exp(Z(t+1,:))*Xexact(t,:)^alf;
end


if do3
    plotdata = [XSSL XCSL XCSL2 Xexact];
    ratiodata = [log(XSSL./Xexact) log(XCSL./Xexact) log(XCSL2./Xexact)];
else
    plotdata = [XSSL XCSL Xexact];
    ratiodata = [log(XSSL./Xexact) log(XCSL./Xexact)];
end

figure;
plot(plotdata) 
figure;
plot(ratiodata)
