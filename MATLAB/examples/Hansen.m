% Hansen's model without labor/leisure decision
clear

%set model parameters
A = 1;
theta = .33;
del = .025;
bet = .995;
gam = 1;
rho = .9;
sig = .02;
D = 2.5;
% set up parameter vector to pass to DSGE function file
param = [A theta del bet D gam rho sig];

%set numerical parameters
nx = 2;
ny = 0;
nz = 1;
options = optimset('Display','iter');
nobs = 250;
logX = true;
DO_QZ = false;
do3 = 1;


Zbar = 0;
% find SS numerically
guessXY = [.1; .33];
XYbar = LinApp_FindSS(@Hansen_dyn,param,guessXY,Zbar,nx,ny);
Xbar = XYbar(1:nx);
Ybar = XYbar(nx+1:nx+ny);
theta0 = [Xbar; Xbar; Xbar; Ybar; Ybar; Zbar; Zbar];

NN = 0;
%find derivatives and coefficients numerically

[AA, BB, CC, DD, FF, GG, HH, JJ, KK, LL, MM, WW, TT] = ...
    LinApp_Deriv(@Hansen_dyn,param,theta0,nx,ny,nz,logX);

[PP, QQ, UU, RR, SS, VV] = ...
    LinApp_Solve(AA,BB,CC,DD,FF,GG,HH,JJ,KK,LL,MM,WW,TT,NN,Zbar);

%find quadratic approximation's steady state\
%[XQtil1, XQtil2] = QuadRoots(.5*HXX,HX-1,H0+sig^2*Hvv/2);

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
tic;
[XSSL, ~] = LinApp_SSL(X0',Z,XYbar',logX,PP,QQ,UU);
toc

tic;
[XCSL, ~] = LinApp_CSL(@Hansen_dyn,param,X0',Z,NN,logX);
toc

if do3
    tic;
    [XCSL2, ~] = LinApp_CSL(@Hansen_dyn,param,X0',Z,NN,logX,[]...
        ,1,PP,QQ);
    toc
end

% plot results
if do3
    plotdatak = [XSSL(:,1) XCSL(:,1) XCSL2(:,1)];
    plotdatah = [XSSL(:,2) XCSL(:,2) XCSL2(:,2)];
else
    plotdatak = [XSSL(:,1) XCSL(:,1)];
    plotdatah = [XSSL(:,2) XCSL(:,2)];
end
figure;
subplot(2,1,1)
plot(plotdatak)
subplot(2,1,2)
plot(plotdatah)







































