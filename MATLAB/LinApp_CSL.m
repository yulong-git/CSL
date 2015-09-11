function [X, Y] = LinApp_CSL(funcname,param,X0,Z,NN,logX,Sylv,Y0)

% Version 1.0, written by Kerk Phillips, April 2014
%  
% Generates a history of X & Y variables by linearizing the policy function
% about the current state as in Evans & Phillips cited below.
%
% This function takes the following inputs:
%  funcname - is the name of the function which generates a column vector 
%          from ny+nx dynamic equations.
%  param - is a vector of parameter values to be passed to funcname.
%  X0    - 1-by-nx vector of X(1) starting values values.
%  Z     - nobs-by-nz matrix of Z values.
%  NN    - nz-by-nz matrix of VAR coefficients from the law of motion for Z
%  logX  - is an indicator that determines if the X & Y variables are
%          log-linearized (true) or simply linearized (false).  Z variables
%          are always simply linearized, default is 1.
%  Sylv  - is an indicator variable telling the program to use the built-in
%          function sylvester() to solve for QQ and SS, if possible.  
%          Default is Sylv=1.
%  Y0    - 1-by-ny vector of Y(1) starting values values.
%
% This function outputs the following:
%  X     - nobs-by-nx matrix containing the value of the endogenous
%          state variables for next period.
%  Y     - nobs-by-ny matix vector containing the value of the endogenous
%          non-state variables for next period.
%
% Source: R. Evans and K. Phillips (2014) "Linearization about the Current
% State: A Computational Method for Approximating Nonlinear Policy 
% Functions during Simulation," mimeo, Brigham Young University Department
% of Economics.
%
% Copyright: K. Phillips.  Feel free to copy, modify and use at your own 
% risk.  However, you are not allowed to sell this software or otherwise 
% impinge on its free distribution.

% Use log-linearized X & Y if no value is specified for logX
if (~exist('logX', 'var'))
    logX = true;
end
% set Y0 to empty vector if not passed.
if (~exist('Y0', 'var'))
    Y0 = [];
end
% Use Sylvester equation solution for QQ and SS by default
if (~exist('Sylv', 'var'))
    Sylv = 1;
end

% get values for nx, ny, nz and nobs
[nobs,nz] = size(Z);
[~,nx] = size(X0);
[~,ny] = size(Y0);

% Generate a history of X's and Y's
X = zeros(nobs,nx);
Y = zeros(nobs,ny);

% set starting values
X(1,:) = X0;
if ny>0
    Y(1,:) = Y0;
else
    Y(1,:) = [];
end

for t=1:nobs-1
    % set the linearization point to the current state
    if ny>0
        theta0 = [X(t,:) X(t,:) X(t,:) Y(t,:) Y(t,:)...
            Z(t+1,:) Z(t+1,:)]';
    else
        theta0 = [X(t,:) X(t,:) X(t,:) Z(t+1,:) Z(t+1,:)]';
    end
    % take derivatives at the linearization point
    [AA, BB, CC, DD, FF, GG, HH, JJ, KK, LL, MM, WW, TT] = ...
    LinApp_Deriv(funcname,param,theta0,nx,ny,nz,logX);
    % solve for coefficient matrices
    [PP, QQ, UU, RR, SS, VV] = LinApp_Solve(AA,BB,CC,DD,FF,GG,HH,JJ,...
        KK,LL,MM,WW,TT,NN,Z(t+1,:),Sylv);
    Xdev = zeros(1,nx);
    Zdev = zeros(1,nz);
    % Since LinApp_Sim uses column vectors and inputs, transpose
    if ny>0
        [Xtil, Ytil] = LinApp_Sim(Xdev',Zdev',PP,QQ,UU,RR,SS,VV);
        Ytil = Ytil';
    else
        [Xtil, ~] = LinApp_Sim(Xdev',Zdev',PP,QQ,UU);
    end
    Xtil = Xtil';
    
    % Convert to levels
    if logX
        X(t+1,:) = X(t,:).*exp(Xtil); 
        if ny> 0
            Y(t+1,:) = Y(t,:).*exp(Ytil);
        else
            Y = [];
        end
    else
        X(t+1,:) = X(t,:) + Xtil;
        if ny>0
            Y(t+1,:) = Y(t,:) + Ytil;
        else
            Y = [];
        end
    end
    
end
