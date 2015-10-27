function [AA, BB, CC, DD, FF, GG, HH, JJ, KK, LL, MM, WW, TT] = ...
    LinApp_Deriv(funcname,param,theta0,nx,ny,nz,logX)

% Version 2.1, written by Kerk Phillips, previously named RBCnumerderiv.m
% Updated April 1, 2014
% This function numerically differentiates a set of functions from a DSGE
% model to yield the matrix of derivatives needed to find the linear
% approximation to the model's policy function.  The names of the subsets
% of this matrix of derivatives are consistent with the notation and code
% from Uhlig's toolkit.
% Implements Heer & Maussner (2007), chapter 11.3.1 numerical
% differentiation.
%
% This function takes the following inputs:
%  funcname - is the name of the function which generates a column vector 
%  from ny+nx dynamic equations. 
%    The ny equations to be linearized into the form below in the first 
%    ny rows.
%     A X(t) + B X(t-1) + C Y(t) + D Z(t) = 0 
%    The function must be written so as to evaluate to zero for all rows
%    in the steady state.
%  param - is a vector of parameter values to be passed to funcname.
%  theta0 - is a matrix of nx+ny values for the linearization point,
%    with the values of X in the first nx rows, and the values of Y in
%    the middle ny rows.
%  nx - is the number of elements in X.
%  ny - is the number of elements in Y.
%  nz - is the number of elements in Z.
%  logX - is an indicator that determines if the X & Y variables are
%    log-linearized (true) or simply linearized (false).  Z variables are
%    always simply linearized.
%
% The function generates matrices AA thru MM from the log-linearization of 
% the equations in the function "funcname".
%  WW is an ny-by-1 vector of constants from the first ny equations.
%  TT is an nx-by-1 vector of constants from the last nx equations.
%   When linearizing about the steady state WW and TT are zeroes.
%
% Copyright: K. Phillips.  Feel free to copy, modify and use at your own 
% risk.  However, you are not allowed to sell this software or otherwise 
% impinge on its free distribution.

% Use log-linearized X & Y if no value is specified for logX
if (~exist('logX', 'var'))
    logX = true;
end

% Calculate the increments to use in calculating numerical derivatives
leng = 3*nx+2*(ny+nz);
eps = 2.2E-16;  % Machine epsilon for double precision 2.2E-16
incr = sqrt(eps)*max(theta0,1000*sqrt(eps)*ones(leng,1));
thetap = theta0 + incr;
thetam = theta0 - incr;
dx = thetap - thetam;

% Value of the function at the linearization point 
%  This should be zero or very close to it when evaluating at the SS
T0 = funcname(theta0,param);

% Create matrices of deviations for input, devp and devm 
% Deviate columns by adding and subtracting incr
devp=repmat(theta0,1,leng);
devm=repmat(theta0,1,leng);
for i=1:leng
  devp(i,i)=devp(i,i)+incr(i);
  devm(i,i)=devm(i,i)-incr(i);
end

% Calculate matrices of deviations for dynamic equations, bigp and bigm
%  rows correspond to equations
%  colums correspond to variables from "theta0" vector being changed
%  note output of funcname must be a column vector

% Initialize the function calls when inputs are deviated
bigp = zeros(nx+ny,leng);
bigm = zeros(nx+ny,leng);

% Calculate matrices of deviations for dynamic equations, bigp and bigm
%  rows correspond to equations
%  colums correspond to variables from "theta0" vector being changed
%  note output of funcname must be a column vector
for i = 1:leng
    if logX
        if i<3*nx+2*ny+1
            bigp(:,i) = theta0(i)*funcname(devp(:,i),param)./(1+T0);
            bigm(:,i) = theta0(i)*funcname(devm(:,i),param)./(1+T0);
        else
            bigp(:,i) = funcname(devp(:,i),param)./(1+T0);
            bigm(:,i) = funcname(devm(:,i),param)./(1+T0);
        end
    else
        bigp(:,i) = funcname(devp(:,i),param);
        bigm(:,i) = funcname(devm(:,i),param);
    end
end
% Calclate the derivatives using the central difference formula
big = (bigp-bigm)./repmat(dx',nx+ny,1);

% Pull out appropriate parts of "big" into Uhlig's matrices
AA = big(1:ny,nx+1:2*nx);
BB = big(1:ny,2*nx+1:3*nx);
CC = big(1:ny,3*nx+ny+1:3*nx+2*ny);
DD = big(1:ny,3*nx+2*ny+nz+1:leng);
FF = big(ny+1:ny+nx,1:nx);
GG = big(ny+1:ny+nx,nx+1:2*nx);
HH = big(ny+1:ny+nx,2*nx+1:3*nx);
JJ = big(ny+1:ny+nx,3*nx+1:3*nx+ny);
KK = big(ny+1:ny+nx,3*nx+ny+1:3*nx+2*ny);
LL = big(ny+1:ny+nx,3*nx+2*ny+1:3*nx+2*ny+nz);
MM = big(ny+1:ny+nx,3*nx+2*ny+nz+1:leng);
if logX
    TT = log(1+T0);
else
    TT = T0;
end
WW = TT(1:ny,:);
TT = TT(ny+1:nx+ny,:);
end
