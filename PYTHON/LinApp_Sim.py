'''
Version 1.0, written by Kerk Phillips, April 2014

Adapted by Yulong Li, November 2015 
'''
from __future__ import division
from numpy import tile, array, empty, ndarray, zeros, log, asarray

def LinApp_Sim(Xm,Z,PP,QQ,UU,RR,SS,VV):
    '''
    Uses the coefficients from a linear approximation to generate data for
    next period given today's state. The set of endogenous state variables 
    known today is Xm and the set of exogenous state variables is Z.
    This program generates X.  

    The input and output values are in deviation from the linearization point 
    (almost always the steady  state, but not necessarily so).  This means 
    you will need to add back the steady state or other values after you have 
    called this function.  How you do this depends on whether you used 
    log-linearization or simple linearization in deriving the values of the 
    input coefficients.

    Parameters
    -----------
    Xm: array, dtype=float
        nx-by-1 vector of X(t-1) values

    Z: array, dtype=float
        nz-by-1 vector of Z(t) values

    PP: 2D-array, dtype=float
        nx-by-nx  matrix of X(t-1) on X(t) coefficients

    QQ: 2D-array, dtype=float
        nx-by-nz  matrix of Z(t) on X(t) coefficients

    UU: array, dtype=float
        nx-by-1 vector of X(t) constants

    RR: 2D-array, dtype=float
        ny-by-nx  matrix of X(t-1) on Y(t) coefficients

    SS: 2D-array, dtype=float
        ny-by-nz  matrix of Z(t) on Y(t) coefficients

    VV: array, dtype=float
        ny-by-1 vector of Y(t) constants

    Returns
    --------
    X: array, dtype=float
        nx-by-1 vector containing the value of the endogenous
        state variables for next period
    
    Y: array, dtype=float
        ny-by-1 vector containing the value of the endogenous
        non-state variables for the current period
    '''
% set RR, SS, and VV to empty matrices if not passed.
if (~exist('RR', 'var'))
    RR = [];
end
if (~exist('SS', 'var'))
    SS = [];
end
if (~exist('VV', 'var'))
    VV = [];
end

% Find the number of each kind of state variable
%  Using Xm find nx
[nrow,ncol] = size(Xm);
nx = nrow;
if nrow > 1 && ncol > 1
    disp('Xm must be a column vector')
    disp('you have input a 2-dimensional matrix')
elseif nrow == 1 && ncol > 1
    disp('Xm must be a column vector')
    disp('you have input row vector, which we will now transpose')
    Xm = Xm';
    nx = ncol;
end

%  Using RR find ny
[nrow,~] = size(RR);
ny = nrow;

%  Using Z find nz
[nrow,ncol]  = size(Z);
nz = nrow;
if nrow > 1 && ncol > 1
    disp('Z must be a column vector')
    disp('you have input a 2-dimensional matrix')
elseif nrow == 1 && ncol > 1
    disp('Z must be a column vector')
    disp('you have input row vector, which we will now transpose')
    Z = Z';
    nz = ncol;
end

% Check conformity of input coefficient matrices
[d1,d2] = size(PP);
if d1 ~= nx || d2 ~= nx
    disp('dimensions of PP incorrect')
end
[d1,d2] = size(QQ);
if d1 ~= nx || d2 ~= nz
    disp('dimensions of QQ incorrect')
end


% Generate data for next period, one equation at a time
X = PP*Xm + QQ*Z + UU;
if ny>0
    Y = RR*Xm + SS*Z + VV;
else
    Y = [];
end

end