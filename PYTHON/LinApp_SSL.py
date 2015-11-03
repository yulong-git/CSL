'''
Version 1.0, written by Kerk Phillips, April 2014

Adapted by Yulong Li, November 2015 
'''
from __future__ import division
from numpy import tile, array, empty, ndarray, zeros, log, asarray

def LinApp_SSL(X0,Z,XYbar,logX=1,PP,QQ,UU,Y0,RR,SS,VV):
    '''
    Generates a history of X & Y variables by linearizing the policy function
    about the steady state as in Uhlig's toolkit.
    
    Parameters
    -----------    
    X0: array, dtype=float
        1-by-nx vector of X(1) starting values values
    
    Z: 2D-array, dtype=float
        nobs-by-nz matrix of Z values
    
    XYbar: array, dtype=float
        1-by-(nx+ny) vector of X and Y steady state values
    
    logX: binary
        an indicator that determines if the X & Y variables are
        log-linearized (true) or simply linearized (false).  Z variables
        are always simply linearized.
    
    PP: 2D-array, dtype=float
        nx-by-nx matrix of X(t-1) on X(t) coefficients
    
    QQ: 2D-array, dtype=float
        nx-by-nz  matrix of Z(t) on X(t) coefficients
    
    UU: array, dtype=float
        nx-by-1 vector of X(t) constants
    
    Y0: array, dtype=float
        1-by-ny vector of Y(1) starting values values.
    
    RR: 2D-array, dtype=float
        ny-by-nx  matrix of X(t-1) on Y(t) coefficients
    
    SS: 2D-array, dtype=float
        ny-by-nz  matrix of Z(t) on Y(t) coefficients
    
    VV: 2D-array, dtype=float
        ny-by-1 vector of Y(t) constants
    
    Returns
    --------
    X: 2D-array, dtype=float
        nobs-by-nx matrix containing the value of the endogenous
        state variables
    
    Y: 2D-array, dtype=float
        nobs-by-ny matrix vector containing the value of the endogenous
        non-state variables
    '''
    % set Y0, RR, SS, and VV to empty matrices if not passed.
    if ()
        Y0 = array([]);
    end
    if (~exist('RR', 'var'))
        RR = [];
    end
    if (~exist('SS', 'var'))
        SS = [];
    end
    if (~exist('VV', 'var'))
        VV = [];
    end

    % get values for nx, ny, nz and nobs
    [nobs,nz] = size(Z);
    [~,nx] = size(X0);
    [~,nxy] = size(XYbar);
    ny = nxy - nx;

    % get Xbar and Ybar
    Xbar = XYbar(:,1:nx);
    Ybar = XYbar(:,nx+1:nx+ny);

    % Generate a history of X's and Y's
    Xtil = zeros(nobs,nx);
    Ytil = zeros(nobs,ny);
    % set starting values
    X(1,:) = X0;
    if ny>0
        Y(1,:) = Y0;
    end
    if logX
        Xtil(1,:) = log(X(1,:)./Xbar);
        if ny>0
            Ytil(1,:) = log(Y(1,:)./Ybar);
        end
    else
        Xtil(1,:) = X(1,:) - Xbar;
        if ny>0
            Ytil(1,:) = Y(1,:) - Ybar;
        end
    end
    for t=1:nobs-1:
        # Since LinApp_Sim uses column vectors and inputs, transpose
        if ny>0
            [Xtemp, Ytemp] =\
                LinApp_Sim(Xtil(t,:)',Z(t+1,:)',PP,QQ,UU,RR,SS,VV)
            Ytil(t+1,:) = Ytemp.T
        else
            [Xtemp, ~] =\
                LinApp_Sim(Xtil(t,:)',Z(t+1,:)',PP,QQ,UU)
        Xtil(t+1,:) = Xtemp.T

    # Convert to levels
    if logX:
        X = repmat(Xbar,nobs,1).*exp(Xtil)
        if ny> 0:
            Y = repmat(Ybar,nobs,1).*exp(Ytil)
        else:
            Y = []
    else:
        X = repmat(Xbar,nobs,1)+Xtil
        if ny>0:
            Y = repmat(Ybar,nobs,1)+Ytil
        else:
            Y = []