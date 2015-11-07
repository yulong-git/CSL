'''
MATLAB version 1.0 written by Kerk Phillips, April 2014

PYTHON version adapted by Yulong Li, November 2015  
'''

from __future__ import division
from numpy import tile, array, concatenate, zeros, log, exp
from LinApp_Deriv import LinApp_Deriv
from LinApp_Solve import LinApp_Solve
from LinApp_Sim import LinApp_Sim

def LinApp_CSL(funcname,param,X0,Z,NN,logX,Sylv,Y0):
    '''    
    Generates a history of X & Y variables by linearizing the policy function
    about the current state as in Evans & Phillips cited below.

    Parameters
    -----------
    funcname: function
        the name of the function which generates a column vector 
        from ny+nx dynamic equations.
    
    param: array, dtype=float
        A vector of parameter values to be passed to funcname.
    
    X0: array, dtype=float
        nx vector of X(1) starting values values.
    
    Z: 2D-array, dtype=float
        nobs-by-nz matrix of Z values.
    
    NN: 2D-array, dtype=float
        nz-by-nz matrix of VAR coefficients from the law of motion for Z
    
    logX: binary, dtype=int
        an indicator that determines if the X & Y variables are
        log-linearized (true) or simply linearized (false).  Z variables
        are always simply linearized, default is 1. Use log-linearized 
        X & Y if no value is specified for logX
    
    Sylv: binary, dtype=int
        an indicator variable telling the program to use the built-in
        function sylvester() to solve for QQ and SS, if possible.  
        Default is Sylv=0.    
    
    Y0: array, dtype=float
        ny vector of Y(1) starting values values.

    Returns
    --------
    X: 2D-array, dtype=float
        nobs-by-nx matrix containing the value of the endogenous
        state variables.

    Y: 2D-array, dtype=float
        nobs-by-ny matrix vector containing the value of the endogenous
        non-state variables.
    Notes
    ------
    Source: R. Evans and K. Phillips (2014) "Linearization about the Current
    State: A Computational Method for Approximating Nonlinear Policy 
    Functions during Simulation," mimeo, Brigham Young University Department
    of Economics.
    '''
    # Formating
    X0 = array(X0)
    Y0 = array(Y0)
    
    # get values for nx, ny, nz and nobs
    nobs,nz = Z.shape
    nx = X0.shape[0]
    ny = Y0.shape[0]

    # Generate a history of X's and Y's
    X = zeros((nobs,nx))
    Y = zeros((nobs,ny))

    # set starting values
    X[0,:] = X0
    if ny>0:
        Y[0,:] = Y0

    for t in xrange(1, nobs):
        # set the linearization point to the current state
        if ny>0:
            theta0 = concatenate( (X[t-1,:], X[t-1,:], X[t-1,:], Y[t-1,:],\
                Y[t-1,:], Z[t,:], Z[t,:]) )
        else:
            theta0 = concatenate( (X[t-1,:], X[t-1,:], X[t-1,:], Z[t,:],\
                Z[t,:]) )
        # take derivatives at the linearization point
        AA, BB, CC, DD, FF, GG, HH, JJ, KK, LL, MM, WW, TT = \
        LinApp_Deriv(funcname,param,theta0,nx,ny,nz,logX)
        # solve for coefficient matrices
        PP, QQ, UU, RR, SS, VV = \
        LinApp_Solve(AA,BB,CC,DD,FF,GG,HH,JJ, KK,LL,MM,WW,TT,NN,Z[t,:],Sylv)
        
        Xdev = zeros((nx))
        Zdev = zeros((nz))
        Xtil, Ytil = LinApp_Sim(Xdev,Zdev,PP,QQ,UU,RR,SS,VV)

        if Ytil.any():
            # used to check whether ny=0 or not
            print('Note ny=0!')
        
        # Convert to levels
        if logX:
            X[t,:] = X[t-1,:]*exp(Xtil)
            if ny> 0:
                Y[t,:] = Y[t-1,:]*exp(Ytil)
        else:
            X[t,:] = X[t-1,:] + Xtil
            if ny>0:
                Y[t+1,:] = Y[t,:] + Ytil

    return array(X), array(Y)
