from __future__ import division
import numpy as np
from Hansen_dyn import Hansen_dyn
from LinApp_CSL import LinApp_CSL
from LinApp_SSL import LinApp_SSL
from LinApp_FindSS import LinApp_FindSS
from LinApp_Deriv import LinApp_Deriv
from LinApp_Solve import LinApp_Solve

# Hansen's model without labor/leisure decision
print("Hansen's model without labor/leisure decision")

#set model parameters
A = 1
theta = .33
delta = .025
bet = .995
gam = 1
rho = .9
sig = .02
D = 2.5
# set up parameter vector to pass to DSGE function file
param = [A, theta, delta, bet, D, gam, rho, sig]

#set numerical parameters
nx = 2
ny = 0
nz = 1
nobs = 250
logX = 1
DO_QZ = 0
do3 = 1


Zbar = [0]
# find SS numerically
guessXY = [.1, .33]
XYbar = LinApp_FindSS(Hansen_dyn,param,guessXY,Zbar,nx,ny)
print 'XYbar', XYbar
Xbar = XYbar[0:nx]
Ybar = XYbar[nx:nx+ny]
theta0 = np.append(np.concatenate((Xbar, Xbar, Xbar)),\
            np.concatenate((Zbar, Zbar)) )

NN = 0
#find derivatives and coefficients numerically
AA, BB, CC, DD, FF, GG, HH, JJ, KK, LL, MM, WW, TT = \
    LinApp_Deriv(Hansen_dyn,param,theta0,nx,ny,nz,logX)

PP, QQ, UU, RR, SS, VV = \
    LinApp_Solve(AA,BB,CC,DD,FF,GG,HH,JJ,KK,LL,MM,WW,TT,NN,Zbar,1)

print "PP\n", PP
print "QQ\n", QQ
print "RR\n", RR
print "SS\n", SS

#find quadratic approximation's steady state\
#[XQtil1, XQtil2] = QuadRoots(.5*HXX,HX-1,H0+sig^2*Hvv/2)

#generate a history of Z's
Z = np.zeros((nobs,nz))
# uncomment for simulation
eps = sig*np.random.randn(nobs,nz)
# uncomment for IRF
# eps = zeros(nobs,nz)
# eps(3,1) = sig
for t in xrange(1,nobs):
    Z[t,:] = Z[t-1,:].dot(NN) + eps[t,:]


# set starting values and simulate
XYbar = Xbar
X0 = Xbar

empty_vec= np.zeros(0)
empty_mat= np.zeros((0,nz))

XSSL, temp_SSL = LinApp_SSL(X0,Z,XYbar,logX,PP,QQ,UU,\
                            [], empty_mat,empty_mat,empty_vec)

XCSL, temp_CSL = LinApp_CSL(Hansen_dyn,param,X0,Z,NN,logX,0,[])


if do3:
    XCSL2, temp_sylv = LinApp_CSL(Hansen_dyn,param,X0,Z,NN,logX,1,[])

print "XSSL", XSSL
print "Difference with CSL",  XSSL - XCSL
print "Difference with CSL2", XSSL - XCSL2

'''
# plot results
if do3
    plotdatak = [XSSL(:,1) XCSL(:,1) XCSL2(:,1)]
    plotdatah = [XSSL(:,2) XCSL(:,2) XCSL2(:,2)]
else
    plotdatak = [XSSL(:,1) XCSL(:,1)]
    plotdatah = [XSSL(:,2) XCSL(:,2)]
end
figure
subplot(2,1,1)
plot(plotdatak)
subplot(2,1,2)
plot(plotdatah)
'''

