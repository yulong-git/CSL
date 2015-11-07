from __future__ import division
import numpy as np
from BrockMirman_dyn import BrockMirman_dyn
from LinApp_CSL import LinApp_CSL
from LinApp_SSL import LinApp_SSL
from LinApp_FindSS import LinApp_FindSS
from LinApp_Deriv import LinApp_Deriv
from LinApp_Solve import LinApp_Solve

##### Borck & Mirman model #####
print("Borck & Mirman model")

#set model parameters
alf = .35
bet = .98
sig = .02
rho = .95
# set up parameter vector to pass to DSGE function file
param = [alf, bet, sig, rho]

#set numerical parameters
nx = 1
ny = 0
nz = 1
nobs = 250
logX = 0
do3 = 1


Zbar = [0]
# find SS numerically
XYbar = LinApp_FindSS(BrockMirman_dyn,param,.1,Zbar,nx,ny)
print 'XYbar', XYbar
Xbar = XYbar[0:nx]
Ybar = XYbar[nx:nx+ny]
theta0 = np.concatenate((Xbar, Xbar, Xbar, Zbar, Zbar))

NN = rho
#find derivatives and coefficients numerically

AA, BB, CC, DD, FF, GG, HH, JJ, KK, LL, MM, WW, TT = \
    LinApp_Deriv(BrockMirman_dyn,param,theta0,nx,ny,nz,logX)

PP, QQ, UU, RR, SS, VV = \
    LinApp_Solve(AA,BB,CC,DD,FF,GG,HH,JJ,KK,LL,MM,WW,TT,NN,Zbar,0)

print "PP\n", PP
print "QQ\n", QQ
print "RR\n", RR
print "SS\n", SS

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

#  steady state linarization
empty_vec= np.zeros(0)
empty_mat= np.zeros((0,0))

XSSL, temp_SSL = LinApp_SSL(X0,Z,XYbar,logX,PP,QQ,UU,\
                            [],empty_mat,empty_mat,empty_vec)

#  current state linarization
XCSL, temp_CSL = LinApp_CSL(BrockMirman_dyn,param,X0,Z,NN,logX,0,[])

#  current state linarization with steady state PP & QQ
if do3:
    XCSL2, temp_slyv = LinApp_CSL(BrockMirman_dyn,param,X0,Z,NN,logX,1,[])

#  exact solution
Xexact = np.zeros((nobs,nx))
Xexact[0,:] = X0
for t in xrange(1,nobs):
    Xexact[t,:] = alf*bet*np.exp(Z[t,:])*Xexact[t-1,:]**alf

print "Exact", Xexact
print "Difference with CSL", Xexact - XCSL
print "Difference with CSL2", Xexact - XCSL2
print "Difference with SSL", Xexact - XSSL

'''
if do3
    plotdata = [XSSL XCSL XCSL2 Xexact]
    ratiodata = [log(XSSL./Xexact) log(XCSL./Xexact) log(XCSL2./Xexact)]
else
    plotdata = [XSSL XCSL Xexact]
    ratiodata = [log(XSSL./Xexact) log(XCSL./Xexact)]

figure
plot(plotdata) 
figure
plot(ratiodata)
'''