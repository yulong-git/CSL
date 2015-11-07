from __future__ import division
import numpy as np
from Hansen_defs import Hansen_defs

def Hansen_dyn(In,param):

	kpp =	In[0]
	hpp =	In[1]
	kp =	In[2]
	hp =	In[3]
	k = 	In[4]
	h = 	In[5]
	zp =	In[6]
	z =		In[7]

	delta = param[2]
	bet = param[3]
	D   = param[4]
	gam = param[5]


	yp, ip, cp, rp, wp = Hansen_defs(kp,hp,zp,kpp,hpp,param)
	y, i, c, r, w = Hansen_defs(k,h,z,kp,hp,param)

	out1 = c**(-gam)*w - D*(1-hp)**(-gam)
	out2 = bet*((c/cp)**gam)*(rp+1-delta)-1
	Out = np.array([out1, out2])
	return  Out