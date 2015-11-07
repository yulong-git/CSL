from __future__ import division
import numpy as np

def BrockMirman_dyn(In,param):

	kpp =  In[0]
	kp =   In[1]
	k =    In[2]
	zp =   In[3]
	z =    In[4]

	alf = param[0]
	bet = param[1]

	c  = np.exp(z)*k**alf - kp
	cp = np.exp(zp)*kp**alf - kpp
	rp = alf*np.exp(zp)*kp**(alf-1)
	Out = bet*(c/cp)*rp - 1
	return Out