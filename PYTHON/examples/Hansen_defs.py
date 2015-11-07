from __future__ import division
import numpy as np

def Hansen_defs(k,h,z,kp,hp,param):

	A = param[0]
	theta = param[1]
	delta = param[2]


	y = A*(k**theta*(np.exp(z)*hp)**(1-theta))
	i = kp - (1-delta)*k
	c = y - i
	r = theta*y/k
	w = (1-theta)*y/hp

	return y, i, c, r, w