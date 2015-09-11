function [y, i, c, r, w] = Hansen_defs(k,h,z,kp,hp,param)

A = param(1);
theta = param(2);
del = param(3);


y = A*(k^theta*(exp(z)*hp)^(1-theta));
i = kp - (1-del)*k;
c = y - i;
r = theta*y/k;
w = (1-theta)*y/hp;