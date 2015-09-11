function [y, i, c, r, w] = simple_defs(X,Z,Xp,param)

A = param(1);
theta = param(2);
del = param(3);
hbar = param(5);

k = X(1);
z = Z(1);
kp = Xp(1);

y = A*(k^theta*(exp(z)*hbar)^(1-theta));
i = kp - (1-del)*k;
c = y - i;
r = theta*y/k;
w = (1-theta)*y/hbar;