function out = simple_dyn(in,param)

kpp =  in(1);
kp =   in(2);
k =    in(3);
zp =   in(4);
z =    in(5);

del = param(3);
bet = param(4);
gam = param(6);

[yp, ip, cp, rp] = simple_defs(kp,zp,kpp,param);
[y, i, c, r] = simple_defs(k,z,kp,param);

out = bet*((c/cp)^gam)*(rp+1-del)-1;