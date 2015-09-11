function out = Hansen_dyn(in,param)

kpp =  in(1);
hpp =  in(2);
kp =   in(3);
hp =   in(4);
k =    in(5);
h =    in(6);
zp =   in(7);
z =    in(8);

del = param(3);
bet = param(4);
D   = param(5);
gam = param(6);


[yp, ip, cp, rp, wp] = Hansen_defs(kp,hp,zp,kpp,hpp,param);
[y, i, c, r, w] = Hansen_defs(k,h,z,kp,hp,param);

out1 = c^(-gam)*w - D*(1-hp)^(-gam);
out2 = bet*((c/cp)^gam)*(rp+1-del)-1;
out = [out1; out2];