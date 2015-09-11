function out = BrockMirman_dyn(in,param)

kpp =  in(1);
kp =   in(2);
k =    in(3);
zp =   in(4);
z =    in(5);

alf = param(1);
bet = param(2);

c  = exp(z)*k^alf - kp;
cp = exp(zp)*kp^alf - kpp;
rp = alf*exp(zp)*kp^(alf-1);
out = bet*(c/cp)*rp - 1;