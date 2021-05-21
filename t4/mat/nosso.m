

beta1 = 178.7;
beta2 = 227.3;
R1 = 80*10^3;
R2 = 20*10^3;
Rc = 1*10^3;
Re = 100;
Rout = 100;
Rl = 8;
Vcc = 12;
Vbeon = 0.7;
VT=25e-3;
VAFN=69.7;

#    Ib1    Ic1 Ie1 Ib2   Ic2 Ie2 IR2 IR1
A = [-beta1  1   0  0     0   0    0   0;
     0       0   0 -beta2 1   0    0   0;
     1       0   0  0     0   0    1  -1;
     1       1  -1  0     0   0    0   0;
     0       0   0  1     1  -1    0   0;
     0       0   0  0     0   0   R2  R1;
     0      Rc  Re  0     0  Rc    0   0;
     0       0   0  0  Rout   0    0   0;]
     
v = [0; 0; 0; 0; 0; -Vcc; Vcc-Vbeon; Vcc-Vbeon]

res = A\v

#logo, v0 = vcc - Ic2*R0

V0 = Vcc-res(5)*Rout

#coiso e tal

gm1=res(2)/VT
rpi1=beta1/gm1
ro1=VAFN/res(2)
