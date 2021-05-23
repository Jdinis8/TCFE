%gain stage

VT=25e-3
BFN=178.7
VAFN=69.7
RE1=100
RC1=1000
RB1=80000
RB2=20000
VBEON=0.7
VCC=12
RS=100

RB=1/(1/RB1+1/RB2)
VEQ=RB2/(RB1+RB2)*VCC
IB1=(VEQ-VBEON)/(RB+(1+BFN)*RE1)
IC1=BFN*IB1
IE1=(1+BFN)*IB1
VE1=RE1*IE1
VO1=VCC-RC1*IC1
VCE=VO1-VE1


gm1=IC1/VT
rpi1=BFN/gm1
ro1=VAFN/IC1

RSB=RB*RS/(RB+RS)

AV1 = RSB/RS * RC1*(RE1-gm1*rpi1*ro1)/((ro1+RC1+RE1)*(RSB+rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)
AVI_DB = 20*log10(abs(AV1))
AV1simple = RB/(RB+RS) * gm1*RC1/(1+gm1*RE1)
AVIsimple_DB = 20*log10(abs(AV1simple))

RE1=0
AV1 = RSB/RS * RC1*(RE1-gm1*rpi1*ro1)/((ro1+RC1+RE1)*(RSB+rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)
AVI_DB = 20*log10(abs(AV1))
AV1simple =  - RSB/RS * gm1*RC1/(1+gm1*RE1)
AVIsimple_DB = 20*log10(abs(AV1simple))

RE1=100
ZI1 = 1/(1/RB+1/(((ro1+RC1+RE1)*(rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)/(ro1+RC1+RE1)))
ZX = ro1*((RSB+rpi1)*RE1/(RSB+rpi1+RE1))/(1/(1/ro1+1/(rpi1+RSB)+1/RE1+gm1*rpi1/(rpi1+RSB)))
ZX = ro1*(   1/RE1+1/(rpi1+RSB)+1/ro1+gm1*rpi1/(rpi1+RSB)  )/(   1/RE1+1/(rpi1+RSB) ) 
ZO1 = 1/(1/ZX+1/RC1)

RE1=0
ZI1 = 1/(1/RB+1/(((ro1+RC1+RE1)*(rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)/(ro1+RC1+RE1)))
ZO1 = 1/(1/ro1+1/RC1)

%ouput stage
BFP = 227.3
VAFP = 37.2
RE2 = 100
VEBON = 0.7
VI2 = VO1
IE2 = (VCC-VEBON-VI2)/RE2
IC2 = BFP/(BFP+1)*IE2
VO2 = VCC - RE2*IE2

gm2 = IC2/VT
rpi2 = beta2/gm2
go2 = IC2/VAFP
gpi2 = gm2/BFP
ge2 = 1/RE2

AV2 = gm2/(gm2+gpi2+go2+ge2)
ZI2 = (gm2+gpi2+go2+ge2)/gpi2/(gpi2+go2+ge2)
ZO2 = 1/(gm2+gpi2+go2+ge2)


%total
gB = 1/(1/gpi2+ZO1)
AV = (gB+gm2/gpi2*gB)/(gB+ge2+go2+gm2/gpi2*gB)*AV1
AV_DB = 20*log10(abs(AV))
ZI=ZI1
ZO=1/(go2+gm2/gpi2*gB+ge2+gB)

pkg load symbolic

syms w

##Zin(w) = RS+1/(j*w*Cin)
##Zed(w) = Re + 1/(j*w*Cb)
##Zl = Rl + 1/(j*w*Co)

syms w
syms Gin(w)
syms Ged(w)
syms G1
syms G2
syms Gpi1
syms Gpi2
syms Go1
syms Go2
syms Gc
syms Gl
syms Gout
syms Gm1
syms Gm2
syms Vin
syms acres(w)
syms Go(w)
syms v5f(w)

beta1 = 178.7;
beta2 = 227.3;
R1 = 95*10^3;
R2 = 10*10^3;
Rc = 2.5*10^3;
Re = 0.13*10^3;
Rout = 0.25*10^3;
Rl = 8;
Vcc = 12;
Vbeon = 0.7;
VT=25e-3;
VAFN=69.7;
VAFP = 37.2
RB=1/(1/R1+1/R2);
RS=100;
Cin = 135*10^(-6);
Cb = 1900*10^(-6);
Co = 1190*10^(-6);
Rl = 8;
Cb1 = 1.61*10^(-11);
Cb2 = 1.4*10^(-11);
Ce1 = 4.388*10^(-12);
Ce2 = 1.11*10^(-11);


Gin(w) = sym(1)/(sym(RS, 'f')+sym(1)/(j*w*sym(Cin, 'f')));
Ged(w) = sym(1/Re, 'f') + j*w*sym(Cb,'f');
Go(w) = j*w*sym(Co,'f');

G1 = sym(1/R1, 'f');
G2 = sym(1/R2, 'f');
Gpi1 = sym(1/rpi1,'f')+j*w*sym(Cb1, 'f');
Gpi2 = sym(1/rpi2) + j*w*sym(Cb2, 'f');
Ge2 = j*w*sym(Ce2);
Ge1 = j*w*sym(Ce1);
Go2 = sym(1/ro2, 'f');
Go1 = sym(1/ro1, 'f');
Gout = sym(1/Rout, 'f');
Gm1 = sym(gm1, 'f');
Gc = sym(1/Rc, 'f');
Gm2 = sym(gm2, 'f');
Gl = sym(1/Rl, 'f');
Vin = sym(0.01, 'f');


     #v1                      #v2                   #v3                   #v4                      #v5
B(w) = [Gin(w)+G2+G1+Gpi1+Ge1 -Gpi1                 -Ge1                   0                        0;
       -Gpi1-Gm1               Gpi1+Ged(w)+Go1+Gm1  -Go1                   0                        0;
        Gm1-Ge1               -Gm1-Go1               Go1+Gc+Gpi2+Ge1+Ge2  -Gpi2-Ge2                 0;
        0                      0                    -Gpi2-Gm2              Gpi2+Go2+Gout+Gm2+Go(w) -Go(w);
        0                      0                     0                    -Go(w)                    Go(w)+Gl]

vb(w) = [Vin*Gin(w); 0; 0; 0; 0];

acres(w) = B(w)\vb(w);

h = matlabFunction(abs(vpa((acres(w)))));


##
fif = @(w) w(5);
v5f = @(w) log10(fif(h(10^(w)))/0.01)

fplot(v5f, [0, 10], 2000)
