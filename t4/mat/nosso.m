clear all
close all

format long

pkg load symbolic

filename1 = "../doc/vce.tex";
filename2 = "../doc/gainstage_gain.tex";
filename3 = "../doc/outputstage_gain.tex";

fid1 = fopen (filename1, "w");
fid2 = fopen (filename2, "w");
fid3 = fopen (filename3, "w");

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
VAFP = 37.2
RB=1/(1/R1+1/R2);
RS=100;

##big system to calculate all currents :)
#gain stage

#    Ib1    Ic1 Ie1 Ib2   Ic2 Ie2 IR2 IR1
A = [-beta1  1   0  0     0   0    0   0;
     0       0   0 -beta2 1   0    0   0;
     1       0   0  0     0   0   -1   1;
     1       1  -1  0     0   0    0   0;
     0       0   0  1     1  -1    0   0;
     0       0   0  0     0   0   R2  R1;
     0      Rc  Re  Rc    0   0    0   0;
     0       0   0  0     0 Rout   0   0;]
     
v = [0; 0; 0; 0; 0; -Vcc; Vcc-Vbeon; Vcc-Vbeon]

res = A\v

#logo, v0 = vcc - Ic2*R0

V0 = Vcc-res(5)*Rout

#coiso e tal

gm1=res(2)/VT
rpi1=beta1/gm1
ro1=VAFN/res(2)

VE1=Re*res(3)
VO1=Vcc-Rc*res(2)
VCE=VO1-VE1 ##we just want this value bigger than Vbeon to be in the FAR

fprintf(fid1, "%f", VCE);

RSB=RB*RS/(RB+RS)
AV1 = RSB/RS * Rc*(Re-gm1*rpi1*ro1)/((ro1+Rc+Re)*(RSB+rpi1+Re)+gm1*Re*ro1*rpi1 - Re^2)
AVI_DB = 20*log10(abs(AV1))
AV1simple = RB/(RB+RS) * gm1*Rc/(1+gm1*Re)
AVIsimple_DB = 20*log10(abs(AV1simple))

ZI1 = 1/(1/RB+1/(((ro1+Rc+Re)*(rpi1+Re)+gm1*Re*ro1*rpi1 - Re^2)/(ro1+Rc+Re)))
ZO1 = 1/(1/ro1+1/Rc)

fprintf(fid2, "$%f$ & $%f$", ZI1, ZO1);

#now time for output stage
IC2 = res(5)
VI2 = VO1
VO2 = Vcc - Rout*res(6)

gm2 = IC2/VT
rpi2 = beta2/gm2
go2 = IC2/VAFP
gpi2 = gm2/beta2
ge2 = 1/Rout
ro2 = VAFN/res(5)

AV2 = gm2/(gm2+gpi2+go2+ge2)
ZI2 = (gm2+gpi2+go2+ge2)/gpi2/(gpi2+go2+ge2)
ZO2 = 1/(gm2+gpi2+go2+ge2)

fprintf(fid3, "$%f$ & $%f$", ZI2, ZO2);


%total
gB = 1/(1/gpi2+ZO1)
AV = (gB+gm2/gpi2*gB)/(gB+ge2+go2+gm2/gpi2*gB)*AV1
AV_DB = 20*log10(abs(AV))
ZI=ZI1
ZO=1/(go2+gm2/gpi2*gB+ge2+gB)

##this concludes point 2 of the theoretical analysis
##for point 3:

Cin = 1*10^(-6);
Cb = 1*10^(-6);
Co = 1*10^(-6);
Rl = 8;

syms w

##Zin(w) = RS+1/(j*w*Cin)
##Zed(w) = Re + 1/(j*w*Cb)
##Zl = Rl + 1/(j*w*Co)

syms Gin(w)
syms Ged(w)
syms Go1
syms Go2
syms Gc
syms Gout
syms Gm1
syms Gm2
syms Vin
syms acres(w)
syms Go(w)
syms v5f(w)

Gin(w) = 1/(RS+1/(j*w*Cin))
Ged(w) = 1/(1/Re + j*w*Cb)
Go(w) = j*w*Co

G1 = 1/R1
G2 = 1/R2
Gpi1 = 1/rpi1
Gpi2 = 1/rpi2
Go2 = 1/ro2
Go1 = 1/ro1
Gout = 1/Rout
Gm1 = gm1
Gc = 1/Rc
#Gm2 = gm2
Vin = 0.1


     #v1             #v2                 #v3           #v4 #v5
B = [Gin+G2+G1+Gpi1  -Gpi1                0             0   0;
     -Gpi1-Gm1       Gpi1+Ged+Go1+Gm1    -Go1+Gm1       0   0;
     Gm1             -Gm1-Go1            Go1+sym(Gc)+Gpi2   -Gpi2 0;
     0                0                  -Gpi2-gm2     Gpi2+Go2+Gout+sym(gm2)+Go -Go;
     0                0                   0            -Go  Go+sym(1/Rl)]

vb = [Vin*Gin; 0; 0; 0; 0]

acres(w) = B\vb

acres(w) = matlabFunction(vpa(abs(acres(w))))

fif = @(w) w(5);

v5f = @(w) fif(acres(w))

vpa(v5f(10))
vpa(v5f(100))

x = linspace(0, 1e9, 90);
y = cell(90);

#fplot(v5f, [0, 0.02], 201)

##for i=1:length(x)
##  y(i) = v5f(x(i));
##endfor


close all

fclose(fid1);
fclose(fid2);
fclose(fid3);