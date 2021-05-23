clear all
close all

format long

pkg load symbolic

filename1 = "../doc/vce.tex";
filename2 = "../doc/gainstage_gain.tex";
filename3 = "../doc/outputstage_gain.tex";
filename4 = "../doc/merit_sim.tex";
filename5 = "../doc/total_impedance.tex";

fid1 = fopen (filename1, "w");
fid2 = fopen (filename2, "w");
fid3 = fopen (filename3, "w");
fid4 = fopen (filename4, "w");
fid5 = fopen (filename5, "w");

beta1 = 178.7;
beta2 = 227.3;
R1 = 90*10^3;
R2 = 10*10^3;
Rc = 2.3*10^3;
Re = 0.12*10^3;
Rout = 0.33*10^3;
Rl = 8;
Vcc = 12;
Vbeon = 0.7;
VT=25e-3;
VAFN=69.7;
VAFP = 37.2;
RB=1/(1/R1+1/R2);
RS=100;
Cin = 135*10^(-6);
Cb = 1840*10^(-6);
Co = 1200*10^(-6);
Rl = 8;

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
     0       0   0  0     0 Rout   0   0;];
     
v = [0; 0; 0; 0; 0; -Vcc; Vcc-Vbeon; Vcc-Vbeon];

res = A\v;

#logo, v0 = vcc - Ic2*R0

V0 = Vcc-res(5)*Rout;

#coiso e tal

gm1=res(2)/VT;
rpi1=beta1/gm1;
ro1=VAFN/res(2);

VE1=Re*res(3);
VO1=Vcc-Rc*res(2);
VCE=VO1-VE1; ##we just want this value bigger than Vbeon to be in the FAR

fprintf(fid1, "%f", VCE);

RSB=RB*RS/(RB+RS);
AV1 = RSB/RS * Rc*(Re-gm1*rpi1*ro1)/((ro1+Rc+Re)*(RSB+rpi1+Re)+gm1*Re*ro1*rpi1 - Re^2);
AVI_DB = 20*log10(abs(AV1));
AV1simple = RB/(RB+RS) * gm1*Rc/(1+gm1*Re);
AVIsimple_DB = 20*log10(abs(AV1simple));

ZI1 = 1/(1/RB+1/(((ro1+Rc+Re)*(rpi1+Re)+gm1*Re*ro1*rpi1 - Re^2)/(ro1+Rc+Re)));
ZO1 = 1/(1/ro1+1/Rc);

fprintf(fid2, "$%f$ & $%f$\\\\", ZI1, ZO1);

#now time for output stage
IC2 = res(5);
VI2 = VO1;
VO2 = Vcc - Rout*res(6);

gm2 = IC2/VT;
rpi2 = beta2/gm2;
go2 = IC2/VAFP;
gpi2 = gm2/beta2;
ge2 = 1/Rout;
ro2 = VAFN/res(5);

AV2 = gm2/(gm2+gpi2+go2+ge2);
ZI2 = (gm2+gpi2+go2+ge2)/gpi2/(gpi2+go2+ge2);
ZO2 = 1/(gm2+gpi2+go2+ge2);

fprintf(fid3, "$%f$ & $%f$\\\\", ZI2, ZO2);


%total
gB = 1/(1/gpi2+ZO1);
AV = (gB+gm2/gpi2*gB)/(gB+ge2+go2+gm2/gpi2*gB)*AV1;
AV_DB = 20*log10(abs(AV));
ZI=ZI1;
ZO=1/(go2+gm2/gpi2*gB+ge2+gB)

fprintf(fid5, "$%f$\\\\", Z0);

##this concludes point 2 of the theoretical analysis
##for point 3:

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
Gpi2 = sym(1/rpi2, 'f') + j*w*sym(Cb2, 'f');
Ge2 = j*w*sym(Ce2, 'f');
Ge1 = j*w*sym(Ce1, 'f');
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
        0                      0                     0                    -Go(w)                    Go(w)+Gl];

vb(w) = [Vin*Gin(w); 0; 0; 0; 0];

acres(w) = B(w)\vb(w);

h = matlabFunction(abs(vpa((acres(w)))));

##
fif = @(w) w(5);
v5f = @(w) 20*log10(fif(h(2*pi*10^(w)))/0.01);

third = @(w) w(3);
v3f = @(w) 20*log10(third(h(2*pi*10^(w)))/0.01);

fplot(v5f, [0, 9], 2000);
legend("Gain_{db}");
xlabel("log10(f) [f] = Hz");
ylabel("log10(y)");
print ("gaindb.png", "-dpng");

fplot(v3f, [0, 9], 2000);
legend("Gain_{db}");
xlabel("log10(f) [f] = Hz");
ylabel("log10(y)");
print ("gaindb_coll.png", "-dpng");

f = linspace(0, 8, 91);
y = linspace(0, 10, 91);


for i=1:length(f)
  y(i) = v5f(f(i));
endfor

v5fmax = max(y);
lowerF = 0;
upperF = 0;
cond1 = true;
cond2 = true;

for i=1:length(f)
  if(abs(y(i)-v5fmax+3) < 2 && f(i) > 0.5)
    if(cond1 == true)
      lowerF = f(i);
      cond1 = false;
    endif
   endif
   
   if(abs(y(i)-v5fmax+3) < 1 && f(i) > 4)
    if(cond2 == true)
      upperF=f(i);
      cond2 = false;
    endif
   endif
endfor

v5fmaxlinear = 10^(v5fmax/20);

cost = (R1+R2+RS+Rc+Re+Rout+Rl)/10^3 + (Cin+Cb+Co)*10^(6);
merit = (10^(upperF)-10^(lowerF))*v5fmaxlinear/(10^(lowerF)*cost);

fprintf(fid4, "Lower Frequency & $%f$\\\\", 10^(lowerF));
fprintf(fid4, "Higher Frequency & $%f$\\\\", 10^(upperF));
fprintf(fid4, "Cost & $%f$\\\\", cost);
fprintf(fid4, "Merit & $%f$\\\\", merit);

close all

fclose(fid1);
fclose(fid2);
fclose(fid3);
fclose(fid4);
fclose(fid5);