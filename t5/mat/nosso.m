clear all
close all

format long

pkg load symbolic

#filename1 = "../doc/vce.tex";

#fid1 = fopen (filename1, "w");

R0 = 5*10^3;
R1 = 1000;
R2 = 100*10^3;
RL = 8;
C1 = 220*10^(-9);
C2 = 220*10^(-9);

syms Vin
syms um
syms G1
syms G2
syms GL
syms GC1(w)
syms acres(w)
syms vout1(w)
syms vout(w)

G0 = sym(1/R0, 'f');
G1 = sym(1/R1, 'f');
G2 = sym(1/R2, 'f');
GL = sym(1/RL, 'f');
G4 = sym(R2/R1,'f');
Vin = sym(1, 'f');
um = sym(1, 'f');
GC1(w) = j*w*sym(C1, 'f');
GC2(w) = j*w*sym(C2, 'f');


 #v1                      #v2                   #v3         
B(w) = [GC1(w)+G0        0                     0;
            0            -GL              GL+GC2(w);
        -um-G4           1                    0;]

vb(w) = [Vin*GC1(w); 0; 0];

acres = B(w)\vb(w);

vout1 = acres(2)-acres(3);

vout(w) = vout1

w=1000
subs(vpa(abs(vout)),w)

abs(vout(w))

#fclose(fid1);