clear all
close all

format long

pkg load symbolic

#filename1 = "../doc/vce.tex";

#fid1 = fopen (filename1, "w");

#ZI = 4*10^6;
#ZO = 40;
R1 = 909.0909;
R2 = 200*10^3;
R3 = 2000;
R4 = 50000/3;
RL = 8;
C1 = 220*10^(-9);
C2 = 0.1061093248*10^(-6);
#g_opamp = 10^5;

syms R1s
syms R2s
syms R3s
syms R4s
syms RLs
syms ZIs
syms ZOs
syms As
syms Vins
syms ums
syms vb(w)
syms GC1(w)
syms GC2(w)
syms Vout(w)
syms vout1(w)
syms VI(w)
syms V2(w)

R1s = sym(R1, 'f');
R2s = sym(R2, 'f');
R3s = sym(R3,'f');
R4s = sym(R4,'f');
RLs = sym(RL, 'f');
#ZIs = sym(ZI,'f');
#ZOs = sym(ZO,'f');
#As  = sym(g_opamp,'f');
Vins = sym(1, 'f');
ums = sym(1, 'f');
GC1(w) = j*w*sym(C1, 'f');
GC2(w) = j*w*sym(C2, 'f');

VI(w) = R1s/(R1s+1/GC1(w))

V2(w) = (1+R2/R3)*VI(w)

#B(w) = [ 1/R4s+GC2(w)+1/RLs;]

#vb(w) = [V2(w)/R4s];

#acres(w) = B(w)\vb(w)

#h = matlabFunction(abs(vpa((acres(10^w)))))

#first = @(w) w(1);


Vout(w) = R1*V2(w)/((1+RLs*GC2(w))*(R4s+R1s/(1+RLs*GC2(w))))


h = matlabFunction(abs(vpa((Vout(10^w)))))

vout1 = @(w) 20*log10(h(w))

#vout1 = @(w) 20*log10(first(h(w)))

fplot(vout1, [1, 8], 750);
legend("Gain_{db}");
xlabel("log10(f) [f] = Hz");
ylabel("Gain (dB)");
print ("gain_banda.png", "-dpng");

w_central = 1000; #vai ser melhor medido

syms Z_in
syms Z_out

Z_in = sym(R1+1/(j*w_central*C1), 'f');
Z_out = sym((j*w_central*C2+1/R4+1/RL)^(-1), 'f');

vpa(Z_in)
vpa(Z_out)

vpa(abs(Z_in))
vpa(abs(Z_out))


#fclose(fid1);