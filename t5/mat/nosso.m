clear all
close all

format long

pkg load symbolic


fid1  = fopen ("../doc/mat_low_cof.tex" , "w");
fid2  = fopen ("../doc/mat_high_cof.tex" , "w");
fid3  = fopen ("../doc/mat_f_central.tex" , "w");
fid4  = fopen ("../doc/mat_gain_central.tex" , "w");
fid5  = fopen ("../doc/mat_Z_in_real.tex" , "w");
fid6  = fopen ("../doc/mat_Z_in_imag.tex" , "w");
fid7  = fopen ("../doc/mat_Z_out_real.tex" , "w");
fid8  = fopen ("../doc/mat_Z_out_imag.tex" , "w");
fid9  = fopen ("../doc/mat_Z_in_abs.tex" , "w");
fid10 = fopen ("../doc/mat_Z_out_abs.tex" , "w");
fid11 = fopen ("../doc/mat_merit.tex" , "w");
fid12 = fopen ("../doc/mat_merit_sem_opAMP.tex" , "w");


#ZI = 4*10^6;
#ZO = 40;
R1 = 10^3;
R2 = 150*10^3;
R3 = 10^3;
R4 = 10^3;
C1 = 220*10^(-9);
C2 = 0.1061093248*10^(-6);
#g_opamp = 10^5;

syms R1s
syms R2s
syms R3s
syms R4s
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
#ZIs = sym(ZI,'f');
#ZOs = sym(ZO,'f');
#As  = sym(g_opamp,'f');
Vins = sym(1, 'f');
ums = sym(1, 'f');
GC1(w) = j*2*sym(pi)*w*sym(C1, 'f');
GC2(w) = j*2*sym(pi)*w*sym(C2, 'f');

VI(w) = R1/(R1+1/GC1(w));

V2(w) = (1+R2/R3)*VI(w);

#B(w) = [ 1/R4s+GC2(w)+1/RLs;]

#vb(w) = [V2(w)/R4s];

#acres(w) = B(w)\vb(w)

#h = matlabFunction(abs(vpa((acres(10^w)))))

#first = @(w) w(1);


Vout(w) = V2(w)/(1+R4*GC2);


h = matlabFunction(abs(vpa((Vout(10^w)))))

vout1 = @(w) 20*log10(h(w))

#vout1 = @(w) 20*log10(first(h(w)))

fplot(vout1, [1, 8], 750);
legend("Gain_{db}");
xlabel("log10(f) [f] = Hz");
ylabel("Gain (dB)");
print ("gain_banda.png", "-dpng");



f = linspace(0, 8, 750);
y = linspace(0, 8, 750);


for i=1:length(f)
  y(i) = vout1(f(i));
endfor

v_max = max(y)

low = -1
high = -1

for i=1:length(y)
	if (y(i) >= v_max-3 && low<0)
	  low = 10^f(i)
	endif
	if (y(i) <= v_max-3 && low>0 && high<0)
	  high = 10^f(i)
	endif
endfor



f_central = sqrt(high*low)
gain_central = vout1(log10(f_central))

syms Z_in
syms Z_out

Z_in = vpa(sym(R1+1/(j*f_central*C1), 'f'));
Z_out = vpa(sym((j*f_central*C2+1/R4)^(-1), 'f'));

Z_in_real = vpa(real(Z_in))
Z_in_imag = vpa(imag(Z_in))

Z_out_real = vpa(real(Z_out))
Z_out_imag = vpa(imag(Z_out))

Z_in_abs = vpa(abs(Z_in))
Z_out_abs = vpa(abs(Z_out))

cost_opAMP = 13322.792;
cost_coisas = 0.22+1+1+300+1+3.44

merit = 1/(10^(-6)+(cost_opAMP+cost_coisas)*abs(40-gain_central)*abs(1000-f_central))
merit_clean = 1/(10^(-6)+(cost_coisas)*abs(40-gain_central)*abs(1000-f_central))

fprintf(fid1, "%s", char(vpa(low,6)));
fprintf(fid2, "%s", char(vpa(high,6)));
fprintf(fid3, "%s", char(vpa(f_central,6)));
fprintf(fid4, "%s", char(vpa(gain_central,6)));
fprintf(fid5, "%s", char(vpa(Z_in_real,6)));
fprintf(fid6, "%s", char(vpa(Z_in_imag,6)));
fprintf(fid7, "%s", char(vpa(Z_out_real,6)));
fprintf(fid8, "%s", char(vpa(Z_out_imag,6)));
fprintf(fid9, "%s", char(vpa(Z_in_abs,6)));
fprintf(fid10, "%s", char(vpa(Z_out_abs,6)));
fprintf(fid11, "%s", char(vpa(merit,6)));
fprintf(fid12, "%s", char(vpa(merit_clean,6)));


fclose(fid1);
fclose(fid2);
fclose(fid3);
fclose(fid4);
fclose(fid5);
fclose(fid6);
fclose(fid7);
fclose(fid8);
fclose(fid9);
fclose(fid10);
fclose(fid11);
fclose(fid12);