clear all
close all

format long

filename1 = "../doc/ripple.tex";
filename2 = "../doc/merit.tex";
filename3 = "../doc/merit_corrected.tex";
filename4 = "../doc/rd.tex";
filename5 = "../doc/media.tex";

fid1       = fopen (filename1, "w");
fid2       = fopen (filename2, "w");
fid3       = fopen (filename3, "w");
fid4       = fopen (filename4, "w");
fid5       = fopen (filename5, "w");

period = 12;
n = 17.78988;
A = 230/n;
f=50;
#t0 = 0
t0 = 27500e-3;
t=linspace(t0, t0+period/f, 10000);
w=2*pi*f;
V_ON = 0.7;
vT = 0.025;
eta = 1;
N = 31;
Is  = 1e-14;
R = 200;
C = 0.00000133;

#Equivelant resistance in each diode
rd = eta*vT/(Is*exp((12.5/N)/(eta*vT)));

R1 = N*rd; #symbolizing the absence of resistor in parallel and the capacitor discharging from the diodes, very high resistance

v_in = 230*cos(w*t);
v = A*cos(w*t);

v_abs = zeros(1, length(t));
v_toff = zeros(1, length(t));
vo = zeros(1, length(t));
vC = zeros(1, length(t));
vL = zeros(1, length(t));
final_vo = zeros(1, length(t));

tOFF = (1/(2*w)) * atan(1/((2*w)*R1*C));

vC = A*cos(2*w*tOFF)*exp(-(t-t0-tOFF)/(R1*C));

##voltage after wave rectified (no need to simulate that part, pretty straightforward)
for i=1:length(t)
    v_abs(i) = abs(v(i));
endfor

umPeriodo = 1/(2*f);

for i=1:length(t)
  if t(i) < tOFF
    v_toff(i) = v_abs(i);
  elseif vC(i) > v_abs(i)
    v_toff(i) = vC(i);
  else 
    tOFF = tOFF + 1/f/2;
    vC = A*abs(cos(2*w*tOFF)) * exp(-(t-t0-tOFF)/R1/C);
    v_toff(i) = v_abs(i);
  endif
endfor


%restriction depending on the number of diodes and V_on
for i = 1:length(t)
  if (v_toff(i)>N*V_ON)
    vL(i) = N*V_ON;
    else
    vL(i) = v_toff(i);
  endif
endfor

%relation between alternating current component and diodes
%incremental resistor
for i = 1:length(t)
    vo(i) = N*rd/(N*rd+R) * (v_toff(i)-mean(v_toff));
endfor 

final_vo = vL + vo;

ripples = max(final_vo) - min(final_vo);

#if theoretical model is as the teacher asked
merit = 1/((R/1000 + C*10^6 + N*0.1)*(ripples + abs(mean(final_vo)-12) + 10^(-6)));

#if the theoretical model is approximately what ngspice gives
merit_corrected = 1/((R/1000 + C*10^6 + N*0.1)*(ripples + 10^(-6)));


fprintf(fid1, "%f", ripples);
fprintf(fid2, "%f", merit);
fprintf(fid3, "%f", merit_corrected);
fprintf(fid4, "%f", rd);
fprintf(fid5, "%f", mean(final_vo));

plot(t, v_toff);
ylim([12.91 12.94]);
legend("Voltage after Envelope Detector");
print ("antes_mat.png", "-dpng");

plot(t, final_vo);
ylim([12.90 12.94]);
legend("Output Voltage");
print ("zauzau_mat.png", "-dpng");

plot(t, final_vo-12);
ylim([0.9 0.94]);
legend("Output Voltage-12");
print ("deviation_mat.png", "-dpng");

plot(t, final_vo,t, v_in);
ylim([-230 230]);
legend("Output Voltage","Input Voltage");
print ("acdc.png", "-dpng"); 

close all

fclose(fid1);
fclose(fid2);
fclose(fid3);
fclose(fid4);
fclose(fid5);