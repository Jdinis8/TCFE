clear all
close all

format long

period = 10;
n = 17.78988;
A = 230/n;
f=50;
#t0 = 0
t0 = 27600e-3
t=linspace(t0, t0+period/f, 10000);
w=2*pi*f;
V_ON = 0.7;
vT = 0.025;
eta = 1;
N = 31;
Is  = 1e-14;
R = 200;
C = 0.00000133;
R1 = 1e8; #symbolizing the absence of resistor in parallel

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

#Regulating our voltage
rd = eta*vT/(Is*exp((12.5/N)/(eta*vT)))

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

ripples = max(final_vo) - min(final_vo)



plot(t, v_toff);
ylim([12.9275 12.929]);
legend("Voltage after Envelope Detector");
print ("antes_mat.png", "-dpng");

plot(t, final_vo);
ylim([12.927 12.9295]);
legend("Output Voltage");
print ("zauzau_mat.png", "-dpng");

plot(t, final_vo-12);
ylim([12.927-12 12.9295-12]);
legend("Output Voltage-12");
print ("deviation_mat.png", "-dpng");

plot(t, final_vo,t, v_in);
ylim([-230 230]);
legend("Output Voltage","Input Voltage");
print ("acdc.png", "-dpng"); 