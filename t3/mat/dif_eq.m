clear all
close all

format long

period = 20;
n = 18.41351;
A = 230/n;
f=50;
t0 = 200e-3
t=linspace(t0, t0+period/f, 2000);
w=2*pi*f;
V_ON = 0.7;
vT = 0.025;
eta = 1;
N = 31;
Is  = 1e-14;
R = 200;
C = 0.00000133;
R1 = 1e8;

v = A*cos(w*t);

v_abs = zeros(1, length(t));
v_toff = zeros(1, length(t));
vo = zeros(1, length(t));
vC = zeros(1, length(t));
vL = zeros(1, length(t));
final_vo = zeros(1, length(t));

#R1->infinity so atan->0
#tOFF = (1/(2*w)) * atan(1/((2*w)*R1*C))
tOFF = (1/(2*w)) * atan(1/((2*w)*R1*C));

#vC = A*cos(2*w*tOFF)*exp(-(t-tOFF)/(R1*C));
vC = A*cos(2*w*tOFF)*exp(-(t-tOFF)/(R1*C));

##voltage after wave rectified (no need to simulate that part, pretty straightforward)
for i=1:length(t)
    v_abs(i) = abs(v(i));
endfor

#tOFF = 0, but we left this for loop for any tOFF changes
for i=1:length(t)
  if t(i) < tOFF
    v_toff(i) = v_abs(i);
  elseif vC(i) > v_abs(i)
    v_toff(i) = vC(i);
  else 
    tOFF = tOFF + 1/f/2;
    vC = A*abs(cos(2*w*tOFF)) * exp(-(t-tOFF)/R1/C);
    v_toff(i) = v_abs(i);
  endif
endfor

ripples = max(v_toff) - min(v_toff)

%%Regulating our voltage
rd = eta*vT/(Is*exp(V_ON/(eta*N*vT)));

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
  if v_toff(i) >= N*V_ON
    vo(i) = N*rd/(N*rd+R) * (v_toff(i)-mean(v_toff));
  else
    vo(i) = v_toff(i) - mean(v_toff);
  endif
endfor 

final_vo = vL + vo;    

final1 = zeros(1, length(t)-400);
final2 = zeros(1, length(t)-400);
final3 = zeros(1, length(t)-400);

j = 1;

for i = 1:length(t)
  if t(i) > 400e-3
    final1(j) = v_abs(i);
    final2(j) = v_toff(i);
    final3(j) = final_vo(i);
    j = j+1;
  endif
endfor

tfinal = linspace(t0 + 400e-3, t0 + period/f + 400e-3, length(t)-400);

plot(tfinal*1000, final1, tfinal*1000, final2, tfinal*1000, final3-12);
xlim([600 725])
ylim([0 17]);
legend("Rectified","Envelope", "Regulator - 12");
print ("varios.png", "-dpng");

plot(tfinal*1000, final3-12);
xlim([600 725]);
ylim([0.3 0.65]);

print ("sozinho.png", "-dpng");