close all
clear all

pkg load symbolic;

global R1;
global R2;
global C;
global w;
global tOFF;

%To be chosen
R1 = 1.e3;
R2 = 1.e3;
C = 0.0005;

%Write values to file
file_chosen_values = fopen("ChosenValues.tex","w");
fprintf(file_chosen_values,"$R_1$ & %f $\\Omega$ \\\\ \\hline\n$R_2$ & %f $\\Omega$ \\\\ \\hline\n$C$ & %f $F$ \\\\ \\hline", R1, R2, C);
fclose (file_chosen_values);

%Given values
f=50;
V1 = 230.;
V2 = 12.;
n=V1/V2;

  
%Envelope detector

%Time vector and angular frequency
t=linspace(0, 5/f, 1000);
w=2*pi*f;

%Voltage in transformer 2
v2 = V2 * cos(w*t);

%Full-wave rectified
vrec = zeros(1, length(t));

for i=1:length(t)
  vrec(i) = abs(v2(i));
endfor

%Final voltage after envelope detector
vOenv = zeros(1, length(t));

%tOFF and exponencial due to capacitor

tOFF = (1/w) * atan(1/(w*R1*C))
vOexp = V2*cos(w*tOFF)*exp(-(t-tOFF)/(R1*C));

%Function
function f_t_on = f_t_on(tON)
  global R1;
  global R2;
  global C;
  global w;
  global tOFF;
  f_t_on=cos(w*tON)-cos(w*tOFF)*exp(-(tON-tOFF)/(R1*C));
endfunction

%Derivative
function f_t_on_derivative = f_t_on_derivative(tON)
  global R1;
  global R2;
  global C;
  global w;
  global tOFF;
  f_t_on_derivative=-w*sin(w*tON)+(1/(R1*C))*cos(w*tOFF)*exp(-(tON-tOFF)/(R1*C));
endfunction

%Newton-Raphson
function f_t_on_solve = f_t_on_solve ()
  global R1;
  global R2;
  global C;
  global w;
  global tOFF;
  delta=1e-6;
  x_next=0.01; 
  do
    x=x_next;
    x_next=x-f_t_on(x)/f_t_on_derivative(x);
  until (abs(x_next-x)<delta)
  f_t_on_solve=x_next;
endfunction

%tON=f_t_on_solve()
tON=1/(2*f)

  %{
%Obtain the final voltage after the envelope detector
for i=1:length(t)
  if t(i) < tOFF
    vOenv(i) = vrec(i);
        elseif ((t(i) >= tOFF) & (t(i)<tOFF+tON))
    vOenv(i) = vOexp(i);
  else
      tOFF = tOFF + tON;
            vOexp = V2*abs(cos(w*tOFF))*exp(-(t-tOFF)/(R1*C));
      vOenv(i) = vrec(i);
        endif
endfor
%}


%Write tON and tOFF values to file
file_toff_ton = fopen("tOFF_tON.tex","w");
fprintf(file_chosen_values,"$t_{OFF}$ & %f $s$ \\\\ \\hline\n$t_{ON}$ & %f $s$ \\\\ \\hline", tOFF, tON);
fclose (file_toff_ton);

  
for i=1:length(t)
  if t(i) < tOFF
    vOenv(i) = vrec(i);
    elseif vOexp(i) > vrec(i)
      vOenv(i) = vOexp(i);
    else
      tOFF = tOFF + tON;
            vOexp = V2*abs(cos(w*tOFF))*exp(-(t-tOFF)/(R1*C));
      vOenv(i) = vrec(i);
   endif
endfor


%Values for voltage regulator
eta=1.;
IS=1.0e-14;
k=1.38064852e-23;
T=298.15;
q=1.60217662e-19;
VT=(k*T)/q;
n_diodes=18;
VD=mean(vOenv)/n_diodes;
rd=(eta*VT)/(IS*exp(VD/(eta*VT)))


vOUT = zeros(1, length(t));
vOUT_variations = zeros(1, length(t));
vo = zeros(1, length(t));

for i=1:length(t)
  vs=V2-vOenv(i);
        vo(i)=((n_diodes*rd)/(n_diodes*rd+R2))*vs;
        vOUT(i)=V2-vo(i);
        vOUT_variations(i)=vo(i);
endfor

%Plot

fig_voltages = figure ("Visible", "off");
title("Voltages in Envelope Detector and Voltage Regulator")
plot (t*1000, vrec, ";v_2 (rectified) (t);", t*1000, vOenv, ";v_O (envelope) (t);", t*1000, vOUT, ";v_{OUT} (t);");
xlabel ("t [ms]")
ylabel ("v [V]")
legend("Location", "southwest");
print (fig_voltages, "envelope.eps", "-depsc");

min_var=min(vOUT_variations);
max_var=max(vOUT_variations);

fig_variation = figure ("Visible", "off");
title("Deviation of final output voltage from 12V")
plot (t*1000, vOUT_variations, "r");
hleg=legend("v_{OUT}-12V (t)","Location","southwest");
ylim([min_var+0.5*(min_var-0.1) max_var+0.5*(max_var+0.1)]);
xlabel ("t [ms]")
ylabel ("v [V]")
print (fig_variation, "variations.eps", "-depsc");

%Merit
cost=R1*0.001+R2*0.001+C*10e-6+0.1*(n_diodes+4);
merit=1/(cost*(mean(vo)+mean(vo)+10e-6));
file_merit = fopen("Merit.tex","w");
fprintf(file_merit,"Cost & %f \\\\ \\hline\nMerit ($M$) & %f\\\\ \\hline", cost, merit);
fclose (file_merit);

close all;