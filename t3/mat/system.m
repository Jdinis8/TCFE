#Asd

#voltage limit = 0.7



n = 14;
A = 230/n;
f=50;
t0 = 200e-3
t=linspace(t0, t0+10/f, 1000);
w=2*pi*f;
vO = A*sin(w*t);

#input voltage
plot(t*1000, vO)
title("Input voltage v_o(t)")
xlabel ("t[ms]")
ylabel ("v_o[V]")
print ("vo.png", "-dpng");

#full wave rectifier
t=linspace(t0, t0+10/f, 1000);
vO = A*abs(sin(w*t));

plot(t*1000, vO)
title("Fullwave rectified voltage v_f(t)")
xlabel("t[ms]")
ylabel("v_f[V]")
print ("vf.png", "-dpng");

close all 
clear all

#envelope and rectifier 
R = 90*10^(3);
C = 110*10^(-6);

n = 14;
A = 230/n;
f=50;
t0 = 200e-3;
t=linspace(t0, t0+10/f, 1000);
w=2*pi*f;
vS = A*abs(sin(w*t));
vOhr = zeros(1, length(t));
vO = zeros(1, length(t));

tOFF = 1/w * atan(1/w/R/C);

vOnexp = A*cos(w*tOFF)*exp(-(t-tOFF)/R/C);

figure
for i=1:length(t)
  if (vS(i) > 0)
    vOhr(i) = vS(i);
  else
    vOhr(i) = 0;
  endif
endfor

plot(t*1000, vOhr)
hold

for i=1:length(t)
  if t(i) < tOFF
    vO(i) = vS(i);
  elseif vOnexp(i) > vOhr(i)
    vO(i) = vOnexp(i);
  else 
    vO(i) = vS(i);
  endif
endfor

plot(t*1000, vO)
title("Output voltage v_o(t)")
xlabel ("t[ms]")
legend("rectified","envelope")
print ("venvlope.png", "-dpng");


#knowing that diode votlage is 0.7, we use the closest limiter

##VON = 0.7;
##nvlim = ceil(12/VON)+1;
##vlim = nvlim*VON


##
##%limit
##for i=1:length(t)
##  if vO(i) > vlim
##    vO(i) = vlim;
##  elseif vO(i) < -vlim
##    vO(i) = -vlim;
##  endif
##endfor
##
##plot(t*1000, vO)
##title("Limited voltage v_f(t)")
##xlabel("t[ms]")
##ylabel("v_limited[V]")
##print ("vlimited.png", "-dpng");