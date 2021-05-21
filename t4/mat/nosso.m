
##to calculate the DC component, we admit that the capacitor cin is blocking the current
##from vin, thus v_base is given by the bias circuit coming from vcc
##thus:

R1 = 80*10^3
R2 = 20*10^3
Vcc= 12


Req = R2*R1/(R1+R2)
Veq = -R2/(R1+R2)*Vcc

##lets now do the rest of the circuit

beta1 = 178.7
beta2 = 227.3
Vbeon = 0.7
Rc = 1*10^3
Re = 100
Rout = 100

##first, do thevenin equivalent for the rload and capacitor
  
   ##Ib1     Ie1  Ic1  Ib2       Ie2 Ic2
A = [1+beta1 -1    0    0        0   0;
     0        0    0    1+beta2  -1  0;
     Req      1    0    0        0   0;
     0        Re   Rc   Rc       0   0;
     1       -1    0    0        0   0;
     0        0    0    1       -1   1;]


v = [0; 0; -Veq-Vbeon; -Vcc-Vbeon; 0; 0]

A\v