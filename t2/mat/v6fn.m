clear all

pkg load symbolic

#files names
filename = "data.txt";

#opening files
fid        = fopen (filename, "r");

#reading data
Kyu = textscan(fid, "%s %f");
[a, b] = Kyu{:};

#assigning data
R1 = b(1);
R2 = b(2);
R3 = b(3);
R4 = b(4);
R5 = b(5);
R6 = b(6);
R7 = b(7);

G1 = 1/R1;
G2 = 1/R2;
G3 = 1/R3;
G4 = 1/R4;
G5 = 1/R5;
G6 = 1/R6;
G7 = 1/R7;

Vs = b(8);
C =  b(9);
Kb = b(10);
Kd = b(11);

#alinea 1
#solving nodal analysis
     #V1    V2      V3    V4      V5       V6     V7      V8 
A = [1      0        0    -1       0        0      0      0;
    -G1  G1+G2+G3  -G2     0      -G3       0      0      0; 
     0    -G2-Kb    G2     0       Kb       0      0      0;
     0      0       0      1        0       0      0      0;
     0    -G3       0     -G4   G3+G4+G5   -G5    -G7    G7;
     0     Kb       0      0     -G5-Kb    G5      0      0;
     0     0        0     -G6       0       0    G6+G7  -G7;
     0     0        0    -Kd*G6     1       0    Kd*G6   -1;];

    
v = [Vs; 0; 0; 0; 0; 0; 0; 0];
v = A\v;

#alinea 2
#aqui o valor de V6-V8 = cte.

     #V1    V2      V3    V4      V5       V6     V7      V8 
A2 = [1      0        0    0       0        0      0      0;
    -G1  G1+G2+G3  -G2     0      -G3       0      0      0; 
     0    -G2-Kb    G2     0       Kb       0      0      0;
     0      0       0      1        0       0      0      0;
     0    -G3+Kb    0    -G4    G4+G3-Kb    0     -G7    G7;
     0     0        0      0        0       1      0     -1;
     0     0        0     -G6       0       0    G6+G7  -G7;
     0     0        0    -Kd*G6     1       0    Kd*G6   -1;];
      
v2 = [0; 0; 0; 0; 0; v(6)-v(8); 0; 0];

v2 = A2\v2;

Ib = Kb*(v2(2)-v2(5));
I5 = (v2(5) - v2(6))/R5;
Ic = -(Ib + I5);
Req = (v(6)-v(8))/Ic;
tau = Req*C;

#alinea 3
% time vector
t = 0: 1e-6: 20e-3;

     #V1    V2      V3    V4      V5       V6     V7      V8 
A3= [1      0        0     -1       0       0      0      0;
    -G1  G1+G2+G3  -G2     0      -G3       0      0      0; 
     0    -G2-Kb    G2     0       Kb       0      0      0;
     0      0       0      1        0       0      0      0;
     0    -G3       0     -G4   G3+G4+G5   -G5    -G7    G7;
     0     Kb       0      0     -G5-Kb    G5      0      0;
     0     0        0     -G6       0       0    G6+G7  -G7;
     0     0        0    -Kd*G6     1       0    Kd*G6   -1;];
    

v3 = [0; 0; 0; 0; 0; 0; 0; 0];

v3 = A3\v3;


############################
G1 = sym(G1, 'r');
G2 = sym(G2, 'r');
G3 = sym(G3, 'r');
G4 = sym(G4, 'r');
G5 = sym(G5, 'r');
G6 = sym(G6, 'r');
G7 = sym(G7, 'r');

C =  sym(C, 'r');
Kb = sym(Kb, 'r');
Kd = sym(Kd, 'r');

f = 1000; #1kHz
omega = pi*2*f;
omega = sym(omega, 'r');

syms vs(x)
syms vl(x)
syms Zc

vs(x) = -j*exp(j*omega*x);

Zc = 1/(j*omega*C);

one = sym(1);

notone = sym(-1);
  #   V1    V2      V3     V5       V6     V7      V8 
Aze = [sym(1) 0        0     0        0      0      0;
    -G1  G1+G2+G3  -G2    -G3       0      0      0; 
     0    -G2-Kb    G2     Kb       0      0      0;
     0    -G3       0   G3+G4+G5 -G5-j*C*omega -G7    G7+j*C*omega;
     0     Kb       0   -G5-Kb    G5+j*C*omega  0    -j*C*omega;
     0     0        0      0        0     G6+G7  -G7;
     0     0        0     sym(1)    0     Kd*G6 sym(-1);];
     
const = 1000;
const = sym(const, 'r');

Aze = Aze*const;

vl(x) = [vs(x)*const; 0; 0; 0; 0; 0; 0];

vl(x) = Aze\vl(x);

ht = matlabFunction(vpa(vl(x)));
ze = ht(0)

############################

omega = pi*2*f;

t = 0: 1e-6: 20e-3;
v6 = v3(6) + (v2(6)-v3(6))*exp(-t./(tau));
t = 0: 1e-6: 20e-3;
v6f = abs(ze(5))*cos(angle(ze(5))+omega*t);

v6tots = v6f + v6;

plot (t*1e3, v6tots)
xlabel("t [ms]");
ylabel("v_{6f}+v_{6n}(t) [V]");
print ("v6fn.eps", "-depsc");

x = -5e-3: 1e-6: 20e-3;
x2 = 1e-6: 1e-6: 20e-3;
I = zeros(size(x));
cond1 = x<0;
cond2 = x == 0;
cond3 = x > 0;
I(cond1) = v(6);
I(cond2) = v2(6);
I(cond3) = v3(6) + (v2(6)-v3(6))*exp(-x2./(tau)) + abs(ze(5))*cos(angle(ze(5))+omega*x2);

x3 = 0: 1e-6: 20e-3;
I2 = zeros(size(x));
cond1 = x<0;
cond2 = x >= 0;
I2(cond1) = Vs;
I2(cond2) = sin(omega*x3);

plot (x*1e3,I), grid
hold on
plot (x*1e3,I2)
xlabel("t [ms]");
ylabel("v_{6f}+v_{6n}(t) [V]");
print ("v6totsze.eps", "-depsc");



#closing files
fclose(fid);