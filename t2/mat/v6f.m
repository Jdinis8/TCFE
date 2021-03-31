
pkg load symbolic

filename = "data.txt";
fid        = fopen (filename, "r");
Kyu = textscan(fid, "%s %f");
[a, b] = Kyu{:};



#alinea 4

##syms vs(t)
##syms va(t)
##syms omega
##syms C
##syms G1
##syms G2
##syms G3
##syms G4
##syms G5
##syms G6
##syms G7
##syms Kd
##syms Kb

#assigning data
f = 1000; #1kHz
omega = pi*2*f;

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
C =  b(9)*0.001;
Kb = b(10);
Kd = b(11);


syms vs(t)
syms v(t)

vs(t) = -j*exp(j*sym('omega')*t);

Zc = 1/(j*omega*C);

    #V1    V2      V3    V4      V5       V6     V7      V8 
A = [1      0        0    -1       0        0      0      0;
    -G1  G1+G2+G3  -G2     0      -G3       0      0      0; 
     0    -G2-Kb    G2     0       Kb       0      0      0;
     0      0       0      1        0       0      0      0;
     0    -G3       0     -G4   G3+G4+G5   -G5-1/Zc    -G7    G7+1/Zc;
     0     Kb       0      0     -G5-Kb    G5+1/Zc      0     -1/Zc;
     0     0        0     -G6       0       0    G6+G7  -G7;
     0     0        0    -Kd*G6     1       0    Kd*G6   -1;]


v(t) = [vs(t); 0; 0; 0; 0; 0; 0; 0];

v(t) = A\v(t)

