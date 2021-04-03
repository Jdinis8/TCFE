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


syms vs(t)
syms vl(t)
syms Zc
syms omega

vs(t) = -j*exp(j*omega*t);

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


vl(t) = [vs(t); 0; 0; 0; 0; 0; 0];

vl(t) = vpa(Aze\vl(t))

ht = matlabFunction(vpa(vl(t)))

Coisoth = @(ze) ze(5) - ze(7);

vc = @(omega) angle(Coisoth(ht(2*pi*omega, 0)))

#~fplot(vc, [1, 1000000], 201)

x = logspace(-8, 6)
y = vc(x)
semilogx(x,y)
grid on

##Coisoth = @(t) t(5) - t(7);
##
##
##angle(vc)
##norm(vc)

##f = 1000; #1kHz
##omega = pi*2*f;
##
##
##
##ht = matlabFunction(vpa(vl(x)));
##ze = ht(0)




#closing files
fclose(fid);