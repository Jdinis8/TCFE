clear all 
close all

pkg load symbolic

#files names
filename = "data.txt";
op       = "../sim/values.inc";
tex_file = "../doc/nodal_an_less_0.tex"

#opening files
fid        = fopen (filename, "r");
output     = fopen(op, "w");
output_tex = fopen(tex_file, "w");

#reading data
Kyu = textscan(fid, "%s %f");
[a, b] = Kyu{:};

for i = 1:rows(a)
  #printf("%s = %d\n",char(a(i)), b(i))
endfor

#printing data for ngspice
#         r1    r2     r3     r4     r5     r6     r7     vs     c    kb            kd
nodes = {"1 2"; "2 3"; "2 5"; "0 5"; "5 6"; "0 7"; "7 8"; "0 1"; "6 8"; "6 3 2 5"; "5 8 Vdumb"};

for i = 1:rows(a)
  fprintf(output, "%s %s %f\n", char(a(i)), nodes{i}, b(i));
endfor

#printing data for .tex
mtx_res = {"V1"; "V2"; "V3"; "V4"; "V5"; "V6"; "V7"; "V8"};

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

     #V2        V3  V5     V7
B = [-G1         0 -G4    -G6;
      G1+G2+G3 -G2 -G3      0;
      0          0   0  G6+G7;
      -G2-Kb    G2  Kb      0;];
      
v2 = [0; 0; G7*v(8); 0];

v2 = B\v2;

Ib = Kb*(v2(1)-v2(3));
I5 = (v(6) - v2(3))/R5;

Ic = Ib + I5

Req = (v(6)-v(8))/Ic


#alinea 3
% time vector
t = 0: 1e-6: 20e-3;

C = C*0.001;

v6 = (v(6)-v(8))*exp(-t./(Req*C));


plot (t*1e3, v6)
xlabel("t [ms]");
ylabel("v_{6n}(t) [V]");
legend("v_{6n}", "vR");

print ("v6n.eps", "-depsc");

#alinea 4

t = 0: 1e-6: 20e-3;

f = 1000; #1kHz

omega = pi*2*f;

vs = -j*exp(j*omega*t);


Zc = 1/(j*omega*C);

##G1 = sym('G1');
##G2 = sym('G2');
##G3 = sym('G3');
##G4 = sym('G4');
##G5 = sym('G5');
##G6 = sym('G6');
##G7 = sym('G7');
##Kd = sym('Kd');
##Kb = sym('Kb');


     #V1    V2      V3    V4      V5       V6     V7      V8 
A2= [1      0        0    -1       0        0      0      0;
    -G1  G1+G2+G3  -G2     0      -G3       0      0      0; 
     0    -G2-Kb    G2     0       Kb       0      0      0;
     0      0       0      1        0       0      0      0;
     0    -G3       0     -G4   G3+G4+G5   -G5-1/Zc    -G7    G7+1/Zc;
     0     Kb       0      0     -G5-Kb    G5+1/Zc      0     -1/Zc;
     0     0        0     -G6       0       0    G6+G7  -G7;
     0     0        0    -Kd*G6     1       0    Kd*G6   -1;];



for i = 1:rows(vs)
  #printf("%s = %d\n",char(a(i)), b(i))
endfor
     
va = [-j*exp(j*omega*t); 0; 0; 0; 0; 0; 0; 0];
##va = A2\va;

#pretty(va);

#real(va(0))

     
#print out the data for the nodal analysis

#variables
fprintf(output_tex, '\\begin{equation} \n\\begin{bmatrix}\n');
for i = 1:rows(mtx_res)
  fprintf(output_tex, "%s \\\\", mtx_res{i});
endfor
fprintf(output_tex, '\\end{bmatrix} = \n');

#results
fprintf(output_tex, '\\begin{bmatrix}\n');
for i = 1:rows(v)
  fprintf(output_tex, "%f \\\\", v(i));
endfor
fprintf(output_tex, '\\end{bmatrix} V\n');

#closing matrix
fprintf(output_tex, '\\label{eqsol}\n \\end{equation}');

#fprintf(output_tex, '\\begin{equation} \\Vec{b} = \\begin{bmatrix} %s \\end{bmatrix} \\label{eqsol} \\end{equation}', strjoin(cellstr(num2str(v(:))),'\\\\ '));

#closing files
fclose(op);
fclose(fid);
fclose(tex_file);