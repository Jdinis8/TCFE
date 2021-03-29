#files names
filename = "data.txt";
op       = "../sim/values.inc";
tex_file = "../doc/nodal_an_less_0.tex"

#opening files
fid        = fopen (filename, "r");
output     = fopen(op, "w");
output_tex = fopen(tex_file, "w");

#reading data
[a, b] = textread (filename, "%s %f");
for i = 1:rows(a)
  printf("%s = %d\n",char(a(i)), b(i))
endfor

#printing data for ngspice
#         r1    r2     r3     r4     r5     r6     r7     vs     c    kb            kd
nodes = {"1 2"; "2 3"; "2 5"; "0 5"; "5 6"; "0 7"; "7 8"; "0 1"; "6 8"; "6 3 2 5"; "5 8 Vdumb"};

for i = 1:rows(a)
  fprintf(output, "%s %s %f\n", char(a(i)), nodes{i}, b(i));
endfor

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

fprintf(output_tex, '\\begin{equation} \\Vec{b} = \\begin{bmatrix} %s \\end{bmatrix} \\label{eqsol} \\end{equation}', strjoin(cellstr(num2str(v(:))),'\\\\ '));

#closing files
fclose(op);
fclose(fid);
fclose(tex_file);