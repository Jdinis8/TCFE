R1 = 1.03994439216 * 1000;
R2 = 2.07923431764 * 1000;
R3 = 3.06168544529 * 1000;
R4 = 4.09516986362 * 1000;
R5 = 3.00136467001 * 1000;
R6 = 2.03324628446 * 1000;
R7 = 1.02216788331 * 1000;

G1 = 1/R1;
G2 = 1/R2;
G3 = 1/R3;
G4 = 1/R4;
G5 = 1/R5;
G6 = 1/R6;
G7 = 1/R7;

Va = 5.03847501972; 
Id = 1.01674167773 * 0.001;
Kb = 7.01505323139 * 0.001;
Kc = 8.37372457746 * 1000;

     #V1    V2       V3   V4     V5        V6   V7     Ia   Ic
A = [G1   -G1       0    0      0         0    0      1     0;
     -G1  G1+G2+G3  -G2  0      -G3       0    0      0     0; 
     0    -G2-Kb    G2   0      Kb        0    0      0     0;
     0    0         0    G4+G6  -G4       0    -G6    -1    0;
     0    -G3       0    -G4    G3+G4+G5  -G5  0      0    -1;
     0    Kb        0    0      -G5-Kb    G5   0      0     0;
     0    0         0    -G6    0         0    G6+G7  0     0;
     1    0         0    -1     0         0    0      0     0;
     0    0         0    0      1         0    0      0   -Kc;];
     

b = [0; 0; 0; 0; -Id; Id; 0; Va; 0];

A = A*0.001;

v = A\b;

v = v*1000;

filename = "nodal_an.tex"

fid = fopen(filename, "w");

fprintf(fid, '\\begin{equation} \\Vec{b} = \\begin{bmatrix} %s \\end{bmatrix} \\label{eqsol} \\end{equation}', strjoin(cellstr(num2str(v(:))),'\\\\ '));

fclose(fid);


vc = v(5);

filename = "vc_nodal.tex"

fid = fopen(filename, "w");

fprintf(fid, '%d', vc);

fclose(fid);


vb = v(2) - v(5);

filename = "vb_nodal.tex"

fid = fopen(filename, "w");

fprintf(fid, '%d', vb);

fclose(fid);