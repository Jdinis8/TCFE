#files names
filename = "data.txt";
op       = "../sim/values.inc";
tex_file = "../doc/nodal_an_less_0.tex";
inc2     = "../sim/values2.inc";
inc3     = "../doc/ic_req.tex";
inc4     = "../sim/values3.inc";
inc5     = "../sim/values4.inc";

#opening files
fid        = fopen (filename, "r");
output     = fopen(op, "w");
output_tex = fopen(tex_file, "w");
output_1   = fopen(inc2, "w");
output_2   = fopen("../doc/ic_req.tex", "w");
output_3   = fopen(inc4, "w");
output_4   = fopen(inc5, "w");

#reading data
Kyu = textscan(fid, "%s %f");
[a, b] = Kyu{:};

rowasd = rows(a) - 2;
rowasb = rows(a) - 4;

for i = 1:rows(a)
  #printf("%s = %d\n",char(a(i)), b(i))
endfor

#printing data for ngspice
#         r1    r2     r3     r4     r5     r6     r7          vs          c    kb            kd
nodes  = {"2 1"; "2 3"; "2 5"; "0 5"; "5 6"; "0 CFP1"; "7 8"; "1 0"; "6 8"; "6 3 2 5"; "5 8 Vdumb"};
nodes2 = {"2 1"; "2 3"; "2 5"; "0 5"; "5 6"; "0 CFP1"; "7 8"; "1 0 0"; "6 8"; "6 3 2 5"; "5 8 Vdumb"};

for i = 1:rowasd
  fprintf(output, "%s %s %f\n", char(a(i)), nodes{i}, b(i));
endfor

fprintf(output, "Gb %s %f\n", nodes{rows(nodes)-1}, b(rows(b)-1));
fprintf(output, "Hd %s %f\n", nodes{rows(nodes)}, b(rows(b)));

for i = 1:rowasb
  fprintf(output_1, "%s %s %f\n", char(a(i)), nodes2{i}, b(i));
  fprintf(output_3, "%s %s %f\n", char(a(i)), nodes{i}, b(i));
  fprintf(output_4, "%s %s %f\n", char(a(i)), nodes{i}, b(i));
endfor

fprintf(output_1, "Vs 1 0 0\n");
fprintf(output_1, "Gb %s %f\n", nodes2{rows(nodes2)-1}, b(rows(b)-1));
fprintf(output_1, "Hd %s %f\n", nodes2{rows(nodes2)}, b(rows(b)));


fprintf(output_3, "Vs 1 0 0\n");
fprintf(output_3, "C1 %s %f\n", nodes{rows(nodes)-2}, b(rows(b)-2));
fprintf(output_3, "Gb %s %f\n", nodes{rows(nodes)-1}, b(rows(b)-1));
fprintf(output_3, "Hd %s %f\n", nodes{rows(nodes)}, b(rows(b)));


fprintf(output_4, "Vs 1 0 SINE(0 1.0 1000 0 0 0)\n");
fprintf(output_4, "C1 %s %f\n", nodes{rows(nodes)-2}, b(rows(b)-2));
fprintf(output_4, "Gb %s %f\n", nodes{rows(nodes)-1}, b(rows(b)-1));
fprintf(output_4, "Hd %s %f\n", nodes{rows(nodes)}, b(rows(b)));

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

fprintf(output_1, "V2 %s %f\n", nodes2{rows(nodes2)-2}, v(6)-v(8));

Ic = -(Ib + I5);

Req = (v(6)-v(8))/Ic;

tau = Req*C;

fprintf(output_2, "$I_c = %f \\rightarrow R_{eq} = %f \\rightarrow \\tau = %f$", Ic, Req, tau);

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

t = 0: 1e-6: 20e-3;

v6 = v3(6) + (v2(6)-v3(6))*exp(-t./(tau));


plot (t*1e3, v6)
xlabel("t [ms]");
ylabel("v_{6n}(t) [V]");
legend("v_{6n}", "vR");

print ("v6n.eps", "-depsc");

     
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
fclose(inc2);
fclose(inc3);
fclose(inc4);
fclose(inc5);