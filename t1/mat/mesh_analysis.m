R1 = 1.03994439216*1000;
R2 = 2.07923431764*1000; 
R3 = 3.06168544529*1000;
R4 = 4.09516986362*1000;
R5 = 3.00136467001*1000;
R6 = 2.03324628446*1000;
R7 = 1.02216788331*1000;
Va = 5.03847501972;
Id = 1.01674167773*0.001;
Kb = 7.01505323139*0.001;
Kc = 8.37372457746*1000;


A = [R1+R4, -1/Kb, -R4; -R4, 0, R4+R6+R7-Kc; R3, 1/Kb-R3, 0]
b = [-Va; 0; 0]

A = A*0.001;
b = b;

v = A\b;

%v = v*1000*1000;

filename = "mesh_an.tex"

fid = fopen(filename, "w");

fprintf(fid, '\\begin{bmatrix}  %s  \\end{bmatrix}', strjoin(cellstr(num2str(v(:))),'\\\\ '));

fclose(fid);


filename = "vc.tex"

vc = Kc*v(3)*0.001;

fid = fopen(filename, "w");

fprintf(fid, '%d', vc);

fclose(fid);


filename = "vb.tex"

vb = (v(2)-v(1))*R3;

fid = fopen(filename, "w");

fprintf(fid, '%d', vb*0.001);

fclose(fid);