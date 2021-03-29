#files names
filename = "data.txt";
op       = "../sim/values.inc";

#opening files
fid    = fopen (filename, "r");
output = fopen(op, "w");

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
Vs = b(8);
C =  b(9);
Kb = b(10);
Kd = b(11);

#solving nodal analysis
A = [R4+R3+R1 -R3 -R4;
     R4 0 R4+R6+R7-Kd;
     Kb*R3 1-Kb*R3 0;];
    
v = [-Vs; 0; 0];

A\v

#closing files
fclose(op);
fclose(fid);