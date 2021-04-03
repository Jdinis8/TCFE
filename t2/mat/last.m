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

vl(t) = vpa(Aze\vl(t));

ht = matlabFunction(vpa(vl(t)));

Coisoth = @(ze) ze(5) - ze(7);
Sixth = @(ze) ze(5);
First = @(ze) ze(1);

function retval = saca_argumento(x)
  if(abs(real(x)) < 1e-9)
    if(imag(x) >= 0) retval = pi/2;
    else retval = -pi/2;
    endif
  else retval = arg(x);
  if(retval > 0) retval = retval - 2*pi;
  endif
  endif
endfunction;

vscenas = @(omega) 180/pi*saca_argumento((First(ht(2*pi*power(10,omega), 0))));

vc = @(omega) 180/pi*saca_argumento((Coisoth(ht(2*pi*power(10,omega), 0))));
v6cenas = @(omega) 180/pi*saca_argumento((Sixth(ht(2*pi*power(10,omega), 0))));

ground = @(omega) 0;

fplot(vc, [-1, 6], 201);
hold on
fplot(v6cenas, [-1, 6], 201);
hold on
fplot(vscenas, [-1, 6], 201);
legend("\\phi(VC)", "\\phi(V6)","\\phi(VS)");
xlabel("log10(f) [f] = Hz");
ylabel("vc_{phi}(f) [V]");
h = findobj(gca, 'type', 'line');
set(h, 'LineWidth', 2)
print ("vcphi.eps", "-depsc");

clf

vscenas = @(omega) 20*log10(abs((First(ht(2*pi*power(10,omega), 0)))));
vc = @(omega) 20*log10(abs((Coisoth(ht(2*pi*power(10,omega), 0)))));
v6cenas = @(omega) 20*log10(abs((Sixth(ht(2*pi*power(10,omega), 0)))));

fplot(vc, [-1, 6], 201);
hold on
fplot(v6cenas, [-1, 6], 201);
hold on
fplot(vscenas, [-1, 6], 201);
legend("MAG(VC)", "MAG(V6)","MAG(VS)");
xlabel("log10(f) [f] = Hz");
ylabel("Mag [dB]");
h = findobj(gca, 'type', 'line');
set(h, 'LineWidth', 2)
print ("vcmag.eps", "-depsc");

fclose(fid);