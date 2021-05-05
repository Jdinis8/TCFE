#using newton raphson method

pkg load symbolic

syms  n
syms  A
syms  f
syms  t0
syms  t
syms  w
syms  V_ON
syms  vT
syms  eta
syms  N
syms  Is
syms  R
syms  nC



n = 18.41351;
A = 230/n;
f=50;
t0 = 200e-3
t=linspace(t0, t0+10/f, 1000);
w=2*pi*f;
V_ON = 0.6;
vT = 0.025;
eta = 1;
N = 31;
Is  = 1e-14;
R = 200;
C = 0.00000133

syms w
syms A
syms v0(t)


##voltage after wave rectified (no need to simulate that part, pretty straightforward)
syms vR = abs(A*sin(w*t));

#R1->infinity so atan->0
#tOFF = (1/(2*w)) * atan(1/((2*w)*R1*C))
syms tOFF = 0

#vC = A*cos(2*w*tOFF)*exp(-(t-tOFF)/(R1*C));
syms vC = A*cos(2*w*tOFF)

syms vL

if (vC>N*V_ON)
	vL = N*V_ON;
	else
	vL = vC;
endif

syms vD = vL/N;

syms rd = eta*vT/(Is*exp(vD/(eta*vT)));

syms oscilation = N*rd/(N*rd+R);

syms vO(t) = oscilation*cos(w*t);