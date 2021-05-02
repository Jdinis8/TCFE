#using newton raphson method

pkg load symbolic

##n = 14;
##A = 230/n;
##f=50;
##t0 = 200e-3
##t=linspace(t0, t0+10/f, 1000);
##w=2*pi*f;

syms w
syms A
syms v0(t)


##voltage after wave rectified (no need to simulate that part, pretty straightforward)
vO = abs(A*sin(w*t))

##now we have the top node v4 = v0, and the bottom node is gnd, so it's zero

##diode equation
syms nd #number of diodes = 31
#eta = 1, so no need to declare it
syms Is #Is = 1
syms VT #Vt = ?
syms Id(v) #current diode

Id(v) = Is*(exp(v/(nd*VT)) - 1)

#current in capacitor is pretty simple: its a phasor and the nodes
#are the same as the of the diodes, so actually we have 

syms ic(t)
syms C

ic(t) = C*A*w*sin(w*t)

#because dv0/dt * C = ic, then we can approximate the derivative by
#delta's and reach a good approximation of the next iteration
#v0+ - v0- = ic*delta(t)/C
#i.e., v0+ = v0- + ic*delta(t)/C,
#where the delta is something of our choosing

syms vnext(vbefore)
syms h

vnext(vbefore) = vbefore + ic(t)*h/C

#then what is left is just the equation for the diode which can be written as

syms f(v5)
syms R

f(v5) = vnext(vbefore) - v5 - R*Id(v5)

#in our calculations, everything will have values, except for v5, which is what
#we want to find. Thus, we use newton-raphsons method to find it.
