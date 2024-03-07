function Ac = getCrestAcceleration(simgrids,simparams)

Q = simgrids.Q; 
k = simgrids.k;
N = simparams.N;
d = simparams.d; 

%% First Step

temp = simparams;
temp.dt = temp.dt/100;
Q = rk4(simgrids,temp);

x = Q(1:N).';  
x(1) = 0; 
x(N/2+1) = pi; 
x(N/2+2:N) = 2*pi - x(N/2:-1:2);

y = Q(N+1:2*N).'; 
y(N/2+2:N) = fliplr(y(2:N/2));

phi = Q(2*N+1:3*N).'; 
phi(N/2+2:N) = fliplr(phi(2:N/2));

Z = x+1i*y;
ZF = fft(Z - (0:N-1)*(2*pi/N));
ZD = 1i*k.*ZF; 
ZD(1) = 2*pi;
ZD = ifft(ZD);

ZDD = -k.^2.*ZF; 
ZDD = ifft(ZDD);

%derivative of phi
dphi_dj = ifft(1i*k.*fft(phi));

%Matrices for vortex strengths
K = matrixD(Z,ZD,ZDD,N,d);
%solve vortex strengths and derivatives
a = real(K\dphi_dj.').';
ad = ifft(1i*k.*fft(a));
%calculate velocities
w = velocity(a,ad,Z,ZD,ZDD,N,d);
vdt = -imag(w(1));

%% Second Step
temp2 = simgrids;
temp2.Q = Q;
Q = rk4(temp2,temp);

x = Q(1:N).'; 
x(1) = 0; 
x(N/2+1) = pi; 
x(N/2+2:N) = 2*pi - x(N/2:-1:2);
y = Q(N+1:2*N).'; 
y(N/2+2:N) = fliplr(y(2:N/2));
phi = Q(2*N+1:3*N).'; 
phi(N/2+2:N) = fliplr(phi(2:N/2)); 

Z = x+1i*y;

ZF = fft(Z - (0:N-1)*(2*pi/N));
ZD = 1i*k.*ZF; 
ZD(1) = 2*pi;
ZD = ifft(ZD);
ZDD = -k.^2.*ZF; 
ZDD = ifft(ZDD);

%derivative of phi
dphi_dj = ifft(1i*k.*fft(phi));

%Matrices for vortex strengths
K = matrixD(Z,ZD,ZDD,N,d);
%solve vortex strengths and derivatives
a = real(K\dphi_dj.').';
ad = ifft(1i*k.*fft(a));
%calculate velocities
w = velocity(a,ad,Z,ZD,ZDD,N,d);
v2dt = -imag(w(1));


Ac = (4*vdt-v2dt)/temp.dt/2;


end