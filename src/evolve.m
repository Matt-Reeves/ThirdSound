function dQdt = evolve(t,Q,k,N,d)
%Derivatives
 
x = Q(1:N).';  
x(1) = 0; 
x(N/2+1) = pi; 
x(N/2+2:N) = 2*pi - x(N/2:-1:2);

y = Q(N+1:2*N).'; 
y(N/2+2:N) = fliplr(y(2:N/2));
%y = project(y);

phi = Q(2*N+1:3*N).'; 
phi(N/2+2:N) = fliplr(phi(2:N/2));
%phi = project(phi);

if (any(diff(x) <0) )
    error('Non-monatonic x. Aborting.')
end

%x(1) = 0; x(N/2+1) = pi; x(N/2+2:end) = 2*pi-fliplr(x(2:N/2));

Z = x+1i*y;
ZF = fft(Z - (0:N-1)*(2*pi/N));

kz = fftshift((-N/2:N/2-1)*2*pi/N);
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
u = real(w); v = -imag(w);

%Equation of motion for phi. 
dphi_dt = d/3*(1./(1+y/d).^3 - 1) + (u.^2 + v.^2)/2;
%dphi_dt = -y + (u.^2 + v.^2)/2;

u(1) = 0; u(N/2+1) = 0;
dQdt = [u.';v.';dphi_dt.'];

%Tony's versions -- verified they are the same 
%  B = matrixdd(Z,ZD,ZDD,N,d,0,0);
%  vx = real(B\dphi_dj.').';
%  vxd = ifft(1i*k.*fft(vx));
%  zq = velocity_old(Z,ZD,ZDD,N,vx,vxd,d,0,0);
%  uu =  real(zq(1,:));
%  vv = -imag(zq(1,:));
 
end

function x = project(x)
 p = ones(size(x));
 %N = length(x);
 %The analysis in Roberts shows that the highest mode needs to be zeroed
 %for the method to be stable. In the VdW problem, we find that sometimes
 %the next highest modes also need to be zeroed for stability ... 
 %k = -N/2:N/2-1;
 %p = exp(-(k/(N/3)).^36);
 p(1) =0; %p(2) = 0; p(end) =0; % p(3) = 0; p(end-1) = 0;
 p = fftshift(p);
 x = real(ifft(p.*fft(x)));
end
