function quantities = conservedQuantities(simgrids,simparams)

x = simgrids.Q(1:simparams.N).';
y = simgrids.Q(simparams.N+1:simparams.N*2).';
phi = simgrids.Q(simparams.N*2+1:simparams.N*3).';
[x,y,phi] = imposeSymmetry(x,y,phi,simparams.N);

Z=x+1i*y;

ZF = fft(Z - (0:simparams.N-1)*(2*pi/simparams.N));

ZD = 1i*simgrids.k.*ZF; 
ZD(1) = 2*pi;
ZD = ifft(ZD);

ZDD = -simgrids.k.^2.*ZF; 
ZDD = ifft(ZDD);

%derivative of phi
dphi_dj = ifft(1i*simgrids.k.*fft(phi));

%Matrices for vortex strengths
K = matrixD(Z,ZD,ZDD,simparams.N,simparams.d);

%solve vortex strengths and derivatives
a = real(K\dphi_dj.').';
ad = ifft(1i*simgrids.k.*fft(a));

%calculate velocities
w = velocity(a,ad,Z,ZD,ZDD,simparams.N,simparams.d);
u = real(w); v = -imag(w);

Xd = real(ZD); Yd = imag(ZD);

%Kinetic Energy
T = 1/4/pi*sum( phi.*(-Yd.*u +Xd.*v));

volumeFlux = 1/2/pi*sum(v.*Xd - u.*Yd);
%Potential Energy
switch simparams.waveType
    case 'VdW'
        V = 1/4/pi*simparams.d^4/3*sum((1./(simparams.d+ imag(Z)).^2 - 1./simparams.d^2).*Xd);
    case 'gravity'
        V = 1/4/pi*sum(imag(Z).^2.*Xd);
end

Es = simparams.kappa/2/pi*(sum(sqrt(Xd.^2 + Yd.^2)) - 2*pi);
Ybar = 1/2/pi*sum(imag(Z).*Xd);

quantities.T = T;
quantities.V = V;
quantities.Es = Es;
quantities.E = T+V+Es;
quantities.Ybar = Ybar;
quantities.volumeFlux = volumeFlux;

return