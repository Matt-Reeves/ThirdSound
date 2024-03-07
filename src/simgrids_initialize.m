function simgrids = simgrids_initialize(simparams)

N = simparams.N; 

J = (0:N-1);
x0 = J*(2*pi/N);
y0 = -simparams.A*cos(x0);
phi0 = zeros(size(x0));

simgrids = struct;
simgrids.Q = [x0(:);y0(:);phi0(:)];
simgrids.v = [x0(2:N/2).'; y0(1:N/2+1).'; simparams.T; simparams.alpha ];
clear x0 y0 phi0

%The highest fourier mode must be supressed, for the reasons outlined by Roberts. 
k = ((-N/2:N/2-1)*2*pi/N);
k = fftshift(k);
k(N/2+1) =0; 
simgrids.k = k; 
clear k; 

return
