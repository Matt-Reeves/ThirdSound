
function [R] = nonlinearProblem(v,simparams,simgrids)
%global figureHandle axesHandle

xx = [0; v(1:simparams.N/2-1); pi].';
yy = v(simparams.N/2:simparams.N).';
phi = zeros(1,simparams.N);
[xx,yy,phi] = imposeSymmetry(xx,yy,phi,simparams.N);


%figure(1); hold on
%plot(xx(1:(simparams.N)/2+1),yy(1:(simparams.N)/2+1),'.-')
%drawnow

T = v(end-1);
alpha = v(end);

simgrids.Q = [xx(:); yy(:); phi(:)];
simgrids.v = v;

simparams.T = T; 
simparams.tfin   = simparams.T*0.25; 
simparams.DT     = simparams.tfin;
simparams.dt     = simparams.DT/300;
simparams.Nt     = (simparams.tfin/simparams.DT);
simparams.nt     = (simparams.DT/simparams.dt);

simparams.alpha = alpha; 

[residuals] = calculateResiduals(simgrids,simparams);

%J = calculateJacobian(simgrids,simparams,residuals).';

R = residuals.all;

end

function residuals = calculateResiduals(simgrids,simparams)

eta_bar = meanHeight(simgrids,simparams);
Ac = getCrestAcceleration(simgrids,simparams);
%x-residuals

tempgrid = simgrids;
for jj = 1:simparams.Nt
    for kk = 1:simparams.nt
    tempgrid.Q = rk4(tempgrid,simparams);
%     figure(999)
%     clf
%     plot(tempgrid.Q(1:simparams.N),tempgrid.Q((1:simparams.N)+simparams.N),'.-')
%     ylim([-2 2]*simparams.A)
%     drawnow
    end

end

%% Residuals at T/4 guess
N = simparams.N; 
Q = tempgrid.Q;
x = Q(1:N).'; y = Q(N+1:2*N).'; phi = Q(2*N+1:3*N).';

resx = x(2:N/2) - 2*pi/N*(1:N/2-1);

%y-residuals
j = 1:N/4; 
resy = y(j) - y(N/2+2-j);

%phi-residuals
j = 1:N/4+1;
resphi = phi(j) + phi(N/2+2-j) - simparams.alpha;

%constructs the residuals vector of the N+2 unknowns
residuals.all = [eta_bar; Ac-simparams.Act; resx.'; resy.'; resphi.' ];
residuals.eta_bar = eta_bar;
residuals.Ac = Ac;
residuals.resx = resx.'; 
residuals.resy = resy.';
residuals.phi = resphi.';
residuals.SSE = sum(residuals.all.^2);


end

function eta_bar = meanHeight(simgrids,simparams)
N = simparams.N;
k = simgrids.k;
x = simgrids.Q(1:N).'; y = simgrids.Q(N+1:2*N).';

xF = 1i*k.*fft(x - (0:N-1)*(2*pi/N)); 
xF(1) = 2*pi; 
dxdj = ifft(xF);

eta_bar = sum(y.*dxdj);
end

% function [x,y,phi] = unpack(Q)
%  N = length(Q)/3;
%  x = Q(1:N).'; 
%  y = Q(N+1:2*N).';
%  phi = Q(2*N+1:3*N).'; 
% end
% 
% function Q = pack(x,y,phi)
% Q = [x(:);y(:);phi(:)]; 
% end

%function Jacobi = calculateJacobian(simgrids,simparams,residuals)              % Internal Jacobian calculation from fsovle seems to perform much
                                                                                %better, so no point using  this... 
%{
N = simparams.N;

[x0,y0,phi0] = unpack(simgrids.Q); 

%% x-part of the Jacobian
%disp('Calculating X-Jacobian')
change = 0.01;
%change = sqrt(eps);
x = x0;
delta = change*2*pi/N;
xvals = x(2:N/2);
Jacobi_x = zeros(N/2-1,N+2);
old_residuals = residuals.all;

parfor kk = 1:length(xvals)
   
    xtemp = xvals;
    xtemp(kk) = xtemp(kk)+delta;
    x = [0, xtemp, pi, 2*pi-fliplr(xtemp)]; y = y0; phi = phi0;
    [x,y,phi] = imposeSymmetry(x,y,phi,N);
    
    temp = simgrids;
    temp.Q = pack(x,y,phi);
   
    residuals_new = calculateResiduals(temp,simparams);
    Jacobi_x(kk,:) = (residuals_new.all - old_residuals)/delta;
end
%% y-part of the Jacobian
%disp('Calculating Y-Jacobian')

Jacobi_y = zeros(N/2+1,N+2);
A = max(y0);
parfor (kk = 1:N/2+1,8)
    delta = change*A;
    x= x0; y = y0; phi = phi0;
    y(kk) = y(kk) + delta; 
    [x,y,phi] = imposeSymmetry(x,y,phi,N);
    temp = simgrids;
    temp.Q = pack(x,y,phi);
    residuals_new = calculateResiduals(temp,simparams);
    Jacobi_y(kk,:) = (residuals_new.all - old_residuals)/delta;
end
%% T-part
x= x0; y =y0; phi = phi0;
simgrids.Q = pack(x,y,phi); 
temp = simparams; 
delta = change*temp.tfin; 
temp.tfin = temp.tfin + delta;
temp.DT = temp.tfin;
temp.dt = temp.DT/500;
temp.nt = (temp.DT/temp.dt);
temp.Nt = (temp.tfin/temp.DT);
residuals_new = calculateResiduals(simgrids,temp);
Jacobi_T(1,1:N+2) = (residuals_new.all - old_residuals)/delta;
Jacobi_T(1,1:2) = 0; %changing T should not change eta_bar or Ac
%% alpha part

%delta = max(sqrt(eps),change*simparams.alpha);
delta = change*simparams.alpha;
temp = simparams;
temp.alpha = temp.alpha + delta; 

residuals_new = calculateResiduals(simgrids,temp);
Jacobi_alpha(1,1:N+2) = (residuals_new.all - old_residuals)/delta;
Jacobi_alpha(1,1:2) = 0; %changing alpha should not change eta_bar or Ac
Jacobi = [Jacobi_x; Jacobi_y; Jacobi_T; Jacobi_alpha];


end
%}


