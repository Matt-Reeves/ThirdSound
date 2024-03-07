
clear 
clc
%load ./Ac-1.1.mat
%load ./data/VdW/d0.2/Ac-0.27.mat
%load ./Ac-0.00155_quadrupleres.mat
load ./Ac0.665.mat

%simparams_initialize;
%simparams.N = 16;
simgrids = simgrids_initialize(simparams);

x = [0; sol(1:simparams.N/2-1); pi]; 
y = sol(simparams.N/2:simparams.N); 
phi = sol(end)*ones(1,simparams.N);
T = sol(end-1);
[x,y,phi] = imposeSymmetry(x.',y.',phi,simparams.N);

Q = [x(:); y(:); phi(:)];
simgrids.Q = Q;

simparams.T = T;
simparams.tfin = T;
simparams.DT     = simparams.tfin;                      %coarse time step
simparams.dt     = simparams.DT/1600/4;                  %fine time step
simparams.Nt     = (simparams.tfin/simparams.DT);       %coarse stepper
simparams.nt     = (simparams.DT/simparams.dt);         %fine stepper
simparams.t      = 0;                                   %current time
simparams.alpha = sol(end);

solution(1:3*simparams.N,1) = Q;
%%
for jj = 1:simparams.nt-1
    simgrids.Q = rk4(simgrids,simparams); 
    solution(1:3*simparams.N,jj+1) = simgrids.Q;
    
    quantities = conservedQuantities(simgrids,simparams);
    energy.E(jj)   = quantities.E;
    energy.T(jj)   = quantities.T;
    energy.V(jj)   = quantities.V;
    volumeFlux(jj) = quantities.volumeFlux;
    eta_bar(jj)    = quantities.Ybar;
end
%%
energyDeviation = max(abs(energy.E-energy.E(1))/energy.E(1))
heightDeviation = max(abs(eta_bar))
fluxDeviation = max(abs(volumeFlux))

%%

%%

filename = 'd0.9_.gif';

for jj = 64:64:simparams.nt-1
    figure(12)
    clf
    
    
    x = solution(1:simparams.N,jj).';
    y = solution(simparams.N+1:simparams.N*2,jj).';
    
    time = (1:jj)*simparams.dt/simparams.T;
    subplot(131)
    plot(time,solution(simparams.N+1,1:jj)/simparams.d)
    xlim([0 simparams.nt]*simparams.dt/simparams.T)
    ylim([-1 1])
    set(gca,'Fontsize',14)
    xlabel('Time, $t/T$','Interpreter','Latex','Fontsize',18)
    ylabel('Signal, $\eta(0,t))/d$','Interpreter','Latex','Fontsize',18)
    
    N2 = 3*simparams.N;
    xx = interpft(x -2*pi/simparams.N*(0:simparams.N-1),N2);
    xx = xx + 2*pi/N2*(0:N2-1);
    yy = interpft(y,N2);
     
    subplot(132)
    fill([xx pi 0],[yy/simparams.d -1 -1],'-b','Linewidth',1,'FaceAlpha',0.3)
    hold on 
    plot(x,y/simparams.d,'.k')
    ylim([-1 1 ])
    xlim([0 pi])
    set(gca,'Fontsize',14)
    xlabel('$2\pi x / L $','Interpreter','Latex','Fontsize',18)
    ylabel('$\eta(x)/d$','Interpreter','Latex','Fontsize',18)
    set(gca,'Xtick',[0 pi/2 pi])
    set(gca,'Xticklabel',{'0' '\pi/2', '\pi'})
    
    subplot(133)
    plot(time,energy.E(1:jj))
    hold on
    plot(time,energy.T(1:jj))
    plot(time, energy.V(1:jj))
    xlim([0 1])
    ylim([0 1.1*max(energy.E)])
    xlabel('Time, $t/T$','Interpreter','Latex','Fontsize',18)
    ylabel('Wave Energy','Interpreter','Latex','Fontsize',18)
    set(gca,'Fontsize',14)
    legend('$T+V$','$T$','$V$')
    set(legend,...
        'Interpreter','Latex',...
        'Fontsize',10,...
        'Location','North',...
        'Orientation','Horizontal')
    
    drawnow
    
    
    
    frame = getframe(12);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if jj == 64
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.0);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.0);
    end
    
end

%%

clear height omega
%Acvals = [-0.005:-0.0025:-0.145 -0.15:-0.0025:-0.7]
Acvals = -0.0025:-0.0025:-0.43
%Acvals = [-0.0025:-0.0025:(-0.03+0.0025) -0.045:-0.0025:-0.45 -0.4605:-0.0025:-0.5005]
%Acvals = 0.005:0.005:0.785
Acvals = -0.0005:-0.0005:-0.0165
Acvals = [0.01:0.005:0.525 0.55:0.005:0.565 0.6:0.005:0.665]
N2 = 2*simparams.N; 
%simparams.N = 16;
for jj = 1:length(Acvals)
    %load(['./data/VdW/d0.3/Ac' num2str(Acvals(jj)) '.mat']);
    load(['./Ac' num2str(Acvals(jj)) '.mat']);
    x = [0; sol(1:simparams.N/2-1); pi]; 
    y = sol(simparams.N/2:simparams.N); 
    phi = sol(end)*ones(1,simparams.N);
    T = sol(end-1);
    [x,y,phi] = imposeSymmetry(x.',y.',phi,simparams.N);
    xx = interpft(x -2*pi/simparams.N*(0:simparams.N-1),N2);
    xx = xx + 2*pi/N2*(0:N2-1);
    yy = interpft(y,N2);
    height(jj) = (max(yy)-min(yy))/2;
    height2(jj) = (y(1)-y(end/2+1))/2;
    omega(jj)  = 2*pi/T;
end

figure
subplot(121)
plot(Acvals,omega/simparams.omega0,'.')
subplot(122)
plot(height,omega/simparams.omega0,'.')

%%
load ./Ac-0.4275.mat
x = [0; sol(1:simparams.N/2-1); pi]; 
y = sol(simparams.N/2:simparams.N); 
phi = sol(end)*ones(1,simparams.N);
T = sol(end-1);
[x,y,phi] = imposeSymmetry(x.',y.',phi,simparams.N);

N2 = 4*simparams.N;
    xx = interpft(x -2*pi/simparams.N*(0:simparams.N-1),N2);
    xx = xx + 2*pi/N2*(0:N2-1);
    yy = interpft(y,N2);
     
    figure
    fill([xx pi 0],[yy/simparams.d -1 -1],'-b','Linewidth',1,'FaceAlpha',0.3)
    hold on 
    plot(x,y/simparams.d,'.k')
    ylim([-1 1 ])
    xlim([0 pi])
    set(gca,'Fontsize',14)
    xlabel('$2\pi x / L $','Interpreter','Latex','Fontsize',18)
    ylabel('$\eta(x)/d$','Interpreter','Latex','Fontsize',18)
    set(gca,'Xtick',[0 pi/2 pi])
    set(gca,'Xticklabel',{'0' '\pi/2', '\pi'})

%%
function eta_bar = meanHeight(simgrids,simparams)
N = simparams.N;
k = simgrids.k;
x = simgrids.Q(1:N).'; y = simgrids.Q(N+1:2*N).';
xF = 1i*k.*fft(x - (0:N-1)*(2*pi/N)); xF(1) = 2*pi; dxdj = real(ifft(xF));
eta_bar = sum(y.*dxdj);
end
