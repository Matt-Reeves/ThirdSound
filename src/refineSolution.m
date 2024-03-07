function [sol,simparams] = refineSolution(loaddir,loadfile,N2)
load([loaddir loadfile])

N2 = str2double(N2);
N = simparams.N;
disp(['Solution has ' num2str(N) 'points...'])
x = [0; sol(1:N/2-1); pi].';
y = sol(N/2:N).';
phi = zeros(1,N);
[x,y,~] = imposeSymmetry(x,y,phi,N);

xx = interpft(x -2*pi/simparams.N*(0:simparams.N-1),N2);
xx = xx + 2*pi/N2*(0:N2-1);
yy = interpft(y,N2);
phi = zeros(1,N2);
[x,y,~] = imposeSymmetry(xx,yy,phi,N2);

sol = [x(2:N2/2).'; y(1:N2/2+1).'; sol(end-1); sol(end)  ];
v0 = sol;

simparams = simparams_initialize(N2,simparams.d,simparams.waveType,simparams.Act);

simparams.T = sol(end-1);
simparams.alpha = sol(end); 

simgrids = simgrids_initialize(simparams);

disp(['Refining to ' num2str(simparams.N) 'points...'])

options = optimoptions('fsolve',...
    'Algorithm','trust-region-dogleg',...
    'Display','iter',...
    'SpecifyObjectiveGradient',false,...
    'UseParallel',true,...
    'OptimalityTolerance',1e-12,...
    'FunctionTolerance',1e-20,...
    'MaxFunctionEvaluations',30*simparams.N);

f = @(v) nonlinearProblem(v,simparams,simgrids);
[sol,FVAL,EXITFLAG,OUTPUT] = fsolve(f,v0,options);
savefile = ['Ac' num2str(simparams.Act) '.mat' ]
save([loaddir '../N' num2str(N2) '/' savefile],'sol','simparams','simgrids','FVAL','EXITFLAG','OUTPUT');

return

