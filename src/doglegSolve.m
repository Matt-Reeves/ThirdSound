function doglegSolve(N,d,Ac,waveType,delta_Ac,Acmax)


%Designed to run on the custer, so input arguments must be given as strings
%and converted to the appropriate types afterwards. 

%Problem setup and initialization
N = str2double(N);          % # of grid points 
d = str2double(d);          % fluid depth
Ac = str2double(Ac);        % crest acceleration: + for VdW, - for gravity
Acmax = str2double(Acmax);
delta_Ac = str2double(delta_Ac);

%waveType = 'VdW';   % must be "VdW" or "gravity"

simparams = simparams_initialize(N,d,waveType,Ac);
simgrids = simgrids_initialize(simparams);

%Setup initial guess for the vector of unknowns, v0. 
x = simgrids.Q(1:simparams.N);
y = simgrids.Q(simparams.N+1:2*simparams.N);
v0 = [x(2:simparams.N/2); y(1:simparams.N/2+1); simparams.T; simparams.alpha];
%% Solve problem for linear initial guess

%options for fsolve()
options = optimoptions('fsolve',...
    'Algorithm','trust-region-dogleg',...
    'Display','iter',...
    'SpecifyObjectiveGradient',false,...
    'UseParallel',true,...
    'OptimalityTolerance',1e-12,...
    'FunctionTolerance',1e-20,...
    'MaxFunctionEvaluations',30*simparams.N);

%% Solve the first iteration using linear initial guess
f = @(v) nonlinearProblem(v,simparams,simgrids);
[sol,FVAL,EXITFLAG,OUTPUT] = fsolve(f,v0,options);

%Make save directory and save
savedir  = ['../data/' simparams.waveType '/d' num2str(simparams.d) '_N' num2str(simparams.N)  '/'];
if exist(savedir,'dir') == 0
    mkdir(savedir)
end
savefile = ['Ac' num2str(simparams.Act) '.mat'];
save([savedir savefile],'sol','simparams','simgrids','FVAL','EXITFLAG','OUTPUT');

%% Map out new solutions incrementing through crest/trough acceleration and using previous solution as initial guess 

for vals = (Ac+delta_Ac):delta_Ac:Acmax  
    
    simparams.Act = vals;
    v0 = sol;
    x = [0; v0(1:N/2-1); pi].';
    y = v0(N/2:N).';
    phi = zeros(1,N);
    %[x,y,phi] = imposeSymmetry(x,y,phi,N);
    simparams.T = v0(end-1);
    simparams.alpha = v0(end);
    
    f = @(v) nonlinearProblem(v,simparams,simgrids);
    [sol,FVAL,EXITFLAG,OUTPUT] = fsolve(f,v0,options);

    savefile = ['Ac' num2str(simparams.Act) '.mat'];
    save([savedir savefile],'sol','simparams','simgrids','FVAL','EXITFLAG','OUTPUT');
end
return
