function simparams = simparams_initialize(N,d,waveType,Ac)
simparams = struct; 

%number of points
simparams.N = N;

%fluid parameters
simparams.waveType = waveType;                             %'gravity' or 'VdW'
simparams.d        = d;                                    %fluid depth
simparams.kappa    = 0;                                    %surface tension (ignoring this for now) 

%residual-related parameters.
if     strcmp(simparams.waveType,'VdW'),     simparams.Act   = +Ac;                               %target crest acceleteration
elseif strcmp(simparams.waveType,'gravity'), simparams.Act   = -Ac;                
end
simparams.alpha = -1e-4;                                %potential offset 

%time-related parameters
simparams.omega0 = sqrt(tanh(simparams.d));             %linear frequency
simparams.T      = 2*pi/simparams.omega0;               %linear period
simparams.tfin   = simparams.T*0.25;                    %evolve T/4 
simparams.DT     = simparams.tfin;                      %coarse timestep
simparams.dt     = simparams.DT/600;                    %fine timestep
simparams.Nt     = (simparams.tfin/simparams.DT);       %coarse stepper
simparams.nt     = (simparams.DT/simparams.dt);         %fine stepper
simparams.t      = 0;                                   %current time 

%Guess wave height for starting point, (linearized approximation)
simparams.A     =  simparams.Act/simparams.omega0.^2;

return




