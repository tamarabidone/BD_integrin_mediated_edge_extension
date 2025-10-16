% rng('shuffle')  % execute this line if you're adding new simulations after MATLAB was turned off and on.
clear;
clc;
SaveDirectory = '/uufs/chpc.utah.edu/common/home/bidone-group3/Remi/New_sims';    

% Setup parallel pool (8 cores)
poolobj = gcp('nocreate'); 
if isempty(poolobj)
    parpool(8); % Start a parallel pool with 8 workers
end

% Setup all combinations of parameters to be varied
values = [];
nRuns = 6;
k_a = [0.0001 ,0.0001, 0.001, 0.001, 0.01, 0.01, 0.1, 0.1];
peak = [1, 2, 1, 2, 1, 2, 1, 2];
nLigands = 400;

for m = 1:length(k_a)
    p = peak(m);
    k = k_a(m);
    for r = 1:nRuns
        values = [values; [k, p, nLigands, r]];
    end
end

% Run simulations in parallel
parfor k = 1:size(values,1) 
    ModelParameters = InitializeModelParameters();
    ModelParameters.TimeStep = 1e-4;
    ModelParameters.TotalSimulationTime = 300;
    ModelParameters.FilamentThermalMotionOn = false; 
    ModelParameters.CytoplasmViscosity = 1e5 * 0.0001;
    ModelParameters.VerticalOffSet = -200;
    ModelParameters.StartingNumberOfFilaments = 32;
    ModelParameters.AdhesionSpringConstant = values(k,1);
    ModelParameters.k_off_pointed = 7;
    ModelParameters.k_branch = 2.2;
    ModelParameters.FAL_connection_Distance = 10*2.75;
    ModelParameters.MaximumFilamentMass = 4000;
    ModelParameters.MolecularClutch_PeakNumber = values(k,2);

    SaveName = ['SIMULATION-001__','Ka_',SimFormat(values(k,1)),...
        '__Peak_',sprintf('%02d',values(k,2)),...
        '__nLigands_',sprintf('%04d',values(k,3)),...
        '_run_',sprintf('%02d',values(k,4)),'.mat'];

    if ~isfile(fullfile(SaveDirectory,SaveName))
        LamellipodiumModel_Bidone01(ModelParameters,SaveDirectory,SaveName);
    end
end

delete(gcp('nocreate')); 
quit;

    