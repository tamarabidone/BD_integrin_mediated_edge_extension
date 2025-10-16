function Reducedflow(task_id)

%% Setup directories
SaveDirectory = '/uufs/chpc.utah.edu/common/home/bidone-group3/Remi/Simulations_V1_DifferentClutch/Reduced_flow';

%% Add current folder and subfolders to path (for parfor workers)
addpath(genpath('/uufs/chpc.utah.edu/common/home/bidone-group3/Remi/Simulations_V1_DifferentClutch'));

%% Parameters
nRuns = 6;
k_a = [0.0001 0.001 0.01];

% Simulation conditions: columns = [ExponentialLifetime, ReducedFlow]
conditions = [
    2, false;
    3, false;
    1, true
];

nLigands = 400;

%% Generate all combinations
values = [];
for m = 1:length(k_a)
    for c = 1:size(conditions,1)
        for r = 1:nRuns
            values = [values; k_a(m), conditions(c,1), conditions(c,2), nLigands, r];
        end
    end
end

%% Start parallel pool if none exists
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool;  % uses all available cores
end

%% Run simulations in parallel
parfor idx = 1:size(values,1)
    rng(values(idx,5), 'twister');  % reproducibility

    %% Initialize model parameters
    ModelParameters = InitializeModelParameters();
    ModelParameters.TimeStep = 1e-4;
    ModelParameters.TotalSimulationTime = 30;
    ModelParameters.FilamentThermalMotionOn = false; 
    ModelParameters.CytoplasmViscosity = 1e5 * 0.0001;
    ModelParameters.VerticalOffSet = -200;
    ModelParameters.StartingNumberOfFilaments = 32;
    ModelParameters.AdhesionSpringConstant = values(idx,1);
    ModelParameters.k_off_pointed = 7;
    ModelParameters.k_branch = 2.2;
    ModelParameters.FAL_connection_Distance = 10*2.75;
    ModelParameters.MaximumFilamentMass = 4000;
    ModelParameters.MolecularClutch_PeakNumber = 1;
    ModelParameters.ExponentialLifetime = values(idx,2);
    ModelParameters.ApplyMyosin = false;
    ModelParameters.ReducedFlow = values(idx,3);

    %% Convert logical ReducedFlow to string for filename
    if values(idx,3)
        reducedFlowStr = 'true';
    else
        reducedFlowStr = 'false';
    end

    %% Save name
    SaveName = ['SIMULATION-001__Ks_', SimFormat(values(idx,1)), ...
                '__ExpLife_', sprintf('%02d', values(idx,2)), ...
                '__ReducedFlow_', reducedFlowStr, ...
                '_run_', sprintf('%02d', values(idx,5)), '.mat'];

    %% Run simulation if not already saved
    if ~isfile(fullfile(SaveDirectory, SaveName))
        LamellipodiumModel_Bidone01(ModelParameters, SaveDirectory, SaveName);
    end
end

%% Clean up
delete(gcp('nocreate')); 
quit;
