function MyosinApplication(task_id)

SaveDirectory = '/uufs/chpc.utah.edu/common/home/bidone-group3/Remi/Simulations_V1_DifferentClutch/Myosin_final';

%% Parameters
nRuns = 20;                      % Number of repeats per stiffness
k_a = [0.0001 0.01];             % Substrate stiffness
peak = [1 1];                     % Molecular clutch peak numbers
nLigands = 400;

%% Create all combinations
values = [];
for m = 1:length(k_a)
    for r = 1:nRuns
        values = [values; k_a(m), peak(m), nLigands, r];
    end
end

%% Determine subset of runs for this array task
totalRuns = size(values,1);
nTasks = str2double(getenv('SLURM_ARRAY_TASK_COUNT'));  % Total array tasks
if isempty(nTasks) || isnan(nTasks)
    nTasks = 8;  % fallback if not available
end

runsPerTask = ceil(totalRuns / nTasks);
startIdx = (task_id-1)*runsPerTask + 1;
endIdx   = min(task_id*runsPerTask, totalRuns);

valuesSubset = values(startIdx:endIdx, :);

%% Start parallel pool (uses up to allocated cores)
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool('local');  % MATLAB uses all allocated cores
end

%% Run simulations in parallel
parfor idx = 1:size(valuesSubset,1)
    v = valuesSubset(idx,:);
    rng(v(4), 'twister');  % reproducibility per run

    %% Initialize model parameters
    ModelParameters = InitializeModelParameters();
    ModelParameters.TimeStep = 1e-4;
    ModelParameters.TotalSimulationTime = 30;
    ModelParameters.FilamentThermalMotionOn = false; 
    ModelParameters.CytoplasmViscosity = 1e5 * 0.0001;
    ModelParameters.VerticalOffSet = -200;
    ModelParameters.StartingNumberOfFilaments = 32;
    ModelParameters.AdhesionSpringConstant = v(1);
    ModelParameters.k_off_pointed = 7;
    ModelParameters.k_branch = 2.2;
    ModelParameters.FAL_connection_Distance = 10*2.75;
    ModelParameters.MaximumFilamentMass = 4000;
    ModelParameters.MolecularClutch_PeakNumber = v(2);
    ModelParameters.ExponentialLifetime = 1;
    ModelParameters.ApplyMyosin = true;
    ModelParameters.ReducedFlow = false;

    %% Save filename
    SaveName = ['SIMULATION-001__Ks_', SimFormat(v(1)), ...
                'Myos_true__Peak_', sprintf('%02d', v(2)), ...
                '_run_', sprintf('%02d', v(4)), '.mat'];

    %% Run simulation if file does not exist
    if ~isfile(fullfile(SaveDirectory, SaveName))
        LamellipodiumModel_Bidone01(ModelParameters, SaveDirectory, SaveName);
    end
end

%% Clean up
delete(gcp('nocreate'));
quit;
