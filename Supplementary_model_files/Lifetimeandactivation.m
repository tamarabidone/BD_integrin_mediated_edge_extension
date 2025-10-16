function Lifetimeandactivation(task_id)

SaveDirectory = '/uufs/chpc.utah.edu/common/home/bidone-group3/Remi/Simulations_V1_DifferentClutch/outputs_final';

%% Parameters
nRuns = 6;
k_a = [0.0001 0.001 0.01];             % substrate stiffness
peakValues = [1 3 4];                   % MolecularClutch_PeakNumber
k_int_values = [0.01 0.1 1];            % integrin activation rate
nLigands = 400;

%% Create all combinations
values = [];
for m = 1:length(k_a)
    for p = 1:length(peakValues)
        for kint = 1:length(k_int_values)
            for r = 1:nRuns
                values = [values; k_a(m), peakValues(p), k_int_values(kint), nLigands, r];
            end
        end
    end
end

%% Determine subset of runs for this array task
totalRuns = size(values,1);
nTasks = str2double(getenv('SLURM_ARRAY_TASK_COUNT'));  % total array tasks
if isempty(nTasks) || isnan(nTasks)
    nTasks = 8;  % fallback
end

runsPerTask = ceil(totalRuns / nTasks);
startIdx = (task_id-1)*runsPerTask + 1;
endIdx   = min(task_id*runsPerTask, totalRuns);
valuesSubset = values(startIdx:endIdx, :);

%% Start parallel pool (uses all allocated cores)
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool('local');  % will use all available cores
end

%% Run simulations in parallel
parfor idx = 1:size(valuesSubset,1)
    v = valuesSubset(idx,:);
    rng(v(5), 'twister');  % seed for reproducibility

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
    ModelParameters.Adhesion_ActivationRate = v(3);
    ModelParameters.ExponentialLifetime = 1;
    ModelParameters.ApplyMyosin = false;
    ModelParameters.ReducedFlow = false;

    %% Save filename
    SaveName = ['SIMULATION-004__Ks_', SimFormat(v(1)), ...
                '__Peak_', sprintf('%02d', v(2)), ...
                '__kInt_', strrep(num2str(v(3)), '.', 'p'), ...
                '_run_', sprintf('%02d', v(5)), '.mat'];

    %% Run simulation if file does not exist
    if ~isfile(fullfile(SaveDirectory, SaveName))
        LamellipodiumModel_Bidone01(ModelParameters, SaveDirectory, SaveName);
    end
end

%% Clean up
delete(gcp('nocreate'));
quit;
