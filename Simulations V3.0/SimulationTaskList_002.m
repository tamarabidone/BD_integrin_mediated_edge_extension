function SimulationTaskList_002(task_id)

SaveDirectory = '/uufs/chpc.utah.edu/common/home/bidone-group3/Remi/3repsvinculin';   

% -----------------------------
% Define all simulations
% -----------------------------
nRuns = 20;
k_a = [0.0001 0.0001 0.001 0.001 0.01 0.01];
peak = [1 2 1 2 1 2];
nLigands = 400;

values = [];
for m = 1:length(k_a)
    p = peak(m);
    k = k_a(m);
    for r = 1:nRuns
        values = [values; [k, p, nLigands, r]];
    end
end

totalSims = size(values,1);

% -----------------------------
% Map array task to simulation subset
% -----------------------------
% You can adjust nArrayTasks to match SLURM --array count
nArrayTasks = 15;  
taskID = task_id;

simIndices = round(linspace(1, totalSims+1, nArrayTasks+1));
idx_start = simIndices(taskID);
idx_end = simIndices(taskID+1) - 1;
simSubset = idx_start:idx_end;

% -----------------------------
% Start parallel pool if none exists
% -----------------------------
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool;  % Will use all available cores
end

% -----------------------------
% Run simulations in parallel
% -----------------------------
parfor k = simSubset
    rng(values(k,4), 'twister');

    ModelParameters = InitializeModelParameters();
    ModelParameters.TimeStep = 1e-4;
    ModelParameters.TotalSimulationTime = 90;
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
    ModelParameters.MyosinFraction = 0;

    SaveName = ['SIMULATION-001__','Ka_',SimFormat(values(k,1)),...
        '__Peak_',sprintf('%02d',values(k,2)),...
        '__nLigands_0400',...
        '_run_',sprintf('%02d',values(k,4)),'.mat'];

    if ~isfile(fullfile(SaveDirectory,SaveName))
        LamellipodiumModel_Bidone01(ModelParameters,SaveDirectory,SaveName);
    end
end

% Close pool when done
delete(gcp('nocreate'));
quit;
