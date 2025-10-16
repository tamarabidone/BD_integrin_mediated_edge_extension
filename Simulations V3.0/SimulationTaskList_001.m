% rng('shuffle')  % execute this line if you're adding new simulations after MATLAB was turned off and on.
clear; 
clc;
SaveDirectory = '/uufs/chpc.utah.edu/common/home/bidone-group3/Remi/Simulations_V3.0/SaveDir';    
    
   % Setup all combinations of parameters to be varied  -------------------------------------------------------------------
    values = [];
    nRuns = 6; % Total number of runs for each condition
    k_a = 0.01; % adhesion spring constant
    peak = 1;   % WT or Mn
    nLigands = 100;

    for m = 1:length(k_a) % Create combinations of all conditions
        for n = 1:length(peak)
            for r = 1:nRuns
                for u = 1:length(nLigands)
                    values = [values; [k_a(m), peak(n), nLigands(u), r] ];
                end
            end
        end
    end
    % ---------------------------------------------------------------------------------------------------------------------
    
    
    for k = 1:size(values,1)
    %parfor k = 1:size(values,1) 
            ModelParameters = InitializeModelParameters;  % Initialize Default Model paramaeters
            ModelParameters.TimeStep = 1e-4;
            % Special Conditions (or edits to ModelParameters) ==============================================================
                ModelParameters.TotalSimulationTime = 90; % seconds
                ModelParameters.FilamentThermalMotionOn = false; 
                ModelParameters.CytoplasmViscosity = 1e5 * 0.0001; % Pascal * seconds
                ModelParameters.VerticalOffSet = -200;
                ModelParameters.StartingNumberOfFilaments = 32; % This is a ballpark value based on MaximumFilamentMass
                ModelParameters.AdhesionSpringConstant = values(k,1); % Change (substrate rigidity)
                ModelParameters.k_off_pointed = 7; % s^-1
                ModelParameters.k_branch = 2.2; % s^-1
                ModelParameters.FAL_connection_Distance = 10*2.75; % nm
                ModelParameters.MaximumFilamentMass = 4000; % monomers
                ModelParameters.MolecularClutch_PeakNumber = values(k,2); % Catch bond (1 = WildType or 2 = Manganese)
            % ===============================================================================================================
    
            SaveName = ['SIMULATION-001__','Ka_',SimFormat(values(k,1)),'__Peak_',sprintf('%02d',values(k,2)),'__nLigands_',sprintf('%04d',values(k,3)),'__run_',sprintf('%02d',values(k,4)),'.mat'];
        
            if ~isfile(fullfile(SaveDirectory,SaveName))  % If code is re-run it will only run values for files not created/completedsimulation yet (so be sure to delete old runs in the same folder)
            %     disp(SaveName);
                LamellipodiumModel_Bidone01(ModelParameters,SaveDirectory,SaveName); % Run simulation
            end
    end
    % disp('SIMULATION-001 complete') 

    quit 
    
    