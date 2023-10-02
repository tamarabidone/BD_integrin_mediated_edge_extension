clc
clear
clc

%This is the folder where we save the outputs from this analysis 
FileSaveDirectory = '/Users/remisondaz/Desktop/MATLAB/Varying_spring_constant/Velocity';

%This folder contains raw model output files after running simulations code (for example: SimulationTaskList_001.m)
RawSaveDirectory = '/Users/remisondaz/Desktop/MATLAB/Varying_spring_constant/Raw_files';


% Setup all combinations of parameters to be varied  -------------------------------------------------------------------
    values = [];
    nRuns = 3;
    nIntegrins = 100;
    k_a =  [0.0001 0.001 0.01 0.1 1]; % adhesion spring constant
    peak = [1 2];   % WT or Mn
    ligands =  400; %Different number of ligands

    Velocity_Raw = cell(nRuns,1);
    Adhesion_data = cell(nIntegrins, 5,10001);
    Membrane_pos = cell(nRuns, 1);
    Velocity =  NaN(nRuns,1);
   
    nTotal = numel(Velocity);
    index  = 0;

    % Create combinations of all conditions

     for m = 1:length(ligands) % Create combinations of all conditions
        for n = 1:length(peak)
            for k = 1:length(k_a)
                runs_with_same_conditions = 0;
                 for r = 1:nRuns
                    SaveName = ['SIMULATION-001__','Ka_',SimFormat(k_a(k)),'__Peak_',sprintf('%02d',peak(n)), '__nLigands_',sprintf('%04d',ligands(m) ), '__run_', sprintf('%02d',r), '.mat'];
                    if exist(fullfile(RawSaveDirectory, SaveName), 'file') == 2
                        runs_with_same_conditions = runs_with_same_conditions + 1;
                        load(fullfile(RawSaveDirectory,SaveName), '-mat')
                        Velocity_Raw{r,1} = SimData.MemVelocity;
                    end
                 end

                 if runs_with_same_conditions == 3
                 index = index + 1;
                 % disp([num2str(index),' of ',num2str(nTotal)])
                 save(fullfile(FileSaveDirectory, [sprintf('%01d',peak(n)), '_', SimFormat(k_a(k)), '_', num2str(ligands(m)),'.mat']), 'Velocity_Raw');
                 end
            end
        end
     end
