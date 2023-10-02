clear all
close all
clc

%This is the folder where we save the outputs from this analysis 
FileSaveDirectory = '/Users/remisondaz/Desktop/MATLAB/Varying_spring_constant/Membrane_spread';

%This folder contains raw model output files after running simulations code (for example: SimulationTaskList_001.m)
RawSaveDirectory = '/Users/remisondaz/Desktop/MATLAB/Varying_spring_constant/Raw_files';



nRuns = 3;
nIntegrins = 100;
k_a = [0.0001 0.001 0.01 0.1 1];
peak = [1 2];
nligands = 400;

dt = 1e-3;
time_max= 29;

  for m = 1:length(nligands) 
        for n = 1:length(peak)
            for k = 1:length(k_a)
                for r = 1:nRuns

                SaveName = ['SIMULATION-001__','Ka_',SimFormat(k_a(k)),'__Peak_',sprintf('%02d',peak(n)), '__nLigands_',sprintf('%04d',nligands ), '__run_', sprintf('%02d',r), '.mat'];
                load(fullfile(RawSaveDirectory, SaveName));
                
               
                Membrane_position = SimData.MembranePosition; 
                
                Membrane_spread = mean(Membrane_position(27001:29001))/std(Membrane_position(27001:29001));% - mean(Membrane_position(15001:19001))/std(Membrane_position(15001:19001));

                STD_position = Membrane_spread/sqrt(30001);

                 Spreadresults.(['k_a_', SimFormat(k_a(k)), '_peak_', sprintf('%01d',peak(n)), '_nligands_', sprintf('%03d', nligands(m)), 'run_', sprintf('%02d',r)]) = struct(...
        'Membrane_spread', Membrane_spread, ...
        'STD', STD_position);

                end
            end
        end
  end
  

    save(fullfile(FileSaveDirectory, 'Membranespreadresults.mat'), 'Spreadresults');