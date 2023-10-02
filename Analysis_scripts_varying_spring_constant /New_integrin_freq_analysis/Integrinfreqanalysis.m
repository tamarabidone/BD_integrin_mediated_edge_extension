clear
close all
clc

%This is the folder where we save the outputs from this analysis 
FileSaveDirectory = '/Users/remisondaz/Desktop/MATLAB/Varying_spring_constant/New_integrin_freq_analysis';

%This folder contains raw model output files after running simulations code (for example: SimulationTaskList_001.m)
RawSaveDirectory = '/Users/remisondaz/Desktop/MATLAB/Varying_spring_constant/Raw_files';

k_a = [0.0001 0.001 0.01 0.1 1];
peak = [1 2];
nligands = 400;
nRuns = 3;

 for m = 1:length(nligands) % Create combinations of all conditions
        for n = 1:length(peak)
            for k = 1:length(k_a)
                for r = 1:nRuns

                    SaveName = ['SIMULATION-001__','Ka_',SimFormat(k_a(k)),'__Peak_',sprintf('%02d',peak(n)), '__nLigands_',sprintf('%04d',nligands(m) ), '__run_', sprintf('%02d',r), '.mat'];
                load(fullfile(RawSaveDirectory, SaveName));


                 for time=1:size(SimData.AdhesionData,3)
                 
                 l = 0;

                 for i = 1:size(SimData.AdhesionData,1)
                 

                  if SimData.AdhesionData(i,4,time) > 0

                  Frequency(time, 1) = l + 1;

                  l = l + 1;

              
                  end
                 end
                 end
                  
                 Frequencyresults.(['k_a_', SimFormat(k_a(k)), '_peak_', sprintf('%01d',peak(n)), '_nligands_', sprintf('%03d', nligands(m)), 'run_', sprintf('%02d',r)]) = struct(...
        'Active_integrin_frequency', Frequency);

                end
            end
        end
 end

        save(fullfile(FileSaveDirectory, 'AdhesionFrequencyresults.mat'), 'Frequencyresults');