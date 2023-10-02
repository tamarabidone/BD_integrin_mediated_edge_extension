clear all
close all
clc


FileSaveDirectory = '/Users/remisondaz/Desktop/MATLAB/Varying_spring_constant/Adhesion_frequency';

RawSaveDirectory = '/Users/remisondaz/Desktop/MATLAB/Varying_spring_constant/Raw_files';



nRuns = 3;
nIntegrins = 100;
k_a = [0.0001 0.001 0.01 0.1 1];
peak = [1 2];
nligands = 400 ;

dt = 1e-3;
time_max= 29;

  for m = 1:length(nligands) % Create combinations of all conditions
        for n = 1:length(peak)
            for k = 1:length(k_a)
                for r = 1:nRuns

                    SaveName = ['SIMULATION-001__','Ka_',SimFormat(k_a(k)),'__Peak_',sprintf('%02d',peak(n)), '__nLigands_',sprintf('%04d',nligands(m) ), '__run_', sprintf('%02d',r), '.mat'];
                load(fullfile(RawSaveDirectory, SaveName));
                
                    Active_integrin_frequency = SimData.nAdhesions;

                     Frequencyresults.(['k_a_', SimFormat(k_a(k)), '_peak_', sprintf('%01d',peak(n)), '_nligands_', sprintf('%03d', nligands(m)), 'run_', sprintf('%02d',r)]) = struct(...
        'Active_integrin_frequency', Active_integrin_frequency);

                end
            end
        end
  end

  save(fullfile(FileSaveDirectory, 'AdhesionFrequencyresults.mat'), 'Frequencyresults');


