clear all
close all
clc


FileSaveDirectory = '/Users/remisondaz/Desktop/MATLAB/Varying_spring_constant/Adhesion_position';

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
                
                memb_pos=SimData.MembranePosition;
                start_time=19001;
                end_time=29000;

                for time=1:size(SimData.AdhesionData,3)
                   Distance(time,:)= memb_pos(time)-SimData.AdhesionData(:,2, time);
                   activation(time,:)=SimData.AdhesionData(:,3, time);
                end
                
                Distance=Distance(start_time:end_time,:);
                activation=activation(start_time:end_time,:);
                active_integrins_distance=(Distance.*activation);

                 Positionresults.(['k_a_', SimFormat(k_a(k)), '_peak_', sprintf('%01d',peak(n)), '_nligands_', sprintf('%03d', nligands(m)), 'run_', sprintf('%02d',r)]) = struct(...
        'Adhesion_position', active_integrins_distance);

                end
            end
        end
  end

    save(fullfile(FileSaveDirectory, 'IntegrinPositionresults.mat'), 'Positionresults');
  

