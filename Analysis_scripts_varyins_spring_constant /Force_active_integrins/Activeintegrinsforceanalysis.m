clear all
close all
clc


FileSaveDirectory = '/Users/remisondaz/Desktop/MATLAB/Varying_spring_constant/Force_active_integrins';

RawSaveDirectory = '/Users/remisondaz/Desktop/MATLAB/Varying_spring_constant/Raw_files';



nRuns = 3;
nIntegrins = 100;
k_a = [0.0001  0.1];
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
                

               Y_speed = NaN(29900,100);

                for time = 1:29900

                    for i = 1:100

                    if SimData.AdhesionData(i,3,time) == 1

                     pos = find(SimData.Data{time, 1}.FilamentName == SimData.AdhesionData(i,5,time));

                 
                          if SimData.MembranePosition(time, 1) - (SimData.Data{time,1}.XYPosition(pos,2)) <= 15

                              if SimData.Data{time,1}.YSpeed(pos,1) > 0 

                                  if SimData.AdhesionData(i, 4, time) > 0
                                   
                                     force(time, i) = SimData.AdhesionData(i, 4, time);
                                     
                                  end

                              end
                          end
                    end
                    end
                end
                

                
                
                

                Activeintegrinsforce.(['k_a_', SimFormat(k_a(k)), '_peak_', sprintf('%01d',peak(n)), '_nligands_', sprintf('%03d', nligands(m)), 'run_', sprintf('%02d',r)]) = struct(...
        'Active_integrins_force', force);
        
            end
        end
        end
  end
 

    save(fullfile(FileSaveDirectory, 'Activeintegrinsforce.mat'), 'Activeintegrinsforce');