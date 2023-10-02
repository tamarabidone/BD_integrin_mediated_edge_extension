clear all
close all
clc


FileSaveDirectory = '/Users/remisondaz/Desktop/MATLAB/Varying_spring_constant/Actin_retrograde_flow';

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
                

                force = NaN(30000,100);
                retrograde_flow = NaN(30000,100);
                
                for i = 1:size(SimData.AdhesionData,1)
                 for time=1:size(SimData.AdhesionData,3)

                  if SimData.AdhesionData(i,3,time) == 1

                  
                  pos = find(SimData.Data{time, 1}.FilamentName == SimData.AdhesionData(i,5,time));

                      if SimData.Data{time,1}.YSpeed(pos,1) < -1

                          if SimData.MembranePosition(time,1) - abs(SimData.Data{time,1}.XYPosition(pos,2)) > -10
    
                              force(time, i) = SimData.AdhesionData(i,4,time);     
            
                              retrograde_flow(time, i) = SimData.Data{time,1}.YSpeed(pos,1);
                    
                          end
                      end
                  end
                 end
                

                 Flowresults.(['k_a_', SimFormat(k_a(k)), '_peak_', sprintf('%01d',peak(n)), '_nligands_', sprintf('%03d', nligands(m)), 'run_', sprintf('%02d',r)]) = struct(...
        'Adhesion_force', force, ...
        'Retrograde_flow', retrograde_flow);

                end
                end
            end
        end
  end
  

    save(fullfile(FileSaveDirectory, 'ActinFlowresults.mat'), 'Flowresults');