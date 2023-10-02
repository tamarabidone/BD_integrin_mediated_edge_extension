clear all
close all
clc

%This is the folder where we save the outputs from this analysis 
FileSaveDirectory = '/Users/remisondaz/Desktop/MATLAB/Varying_spring_constant/Fraction_active_filaments';

%This folder contains raw model output files after running simulations code (for example: SimulationTaskList_001.m)
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
                

               Active_fil = NaN(29900,1);

                for time = 1:29900

                     l = 0;

                    for i = 1:100

                    if SimData.AdhesionData(i,3,time) == 1

                     pos = find(SimData.Data{time, 1}.FilamentName == SimData.AdhesionData(i,5,time));

                 
                          if SimData.MembranePosition(time, 1) - (SimData.Data{time,1}.XYPosition(pos,2)) <= 15

                              if SimData.Data{time,1}.YSpeed(pos,1) > 0 
                                   
                                  Active_fil(time, 1) = l + 1;

                                   l = l + 1;
    
                                  
                   

                              end
                          end
                    end
                    end

                    Active_fil(time, 1) = l/size(SimData.Data{time,1}.FilamentName, 1);

                end
                
                
                
                

                 Filamentresults.(['k_a_', SimFormat(k_a(k)), '_peak_', sprintf('%01d',peak(n)), '_nligands_', sprintf('%03d', nligands(m)), 'run_', sprintf('%02d',r)]) = struct(...
        'Fraction_active_filaments', Active_fil);
        
            end
        end
        end
  end

    save(fullfile(FileSaveDirectory, 'Filamentresults.mat'), 'Filamentresults');