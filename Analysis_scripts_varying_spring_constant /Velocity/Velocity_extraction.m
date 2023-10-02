clc
clear
close all
clc

%This is the folder where we save the outputs from this analysis 
FileSaveDirectory = '/Users/remisondaz/Desktop/MATLAB/Varying_spring_constant/Velocity';

nRuns = 3;
 k_a =  [0.0001 0.001 0.01 0.1 1]; 
peak = [1 2];
ligands = 400;

dt = 1e-3;
time_max= 29;

fc = 10;    % cuttoff freqneucy (Hz). Cuttoff Period = 1/fc
fs = 1/dt;  % sampling freqneucy (Hz)
[b,a] = butter(1,fc/(fs/2));

Time = [0:1e-3:time_max, 1];

Time1 = Time(10/dt:end);


        for n = 1:length(peak)
                for k = 1:length(ligands)
                    for m = 1:length(k_a)

                    SaveName = ([sprintf('%01d',peak(n)), '_', SimFormat(k_a(m)), '_', num2str(ligands(k)),'.mat']);
                    load(fullfile(FileSaveDirectory, SaveName));
                    V = Velocity_Raw;

                    vector_dim = size(V);

                    Velocity = NaN(nRuns,1);
                    Width    = NaN(nRuns,1);
                    Periods  = NaN(nRuns,1);
                    Prominence = NaN(nRuns,1);
                    Frequency = NaN(nRuns, 1);
                   

                  
                    for I = 1:vector_dim(1)
                        for J = 1:vector_dim(2)

                            V = Velocity_Raw{I,J};

                    Vfilt1 = filtfilt(b,a,V(10000:1000*time_max)); % Use low pass filter with high frequency cuttoff to tame the high peaks
                    IMF = emd(Vfilt1);
                    nIMFs = 2; % nummber of IMF's to remove
                    Vfilt = Vfilt1 - sum(IMF(:,1:nIMFs),2); %sum( IMF(:,2:end)' ); % Remove first 
                    % Measure peaks
                    [pks,locs,w,p] = findpeaks( Vfilt );
                    prc  = prctile(p,25); % grab upper percentage of values
                    idx  = find( p >= prc );
                    pks  = pks(idx);
                    locs = locs(idx);
                    w    = (w(idx)/1000);
                    p    = p(idx);

                    Velocity(I, J)   = median( Vfilt );
                    %figure(1);  histogram(Vfilt1,[0:1:60]); 
                    Periods(I, J)    = median( diff(Time1(locs)) );
                    Frequency(I, J) = 1/Periods(I, J);
                    Width(I, J)      = median( w ); % multiply by time step to get seconds
                    Prominence(I, J) = median( p ); %nm/s
                   

                  
                        end
                    end

           % Store results in the struct
            Results.(['k_a_', SimFormat(k_a(m)), '_peak_', num2str(peak(n)), '_nligands_', num2str(ligands(k))]) = struct(...
                'Velocity', mean(Velocity), ...
                'STD_Velocity', std(Vfilt)/10000, ...
                'Periods', mean(Periods), ...
                'STD_Periods', std(Periods), ...
                'Width', mean(Width), ...
                'STD_Width', std(Width), ...
                'Prominence', mean(Prominence), ...
                'STD_Prominence', std(p)/13, ...
                'Frequency', mean(Frequency), ...
                'STD_Frequency', std(Frequency) ...
            );
                    end
                end
        end
       
% Save the struct to a file
save(fullfile(FileSaveDirectory, 'Velocityresults.mat'), 'Results');
