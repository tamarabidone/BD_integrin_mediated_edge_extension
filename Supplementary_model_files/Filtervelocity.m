%% Signal processing

function f = filter(n)

   % Set up low pass filter parameters
                        fc = 10;    % cuttoff freqneucy (Hz). Cuttoff Period = 1/fc
                        fs = 1/dt;  % sampling freqneucy (Hz)
                        [b,a] = butter(1,fc/(fs/2));  % Use a 1st order Butterworth filter (it's a pretty standard filter for smoothing)  
                    % Filter velocity data
                        Vfilt1 = filtfilt(b,a,n); % Use low pass filter with high frequency cuttoff to tame the high peaks
                        IMF = emd(Vfilt1);
                        nIMFs = 2; % nummber of IMF's to remove
                        Vfilt = Vfilt1 - sum(IMF(:,1:nIMFs),2); %sum( IMF(:,2:end)' ); % Remove first 

end
