Directory = '/Users/keithcarney/Dropbox/ENM/LP/outputs';
figure(1); clf
figure(4); clf
figure(3); clf
for n = 1:3
    switch n
        case 1
            FileName = 'SimData_1pN_wt_1.mat';
            SimData  = load( fullfile( Directory, FileName ) );
        case 2
            FileName = 'SimData_01pN_wt_1.mat';
            SimData  = load( fullfile( Directory, FileName ) );   
        case 3
            FileName = 'SimData_001pN_wt_1.mat';
            SimData  = load( fullfile( Directory, FileName ) );
    end
    
    V = SimData.MemVelocity;
    dt = SimData.ModelParameters.TimeStep;
    V = V( 10/dt:end);
    Time = SimData.TimeVector';
    Time = Time( 10/dt:end );


    % Set up low pass filter parameters
        fc = 10;    % cuttoff freqneucy (Hz). Cuttoff Period = 1/fc
        fs = 1/dt; % sampling freqneucy (Hz)
        [b,a] = butter(3,fc/(fs/2));  % Use a 3rd order Butterworth filter (it's a pretty standard filter for smoothing)  
    % Filter velocity data
        Vfilt = filtfilt(b,a,V); 

    figure(1)
    if n == 1
        set(gcf,'Color',[1,1,1])
        %subplot(2,1,1)
        plot(Time,V    ,'-b','LineWidth',1); hold(gca,'on')
        plot(Time,Vfilt,'-r','LineWidth',2); hold(gca,'off')
        set(gca,'FontSize',28,'LineWidth',2)
        xlim([30,40])
        ylim([0,100])
        xlabel('Time (s)')
        ylabel('Velocity (nm/s)')
        title(FileName,'FontSize',18,'Interpreter','none')
    end
    % Measure peaks
        [pks,locs,w,p] = findpeaks(Vfilt);
        % find the 50% most prominent peaks 
        prc = prctile(p,75); 
        idx  = find( p >= prc );
        pks  = pks(idx);
        locs = locs(idx);
        w    = w(idx);
        p    = p(idx);
        % Add points of measured peaks to plot
    if n == 1
        figure(1)
        hold(gca,'on')
        plot(Time(locs),Vfilt(locs),'or','MarkerSize',10,'LineWidth',2,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,0,0])
        set(gca,'FontSize',22,'LineWidth',2)
    end
    % Calculate periods between peaks
        Periods = diff( Time(locs) );
        figure(4)
        set(gcf,'Color',[1,1,1])
            subplot(3,2,(n-1)*2+1)
                histogram(Periods,0:0.1:1.5,'Normalization','count')
                xlabel('Period between peaks (s)')
                ylabel('Counts')
                set(gca,'FontSize',22,'LineWidth',2)
                title(FileName,'FontSize',18,'Interpreter','none')
                ylim([0,20])
            subplot(3,2,(n-1)*2+2)
                histogram(w*dt,0:0.05:0.5,'Normalization','count')
                xlabel('Base width of peaks (s)')
                ylabel('Counts')
                set(gca,'FontSize',22,'LineWidth',2)
                title(FileName,'FontSize',18,'Interpreter','none')
                ylim([0,100])
    % Apply Fourier Analysis
            Vf = V(10/dt:end);
            ReSamp = 50;
            Vf = resample(Vf,1,ReSamp);
            
            figure(3)
            subplot(3,1,n)
            set(gcf,'Color',[1,1,1])
            [pxx,f,pxxc] = pwelch(Vf,[],[],[],(1/dt)/ReSamp,'ConfidenceLevel',0.95);
            plot(f,10*log10(pxx),'LineWidth',2)
            hold on
            plot(f,10*log10(pxxc),'-.','LineWidth',2)
            hold off
            xticks([0:10])
            xlabel('Frequency (Hz)')
            ylabel('PSD (dB/Hz)')
             set(gca,'FontSize',22,'LineWidth',2)
            title([FileName, '    Power Spectrum +/- 95%-Confidence Bounds'],'FontSize',18,'Interpreter','none')
           
end


    