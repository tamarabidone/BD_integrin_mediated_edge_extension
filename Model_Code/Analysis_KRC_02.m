SaveDirectory = '/Users/keithcarney/Dropbox/ENM/LP/outputs_1212';    
    
   % Setup all combinations of parameters to be varied  -------------------------------------------------------------------
    values = [];
    nRuns = 10; % Total number of runs for each condition
    k_a = [10, 1, 0.1, 0.01, 0.001, 1E-04]; % adhesion spring constant
    peak = [1,2];   % WT or Mn
    
    %Velocity_Raw = cell(nRuns,length(k_a),length(peak));
    Velocity = NaN(nRuns,length(k_a),length(peak));
    Width    = NaN(nRuns,length(k_a),length(peak));
    Periods  = NaN(nRuns,length(k_a),length(peak));
    Prominence = NaN(nRuns,length(k_a),length(peak));
    
    nTotal = numel(Velocity);
    index  = 0;
    
    load( fullfile(SaveDirectory,'VelocityData.mat') ) % Velocities are presaved in a single cell array to speed up processing
    
    for m = 1:length(k_a) % Create combinations of all conditions
        for n = 1:length(peak)
            for r = 1:nRuns
                    index = index + 1;
                    disp([num2str(index),' of ',num2str(nTotal)])
                    % SaveName = ['SIMULATION-001__','Ka_',SimFormat(k_a(m)),'__Peak_',sprintf('%02d',peak(n)) '__run_', sprintf('%02d',r), '.mat'];
                    % load( fullfile(SaveDirectory,SaveName) )
                    % V = SimData.MemVelocity;
                    V = Velocity_Raw{r,m,n};
                    
                   %Velocity_Raw{r,m,n} = V;
                   %  dt = SimData.ModelParameters.TimeStep;
                   %  Time = SimData.TimeVector';
                    
                    % Crop first 10 seconds from Time and Velocity
                    V = V( 10/dt:end );
                    Time1 = Time( 10/dt:end );
                    
                    
                    
                    % Set up low pass filter parameters
                        fc = 10;    % cuttoff freqneucy (Hz). Cuttoff Period = 1/fc
                        fs = 1/dt;  % sampling freqneucy (Hz)
                        [b,a] = butter(1,fc/(fs/2));  % Use a 3rd order Butterworth filter (it's a pretty standard filter for smoothing)  
                    % Filter velocity data
                        Vfilt1 = filtfilt(b,a,V); % Use low pass filter with high frequency cuttoff to tame the high peaks
                        IMF = emd(Vfilt1);
                        nIMFs = 2; % nummber of IMF's to remove
                        Vfilt = Vfilt1 - sum(IMF(:,1:nIMFs),2); %sum( IMF(:,2:end)' ); % Remove first 
                    % Measure peaks
                        [pks,locs,w,p] = findpeaks( Vfilt );
%                         prc  = prctile(p,25); % grab upper percentage of values
%                         idx  = find( p >= prc );
%                         pks  = pks(idx);
%                         locs = locs(idx);
%                         w    = w(idx);
%                         p    = p(idx);
                    
                    % Calculate periods between peaks
                        Velocity(r,m,n)   = median( Vfilt1 );
                        %figure(1);  histogram(Vfilt1,[0:1:60]); 
                        Periods(r,m,n)    = median( diff(Time1(locs)) );
                        Width(r,m,n)      = median( w );
                        Prominence(r,m,n) = median( p );

            end
        end
    end
    
%     figure(2); clf
%     set(gcf,'Color',[1,1,1])
%     T = tiledlayout(2,4);
%     
% % VELOCITY PLOTS
%     nexttile(T,1)
%         B1 = boxplot(Velocity(:,:,1),k_a,'Notch','on','Symbol','+r'); % ,'whisker',0
%         set(B1(:),'LineWidth',2)
%         ylim([  floor(min(Velocity(:))),ceil(max(Velocity(:))) ])
%         %xlabel('Integrin spring constant (pN/nm)')
%         ylabel('Median velocity (nm/s)')
%         title('WildType','FontWeight','normal')
%         xtickangle(30)
%         set(gca,'FontSize',25,'Linewidth',2)
%     nexttile(T,5)
%         B2 = boxplot(Velocity(:,:,2),k_a,'Notch','on','Symbol','+r');
%         set(B2(:),'LineWidth',2)
%         ylim([  floor(min(Velocity(:))),ceil(max(Velocity(:))) ])
%         xlabel('Integrin spring constant (pN/nm)')
%         ylabel('Median velocity (nm/s)')
%         title('Manganese','FontWeight','normal')
%         xtickangle(30)
%         set(gca,'FontSize',25,'Linewidth',2)
%         
% % PERIOD plots
%    nexttile(T,2)
%         B1 = boxplot(Periods(:,:,1),k_a,'Notch','on','Symbol','+r'); % ,'whisker',0
%         set(B1(:),'LineWidth',2)
%         ylim([  floor(min(Periods(:))),ceil(max(Periods(:))) ])
%         ylim([0.2,0.6])
% %         xlabel('Integrin spring constant (pN/nm)')
%         ylabel('Median period (s)')
%         title('WildType','FontWeight','normal')
%         xtickangle(30)
%         set(gca,'FontSize',25,'Linewidth',2)
%     nexttile(T,6)
%         B2 = boxplot(Periods(:,:,2),k_a,'Notch','on','Symbol','+r');
%         set(B2(:),'LineWidth',2)
%         ylim([  floor(min(Periods(:))),ceil(max(Periods(:))) ])
%         ylim([0.2,0.6])
%         xlabel('Integrin spring constant (pN/nm)')
%         ylabel('Median period (s)')
%         title('Manganese','FontWeight','normal')
%         xtickangle(30)
%         set(gca,'FontSize',25,'Linewidth',2)                
%                 
% % WIDTH PLOTS
%     nexttile(T,3)
%         B1 = boxplot(dt*Width(:,:,1),k_a,'Notch','on','Symbol','+r'); % ,'whisker',0
%         set(B1(:),'LineWidth',2)
%         ylim([  floor(min(Width(:))),ceil(max(Width(:))) ])
%         ylim([0,0.3])
% %         xlabel('Integrin spring constant (pN/nm)')
%         ylabel('Median peak width (s)')
%         title('WildType','FontWeight','normal')
%         xtickangle(30)
%         set(gca,'FontSize',25,'Linewidth',2)
%     nexttile(T,7)
%         B2 = boxplot(dt*Width(:,:,2),k_a,'Notch','on','Symbol','+r');
%         set(B2(:),'LineWidth',2)
%         ylim([  floor(min(Width(:))),ceil(max(Width(:))) ])
%         ylim([0,0.3])
%         xlabel('Integrin spring constant (pN/nm)')
%         ylabel('Median peak width (s)')
%         title('Manganese','FontWeight','normal')
%         xtickangle(30)
%         set(gca,'FontSize',25,'Linewidth',2)     
%         
%         
% % PROMINENCE PLOTS
%     nexttile(T,4)
%         B1 = boxplot(Prominence(:,:,1),k_a,'Notch','on','Symbol','+r'); % ,'whisker',0
%         set(B1(:),'LineWidth',2)
%         ylim([  floor(min(Prominence(:))),ceil(max(Prominence(:))) ])
% %         xlabel('Integrin spring constant (pN/nm)')
%         ylabel('Median peak prominence (s)')
%         title('WildType','FontWeight','normal')
%         xtickangle(30)
%         set(gca,'FontSize',25,'Linewidth',2)
%     nexttile(T,8)
%         B2 = boxplot(Prominence(:,:,2),k_a,'Notch','on','Symbol','+r');
%         set(B2(:),'LineWidth',2)
%         ylim([  floor(min(Prominence(:))),ceil(max(Prominence(:))) ])
%         xlabel('Integrin spring constant (pN/nm)')
%         ylabel('Median peak prominence (s)')
%         title('Manganese','FontWeight','normal')
%         xtickangle(30)
%         set(gca,'FontSize',25,'Linewidth',2)   
%         
%         
% Sample Velocity plot        
    figure(3); clf
        set(gcf,'Color',[1,1,1])
        plot(Time1,V,'-b',Time1,Vfilt,'-r',Time1,Vfilt1,'-k')
        xlim([20,40])
        ylim([-25,100])
        xlabel('Time (s)')
        ylabel('Velocity (nm/s)')
        set(gca,'FontSize',25,'Linewidth',2)   
        
        
        
 
        
% Overlaying PLOTS with error bars        
    figure(4); clf
    set(gcf,'Color',[1,1,1])
    T = tiledlayout(1,4);
    LW = 4;
    CS = 10;
    Falpha = 0.15;
    TitleText = 'Mean +/- SEM';
    
    AFontSize = 25;
    XFontSize = 30;
    YFontSize = 32;
    LFontSize = 26;
    TFontSize = 16;
    
% VELOCITY PLOTS
    nexttile(T,1)
        Er1 =  std(Velocity(:,:,1))/sqrt(10); % sort(Velocity(:,:,1));
        M1 = mean(Velocity(:,:,1));
        E1 = plot(k_a,M1,'s-','LineWidth',LW,'Color',[0,0,0.7]);
        hold on
        F1 = fill([k_a,fliplr(k_a)],[M1+Er1, fliplr(M1-Er1)],E1.Color,'FaceAlpha',Falpha,'EdgeColor','none');
        %E1 = errorbar(k_a,M1,Er1,'s-','LineWidth',LW,'CapSize',CS);
       
        Er2 = std(Velocity(:,:,2))/sqrt(10); % sort(Velocity(:,:,2)); 
        M2 = mean(Velocity(:,:,2));
        E2 = plot(k_a,M2,'s-','LineWidth',LW,'Color',[0.5,0,0]);
        F2 = fill([k_a,fliplr(k_a)],[M2+Er2, fliplr(M2-Er2)],E2.Color,'FaceAlpha',Falpha,'EdgeColor','none');
        %E2 = errorbar(k_a, M2, Er2,'s-','LineWidth',LW,'CapSize',CS);
        hold off
        
        set(gca,'LineWidth',LW,'XScale','log','XMinorTick','off','FontSize',AFontSize)
        xlim([2E-5,50])
        ylim([10,15]) 
       
            
        xlabel('\it{k} (pN/nm)','FontWeight','bold','FontSize',XFontSize)
        ylabel('Velocity (nm/s)','FontWeight','bold','FontSize',YFontSize)
        xticks(fliplr(k_a))
        xtickangle(30)
        title(gca,TitleText,'FontSize',TFontSize)
        axis square
        legend([E1,E2],'Wild type','Mn^{2+}','Location','northwest','FontSize',LFontSize)
        legend boxoff
        
% PERIOD plots
   nexttile(T,2)
        Er1 = std(Periods(:,:,1))/sqrt(10); % sort(Velocity(:,:,1));
        M1 = mean(Periods(:,:,1));
        %E1 = errorbar(k_a,M1,Er1,'s-','LineWidth',LW,'CapSize',CS);
        E1 = plot(k_a,M1,'s-','LineWidth',LW,'Color',[0,0,0.7]); hold on
        F1 = fill([k_a,fliplr(k_a)],[M1+Er1, fliplr(M1-Er1)],E1.Color,'FaceAlpha',Falpha,'EdgeColor','none');
        
        Er2 = std(Periods(:,:,2))/sqrt(10); % sort(Velocity(:,:,2)); 
        M2 = mean(Periods(:,:,2));
        E2 = plot(k_a,M2,'s-','LineWidth',LW,'Color',[0.5,0,0]);
        F2 = fill([k_a,fliplr(k_a)],[M2+Er2, fliplr(M2-Er2)],E2.Color,'FaceAlpha',Falpha,'EdgeColor','none');
        %E2 = errorbar(k_a, M2, Er2,'s-','LineWidth',LW,'CapSize',CS);
        hold off
        
        set(gca,'LineWidth',LW,'XScale','log','XMinorTick','off','FontSize',AFontSize)
        
        xlim([2E-5,50])
        if nIMFs == 1
            ylim([0.25,0.45]) 
            set(gca,'YTickLabel',compose('%0.2f',yticks'))
        elseif nIMFs == 2
            ylim([0.7,1.2]) 
            set(gca,'YTickLabel',compose('%0.1f',yticks'))
        end
        
        xlabel('\it{k} (pN/nm)','FontWeight','bold','FontSize',XFontSize)
        ylabel('Period (s)','FontWeight','bold','FontSize',YFontSize)
        xticks(fliplr(k_a))
        xtickangle(30)
        title(gca,TitleText,'FontSize',TFontSize)
        axis square
        legend([E1,E2],'Wild type','Mn^{2+}','Location','northwest','FontSize',LFontSize)
        legend boxoff
   
      
                
% WIDTH PLOTS
    nexttile(T,3)
        Er1 = std(Width(:,:,1))/sqrt(10); % sort(Velocity(:,:,1));
        M1 = mean(Width(:,:,1));
        %E1 = errorbar(k_a,M1,Er1,'s-','LineWidth',LW,'CapSize',CS);
        E1 = plot(k_a,M1,'s-','LineWidth',LW,'Color',[0,0,0.7]); hold on
        F1 = fill([k_a,fliplr(k_a)],[M1+Er1, fliplr(M1-Er1)],E1.Color,'FaceAlpha',Falpha,'EdgeColor','none');
        
        Er2 = std(Width(:,:,2))/sqrt(10); % sort(Velocity(:,:,2)); 
        M2 = mean(Width(:,:,2));
        %E2 = errorbar(k_a, M2, Er2,'s-','LineWidth',LW,'CapSize',CS);
        E2 = plot(k_a,M2,'s-','LineWidth',LW,'Color',[0.5,0,0]);
        F2 = fill([k_a,fliplr(k_a)],[M2+Er2, fliplr(M2-Er2)],E2.Color,'FaceAlpha',Falpha,'EdgeColor','none');

        set(gca,'LineWidth',LW,'XScale','log','XMinorTick','off','FontSize',AFontSize)
        
        xlim([2E-5,50])
        if nIMFs == 1
            ylim([100,180]) 
        elseif nIMFs == 2
            ylim([300,500])
        end
        
        xlabel('\it{k} (pN/nm)','FontWeight','bold','FontSize',XFontSize)
        ylabel('Peak width (s)','FontWeight','bold','FontSize',YFontSize)
        xticks(fliplr(k_a))
        xtickangle(30)
        title(gca,TitleText,'FontSize',TFontSize)
        axis square
        legend([E1,E2],'Wild type','Mn^{2+}','Location','northwest','FontSize',LFontSize)
        legend boxoff
        
% PROMINENCE PLOTS
    nexttile(T,4)
        Er1 = std(Prominence(:,:,1))/sqrt(10); % sort(Velocity(:,:,1));
        M1 = mean(Prominence(:,:,1));
        %E1 = errorbar(k_a,M1,Er1,'s-','LineWidth',LW,'CapSize',CS);
        E1 = plot(k_a,M1,'s-','LineWidth',LW,'Color',[0,0,0.7]); hold on
        F1 = fill([k_a,fliplr(k_a)],[M1+Er1, fliplr(M1-Er1)],E1.Color,'FaceAlpha',Falpha,'EdgeColor','none');
        
        Er2 = std(Prominence(:,:,2))/sqrt(10); % sort(Velocity(:,:,2)); 
        M2 = mean(Prominence(:,:,2));
        %E2 = errorbar(k_a, M2, Er2,'s-','LineWidth',LW,'CapSize',CS);
        E2 = plot(k_a,M2,'s-','LineWidth',LW,'Color',[0.5,0,0]);
        F2 = fill([k_a,fliplr(k_a)],[M2+Er2, fliplr(M2-Er2)],E2.Color,'FaceAlpha',Falpha,'EdgeColor','none');
        
        set(gca,'LineWidth',LW,'XScale','log','XMinorTick','off','FontSize',AFontSize)
        
        xlim([2E-5,50])
        if nIMFs == 1
            ylim([2,10]) 
        elseif nIMFs == 2
            ylim([0,20])
        end
        xlabel('\it{k} (pN/nm)','FontWeight','bold','FontSize',XFontSize)
        ylabel('Peak prominence (nm/s)','FontWeight','bold','FontSize',YFontSize)
        xticks(fliplr(k_a))
        xtickangle(30)
        title(gca,TitleText,'FontSize',TFontSize)
        axis square
        legend([E1,E2],'Wild type','Mn^{2+}','Location','northwest','FontSize',LFontSize)
        legend boxoff
        