SaveDirectory = '/Users/keithcarney/Dropbox/ENM/LP/outputs_1212';  

    values = [];
    nRuns = 10; % Total number of runs for each condition
    k_a = [10, 1, 0.1, 0.01, 0.001, 1E-04]; % adhesion spring constant
    peak = [1,2];   % WT or Mn
    
    PercentLigated = NaN(nRuns,length(k_a),length(peak));
    
    nTotal = numel(PercentLigated);
    index  = 0;
    
    load( fullfile(SaveDirectory,'PercentLigatedData.mat') ) % Velocities are presaved in a single cell array to speed up processing
    
%     for m = 1:length(k_a) % Create combinations of all conditions
%         for n = 1:length(peak)
%             for r = 1:nRuns
%                     index = index + 1;
%                     disp([num2str(index),' of ',num2str(nTotal)])
%                     SaveName = ['SIMULATION-001__','Ka_',SimFormat(k_a(m)),'__Peak_',sprintf('%02d',peak(n)) '__run_', sprintf('%02d',r), '.mat'];
%                     load( fullfile(SaveDirectory,SaveName) )
%                     
%                     
%                     nL = numel(SimData.AdhesionData(:,3,:)); % Total available ligands over all timepoints
%                     LigandActivation = squeeze(SimData.AdhesionData(:,3,:)); 
%                     nLigated = sum(LigandActivation(:));     % Total activated ligands over all timepoints
%                     
%                     PercentLigated(r,m,n) = 100 * nLigated / nL; % Percent activated ligands all timepoints
% 
%             end
%         end
%     end
    
    
    
    figure(1); clf
    set(gcf,'Color',[1,1,1])
    T = tiledlayout(1,1);
    LW = 4;
    CS = 10;
    Falpha = 0.2;
    TitleText = 'Mean +/- SEM';
    
    AFontSize = 25;
    XFontSize = 30;
    YFontSize = 32;
    LFontSize = 26;
    TFontSize = 16;
    
% VELOCITY PLOTS
    nexttile(T,1)
        Er1 =  std(PercentLigated(:,:,1))/sqrt(10); % sort(Velocity(:,:,1));
        M1 = mean(PercentLigated(:,:,1));
        E1 = plot(k_a,M1,'s-','LineWidth',LW,'Color',[0,0,0.7]);
        hold on
        F1 = fill([k_a,fliplr(k_a)],[M1+Er1, fliplr(M1-Er1)],E1.Color,'FaceAlpha',Falpha,'EdgeColor','none');
       
        Er2 = std(PercentLigated(:,:,2))/sqrt(10); % sort(Velocity(:,:,2)); 
        M2 = mean(PercentLigated(:,:,2));
        E2 = plot(k_a,M2,'s-','LineWidth',LW,'Color',[0.5,0,0]);
        F2 = fill([k_a,fliplr(k_a)],[M2+Er2, fliplr(M2-Er2)],E2.Color,'FaceAlpha',Falpha,'EdgeColor','none');
        hold off
        
        set(gca,'LineWidth',LW,'XScale','log','XMinorTick','off','FontSize',AFontSize)
        xlim([2E-5,50])
        ylim([0,100]) 
       
            
        xlabel('\it{k} (pN/nm)','FontWeight','bold','FontSize',XFontSize)
        ylabel('Percent ligated (%)','FontWeight','bold','FontSize',YFontSize)
        xticks(fliplr(k_a))
        xtickangle(30)
        title(gca,TitleText,'FontSize',TFontSize)
        axis square
        legend([E1,E2],'Wild type','Mn^{2+}','Location','northwest','FontSize',LFontSize)
        legend boxoff
        