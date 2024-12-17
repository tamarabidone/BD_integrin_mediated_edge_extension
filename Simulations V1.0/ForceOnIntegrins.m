SaveDirectory = '/Users/keithcarney/Dropbox/ENM/LP/outputs_1212';  

    values = [];
    nRuns = 10; % Total number of runs for each condition
    k_a = [10, 1, 0.1, 0.01, 0.001, 1E-04]; % adhesion spring constant
    peak = [1,2];   % WT or Mn
    
    IntegrinForceValues = cell(nRuns,length(k_a),length(peak));
    
    nTotal = numel(IntegrinForceValues);
    index  = 0;
    
    load( fullfile(SaveDirectory,'IntegrinForceData.mat') ) % Velocities are presaved in a single cell array to speed up processing
    
%     for m = 1:length(k_a) % Create combinations of all conditions
%         for n = 1:length(peak)
%             for r = 1:nRuns
%                     index = index + 1;
%                     disp([num2str(index),' of ',num2str(nTotal)])
%                     SaveName = ['SIMULATION-001__','Ka_',SimFormat(k_a(m)),'__Peak_',sprintf('%02d',peak(n)) '__run_', sprintf('%02d',r), '.mat'];
%                     load( fullfile(SaveDirectory,SaveName) )
%                    
%                     ForceValues = squeeze(SimData.AdhesionData(:,4,:));
%                     IntegrinForceValues{r,m,n} = ForceValues( ForceValues > 0 ); % Grab all integrin force values greater than zero
% 
%             end
%         end
%     end
    


% Combine data from all 10 runs --------------
ForceValues = cell(2,2);
AdhesionIndex = [2,4]; % k_a = 1.00 and 0.01
for n = 1:2
    for peak = 1:2
        for r = 1:10
            a = AdhesionIndex(n);
            ForceValues{n,peak} = [ForceValues{n,peak}; IntegrinForceValues{r,a,peak}];
        end
    end
end
%---------------------------------------------

figure(1); clf
    set(gcf,'Color',[1,1,1])
    T = tiledlayout(1,1);
    LW = 4;
    CS = 10;
    Falpha = 0.15;
    
    AFontSize = 25;
    XFontSize = 30;
    YFontSize = 32;
    LFontSize = 26;
    TFontSize = 16;
    
% Histograms
    nexttile(T,1)
    
        % Wild Type Forces
        Values = ForceValues{1,1};
        ForceBins1 = 0:1:100;
        H1 = histogram( Values, ForceBins1); hold on
        H1.FaceAlpha = 1;
        
        Values = ForceValues{2,1};
        ForceBins2 = 0:0.01:2;
        H2 = histogram( Values, ForceBins2);
        H2.FaceAlpha = 1;
        
        % Mn2+ Forces
        Values = ForceValues{1,2};
        H3 = histogram( Values, ForceBins1);
        H3.FaceAlpha = 1;
        
        Values = ForceValues{2,2};
        H4 = histogram( Values, ForceBins2);
        H4.FaceAlpha = 1;
        
        % Reverse the plotting order so that the histogram colors correspond to the legend colors
        AH = gca;
        set(gca, 'Children',flipud(AH.Children))
        
        
        set(gca,'LineWidth',LW,'XScale','log','XMinorTick','off','FontSize',20)
        xlim([5E-3,50])
        ylim([0,3E6]) 
       
            
        xlabel('Force on integrins (pN)','FontWeight','bold','FontSize',XFontSize)
        ylabel('Counts','FontWeight','bold','FontSize',YFontSize)
        xticks(fliplr(k_a))
       
        legend([H1,H2,H3,H4],'Wild type, \it{k} = 10^{0} pN/nm',...
                              'Wild type, \it{k} = 10^{-2} pN/nm',...
                              'Mn^{2+}, \it{k} = 10^{0} pN/nm',...
                              'Mn^{2+}, \it{k} = 10^{-2} pN/nm',...
                              'Location','north','FontSize',LFontSize,'NumColumns',2)
        legend boxoff


