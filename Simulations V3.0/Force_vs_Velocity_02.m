SaveDirectory = '/Users/remisondaz/Desktop/MATLAB/LP code';  

nRuns = 10;
k_a    = [10, 1, 0.1, 0.01, 0.001, 1E-04]; % adhesion spring constant
peak   = [1,2];   % WT or Mn
Velocities   = cell(nRuns,length(k_a),length(peak));
ForceSum     = cell(nRuns,length(k_a),length(peak));
ForceMean    = cell(nRuns,length(k_a),length(peak));
LigatedIntegrins = cell(nRuns,length(k_a),length(peak));

nTotal = numel(Velocities);
index  = 0;

load( fullfile(SaveDirectory,'ForceVelocityIntegrinData.mat') ) % Velocities are presaved in a single cell array to speed up processing
    
% for m = 1:length(k_a) % Create combinations of all conditions
%     for n = 1:length(peak)
%         for r = 1:nRuns
%                 index = index + 1;
%                 disp([num2str(index),' of ',num2str(nTotal)])
%                 SaveName = ['SIMULATION-001__','Ka_',SimFormat(k_a(m)),'__Peak_',sprintf('%02d',peak(n)) '__run_', sprintf('%02d',r), '.mat'];
%                 load( fullfile(SaveDirectory,SaveName) )
% 
%                 dt = SimData.ModelParameters.TimeStep;
%                 
%                 fc = 10;    % cuttoff freqneucy (Hz). Cuttoff Period = 1/fc
%                 fs = 1/dt;  % sampling freqneucy (Hz)
%                 [b,a] = butter(1,fc/(fs/2));  % Use a 3rd order Butterworth filter (it's a pretty standard filter for smoothing)  
%                 % Filter velocity data
%                 V = SimData.MemVelocity;
%                 %V = filtfilt(b,a,V); 
%                 
%                 % Grab all the force values in 10*dt time steps
%                 AdhesionForceSum  = nan( length(V),1 );
%                 AdhesionForceMean = nan( length(V),1 );
%                 VelocityValues = nan( length(V),1 );
%                 PercentActive  = nan( length(V),1 );
%                 index2 = 0;
%                 interval = 1; % 10
%                 for t = 1:interval:length(V) %10
%                     %index2 = index2 + 1;
%                     ActiveIntegrins = squeeze(SimData.AdhesionData(:,3,t)); %(:,3,t:t+9));
%                     ForceValues = squeeze(SimData.AdhesionData(:,4,t)); %(:,4,t:t+9));
%                     % ForceValues = ForceValues( ForceValues > 0 ); % Grab all integrin force values greater than zero
%                     
%                     AdhesionForceSum(t,1)  = sum(ForceValues);
%                     AdhesionForceMean(t,1) = mean(ForceValues);
%                     VelocityValues(t,1) = V(t); %mean(V(t:t+9,1));
%                     PercentActive(t,1)  = 100*sum(ActiveIntegrins(:))/numel(ActiveIntegrins(:));
%                 end
%                 
%                 Velocities{r,m,n} = VelocityValues;
%                 ForceSum{r,m,n} = AdhesionForceSum;
%                 ForceMean{r,m,n} = AdhesionForceMean;
%                 LigatedIntegrins {r,m,n} = PercentActive;
%         end
%     end
% end


figure(1); clf
    set(gcf,'Color',[1,1,1])
    T = tiledlayout(1,2);
    LW = 4;
    CS = 10;
    Falpha = 0.15;
    
    AFontSize = 25;
    XFontSize = 30;
    YFontSize = 32;
    LFontSize = 20;
    TFontSize = 16;
    
    AdhesionIndex = [2,4]; % k_a = 1.00 and 0.01
% LEFT PLOT -------------------------------------------------------
    nexttile(T,1)
    
   % Wild Type
       VelocityValues = Velocities(:,AdhesionIndex(1),1);
       VelocityValues = [VelocityValues{:}];
       VelocityValues = VelocityValues(:);
       VelocityValues(find(isnan(VelocityValues))) = 0;
       %VelocityValues = smoothdata(VelocityValues,'gaussian',10);
       ForceValues = ForceMean(:,AdhesionIndex(1),1);
       ForceValues = [ForceValues{:}];
       ForceValues = ForceValues(:);
       ForceValues(find(isnan(ForceValues))) = 0;
       %ForceValues = smoothdata(ForceValues,'gaussian',10);
       [N1,xedges1,yedges1] = histcounts2(VelocityValues, ForceValues, 0:0.2:100, logspace(0,1.2,300));
        N1 = imboxfilt(N1,5); % Smooth
        N1 = N1./max(N1(:)); % normalize counts
        N1(N1<0.02) = nan;  % remove zeros and small values
       [Xgrid1, Ygrid1] = meshgrid( mean([ xedges1(1,1:end-1); xedges1(1,2:end)]',2),...
                                    mean([ yedges1(1,1:end-1); yedges1(1,2:end)]',2));
       
       P1 = surface(Xgrid1,Ygrid1,N1');
       P1.EdgeColor = 'none';
       colormap(turbo)
       hold on
       %P1 = plot( VelocityValues, ForceValues,'.'); hold on
       
       VelocityValues = Velocities(:,AdhesionIndex(2),1);
       VelocityValues = [VelocityValues{:}];
       VelocityValues = VelocityValues(:);
       VelocityValues(find(isnan(VelocityValues))) = 0;
       %VelocityValues = smoothdata(VelocityValues,'gaussian',10);
       ForceValues = ForceMean(:,AdhesionIndex(2),1);
       ForceValues = [ForceValues{:}];
       ForceValues = ForceValues(:);
       ForceValues(find(isnan(ForceValues))) = 0;
       %ForceValues = smoothdata(ForceValues,'gaussian',10);
       [N1,xedges1,yedges1] = histcounts2(VelocityValues, ForceValues, 0:0.2:100, logspace(-2,0.4,300));
        N1 = imboxfilt(N1,5);
        N1 = N1./max(N1(:));
        N1(N1<0.02) = nan;
       [Xgrid1, Ygrid1] = meshgrid( mean([ xedges1(1,1:end-1); xedges1(1,2:end)]',2),...
                                    mean([ yedges1(1,1:end-1); yedges1(1,2:end)]',2));
       
       P2 = surface(Xgrid1,Ygrid1,N1');
       P2.EdgeColor = 'none';  
       
        set(gca,'LineWidth',LW,'YScale','log','XMinorTick','off','FontSize',20)
        xlim([0,60])
        ylim([1.2E-2,50]) 
       
            
       xlabel('Velocity (nm/s)','FontWeight','bold','FontSize',XFontSize)
       ylabel('Mean integrin tension (pN)','FontWeight','bold','FontSize',YFontSize)
       
       text(10,20,'Wild type, \it{k} = 1.00 pN/nm','FontSize',LFontSize)
       text(10,0.2,'Wild type, \it{k} = 0.01 pN/nm','FontSize',LFontSize)
       
%        
%         legend([P1,P2],'Wild type, \it{k} = 1.00 pN/nm',...
%                        'Wild type, \it{k} = 0.01 pN/nm',...
%                        'Location','southeast','FontSize',LFontSize)
       % legend boxoff

% MG 2+
nexttile(T,2)
       VelocityValues = Velocities(:,AdhesionIndex(1),2);
       VelocityValues = [VelocityValues{:}];
       VelocityValues = VelocityValues(:);
       VelocityValues(find(isnan(VelocityValues))) = 0;
       %VelocityValues = smoothdata(VelocityValues,'gaussian',10);
       ForceValues = ForceMean(:,AdhesionIndex(1),2);
       ForceValues = [ForceValues{:}];
       ForceValues = ForceValues(:);
       ForceValues(find(isnan(ForceValues))) = 0;
       %ForceValues = smoothdata(ForceValues,'gaussian',10);
       [N1,xedges1,yedges1] = histcounts2(VelocityValues, ForceValues, 0:0.2:100, logspace(0,1.2,300));
        N1 = imboxfilt(N1,5);
        N1 = N1./max(N1(:)); % normalize counts
        N1(N1<0.02) = nan;  % remove zeros and small values
       [Xgrid1, Ygrid1] = meshgrid( mean([ xedges1(1,1:end-1); xedges1(1,2:end)]',2),...
                                    mean([ yedges1(1,1:end-1); yedges1(1,2:end)]',2));
       
       P1 = surface(Xgrid1,Ygrid1,N1');
       P1.EdgeColor = 'none';
       colormap(turbo)
       hold on
       
       VelocityValues = Velocities(:,AdhesionIndex(2),2);
       VelocityValues = [VelocityValues{:}];
       VelocityValues = VelocityValues(:);
       VelocityValues(find(isnan(VelocityValues))) = 0;
       %VelocityValues = smoothdata(VelocityValues,'gaussian',10);
       ForceValues = ForceMean(:,AdhesionIndex(2),2);
       ForceValues = [ForceValues{:}];
       ForceValues = ForceValues(:);
       ForceValues(find(isnan(ForceValues))) = 0;
       %ForceValues = smoothdata(ForceValues,'gaussian',10);
       [N1,xedges1,yedges1] = histcounts2(VelocityValues, ForceValues, 0:0.2:100, logspace(-2,0.4,300));
        N1 = imboxfilt(N1,5);
        N1 = N1./max(N1(:));
        N1(N1<0.02) = nan;
       [Xgrid1, Ygrid1] = meshgrid( mean([ xedges1(1,1:end-1); xedges1(1,2:end)]',2),...
                                    mean([ yedges1(1,1:end-1); yedges1(1,2:end)]',2) );
       
       P2 = surface(Xgrid1,Ygrid1,N1');
       P2.EdgeColor = 'none';  
       
        set(gca,'LineWidth',LW,'YScale','log','XMinorTick','off','FontSize',20)
        xlim([0,60])
        ylim([1.2E-2,50]) 
        
        text(10,20, 'Mn^{2+}, \it{k} = 1.00 pN/nm','FontSize',LFontSize)
        text(10,0.2,'Mn^{2+}, \it{k} = 0.01 pN/nm','FontSize',LFontSize)
            
       xlabel('Velocity (nm/s)','FontWeight','bold','FontSize',XFontSize)
       %ylabel('Mean integrin tension (pN)','FontWeight','bold','FontSize',YFontSize)
%        
     


