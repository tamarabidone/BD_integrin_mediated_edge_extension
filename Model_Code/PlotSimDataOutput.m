Directory = '/Users/keithcarney/Dropbox/ENM/LP/outputs';
SimData_1d00  = load( fullfile( Directory, 'SimData_1pN_wt_1.mat') );
SimData_0d10  = load( fullfile( Directory, 'SimData_01pN_wt_1.mat') );                        
SimData_0d01  = load( fullfile( Directory, 'SimData_001pN_wt_1.mat') );

TimeVec = SimData_1d00.TimeVector;

figure(1); clf
T = tiledlayout(2,1);
set(gcf,'Color',[1,1,1])

nexttile
    plot(TimeVec,SimData_1d00.MemVelocity,'b-',...
         TimeVec,SimData_0d10.MemVelocity,'r-',...
         TimeVec,SimData_0d01.MemVelocity,'-k','LineWidth',2);  
    %xticks(0:5)
    ylim([0,50])
    xlim([10,60])
    set(gca,'FontSize',28,'LineWidth',2)
    xlabel('Time (s)')
    ylabel('Velocity (nm/s)')
    legend({'k_{A} = 0.1 pN/nm','k_{A} = 10 pN/nm','k_{A} = 1000 pN/nm'})
    grid on
    
nexttile
    plot(TimeVec,SimData_1d00.nMonomers,'b-',...
         TimeVec,SimData_0d10.nMonomers,'r-',...
         TimeVec,SimData_0d01.nMonomers,'-k','LineWidth',2);  
    %xticks(0:5)
    %ylim([1e4,1.8e4])
    xlim([10,60])
    set(gca,'FontSize',28,'LineWidth',2)
    xlabel('Time (s)')
    ylabel('Filament mass (monomers)')
    legend({'k_{A} = 0.1 pN/nm','k_{A} = 10 pN/nm','k_{A} = 1000 pN/nm'},'Location','southeast')
    grid on