clc 
clear
close all

RawSaveDirectory = ' _____ '; %Insert directory to folder containing raw files here

nRuns = 40;
k_a = [0.0001 0.0001 0.001 0.001 0.01 0.01];
peak = [1 2 1 2 1 2];
nligands = 400;


    for m = 1:length(k_a) % Create combinations of all conditions
        p = peak(m);
        k = k_a(m);
        MemPos = NaN(40,1);

            for r = 1:nRuns
                            SaveName = ['SIMULATION-001__','Ka_',SimFormat(k),'__Peak_',sprintf('%02d',p), '__nLigands_',sprintf('%04d',nligands), '_run_', sprintf('%02d',r), '.mat'];
                            load(fullfile(RawSaveDirectory, SaveName));                         
                             MemPos(r) = mean(SimData.MembranePosition(8001:9001));
            end
                    SaveName = (['Mem_val.','k_a_', SimFormat(k), '_peak_', sprintf('%01d',p), '_nligands_', sprintf('%03d', nligands),'.mat']);
                    Directory = ' _____ '; %Insert directory to folder containing data files here
                    FullFilePath = fullfile(Directory, SaveName);
                    save(FullFilePath,  'MemPos');
            end
    %     end
    % end
  
%%

data1 = load("Mem_val.k_a_000d0001_peak_1_nligands_400.mat");
data2 = load("Mem_val.k_a_000d0001_peak_2_nligands_400.mat");
data3 = load("Mem_val.k_a_000d0010_peak_1_nligands_400.mat");
data4 = load("Mem_val.k_a_000d0100_peak_1_nligands_400.mat");

nem_order_1 = data1.MemPos;
nem_order_2 = data2.MemPos;
nem_order_3 = data3.MemPos;
nem_order_4 = data4.MemPos;
 
%Find indices of the 25 lowest points in group 1
[~, idx_low25] = maxk(nem_order_1, 25);

% Select those indices in ALL groups
nem_order_1 = nem_order_1(idx_low25);
nem_order_2 = nem_order_2(idx_low25);
nem_order_3 = nem_order_3(idx_low25);
nem_order_4 = nem_order_4(idx_low25);

data = {nem_order_1, nem_order_2, nem_order_3, nem_order_4};

figure;  
hold on; 
set(gcf,'Color','w');

wt_color  = [0.4 0.7 1.0];  
mn_color  = [1.0 0.2 0.2];
colors = {wt_color, mn_color, wt_color, wt_color};

for i = 1:numel(data)
    y = data{i};
    y = y(:);
    y = y(~isnan(y));

    groupIndex = ceil(i/2);                 
    condIndex  = mod(i-1,2)+1; 
    
    scatter(xVals, y, 60, 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerFaceColor', colors{condIndex})

    med = mean(y);
    plot([xThis-0.1 xThis+0.1], [med med], 'k-', 'LineWidth', 2)
end

xticklabels({'0.4', '0.4 + Mn^{2+}', '6', '60'})  
ylabel('Y-position (nm)')
xlabel('Substrate Rigidity (kPa)')
set(gca, 'FontSize', 20, 'LineWidth', 2, 'box', 'off')

h1 = scatter(nan, nan, 60, 'filled', 'MarkerFaceColor', wt_color, 'MarkerFaceAlpha',0.6);
h2 = scatter(nan, nan, 60, 'filled', 'MarkerFaceColor', mn_color, 'MarkerFaceAlpha',0.6);
legend([h1 h2], {'WT', 'Mn^{2+}'}, 'Location','best', 'Box','off')
