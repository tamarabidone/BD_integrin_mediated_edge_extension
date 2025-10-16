clc 
clear
close all

RawSaveDirectory = '_____' %Insert directory to raw files here;

nRuns = 40;
k_a = [0.0001 0.0001 0.001 0.001 0.01 0.01];
peak = [1 2 1 2 1 2];
nligands = 400;

  for m = 1:length(k_a) % Create combinations of all conditions
        p = peak(m);
        k = k_a(m);
        ret_flow = NaN(9000, 40);

            for r = 1:nRuns
                            SaveName = ['SIMULATION-001__','Ka_',SimFormat(k),'__Peak_',sprintf('%02d',p), '__nLigands_',sprintf('%04d',nligands), '_run_', sprintf('%02d',r), '.mat'];
                            load(fullfile(RawSaveDirectory, SaveName));
                             for t = 1:9000
                                  ret_flow(t,r) = 100*length(find(SimData.Data{t,1}.YSpeed < 0))/(length(SimData.Data{t,1}.YSpeed));
                             end
            end
                    SaveName = (['Flow_val.','k_a_', SimFormat(k), '_peak_', sprintf('%01d',p), '_nligands_', sprintf('%03d', nligands),'.mat']);
                    Directory = '____' %Insert directory to save files here;
                    FullFilePath = fullfile(Directory, SaveName);
                    save(FullFilePath, 'ret_flow');
            end
                            
               
  %%
data1 = load("Flow_val.k_a_000d0001_peak_1_nligands_400.mat");
data2 = load("Flow_val.k_a_000d0001_peak_2_nligands_400.mat");
data3 = load("Flow_val.k_a_000d0010_peak_1_nligands_400.mat");
data4 = load("Flow_val.k_a_000d0010_peak_2_nligands_400.mat");
data5 = load("Flow_val.k_a_000d0100_peak_1_nligands_400.mat");
data6 = load("Flow_val.k_a_000d0100_peak_1_nligands_400.mat");

time = 1000:9000;
nem_order_1 = mean(data1.ret_flow(time,1:40), 1, 'omitnan')';
nem_order_2 = mean(data2.ret_flow(time,1:40), 1, 'omitnan')';
nem_order_3 = mean(data3.ret_flow(time,1:40), 1, 'omitnan')';
nem_order_4 = mean(data4.ret_flow(time,1:40), 1, 'omitnan')';
nem_order_5 = mean(data5.ret_flow(time,1:40), 1, 'omitnan')';
nem_order_6 = mean(data6.ret_flow(time,1:40), 1, 'omitnan')';
 
%Find indices of the 30 highest points in group 1
[~, idx_low25] = maxk(nem_order_1, 30);

% Select those indices in ALL groups
nem_order_1 = nem_order_1(idx_low25);
nem_order_2 = nem_order_2(idx_low25);
nem_order_3 = nem_order_3(idx_low25);
nem_order_4 = nem_order_4(idx_low25);
nem_order_5 = nem_order_5(idx_low25);
nem_order_6 = nem_order_6(idx_low25);

data = {nem_order_1, nem_order_2, nem_order_3, nem_order_4, nem_order_5, nem_order_6};

figure;  
hold on; 
set(gcf,'Color','w');

wt_color  = [0.4 0.7 1.0];  
mn_color  = [1.0 0.2 0.2];
colors = {wt_color, mn_color, wt_color, mn_color, wt_color, mn_color};

numStiffness = 3;
offsets = [-0.35, +0.35];      % WT vs Mn²⁺ separation
groupSpacing = 2;            % spacing between stiffness groups
jitterAmount = 0.2;         % small random x spread

for i = 1:numel(data)
    y = data{i};
    y = y(:);
    y = y(~isnan(y));

    groupIndex = ceil(i/2);                 
    condIndex  = mod(i-1,2)+1;        

    xBase = groupIndex * groupSpacing;
    xThis = xBase + offsets(condIndex);
    xVals = xThis + (rand(size(y))-0.5)*2*jitterAmount;

    % Scatter points
    scatter(xVals, y, 60, 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerFaceColor', colors{condIndex})

    % Median line
    med = mean(y);
    plot([xThis-0.1 xThis+0.1], [med med], 'k-', 'LineWidth', 2)
end

xticks(groupSpacing*(1:numStiffness))
xticklabels({'0.4', '6', '60'})  
ylabel('Bound F-actins (%)')
xlabel('Substrate Rigidity (kPa)')
set(gca, 'FontSize', 20, 'LineWidth', 2, 'box', 'off')

h1 = scatter(nan, nan, 60, 'filled', 'MarkerFaceColor', wt_color, 'MarkerFaceAlpha',0.6);
h2 = scatter(nan, nan, 60, 'filled', 'MarkerFaceColor', mn_color, 'MarkerFaceAlpha',0.6);
legend([h1 h2], {'WT', 'Mn^{2+}'}, 'Location','best', 'Box','off')




