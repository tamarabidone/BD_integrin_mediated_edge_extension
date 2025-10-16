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
        branching_region = NaN(9000, 40);

            for r = 1:nRuns
                            SaveName = ['SIMULATION-001__','Ka_',SimFormat(k),'__Peak_',sprintf('%02d',p), '__nLigands_',sprintf('%04d',nligands), '_run_', sprintf('%02d',r), '.mat'];
                            load(fullfile(RawSaveDirectory, SaveName));
                             for t = 1:9000
                                branching_region(t,r) = length(find(SimData.Data{t,1}.DistanceMembrane(:) < 15))/(length(SimData.Data{t,1}.Parent(:)));
                             end
            end
                    SaveName = (['Reg_val.','k_a_', SimFormat(k), '_peak_', sprintf('%01d',p), '_nligands_', sprintf('%03d', nligands),'.mat']);
                    Directory = '____' %Insert directory to save files here;
                    FullFilePath = fullfile(Directory, SaveName);
                    save(FullFilePath, 'branching_region');
            end
                           

  %%
% This script generates violin plots using the Violinplot package.
% The package is available online at: https://github.com/bastibe/Violinplot-Matlab
% Make sure that Violinplot.m (and related files) are either in the same folder 
% as this script or added to your MATLAB path.
data1 = load("Reg_val.k_a_000d0001_peak_1_nligands_400.mat");
data2 = load("Reg_val.k_a_000d0001_peak_2_nligands_400.mat");
data3 = load("Reg_val.k_a_000d0010_peak_1_nligands_400.mat");
data4 = load("Reg_val.k_a_000d0100_peak_1_nligands_400.mat");

time = 7000:9000;
nem_order_1 = mean(data1.branching_region(time,1:40), 1, 'omitnan')';
nem_order_2 = mean(data2.branching_region(time,1:40), 1, 'omitnan')';
nem_order_3 = mean(data3.branching_region(time,1:40), 1, 'omitnan')';
nem_order_4 = mean(data4.branching_region(time,1:40), 1, 'omitnan')';
 
%Find indices of the 25 highest points in group 1
[~, idx_low25] = maxk(nem_order_1, 25);

% Select those indices in ALL groups
nem_order_1 = nem_order_1(idx_low25);
nem_order_2 = nem_order_2(idx_low25);
nem_order_3 = nem_order_3(idx_low25);
nem_order_4 = nem_order_4(idx_low25);

data = {nem_order_1, nem_order_2, nem_order_3, nem_order_4};

colors = {[0.50 0.50 0.50] [0.90 0.40 0.60] [0.25 0.70 0.70] [0.40 0.20 0.60]};

figure; 
hold on;
set(gcf,'Color','w')
 
bh = violinplot(data);
 
for f = 1:length(bh)
     bh(f).ViolinColor = colors(f);
     bh(f).ShowData = true;
     bh(f).ShowMedian = true;
     bh(f).MedianColor = [1 0.5 0];
     bh(f).ShowBox = true;
     bh(f).ShowMean = false;
     bh(f).ShowNotches = false;
     bh(f).EdgeColor = [0 0 0];
     bh(f).ShowWhiskers = true;
     %bh(i).QuartileStyle = 'shadow';
     %bh(i).MedianMarkerSize = 48;
     %bh(i).ViolinAlpha = 1;
 end
 set(gca, 'box', 'off', 'fontsize', 25, 'LineWidth', 3)
 set(gca, 'TickLabelInterpreter', 'latex','FontName', 'Times', 'FontSize', 25);
 xticks([1 2 3 4])
 xticklabels({'0.4', '0.4 + Mn^{2+}', '6','60'})
 ylabel('F-actin ratio in d_{branch} (s^{-1}')
 xlabel('Substrate Rigidity (kPa)')
 box off
 hold off
