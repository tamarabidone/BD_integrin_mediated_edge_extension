clc 
clear
close all

RawSaveDirectory = '/Users/remisondaz/Desktop/MATLAB/Extended_Sims';

nRuns = 6;
nIntegrins = 100;
k_a =  [0.0001, 0.0001, 0.001, 0.01];
peak = [1, 2, 1, 1];
nligands = 400;

for n = 1:length(k_a)
             p = peak(n);
             k = k_a(n);
             Branch = [];
             for r = 1:nRuns
                            tic
                            SaveName = ['SIMULATION-001__','Ka_',SimFormat(k),'__Peak_',sprintf('%02d',p), '__nLigands_',sprintf('%04d',nligands ), '_run_', sprintf('%02d',r), '.mat'];
                            load(fullfile(RawSaveDirectory, SaveName));
                            for t = 300:400
                                 % structNames = SimData.Data{t,1}.StructName;
                                 % bins = 0.5:1:length(structNames)+0.5;
                                 % BinCenters = floor(bins(2:end));
                                 % counts = histcounts(structNames,bins)';
                                 % cIdx = find(counts >= 2);
                                 % Lia = ismember(structNames,BinCenters(cIdx));
                                 % Struct = [Struct;length(structNames(Lia))/length(structNames)];
                                 branches = find(SimData.Data{t,1}.Parent ~= 0);
                                 Branch = [Branch;length(branches)/length(SimData.Data{t,1}.Parent)];
                           end       
                     end   

                    SaveName = (['Branch_vals.','k_a_', SimFormat(k), '_peak_', num2str(p), '_nligands_', sprintf('%03d', nligands),'.mat']);
                    Directory = '/Users/remisondaz/Desktop/MATLAB/Histograms';
                    FullFilePath = fullfile(Directory, SaveName);
                    save(FullFilePath, 'Branch');
end

  %%
  data1 = load("Branch_vals.k_a_000d0000_peak_1_nligands_400.mat");
  data2 = load("Branch_vals.k_a_000d0000_peak_2_nligands_400.mat");
  data3 = load("Branch_vals.k_a_000d0001_peak_1_nligands_400.mat");
  data4 = load("Branch_vals.k_a_000d0001_peak_2_nligands_400.mat");
  data5 = load("Branch_vals.k_a_000d0010_peak_1_nligands_400.mat");
  data6 = load("Branch_vals.k_a_000d0010_peak_2_nligands_400.mat");
  data7 = load("Branch_vals.k_a_000d0100_peak_1_nligands_400.mat");
  data8 = load("Branch_vals.k_a_000d0100_peak_2_nligands_400.mat");
  data9 = load("Branch_vals.k_a_000d1000_peak_1_nligands_400.mat");
  data10 = load("Branch_vals.k_a_000d1000_peak_2_nligands_400.mat");
  data11 = load("Branch_vals.k_a_001d0000_peak_1_nligands_400.mat");
  data12 = load("Branch_vals.k_a_001d0000_peak_2_nligands_400.mat");

  nem = {data3.Branch,data4.Branch, data5.Branch, data7.Branch};
% average_over_seconds = @(x) arrayfun(@(i) mean(x(i:min(i+999, numel(x)))), 1:1000:numel(x));
%   nem = cellfun(average_over_seconds, nem, 'UniformOutput', false);
  [Fout,nCells,Ns] = CellArray2PaddedNanArray(nem');
FH = figure(1); clf
FH.Color = 'w';
MarkerSize = 20;
 bh = violinplot(Fout(:, [1,2,3,4]));
colors = {[0.1, 0.1, 0.1], [1, 0.75, 0.796], [0.678, 0.847, 0.902], [0.502, 0, 0.502]};  % Example colors
% Loop through each violin and set the color
for i = 1:length(bh)
    bh(i).ViolinColor = colors(i);
    bh(i).ShowData = false;
    bh(i).ShowMedian = true;
    bh(i).MedianColor = [1 0.5 0];
    bh(i).ShowBox = true;
    bh(i).ShowMean = false;
    bh(i).ShowNotches = false;
    bh(i).EdgeColor = [0 0 0];
    bh(i).ShowWhiskers = true;
    %bh(i).QuartileStyle = 'shadow';
    %bh(i).MedianMarkerSize = 48;
    %bh(i).ViolinAlpha = 1;
end
 %bh.LineWidth = 3;
  set(gca, 'box', 'off', 'fontsize', 42, 'LineWidth', 6)
  xticks([1 2 3 4])
  %ylim([0 1])
  set(gca, 'TickLabelInterpreter', 'latex','FontName', 'Arial','FontSize', 50);
xticklabels({'0.4', '$0.4 + Mn^{2+}$', '6', '60'});
xlabel('Substrate rigidity (kPa)', 'FontSize', 55);
  ylabel({'F-actin Branching', 'Ratio'}, 'FontSize', 55);
  set(gcf, 'Position',[700 100 700 1000]);
  exportgraphics(gca,'Branchingratio4cases50-60.png','Resolution',350)
