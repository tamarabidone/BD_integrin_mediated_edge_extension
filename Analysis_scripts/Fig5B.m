clc 
clear
close all

RawSaveDirectory = ' _____ '; %Insert directory to folder containing raw files here


nRuns = 6;
nIntegrins = 100;
k_a =  [0.0001, 0.0001, 0.001, 0.01];
peak = [1, 2, 1, 1];
nligands = 400;

for n = 1:length(k_a)
             p = peak(n);
             k = k_a(n);
             Distance = [];
             for r = 1:nRuns
                            tic
                            SaveName = ['SIMULATION-001__','Ka_',SimFormat(k),'__Peak_',sprintf('%02d',p), '__nLigands_',sprintf('%04d',nligands ), '_run_', sprintf('%02d',r), '.mat'];
                            load(fullfile(RawSaveDirectory, SaveName));
                            for t = 300:501
                                idx = find(SimData.Data{t,1}.DistanceMembrane(:,:) < 15);
                            Distance = [Distance;SimData.Data{t,1}.DistanceMembrane(idx)];
                            %Distance = [Distance;SimData.Data{t,1}.DistanceMembrane(:,:)];
                            end

                        end

                    SaveName = (['Distance.','k_a_', SimFormat(k), '_peak_', num2str(p), '_nligands_', sprintf('%03d', nligands),'.mat']);
                    Directory = '/Users/remisondaz/Desktop/MATLAB/Histograms';
                    FullFilePath = fullfile(Directory, SaveName);
                    save(FullFilePath,  'Distance');
end
  
%%

  data1 = load("Distance.k_a_000d0000_peak_1_nligands_400.mat");
  data2 = load("Distance.k_a_000d0000_peak_2_nligands_400.mat");
  data3 = load("Distance.k_a_000d0001_peak_1_nligands_400.mat");
  data4 = load("Distance.k_a_000d0001_peak_2_nligands_400.mat");
  data5 = load("Distance.k_a_000d0010_peak_1_nligands_400.mat");
  data6 = load("Distance.k_a_000d0010_peak_2_nligands_400.mat");
  data7 = load("Distance.k_a_000d0100_peak_1_nligands_400.mat");
  data8 = load("Distance.k_a_000d0100_peak_2_nligands_400.mat");
  data9 = load("Distance.k_a_000d1000_peak_1_nligands_400.mat");
  data10 = load("Distance.k_a_000d1000_peak_2_nligands_400.mat");

nem = {data3.Distance,data4.Distance, data5.Distance,data7.Distance};
% average_over_seconds = @(x) arrayfun(@(i) mean(x(i:min(i+9, numel(x)))), 1:10:numel(x));
% nem = cellfun(average_over_seconds, nem, 'UniformOutput', false);
[Fout, nCells, Ns] = CellArray2PaddedNanArray(nem');
Fout = Fout(:, [1,2,3,4]);
FH = figure(1); clf
FH.Color = 'w';
MarkerSize = 20;
 bh = violinplot(Fout);
colors = {[0, 0, 0], [1, 0.75, 0.796], [0.678, 0.847, 0.902], [0.502, 0, 0.502]};  % Example colors
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
%ylim([0 200])
set(gca, 'box', 'off', 'fontsize', 42, 'LineWidth', 6)
xticks([1 2 3 4])
xticklabels({' ', ' ', ' ', ' '});
xlabel('Substrate rigidity (kPa)','FontSize', 55);
ylabel('y_{0} in d_{branch}', 'FontSize',55);
set(gca, 'TickLabelInterpreter', 'latex','FontName', 'Arial', 'FontSize', 50);
set(gcf, 'Position',[700 100 700 1000]);
exportgraphics(gca,'Ratiodbranchextended.png','Resolution',450)
