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
             Force = [];
             for r = 1:nRuns
                            tic
                            SaveName = ['SIMULATION-001__','Ka_',SimFormat(k),'__Peak_',sprintf('%02d',p), '__nLigands_',sprintf('%04d',nligands ), '_run_', sprintf('%02d',r), '.mat'];
                            load(fullfile(RawSaveDirectory, SaveName));
                            for t = 200:351
                            Adhesion_ratio = SimData.Data{t,1}.ForceonMembrane(:);
                            Force = [Force;Adhesion_ratio];
                        end  
                        end

                    SaveName = (['ForceComp.','k_a_', SimFormat(k), '_peak_', num2str(p), '_nligands_', sprintf('%03d', nligands),'.mat']);
                    Directory = '/Users/remisondaz/Desktop/MATLAB/Histograms';
                    FullFilePath = fullfile(Directory, SaveName);
                    save(FullFilePath,  'Force');
end
  
%%

  data1 = load("ForceComp.k_a_000d0000_peak_1_nligands_400.mat");
  data2 = load("ForceComp.k_a_000d0000_peak_2_nligands_400.mat");
  data3 = load("ForceComp.k_a_000d0001_peak_1_nligands_400.mat");
  data4 = load("ForceComp.k_a_000d0001_peak_2_nligands_400.mat");
  data5 = load("ForceComp.k_a_000d0010_peak_1_nligands_400.mat");
  data6 = load("ForceComp.k_a_000d0010_peak_2_nligands_400.mat");
  data7 = load("ForceComp.k_a_000d0100_peak_1_nligands_400.mat");
  data8 = load("ForceComp.k_a_000d0100_peak_2_nligands_400.mat");
  data9 = load("ForceComp.k_a_000d1000_peak_1_nligands_400.mat");
  data10 = load("ForceComp.k_a_000d1000_peak_2_nligands_400.mat");

  idx3 = data3.Force > 0.0257;
  Comp3 = data3.Force(idx3);
  idx4 = data4.Force > 0.0257; 
  Comp4 = data4.Force(idx4);
  idx5 = data5.Force > 0.0257; 
  Comp5 = data5.Force(idx5);
  idx7 = data7.Force > 0.0257; 
  Comp7 = data7.Force(idx7);

nem = {Comp3,Comp4, Comp5,Comp7};
average_over_seconds = @(x) arrayfun(@(i) mean(x(i:min(i+99, numel(x)))), 1:10:numel(x));
nem = cellfun(average_over_seconds, nem, 'UniformOutput', false);
[Fout, nCells, Ns] = CellArray2PaddedNanArray(nem');
Fout = Fout(:, [1 2 3 4]);
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
%ylim([0.025 0.2])
set(gca, 'box', 'off', 'fontsize', 42, 'LineWidth', 6)
xticks([1 2 3 4])
xticklabels({'0.4', '$0.4 + Mn^{2+}$ ', '6 ', '60'});
xlabel('Substrate rigidity (kPa)','FontSize', 55);
ylabel('{\itF_{BR}} (pN)', 'FontSize',55);
set(gca, 'TickLabelInterpreter', 'latex','FontName', 'Arial', 'FontSize', 50, 'YScale', 'linear');
set(gcf, 'Position',[700 100 700 1000]);
exportgraphics(gca,'FBR4casesextended.png','Resolution',450)
