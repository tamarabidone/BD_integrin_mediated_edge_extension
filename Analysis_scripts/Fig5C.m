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
             struct = [];
             for r = 1:nRuns
                            tic
                            SaveName = ['SIMULATION-001__','Ka_',SimFormat(k),'__Peak_',sprintf('%02d',p), '__nLigands_',sprintf('%04d',nligands ), '_run_', sprintf('%02d',r), '.mat'];
                            load(fullfile(RawSaveDirectory, SaveName));
                            for t = 300:501
                                 structNames = SimData.Data{t,1}.StructName;
                                 bins = 0.5:1:length(structNames)+0.5;
                                 BinCenters = floor(bins(2:end));
                                 counts = histcounts(structNames,bins)';
                                 cIdx = find(counts > 1);
                                 Lia = ismember(structNames,BinCenters(cIdx));
                                 struct=[struct;length(SimData.Data{t,1}.StructName(Lia))/length(SimData.Data{t,1}.StructName)];
                            end 
                     end
                    SaveName = (['Struct_vals.','k_a_', SimFormat(k), '_peak_', num2str(p), '_nligands_', sprintf('%03d', nligands),'.mat']);
                    Directory = '/Users/remisondaz/Desktop/MATLAB/Histograms';
                    FullFilePath = fullfile(Directory, SaveName);
                    save(FullFilePath, 'struct');
                    toc
            end

%%

  data1 = load("Struct_vals.k_a_000d0000_peak_1_nligands_400.mat");
  data2 = load("Struct_vals.k_a_000d0000_peak_2_nligands_400.mat");
  data3 = load("Struct_vals.k_a_000d0001_peak_1_nligands_400.mat");
  data4 = load("Struct_vals.k_a_000d0001_peak_2_nligands_400.mat");
  data5 = load("Struct_vals.k_a_000d0010_peak_1_nligands_400.mat");
  data6 = load("Struct_vals.k_a_000d0010_peak_2_nligands_400.mat");
  data7 = load("Struct_vals.k_a_000d0100_peak_1_nligands_400.mat");
  data8 = load("Struct_vals.k_a_000d0100_peak_2_nligands_400.mat");
  data9 = load("Struct_vals.k_a_000d1000_peak_1_nligands_400.mat");
  data10 = load("Struct_vals.k_a_000d1000_peak_2_nligands_400.mat");
  data11 = load("Struct_vals.k_a_001d0000_peak_1_nligands_400.mat");
  data12 = load("Struct_vals.k_a_001d0000_peak_2_nligands_400.mat");

  idx3 = data3.struct>0;
  struct3 = data3.struct(idx3)-1;
  idx4 = data4.struct>0;
  struct4 = data4.struct(idx4)-1;
  idx5 = data5.struct > 0;
  struct5 = data5.struct(idx5)-1;
  idx7 = data7.struct > 0;
  struct7 = data7.struct(idx7)-1;


  nem = {data1.struct, data2.struct,data3.struct,data4.struct,...
         data5.struct,data6.struct,data7.struct,data8.struct,...
         unique(data9.struct),unique(data10.struct),unique(data11.struct),unique(data12.struct)};

[Fout,nCells,Ns] = CellArray2PaddedNanArray(nem');
FH = figure(1); clf
FH.Color = 'w';
MarkerSize = 20;
 bh = violinplot(Fout(:, [3,4,5,7]));
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
  set(gca, 'box', 'off', 'fontsize', 42, 'LineWidth', 6)
  xticks([1 2 3 4])
  %ylim([0 1])
  set(gca, 'TickLabelInterpreter', 'latex','FontName', 'Arial','FontSize', 45);
xticklabels({' ', ' ', ' ', ' '});
xlabel('Substrate rigidity (kPa)', 'FontSize', 55);
  ylabel({'Ratio of branched F-actin'}, 'FontSize', 45);
  set(gcf, 'Position',[700 100 800 1000]);
  exportgraphics(gca,'RatioBranched.png','Resolution',450)
