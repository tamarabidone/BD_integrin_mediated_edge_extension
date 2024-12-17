clc 
clear
close all

RawSaveDirectory = '/Users/remisondaz/Desktop/MATLAB/Extended_Sims';

% nRuns = 3;
% k_a = [0.0001 0.001 0.01];
% peak = [1 2];
% nligands = 400;
% f = [1 2];


nRuns = 6; % Total number of runs for each condition
k_a = [0.0001 ,0.0001,0.001, 0.01, 0.1]; % adhesion spring constant
peak = [1, 2, 1, 1,1];   % WT or Mn
nLigands = 400;

    for m = 1:length(k_a) % Create combinations of all conditions
        p = peak(m);
        k = k_a(m);
        Nem_order = [];
            for r = 1:nRuns
  % for m = 1:length(nligands) % Create combinations of all conditions
  %       for n = 1:length(peak)
  %           for k = 1:length(k_a)
  %                    for r = 1:nRuns
                            SaveName = ['SIMULATION-001__','Ka_',SimFormat(k),'__Peak_',sprintf('%02d',p), '__nLigands_',sprintf('%04d',nLigands), '_run_', sprintf('%02d',r), '.mat'];
                            load(fullfile(RawSaveDirectory, SaveName));
                            nematic_order = [];
                            tic
                             for t = 300:501
                                 Orientation = atan2d(SimData.Data{t,1}.Orientation(:,2),SimData.Data{t,1}.Orientation(:,1));
                                 Length = SimData.Data{t,1}.FilamentLength(:,1);
                                 Orient = Orientation;
                                 W_O = Orient.*Length;
                                 Tot_L = sum(Length);
                                 Orientation_per_run = sum(W_O)/((Tot_L));
                                 theta = abs(Orientation_per_run-Orient);
                                 nematic_order = (3 * cosd(theta).^2 - 1 )/ 2;
                                 W_N = nematic_order .* Length;
                                 W_N_T = sum(W_N)/Tot_L;
                                 Nem_order = [Nem_order;mean(nematic_order)];
                             end
                     end

                    SaveName = (['Nematic_vals.','k_a_', SimFormat(k), '_peak_', sprintf('%01d',p), '_nligands_', sprintf('%03d', nLigands),'.mat']);
                    Directory = '/Users/remisondaz/Desktop/MATLAB/Histograms';
                    FullFilePath = fullfile(Directory, SaveName);
                    save(FullFilePath,  'Nem_order');
                        
            end
  %       end
  % end

  %%
  data1 = load("Nematic_vals.k_a_000d0000_peak_1_nligands_400.mat");
  data2 = load("Nematic_vals.k_a_000d0000_peak_2_nligands_400.mat");
  data3 = load("Nematic_vals.k_a_000d0001_peak_1_nligands_400.mat");
  data4 = load("Nematic_vals.k_a_000d0001_peak_2_nligands_400.mat");
  data5 = load("Nematic_vals.k_a_000d0010_peak_1_nligands_400.mat");
  data6 = load("Nematic_vals.k_a_000d0010_peak_2_nligands_400.mat");
  data7 = load("Nematic_vals.k_a_000d0100_peak_1_nligands_400.mat");
  data8 = load("Nematic_vals.k_a_000d0100_peak_2_nligands_400.mat");
  data9 = load("Nematic_vals.k_a_000d1000_peak_1_nligands_400.mat");
  data10 = load("Nematic_vals.k_a_000d1000_peak_2_nligands_400.mat");
  data11 = load("Nematic_vals.k_a_001d0000_peak_1_nligands_400.mat");
  data12 = load("Nematic_vals.k_a_001d0000_peak_2_nligands_400.mat");

  
  nem = {unique(data1.Nem_order), unique(data2.Nem_order), data3.Nem_order, ...,
  data4.Nem_order, data5.Nem_order, unique(data6.Nem_order),  ...
  data7.Nem_order, unique(data8.Nem_order),data9.Nem_order,unique(data10.Nem_order),...
  unique(data11.Nem_order),unique(data12.Nem_order)};
  [Fout,nCells,Ns] = CellArray2PaddedNanArray(nem');
Fout = Fout(:, [3,4,5,7]);
FH = figure(1); clf

% subplot(1,3,3)
FH.Color = 'w';
MarkerSize = 20;
 bh = violinplot(Fout);
colors = {[0, 0, 0], [1, 0.75, 0.796], [0.678, 0.847, 0.902], [0.502, 0, 0.502]};  % Example colors
% Loop through each violin and set the color
for f = 1:length(bh)
    bh(f).ViolinColor = colors(f);
    bh(f).ShowData = false;
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
FS = 45;
 %bh.LineWidth = 3;
%ylim([0 1])
set(gca, 'box', 'off', 'fontsize', FS, 'LineWidth', 3)
set(gca, 'TickLabelInterpreter', 'latex','FontName', 'Arial', 'FontSize', FS);
xticks([1 2 3 4])
%xticklabels({'0.4', '0.4\\newline$+ Mn^{2+}$', '6', '60'});
xticklabels({' ', ' ', ' ', ' '});
%set(gca, 'TickLabelInterpreter', 'none');  % Disable LaTeX interpretation for xticklabels
xlabel('Substrate rigidity (kPa)','FontSize', FS);
ylabel('Order Parameter', 'FontSize',FS);
set(gcf, 'Position',[700 100 700 1000]);
exportgraphics(gca,'ExtendedSims2.png','Resolution',450)





