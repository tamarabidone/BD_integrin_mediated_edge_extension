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
             MemPosi = NaN(501,6);
             Avg = [];
             for r = 1:nRuns
                            tic
                            SaveName = ['SIMULATION-001__','Ka_',SimFormat(k),'__Peak_',sprintf('%02d',p), '__nLigands_',sprintf('%04d',nligands ), '_run_', sprintf('%02d',r), '.mat'];
                            load(fullfile(RawSaveDirectory, SaveName));
                            for t = 451:501
                            MemPosi(t,r) = SimData.MembranePosition(t,1);
                            end
             end
             Avg = mean(MemPosi, 2);

                    SaveName = (['MemPosition.','k_a_', SimFormat(k), '_peak_', num2str(p), '_nligands_', sprintf('%03d', nligands),'.mat']);
                    Directory = '/Users/remisondaz/Desktop/MATLAB/Histograms';
                    FullFilePath = fullfile(Directory, SaveName);
                    save(FullFilePath,  'Avg');
end
  
%%

  data1 = load("MemPosition.k_a_000d0001_peak_1_nligands_400.mat");
  data2 = load("MemPosition.k_a_000d0001_peak_2_nligands_400.mat");
  data3 = load("MemPosition.k_a_000d0010_peak_1_nligands_400.mat");
  data5 = load("MemPosition.k_a_000d0100_peak_1_nligands_400.mat");

nem = {data1.Avg,data2.Avg, data3.Avg,data5.Avg};
% average_over_seconds = @(x) arrayfun(@(i) mean(x(i:min(i+9, numel(x)))), 1:10:numel(x));
% nem = cellfun(average_over_seconds, nem, 'UniformOutput', false);
[Fout, nCells, Ns] = CellArray2PaddedNanArray(nem');
Fout = Fout(:, [1,2,3,4]);
 X = [1 2 3 4];
 Xvals = repmat(X, size(Fout,1),1);

FH = figure(1); clf
FH.Color = 'w';
MarkerSize = 20;
%Xlabels = [10e-2,10e-2,10e0,10e0];
SH1 = swarmchart(Xvals,abs(Fout(:, [1,2,3,4])), '.k', 'SizeData',1,'XJitterWidth',0.5,'CData',0.4*[1,1,1]); hold on
  bh = boxplot(abs(Fout(:, [1,2,3,4])),'symbol', ' ', 'Notch', 'on');
  %ylabel('F-actin Flow (nm/s)','FontSize', 40)
  set(bh,'LineWidth',3)
  medians = findobj(gca, 'Tag', 'Median');
  boxes = findobj(gca, 'Tag', 'Box');
  colors = {[0.753, 0.753, 0.753],[1 0 0],[0.502, 0.502, 0.502],[0.251, 0.251, 0.251]};
  boxes = findobj(gca, 'Tag', 'Box');
  boxes = flipud(boxes);
for i = 1:length(boxes)
    set(boxes(i), 'LineWidth', 3, 'Color', colors{i}); % Set the outline color
end
for j = 1:length(boxes)
    % Use patch to set the face to be transparent with colored edges
    xData = get(boxes(j), 'XData');
    yData = get(boxes(j), 'YData');
    patch(xData, yData, colors{j}, 'FaceAlpha', 0.2, 'EdgeColor', colors{j}, 'LineWidth', 3); 
end
  % ylim([0 25])
  %xticks([1 2 3 4 5 6 7 8 9 10])
  set(gca, 'TickLabelInterpreter', 'latex','FontName', 'Arial');
xticklabels({'0.4', '$0.4 + Mn^{2+}$','6','60'});
xlabel('Substrate Rigidity (kPa)','FontSize', 50);
ylabel('Y-position (nm)', 'FontSize',50);
set(gca, 'box', 'off', 'fontsize', 42, 'LineWidth', 3)
set(gcf, 'Position',[700 100 700 1000]);
hold off
print('FlowMag.svg', '-dsvg', '-r400');
exportgraphics(gca,'FlowMag.svg','Resolution',350)
