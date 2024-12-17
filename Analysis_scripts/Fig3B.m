clc 
clear
close all

RawSaveDirectory = '/Users/remisondaz/Desktop/MATLAB/Extended_Sims';

nRuns = 6;
nIntegrins = 100;
k_a =  [0.0001, 0.0001, 0.001, 0.01,0.001, 0.01];
peak = [1, 2, 1, 1,2,2];
nligands = 400;

for n = 1:length(k_a)
             p = peak(n);
             k = k_a(n);
             force_adhesions = [];
             for r = 1:nRuns
                            tic
                            SaveName = ['SIMULATION-001__','Ka_',SimFormat(k),'__Peak_',sprintf('%02d',p), '__nLigands_',sprintf('%04d',nligands ), '_run_', sprintf('%02d',r), '.mat'];
                            load(fullfile(RawSaveDirectory, SaveName));
                            for t = 301:501          
                                 force_index = find(SimData.AdhesionData(:,4,t)~=0);
                                 force_adhesions = [force_adhesions;SimData.AdhesionData(force_index,4,t)];
                             end  
                            toc
                     end

                    SaveName = (['Adhesion_vals.','k_a_', SimFormat(k), '_peak_', num2str(p), '_nligands_', sprintf('%03d', nligands),'.mat']);
                    Directory = '/Users/remisondaz/Desktop/MATLAB/Histograms';
                    FullFilePath = fullfile(Directory, SaveName);
                    save(FullFilePath,  'force_adhesions');
                        
end
                            
               
  %%
  data1 = load("Adhesion_vals.k_a_000d0001_peak_1_nligands_400.mat");
  data2 = load("Adhesion_vals.k_a_000d0001_peak_2_nligands_400.mat");
  data3 = load("Adhesion_vals.k_a_000d0010_peak_1_nligands_400.mat");
  data4 = load("Adhesion_vals.k_a_000d0010_peak_2_nligands_400.mat");
  data5 = load("Adhesion_vals.k_a_000d0100_peak_1_nligands_400.mat");
  data6 = load("Adhesion_vals.k_a_000d0100_peak_2_nligands_400.mat");

  X = [1 2 3 4 5 6];
  Flow = {data1.force_adhesions; data2.force_adhesions;data3.force_adhesions;...
          data4.force_adhesions;data5.force_adhesions;data6.force_adhesions};
  [Fout,nCells,Ns] = CellArray2PaddedNanArray(Flow);
  Medians = median(Fout, "omitnan");
  Xvals = repmat(X, size(Fout,1),1);


FH = figure(1); clf
FH.Color = 'w';
MarkerSize = 20;
%Xlabels = [10e-2,10e-2,10e0,10e0];
SH1 = swarmchart(Xvals,Fout, '.k', 'SizeData',1,'XJitterWidth',0.5,'CData',0.4*[1,1,1]); hold on
  bh = boxplot(Fout(:, [1,2,3,4,5,6]),'symbol', ' ', 'Notch', 'on');
  set(bh,'LineWidth',3)
  medians = findobj(gca, 'Tag', 'Median');
  set(gca, 'box', 'off', 'fontsize', 42, 'LineWidth', 3)
  boxes = findobj(gca, 'Tag', 'Box');
  % Set individual colors for each boxplot
  colors = {[0.753, 0.753, 0.753],[1 0 0],[0.502, 0.502, 0.502], [1 0 0], [0.251, 0.251, 0.251], [1 0 0]};
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
  %ylim([0 25])
  %xticks([1 2 3 4 5 6 7 8 9 10])
  set(gca, 'TickLabelInterpreter', 'latex','FontName', 'Arial');
 xticklabels({'0.4', '$0.4 + Mn^{2+}$','6', '$6 + Mn^{2+}$','60','$60 + Mn^{2+}$'});
xlabel('Substrate Rigidity (kPa)','FontSize', 50);
ylabel('F_{sub} (pN)','FontSize', 50);
set(gcf, 'Position',[700 100 700 1000]);
hold off
print('AdhesionForce.svg', '-dsvg', '-r200');
%exportgraphics(gca,'Flow4caseslastsec.png','Resolution',350)
