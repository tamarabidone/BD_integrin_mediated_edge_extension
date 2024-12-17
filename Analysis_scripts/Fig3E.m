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
             Flow_filaments = [];
             Negative_flow = [];
             for r = 1:nRuns
                            tic
                            SaveName = ['SIMULATION-001__','Ka_',SimFormat(k),'__Peak_',sprintf('%02d',p), '__nLigands_',sprintf('%04d',nligands ), '_run_', sprintf('%02d',r), '.mat'];
                            load(fullfile(RawSaveDirectory, SaveName));
                            for t = 301:501
                                v= SimData.Data{t,1}.YSpeed;
                                idx = find(v > -25 & v < 0);
                                percentage_negative_flow = (length((v(idx)))/length(SimData.Data{t,1}.FilamentName))*100;
                                Negative_flow = [Negative_flow; percentage_negative_flow(:)];
                                % idx2 = find(v>=0);
                                % v(idx2) = 5;
                                flow_filaments = mean(v).*(length(v(idx))/length(v));
                                Flow_filaments = [Flow_filaments;mean(v(:))];
                             end  
                            toc
                     end

                    SaveName = (['Flow_vals.','k_a_', SimFormat(k), '_peak_', num2str(p), '_nligands_', sprintf('%03d', nligands),'.mat']);
                    Directory = '/Users/remisondaz/Desktop/MATLAB/Histograms';
                    FullFilePath = fullfile(Directory, SaveName);
                    save(FullFilePath,  'Flow_filaments', 'Negative_flow');
                        
end
                            
               
  %%
  data1 = load("Flow_vals.k_a_000d0001_peak_1_nligands_400.mat");
  data2 = load("Flow_vals.k_a_000d0001_peak_2_nligands_400.mat");
  data3 = load("Flow_vals.k_a_000d0010_peak_1_nligands_400.mat");
  data4 = load("Flow_vals.k_a_000d0010_peak_2_nligands_400.mat");
  data5 = load("Flow_vals.k_a_000d0100_peak_1_nligands_400.mat");
  data6 = load("Flow_vals.k_a_000d0100_peak_2_nligands_400.mat");
  data7 = load("Flow_vals.k_a_000d1000_peak_1_nligands_400.mat");
  data8 = load("Flow_vals.k_a_000d1000_peak_2_nligands_400.mat");
  data9 = load("Flow_vals.k_a_001d0000_peak_1_nligands_400.mat");
  data10 = load("Flow_vals.k_a_001d0000_peak_2_nligands_400.mat");

  % flow1 = data1.Flow_filaments(:) > -20 ; 
  % flow2 = data2.Flow_filaments(:)> -20 ;
  % flow3 = data3.Flow_filaments(:)> -20 ;
  % flow4 = data4.Flow_filaments(:)> -20 ;
  % flow5 = data5.Flow_filaments(:)> -20 ;
  % flow6 = data6.Flow_filaments(:)> -20 ;
  

  X = [1 2 3 4 5 6];
  Flow = {data1.Flow_filaments; data2.Flow_filaments;data3.Flow_filaments;...
          data4.Flow_filaments;data5.Flow_filaments;data6.Flow_filaments};
  [Fout,nCells,Ns] = CellArray2PaddedNanArray(Flow);
  Medians = median(Fout, "omitnan");
  Xvals = repmat(X, size(Fout,1),1);
  Percentage = {data1.Negative_flow; data2.Negative_flow; data3.Negative_flow; data4.Negative_flow;...
                data5.Negative_flow;data6.Negative_flow};
  [Pout,nCells,Ns] = CellArray2PaddedNanArray(Percentage);
  Medians2 = median(Pout, "omitnan");
  Xvals2 = repmat(X, size(Pout,1),1);


FH = figure(1); clf
FH.Color = 'w';
MarkerSize = 20;
%Xlabels = [10e-2,10e-2,10e0,10e0];
SH1 = swarmchart(Xvals,abs(Fout(:, [1,2,3,4,5,6])), '.k', 'SizeData',1,'XJitterWidth',0.5,'CData',0.4*[1,1,1]); hold on
  bh = boxplot(abs(Fout(:, [1,2,3,4,5,6])),'symbol', ' ', 'Notch', 'on');
  %ylabel('F-actin Flow (nm/s)','FontSize', 40)
  set(bh,'LineWidth',3)
  medians = findobj(gca, 'Tag', 'Median');
  boxes = findobj(gca, 'Tag', 'Box');
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
  ylim([0 3])
  %xticks([1 2 3 4 5 6 7 8 9 10])
  set(gca, 'TickLabelInterpreter', 'latex','FontName', 'Arial');
%xticklabels({'0.4', '$0.4 + Mn^{2+}$','6', '$6 + Mn^{2+}$','60','$60 + Mn^{2+}$'});
%xlabel('Substrate Rigidity (kPa)','FontSize', 50);
%ylabel('v_{flow} (nm/s)', 'FontSize',50);
set(gca, 'box', 'off', 'fontsize', 42, 'LineWidth', 3)
set(gcf, 'Position',[700 100 700 1000]);
hold off
%print('FlowMag.svg', '-dsvg', '-r400');
%exportgraphics(gca,'FlowMag.svg','Resolution',350)

  BH = figure(2); 
  BH.Color = 'w';
MarkerSize = 20;
SH2 = swarmchart(Xvals2,Pout, '.k', 'SizeData',1,'XJitterWidth',0.5,'CData',0.4*[1,1,1]); hold on
  ph = boxplot(Pout(:, [1,2,3,4,5,6]), 'Symbol', ' ',  'Notch', 'on');
  set(ph,'LineWidth',3)
  set(gca, 'box', 'off', 'fontsize', 42, 'LineWidth', 3)
  ylabel('F-actin at v_{flow} > 0 (%)', 'FontSize', 40)
  boxes2 = findobj(gca, 'Tag', 'Box');
  boxes2 = flipud(boxes2);
colors = {[0.753, 0.753, 0.753],[1 0 0],[0.502, 0.502, 0.502], [1 0 0], [0.251, 0.251, 0.251], [1 0 0]};
for i = 1:length(boxes2)
    set(boxes2(i), 'LineWidth', 3, 'Color', colors{i});
end
for j = 1:length(boxes2)
    xData = get(boxes2(j), 'XData');
    yData = get(boxes2(j), 'YData');
    patch(xData, yData, colors{j}, 'FaceAlpha', 0.2, 'EdgeColor', colors{j}, 'LineWidth', 3); 
end
  set(gca, 'TickLabelInterpreter', 'latex','FontName', 'Arial');
 xticklabels({'0.4', '$0.4 + Mn^{2+}$','6', '$6 + Mn^{2+}$','60','$60 + Mn^{2+}$'});
xlabel('Substrate Rigidity (kPa)','FontSize', 50);
  set(gca, 'TickLabelInterpreter', 'latex','FontName', 'Arial');
  set(gcf, 'Position',[700 100 700 1000]);





