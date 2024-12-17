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
             Bound_filaments = [];
             Bound_integrins = [];
             for r = 1:nRuns
                            tic
                            SaveName = ['SIMULATION-001__','Ka_',SimFormat(k),'__Peak_',sprintf('%02d',p), '__nLigands_',sprintf('%04d',nligands ), '_run_', sprintf('%02d',r), '.mat'];
                            load(fullfile(RawSaveDirectory, SaveName));
                            for t = 201:501          
                                 percentage_bound_integrins = SimData.nAdhesions(t,1);
                                 v= SimData.AdhesionData(:,5,t);
                                 idx = ~isnan(v);
                                 percentage_bound_filaments= (length(unique(v(idx)))/length(SimData.Data{t,1}.FilamentName))*100;
                                 Bound_filaments = [Bound_filaments;percentage_bound_filaments(:)];
                                 Bound_integrins = [Bound_integrins; percentage_bound_integrins(:)];
 
                             end

                     end

                    SaveName = (['Ligation_vals.','k_a_', SimFormat(k), '_peak_', num2str(p), '_nligands_', sprintf('%03d', nligands),'.mat']);
                    Directory = '/Users/remisondaz/Desktop/MATLAB/Histograms';
                    FullFilePath = fullfile(Directory, SaveName);
                    save(FullFilePath,  'Bound_integrins','Bound_filaments');
                        
end                 
               
  %%
  data1 = load("Ligation_vals.k_a_000d0001_peak_1_nligands_400.mat");
  data2 = load("Ligation_vals.k_a_000d0001_peak_2_nligands_400.mat");
  data3 = load("Ligation_vals.k_a_000d0010_peak_1_nligands_400.mat");
  data4 = load("Ligation_vals.k_a_000d0010_peak_2_nligands_400.mat");
  data5 = load("Ligation_vals.k_a_000d0100_peak_1_nligands_400.mat");
  data6 = load("Ligation_vals.k_a_000d0100_peak_2_nligands_400.mat");


  X = [1 2 3 4 5 6];

  bound_integrins = [data1.Bound_integrins, data2.Bound_integrins, data3.Bound_integrins, ...,
                     data4.Bound_integrins, data5.Bound_integrins, data6.Bound_integrins];
  Medians = median(bound_integrins, "omitnan");
  Xvals = repmat(X, size(bound_integrins,1),1);
  bound_filaments = [data1.Bound_filaments, data2.Bound_filaments, data3.Bound_filaments, data4.Bound_filaments,...
                     data5.Bound_filaments, data6.Bound_filaments];
  Medians2 = median(bound_filaments, "omitnan");
  Xvals2 = repmat(X, size(bound_filaments,1),1);

  FH = figure(1); clf
  FH.Color = 'w';
  MarkerSize = 20;
  SH1 = swarmchart(Xvals,abs(bound_integrins), '.k', 'SizeData',1,'XJitterWidth',0.5,'CData',0.4*[1,1,1]); hold on
  bh = boxplot(bound_integrins(:, [1,2,3,4,5,6]),'Symbol', ' ',  'Notch', 'on');
  %Xlabels = [10e-2,10e-2,10e0,10e0];
  ylabel('Bound Integrins (%)','FontSize', 40)
  set(bh,'LineWidth',3)
  set(gca, 'box', 'off', 'fontsize', 42, 'LineWidth', 3)
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
set(gca, 'TickLabelInterpreter', 'latex','FontName', 'Arial');
xticklabels({'0.4', '$0.4 + Mn^{2+}$','6', '$6 + Mn^{2+}$','60','$60 + Mn^{2+}$'});
xlabel('Substrate Rigidity (kPa)','FontSize', 50);
set(gcf, 'Position',[700 100 700 1000]);

  BH = figure(2);
  BH.Color = 'w';
  MarkerSize = 20;
  SH2 = swarmchart(Xvals,abs(bound_filaments), '.k', 'SizeData',1,'XJitterWidth',0.5,'CData',0.4*[1,1,1]); hold on
  lh = boxplot(bound_filaments(:, [1,2,3,4,5,6]),'Symbol', ' ',  'Notch', 'on');
  set(lh,'LineWidth',3)
  ylabel('Bound F-actin (%)', 'FontSize', 40)
  %Xlabels = [10e-2,10e-2,10e0,10e0];
  set(bh,'LineWidth',3)
  set(gca, 'box', 'off', 'fontsize', 42, 'LineWidth', 3)
  colors = {[0.753, 0.753, 0.753],[1 0 0],[0.502, 0.502, 0.502], [1 0 0], [0.251, 0.251, 0.251], [1 0 0]};
  boxes2 = findobj(gca, 'Tag', 'Box');
  boxes2 = flipud(boxes2);
for i = 1:length(boxes2)
    set(boxes2(i), 'LineWidth', 3, 'Color', colors{i}); % Set the outline color
end
for j = 1:length(boxes2)
    % Use patch to set the face to be transparent with colored edges
    xData = get(boxes2(j), 'XData');
    yData = get(boxes2(j), 'YData');
    patch(xData, yData, colors{j}, 'FaceAlpha', 0.2, 'EdgeColor', colors{j}, 'LineWidth', 3); 
end
 set(gca, 'TickLabelInterpreter', 'latex','FontName', 'Arial');
xticklabels({'0.4', '$0.4 + Mn^{2+}$','6', '$6 + Mn^{2+}$','60','$60 + Mn^{2+}$'});
xlabel('Substrate Rigidity (kPa)','FontSize', 50);
set(gcf, 'Position',[700 100 700 1000]);




           

     
