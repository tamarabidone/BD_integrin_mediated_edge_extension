clc 
clear
close all

RawSaveDirectory = '/Users/remisondaz/Desktop/MATLAB/Review_sims/Review_sims';

nRuns = 40;
k_a = [0.0001 0.0001 0.01];
peak = [1 1 1];
flow = [true false false];
k_int = [0.01 0.1 1];
nligands = 400;
% f = [1 2];


    for m = 1:length(k_a) % Create combinations of all conditions
        % for p = 1:length(peak)
        %     for v = 1:length(k_int)
        p = peak(m);
        k = k_a(m);
        f = flow(m);
        flow_str = char(string(f));
        MemPos = NaN(14,1);
        %Nem_order = NaN(9000, 39);
        Mean_nem_order = NaN(9000, 40);
        Median_nem_order = NaN(9000, 40);
        Mode_nem_order = NaN(9000, 40);
        Min_nem_order = NaN(9000, 40);
        Max_nem_order = NaN(9000, 40);
        mean_traction_stress = NaN(9000, 40);
        mean_integrins = NaN(9000, 40);
        mean_actins = NaN(9000, 40);
        ret_flow = NaN(9000, 40);
        branches = NaN(9000, 40);
        force_mem = NaN(9000, 40);

            for r = 1:nRuns
                            SaveName = ['SIMULATION-001__','Ka_',SimFormat(k),'__Peak_',sprintf('%02d',p), '__nLigands_',sprintf('%04d',nligands), '_run_', sprintf('%02d',r), '.mat'];
                            %SaveName = ['SIMULATION-001__','Ks_',SimFormat(k),'Myos_true','__Peak_',sprintf('%02d',p), '_run_', sprintf('%02d',r), '.mat'];
                            %SaveName = ['SIMULATION-004__','Ks_',SimFormat(k_a(m)),'__Peak_',sprintf('%02d',peak(p)),'__kInt_', strrep(num2str(k_int(v)), '.', 'p') '_run_', sprintf('%02d',r), '.mat'];
                            %SaveName = ['SIMULATION-001__','Ks_',SimFormat(k_a(m)),'__ExpLife_',sprintf('%02d',peak(p)),'__ReducedFlow_', flow_str, '_run_', sprintf('%02d',r), '.mat'];
                            load(fullfile(RawSaveDirectory, SaveName));
                            nematic_order = [];
                             for t = 1:3000
                                 % Orientation = atan2d(SimData.Data{t,1}.Orientation(:,2),SimData.Data{t,1}.Orientation(:,1));
                                 % Length = SimData.Data{t,1}.FilamentLength(:,1);
                                 % Orient = Orientation;
                                 % W_O = Orient.*Length;
                                 % Tot_L = sum(Length);
                                 % Orientation_per_run = sum(W_O)/((Tot_L));
                                 % theta = abs(Orientation_per_run-Orient);
                                 % nematic_order = (3 * cosd(theta).^2 - 1 )/ 2;
                                 % W_N = nematic_order .* Length;
                                 % W_N_T = sum(W_N)/Tot_L;
                                 % Mean_nem_order(t,r) = mean(W_N_T);
                                 % Median_nem_order(t,r) = median(W_N_T);
                                 % Mode_nem_order(t,r) = mode(W_N_T);
                                 % Min_nem_order(t,r) = min(W_N_T);
                                 % Max_nem_order(t,r) = max(W_N_T);
                                 % 
                                 % mean_traction_stress(t,r) = (sum(SimData.AdhesionData(:,4,t))*1e-12)/(2.5e-13);
                                 mean_integrins(t,r) = length(find(SimData.AdhesionData(:,3,t) == 1));
                                 % actins_bound = SimData.AdhesionData(:,5,t);      
                                 % actins_bound = actins_bound(~isnan(actins_bound)); 
                                 % mean_actins(t,r) = 100 * length(unique(actins_bound)) / ...
                                 %                          length(SimData.Data{t,1}.Parent);
                                 % ret_flow(t,r) = 100*length(find(SimData.Data{t,1}.YSpeed < 0))/(length(SimData.Data{t,1}.YSpeed));
                                 % branches(t,r) = length(find(SimData.Data{t,1}.Parent(:)~=0))/(length(SimData.Data{t,1}.Parent(:)));
                                 % force_mem(t,r) = sum(SimData.Data{t,1}.ForceonMembrane);
                             end
                             MemPos(r) = mean(SimData.MembranePosition(2901:3001));
            end
                    %Mean_nem_order = (mean(Nem_order,1,'omitnan'))';
                    SaveName = (['FinalMem_val.','k_a_', SimFormat(k), '_peak_', sprintf('%01d',p), '_nligands_', sprintf('%03d', nligands),'.mat']);
                    %SaveName = (['Activ_val.','Ks_',SimFormat(k_a(m)),'__Peak_',sprintf('%02d',peak(p)),'__kInt_', strrep(num2str(k_int(v)), '.', 'p'), '.mat']);
                    %SaveName = (['Life_val.','Ks_',SimFormat(k_a(m)),'__ExpLife_',sprintf('%02d',peak(p)),'__ReducedFlow_', flow_str,'.mat']);
                    Directory = '/Users/remisondaz/Desktop/MATLAB/Review';
                    FullFilePath = fullfile(Directory, SaveName);
                    save(FullFilePath,  'MemPos', 'mean_integrins');
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
