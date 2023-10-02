clear
close all
clc


FileSaveDirectory = '/Users/remisondaz/Desktop/MATLAB/Varying_spring_constant/Fraction_active_filaments';

load(fullfile(FileSaveDirectory, 'Filamentresults.mat'));


k_a = [0.0001 0.001 0.01 0.1 1];
peak = [1 2];
nligands = 400;
runs = 3;

start_time=20000;
end_time=29000;


for i = 1:length(peak)
    for j = 1:length(k_a)
        
         
        for r = 1:runs

            Actin_frac_per_run = Filamentresults.(['k_a_', SimFormat(k_a(j)), '_peak_', num2str(peak(i)), '_nligands_', num2str(nligands), 'run_0', num2str(runs)]).Fraction_active_filaments;

            Actin_frac_per_run = Actin_frac_per_run(~isnan(Actin_frac_per_run));

            if r==1

      Active_filaments = Actin_frac_per_run;
      
            else

                Active_filaments(end+1:end+length(Actin_frac_per_run)) = Actin_frac_per_run;
           
            end

            clear Actin_frac_per_run
            
        end
       
        Average_fraction_of_filaments(i,j) = mean(Active_filaments);

    end
end


% % data is a matrix:
% % 1st line Mn2+ conditions, in the order of increasing stiffness (form 
% % 0.0001 to 0.1 pN/nm).
% % 2nd line is WT, in the same order
% 
% data = flipud(Average_fraction_of_filaments(:, :));
% 
% 
% figure(1); clf
% set(gcf,'Color','w')
% 
% 
% 
%          [X,Y] = meshgrid(1:size(data,2), 1:size(data,1));
%          %// Define a finer grid of points
%          [X2,Y2] = meshgrid(1:0.01:size(data,2), 1:0.01:size(data,1));
% 
% 
% 
%         %// Interpolate the data and show the output
%          outData = interp2(X, Y, data, X2, Y2, 'linear');
%          imagesc(outData*100);
%          %// Cosmetic changes for the axes
%          set(gca, 'XTick', linspace(1,size(X2,2),size(X,2)));
%          set(gca, 'YTick', linspace(1,size(X2,1),size(X,1)));
%          set(gca, 'XTickLabel', 1:size(X,2));
%          set(gca, 'YTickLabel', 1:size(X,1));
% 
% 
% %// Add colour bar
% colorbar;
% ylabel('[Mn^2^+]');
%          xlabel('\it{k} (pN/nm)');
%          set (gca,'linewidth', 8, 'fontsize', 50, 'color', [0.5 0.5 0.5])
% 
% xticklabels({'10^-^4', '10^-^3', '10^-^2', '10^-^1', '10^0'});
% yticklabels({' ', ' '});
%          xtickangle(60)
%          axis square
% 
%          set(gcf,'position',[100,100,900,600])
%          CH = colorbar;
%          CH.Label.String = 'Edge filaments (%)';

           r = figure(1);

            r.Color = 'white';

            h1 = semilogx(k_a, Average_fraction_of_filaments(1,:), '-k', 'LineWidth', 4);
            hold on
            h2 = semilogx(k_a, Average_fraction_of_filaments(2,:), '-r', 'LineWidth', 4);
            title('Average Fraction of Active Filaments', 'fontsize', 10);
            xlabel('Spring Constant (pN/nm')
            ylabel('Fraction of active Filaments')
            legend([h1, h2], 'Wild Type', 'Manganese', 'Location','northwest')
            xticks([10^-5 10^-4 10^-3 10^-2 10^-1 10^0 10^1])
            xtickangle(60)
            set (gca,'linewidth', 5, 'fontsize', 15)
            % f.Position = [655 234 650 530];
            axis square
            print('Fraction_of_active_filaments_considering_positive_flow_threshold_15_nm', '-dpng', '-r400')

