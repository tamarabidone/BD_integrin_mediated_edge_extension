clear
close all
clc


FileSaveDirectory = '/Users/remisondaz/Desktop/MATLAB/Varying_spring_constant/Adhesion_force';

load(fullfile(FileSaveDirectory, 'IntegrinForceresults.mat'));


k_a = [0.0001 0.001 0.01 0.1 1];
peak = [1 2];
nligands = 400;
runs = 3;

% for i = 1:length(peak)
%     for j = 1:3
% 
%         for r = 1:runs
% 
%         ka_small_Adhesion_force = Forceresults.(['k_a_', SimFormat(k_a(j)), '_peak_', num2str(peak(i)), '_nligands_', num2str(nligands), 'run_0', num2str(r)]).Adhesion_force;
%         ka_medium_Adhesion_force = Forceresults.(['k_a_', SimFormat(k_a(j+1)), '_peak_', num2str(peak(i)), '_nligands_', num2str(nligands), 'run_0', num2str(r)]).Adhesion_force;
%         ka_large_Adhesion_force = Forceresults.(['k_a_', SimFormat(k_a(j+2)), '_peak_', num2str(peak(i)), '_nligands_', num2str(nligands), 'run_0', num2str(r)]).Adhesion_force;
% 
%         ka_small_Adhesion_force =   ka_small_Adhesion_force(ka_small_Adhesion_force ~= 0);
%         ka_medium_Adhesion_force =   ka_medium_Adhesion_force(ka_medium_Adhesion_force ~= 0);
%         ka_large_Adhesion_force = ka_large_Adhesion_force(ka_large_Adhesion_force ~= 0);
% 
%         if r == 1 
%             Small_force_per_run = ka_small_Adhesion_force;
%             Medium_force_per_run = ka_medium_Adhesion_force; 
%             Large_force_per_run = ka_large_Adhesion_force;
%         else
%             Small_force_per_run(end+1:end+length(ka_small_Adhesion_force)) = ka_small_Adhesion_force;
%             Medium_force_per_run(end+1:end+length(ka_medium_Adhesion_force)) = ka_medium_Adhesion_force;
%             Large_force_per_run(end+1:end+length(ka_large_Adhesion_force)) = ka_large_Adhesion_force;
%         end
% 
%         clear ka_small_Adhesion_force
%         clear ka_medium_Adhesion_force
%         clear ka_large_Adhesion_force
% 
%         end
% 
%         [counts1, binCenters1] = histcounts(Small_force_per_run, 20);
%         [counts2, binCenters2] = histcounts(Medium_force_per_run, 20);
%         [counts3, binCenters3] = histcounts(Large_force_per_run, 20);
% 
%         figure;
%         plot(binCenters1(1:end-1), counts1, 'r-', 'LineWidth',3);
%         hold on;
%         plot(binCenters2(1:end-1), counts2, 'g-', 'LineWidth',3);
%         plot(binCenters3(1:end-1), counts3, 'b-','LineWidth',3);
% 
%         xlabel('Adhesion Force (pN)');
%         ylabel('Frequency');
%         legend(['k_a = ' num2str(k_a(j))], ['k_a = ' num2str(k_a(j+1))], ['k_a = ' num2str(k_a(j+2))]);
%         title(['Histogram of Adhesion Force for peak = ' num2str(peak(i))]);
%         % xtickangle(60)
%         set (gca,'linewidth', 5, 'fontsize', 15)
%         % axis([0 1200 0 200])
%         % f.Position = [655 234 650 530];
%         %axis square
%         print(['Integrinforce_', 'k_a_', SimFormat(k_a(j)), '_vs_' 'k_a_', SimFormat(k_a(j+1)), '_vs_' 'k_a_', SimFormat(k_a(j+2)) '_peak_', num2str(peak(i))], '-dpng', '-r400')
% 
%     end
% end

for i = 1:length(k_a)
    for j = 1:length(peak)

          for r = 1:runs

        Adhesion_f = Forceresults.(['k_a_', SimFormat(k_a(i)), '_peak_', num2str(peak(j)), '_nligands_', num2str(nligands), 'run_0', num2str(r)]).Adhesion_force;
        

        Adhesion_f =   Adhesion_f(Adhesion_f ~= 0);
        
        if r == 1 

            Adhesion_force_per_run = Adhesion_f;
           
        else
            Adhesion_force_per_run(end+1:end+length(Adhesion_f)) = Adhesion_f;
           
        end

        clear Adhesion_f
        
          end

        Mean_Force(i,j) = mean(Adhesion_force_per_run);
        STD(i,j) = std(Adhesion_force_per_run)/sqrt(length(Adhesion_force_per_run));

    end
end

 % f = figure(2);
 %        f.Color = 'white';
 %        h1 = semilogx(k_a, Mean_Force(:,1), 'o-k', 'linewidth', 4);
 %        hold on
 %        h2 = semilogx(k_a, Mean_Force(:,2), '-r', 'linewidth', 4);
 % 
 % 
 % 
 %        xlabel('Spring Constant (pN/nm)');
 %        ylabel('Average Force (pN)');
 %        legend([h1, h2], 'Wild Type', 'Manganese', 'Location','northwest')
 %        set (gca,'linewidth', 5, 'fontsize', 30)
 %        xticks([10^-5 10^-4 10^-3 10^-2 10^-1 10^0 10])
 %        axis([10^-5 10 0 15])
 %        set (gca,'linewidth', 5, 'fontsize', 30)
 %        xtickangle(60)
 %        title('Adhesion Force vs Spring Constant');
 %        print('MeanIntegrinforce_vs_spring_constant', '-dpng', '-r400')



 r = figure(3);
    r.Color = 'white';
         f1 = semilogx(k_a, STD(:,1), ' -k', 'linewidth', 4);
         hold on 
         f2 = semilogx(k_a, STD(:,2), ' -r', 'linewidth', 4);

        xlabel('Spring Constant (pN/nm)');
        ylabel('Standard error of force (pN)');
        legend([f1, f2], 'Wild Type', 'Manganese', 'Location','northwest')
        set (gca,'linewidth', 5, 'fontsize', 30)
        xticks([10^-5 10^-4 10^-3 10^-2 10^-1 10^0 10])
        axis([10^-5 10 0 0.02])
        set (gca,'linewidth', 5, 'fontsize', 30)
        xtickangle(60)
        title('Force Standard Error vs Spring Constant');
        print('IntegrinSTD_vs_spring_constant', '-dpng', '-r400')


