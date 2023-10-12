clear
close all
clc

%This is the folder where we save the outputs from this analysis 
FileSaveDirectory = '/Users/remisondaz/Desktop/MATLAB/Varying_spring_constant/Adhesion_position';

%This folder contains model output after running simulations 
load(fullfile(FileSaveDirectory, 'IntegrinPositionresults.mat'));


k_a = [0.0001 0.001 0.01 0.1 1];
peak = [1 2];
nligands = 400;
runs = 3;

%Comparing three successive spring constant to see how the position of
%adhesion varies. This segment will generate figures plotting histograms for the dirstibution of adhesions positions in the model over the 30 seconds for three successive spring constants.

% for i = 1:length(peak)
%     for j = 1:3
% 
% 
%         for r = 1:length(runs)
% 
%         ka_small_Adhesion_position = Positionresults.(['k_a_', SimFormat(k_a(j)), '_peak_', num2str(peak(i)) '_nligands_', num2str(nligands), 'run_0', num2str(r)]).Adhesion_position;
%         ka_medium_Adhesion_position = Positionresults.(['k_a_', SimFormat(k_a(j+1)), '_peak_', num2str(peak(i)), '_nligands_', num2str(nligands), 'run_0', num2str(r)]).Adhesion_position;
%         ka_large_Adhesion_position = Positionresults.(['k_a_', SimFormat(k_a(j+2)), '_peak_', num2str(peak(i)), '_nligands_', num2str(nligands), 'run_0', num2str(r)]).Adhesion_position;
% 
%         ka_small_Adhesion_position =   ka_small_Adhesion_position(ka_small_Adhesion_position ~= 0);
%         ka_medium_Adhesion_position =   ka_medium_Adhesion_position(ka_medium_Adhesion_position ~= 0);
%         ka_large_Adhesion_position = ka_large_Adhesion_position(ka_large_Adhesion_position ~= 0);
% 
%         if r == 1
%             Small_position = ka_small_Adhesion_position;
%             Medium_position = ka_medium_Adhesion_position;
%             Large_position = ka_large_Adhesion_position;
% 
%         else
% 
%             Small_position(end+1:end+length(ka_small_Adhesion_position)) = ka_small_Adhesion_position;
%             Medium_position(end+1:end+length(ka_medium_Adhesion_position)) = ka_medium_Adhesion_position;
%             Large_position(end+1:end+length(ka_large_Adhesion_position)) = ka_large_Adhesion_position;
%         end
% 
%         clear ka_small_Adhesion_position
%         clear ka_medium_Adhesion_position
%         clear ka_large_Adhesion_position
% 
%         end
% 
% 
%         [counts1, binCenters1] = histcounts(Small_position, 20);
%         [counts2, binCenters2] = histcounts(Medium_position, 20);
%         [counts3, binCenters3] = histcounts(Large_position, 20);
% 
%         figure;
%         plot(binCenters1(1:end-1), counts1, 'r-', 'LineWidth',3);
%         hold on;
%         plot(binCenters2(1:end-1), counts2, 'g-', 'LineWidth',3);
%         plot(binCenters3(1:end-1), counts3, 'b-','LineWidth',3);
% 
% 
%         xlabel('Adhesion Distance from Membrane (nm)');
%         ylabel('Frequency');
%         set (gca,'linewidth', 5, 'fontsize', 15)
%         legend(['k_a = ' num2str(k_a(j))], ['k_a = ' num2str(k_a(j+1))], ['k_a = ' num2str(k_a(j+2))]);
%         title(['Histogram of Adhesion Position for peak = ' num2str(peak(i))]);
%         print(['Integrinposition_', 'k_a_', SimFormat(k_a(j)), '_vs_' 'k_a_', SimFormat(k_a(j+1)), '_vs_' , SimFormat(k_a(j+2)) '_peak_', num2str(peak(i)) ], '-dpng', '-r400')
% 
%     end
% end

%Comparing two successive spring constants to see how the position difffers. This segment of code will generate histograms of the distribution of adhesions positions from the leading edge for two successive spring constants.

% for i = 1:length(peak)
%     for j = 1:3
% 
% 
%         for r = 1:length(runs)
% 
%         ka_small_Adhesion_position = Positionresults.(['k_a_', SimFormat(k_a(j)), '_peak_', num2str(peak(i)) '_nligands_', num2str(nligands), 'run_0', num2str(r)]).Adhesion_position;
%         ka_large_Adhesion_position = Positionresults.(['k_a_', SimFormat(k_a(j+2)), '_peak_', num2str(peak(i)), '_nligands_', num2str(nligands), 'run_0', num2str(r)]).Adhesion_position;
% 
%         ka_small_Adhesion_position =   ka_small_Adhesion_position(ka_small_Adhesion_position ~= 0);
%         ka_large_Adhesion_position = ka_large_Adhesion_position(ka_large_Adhesion_position ~= 0);
% 
%         if r == 1
%             Small_position = ka_small_Adhesion_position;
%             Large_position = ka_large_Adhesion_position;
% 
%         else
% 
%             Small_position(end+1:end+length(ka_small_Adhesion_position)) = ka_small_Adhesion_position;
%             Large_position(end+1:end+length(ka_large_Adhesion_position)) = ka_large_Adhesion_position;
%         end
% 
%         clear ka_small_Adhesion_position
%         clear ka_large_Adhesion_position
% 
%         end
% 
% 
%         [counts1, binCenters1] = histcounts(Small_position, 20);
%         [counts3, binCenters3] = histcounts(Large_position, 20);
% 
%         figure;
%         plot(binCenters1(1:end-1), counts1, 'r-', 'LineWidth',3);
%         hold on;
%         plot(binCenters3(1:end-1), counts3, 'b-','LineWidth',3);
% 
% 
%         xlabel('Adhesion Distance from Membrane (nm)');
%         ylabel('Frequency');
%         set (gca,'linewidth', 5, 'fontsize', 15)
%         legend(['k_a = ' num2str(k_a(j))], ['k_a = ' num2str(k_a(j+2))]);
%         title(['Histogram of Adhesion Position for peak = ' num2str(peak(i))]);
%         print(['Integrinposition_', 'k_a_', SimFormat(k_a(j)), '_vs_' 'k_a_', SimFormat(k_a(j+2)), '_peak_', num2str(peak(i)) ], '-dpng', '-r400')
% 
%     end
% end

% Similarly, this segment of code will generate figures that plot the histograms of adhesion position distribution for bigger differences in spring constants. For instance here we will be comparing differences of 3 in the list of 5 spring constants

% for i = 1:length(peak)
%     for j = 1:2
% 
% 
%         for r = 1:length(runs)
% 
%         ka_small_Adhesion_position = Positionresults.(['k_a_', SimFormat(k_a(j)), '_peak_', num2str(peak(i)) '_nligands_', num2str(nligands), 'run_0', num2str(r)]).Adhesion_position;
%         ka_large_Adhesion_position = Positionresults.(['k_a_', SimFormat(k_a(j+3)), '_peak_', num2str(peak(i)), '_nligands_', num2str(nligands), 'run_0', num2str(r)]).Adhesion_position;
% 
%         ka_small_Adhesion_position =   ka_small_Adhesion_position(ka_small_Adhesion_position ~= 0);
%         ka_large_Adhesion_position = ka_large_Adhesion_position(ka_large_Adhesion_position ~= 0);
% 
%         if r == 1
%             Small_position = ka_small_Adhesion_position;
%             Large_position = ka_large_Adhesion_position;
% 
%         else
% 
%             Small_position(end+1:end+length(ka_small_Adhesion_position)) = ka_small_Adhesion_position;
%             Large_position(end+1:end+length(ka_large_Adhesion_position)) = ka_large_Adhesion_position;
%         end
% 
%         clear ka_small_Adhesion_position
%         clear ka_large_Adhesion_position
% 
%         end
% 
% 
%         [counts1, binCenters1] = histcounts(Small_position, 20);
%         [counts3, binCenters3] = histcounts(Large_position, 20);
% 
%         figure;
%         plot(binCenters1(1:end-1), counts1, 'r-', 'LineWidth',3);
%         hold on;
%         plot(binCenters3(1:end-1), counts3, 'b-','LineWidth',3);
% 
% 
%         xlabel('Adhesion Distance from Membrane (nm)');
%         ylabel('Frequency');
%         set (gca,'linewidth', 5, 'fontsize', 15)
%         legend(['k_a = ' num2str(k_a(j))], ['k_a = ' num2str(k_a(j+3))]);
%         title(['Histogram of Adhesion Position for peak = ' num2str(peak(i))]);
%         print(['Integrinposition_', 'k_a_', SimFormat(k_a(j)), '_vs_' 'k_a_', SimFormat(k_a(j+3)), '_peak_', num2str(peak(i)) ], '-dpng', '-r400')
% 
%     end
% end


 % Comparing a difference of 4 in spring constants to see how the position
%  varies with extreme differences in spring constant. This will generate figures containing plots of adhesion position distribution in the more extreme case, the softest substrate and the stiffest substrate

% for i = 1:length(peak)
%     for j = 1
% 
% 
%         for r = 1:length(runs)
% 
%         ka_small_Adhesion_position = Positionresults.(['k_a_', SimFormat(k_a(j)), '_peak_', num2str(peak(i)) '_nligands_', num2str(nligands), 'run_0', num2str(r)]).Adhesion_position;
%         ka_large_Adhesion_position = Positionresults.(['k_a_', SimFormat(k_a(j+4)), '_peak_', num2str(peak(i)), '_nligands_', num2str(nligands), 'run_0', num2str(r)]).Adhesion_position;
% 
%         ka_small_Adhesion_position =   ka_small_Adhesion_position(ka_small_Adhesion_position ~= 0);
%         ka_large_Adhesion_position = ka_large_Adhesion_position(ka_large_Adhesion_position ~= 0);
% 
%         if r == 1
%             Small_position = ka_small_Adhesion_position;
%             Large_position = ka_large_Adhesion_position;
% 
%         else
% 
%             Small_position(end+1:end+length(ka_small_Adhesion_position)) = ka_small_Adhesion_position;
%             Large_position(end+1:end+length(ka_large_Adhesion_position)) = ka_large_Adhesion_position;
%         end
% 
%         clear ka_small_Adhesion_position
%         clear ka_large_Adhesion_position
% 
%         end
% 
% 
%         [counts1, binCenters1] = histcounts(Small_position, 20);
%         [counts3, binCenters3] = histcounts(Large_position, 20);
% 
%         figure;
%         plot(binCenters1(1:end-1), counts1, 'r-', 'LineWidth',3);
%         hold on;
%         plot(binCenters3(1:end-1), counts3, 'b-','LineWidth',3);
% 
% 
%         xlabel('Adhesion Distance from Membrane (nm)');
%         ylabel('Frequency');
%         set (gca,'linewidth', 5, 'fontsize', 15)
%         legend(['k_a = ' num2str(k_a(j))], ['k_a = ' num2str(k_a(j+4))]);
%         title(['Histogram of Adhesion Position for peak = ' num2str(peak(i))]);
%         print(['Integrinposition_', 'k_a_', SimFormat(k_a(j)), '_vs_' 'k_a_', SimFormat(k_a(j+4)), '_peak_', num2str(peak(i)) ], '-dpng', '-r400')
% 
%     end
% end


%This segment of code will concatenate the position from the 3 runs and then compute the average position of the adhesions for each spring constant, and then generate a figure plotting the average adhesion position for the different spring constants. It will also generate a plot comparing the STD of all adhesion positions for the 5 spring constants

for i = 1:length(k_a)
    for j = 1:length(peak)

        for r = 1:runs

        Adhesion_p = Positionresults.(['k_a_', SimFormat(k_a(i)), '_peak_', sprintf('%01d',peak(j)), '_nligands_', sprintf('%03d',nligands), 'run_0', num2str(r)]).Adhesion_position;

        Adhesion_p = Adhesion_p(Adhesion_p ~= 0);

        if r == 1
            Adhesion_position = Adhesion_p;

        else

            Adhesion_position(end+1:end+length(Adhesion_p)) = Adhesion_p;

        end

        clear Adhesion_p

        end

        Mean_position(i,j) = mean(Adhesion_position);
        Median_position(i,j) = median(Adhesion_position);
        STD(i,j) = std(Adhesion_position)/sqrt(length(Adhesion_position));

    end
end

 % f = figure(2);
 %        f.Color = 'white';
 %        h1 = semilogx(k_a, Mean_position(:,1), 'o-k', 'linewidth', 4);
 %        hold on
 %        h2 = semilogx(k_a, Mean_position(:,2), '-r', 'linewidth', 4);
 % 
 % 
 %        xlabel('Spring Constant (pN/nm)');
 %        ylabel('Average Position (nm)');
 %        legend([h1, h2], 'Wild Type', 'Manganese', 'Location','northwest')
 %        set (gca,'linewidth', 5, 'fontsize', 30)
 %        xticks([10^-5 10^-4 10^-3 10^-2 10^-1 10^0 10])
 %        axis([10^-5 10 0 600])
 %        set (gca,'linewidth', 5, 'fontsize', 30)
 %        xtickangle(60)
 %        title('Adhesion Position vs Spring Constant');
 %        print('MeanIntegrinposition_vs_spring_constant', '-dpng', '-r400')


   f = figure(3);
        f.Color = 'white';
        h1 = semilogx(k_a, Median_position(:,1), 'o-k', 'linewidth', 4);
        hold on
        h2 = semilogx(k_a, Median_position(:,2), '-r', 'linewidth', 4);


        xlabel('Spring Constant (pN/nm)');
        ylabel('Median Position (nm)');
        legend([h1, h2], 'Wild Type', 'Manganese', 'Location','northwest')
        set (gca,'linewidth', 5, 'fontsize', 30)
        xticks([10^-5 10^-4 10^-3 10^-2 10^-1 10^0 10])
        axis([10^-5 10 0 600])
        set (gca,'linewidth', 5, 'fontsize', 30)
        xtickangle(60)
        title('Adhesion Position vs Spring Constant');
        print('MedianIntegrinposition_vs_spring_constant', '-dpng', '-r400')


% r = figure(4);
%     r.Color = 'white';
%          f1 = semilogx(k_a, STD(:,1), ' -k', 'linewidth', 4);
%          hold on 
%          f2 = semilogx(k_a, STD(:,2), ' -r', 'linewidth', 4);
% 
%         xlabel('Spring Constant (pN/nm)');
%         ylabel('Standard error of position (nm)');
%         legend([f1, f2], 'Wild Type', 'Manganese', 'Location','northwest')
%         set (gca,'linewidth', 5, 'fontsize', 30)
%         xticks([10^-5 10^-4 10^-3 10^-2 10^-1 10^0 10])
%         set (gca,'linewidth', 5, 'fontsize', 30)
%         xtickangle(60)
%         title('Position Standard Error vs Spring Constant');
%         print('PositionSTD_vs_spring_constant', '-dpng', '-r400')

