clear
close all
clc


FileSaveDirectory = '/Users/remisondaz/Desktop/MATLAB/Varying_spring_constant/New_integrin_freq_analysis';

load(fullfile(FileSaveDirectory, 'AdhesionFrequencyresults.mat'));

Mean_Frequency = zeros(5,2);
Median_Frequency = zeros(5,2);
STD = zeros(5,2);

k_a = [0.0001 0.001 0.01 0.1 1];
peak = [1 2];
nligands = 400;
runs = 3;

for i = 1:length(k_a)
    for j = 1:length(peak)

        for r =1:runs

        Actin_freq = Frequencyresults.(['k_a_', SimFormat(k_a(i)), '_peak_', sprintf('%01d',peak(j)), '_nligands_', sprintf('%03d',nligands), 'run_0', num2str(r)]).Active_integrin_frequency(10001:29000);

        if r == 1

            Frequency = Actin_freq;

        else

            Frequency(end+1:end+length(Actin_freq)) = Actin_freq;

        end

        clear Actin_freq

        end


        Mean_Frequency(i,j) = mean(Frequency);
        STD(i,j) = std(Frequency);
    end
end

        
       f = figure(1);
        f.Color = 'white';
        h1 = semilogx(k_a, Mean_Frequency(:,1), 'o-k', 'linewidth', 4);
        hold on
        errorbar(k_a, Mean_Frequency(:,1),STD(:,1)./sqrt(3), ' -k', 'linewidth', 4)
        hold on
        h2 = semilogx(k_a, Mean_Frequency(:,2), '-r', 'linewidth', 4);
        hold on 
        errorbar(k_a, Mean_Frequency(:,2),STD(:,2)./sqrt(3), ' -r', 'linewidth', 4)


        xlabel('Spring Constant');
        ylabel('Average Frequency');
        legend([h1, h2], 'Wild Type', 'Manganese', 'Location','northwest')
        set (gca,'linewidth', 5, 'fontsize', 30)
        axis([10^-5 10 40 100])
        xticks([10^-5 10^-4 10^-3 10^-2 10^-1 10^0 10])
        set (gca,'linewidth', 5, 'fontsize', 30)
        xtickangle(60)
        title('Adhesion Frequency vs Spring Constant');
        print('MeanIntegrinfrequency_vs_spring_constant', '-dpng', '-r400')