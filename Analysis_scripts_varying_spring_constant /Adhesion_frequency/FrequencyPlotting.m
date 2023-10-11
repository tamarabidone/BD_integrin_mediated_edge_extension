clear
close all
clc

%This is the folder where we save the outputs from this analysis 
FileSaveDirectory = '/Users/remisondaz/Desktop/MATLAB/Varying_spring_constant/Adhesion_frequency';

%This folder contains raw model output files after running simulations code (for example: SimulationTaskList_001.m)
load(fullfile(FileSaveDirectory, 'AdhesionFrequencyresults.mat'));


k_a = [0.0001 0.001 0.01 0.1 1];
peak = [1 2];
nligands = 400;
runs = 3;

%Concatenating the three runs using the same substrate stiffness
for i = 1:length(k_a)
    for j = 1:length(peak)

        for r = 1:runs

        Actin_f_per_run = Frequencyresults.(['k_a_', SimFormat(k_a(i)), '_peak_', sprintf('%01d',peak(j)), '_nligands_', sprintf('%03d',nligands), 'run_0', num2str(r)]).Active_integrin_frequency(26001:29000);
        
           if r==1

      Actin_freq = Actin_f_per_run;
      

            else
                Actin_freq(end+1:end+length(Actin_f_per_run))=Actin_f_per_run;
               
            end

            clear Actin_f_per_run

        end

        Frequency(i,j) = mean(Actin_freq);
        STD(i,j) = std(Actin_freq);

    end
end

%This segment of code will generate a figure plotting the average frequency of active integrins over the 30 simulation time for the different substrate stiffnesses
        
       f = figure(1);
        f.Color = 'white';
        h1 = semilogx(k_a, Frequency(:,1), 'o-k', 'linewidth', 4);
        hold on
        errorbar(k_a, Frequency(:,1),STD(:,1)./sqrt(3), ' -k', 'linewidth', 4)
        hold on
        h2 = semilogx(k_a, Frequency(:,2), '-r', 'linewidth', 4);
        hold on 
        errorbar(k_a, Frequency(:,2),STD(:,2)./sqrt(3), ' -r', 'linewidth', 4)


        xlabel('Spring Constant (pN/nm)');
        ylabel('Average Frequency');
        legend([h1, h2], 'Wild Type', 'Manganese', 'Location','northwest')
        set (gca,'linewidth', 5, 'fontsize', 30)
        axis([10^-5 10 40 100])
        xticks([10^-5 10^-4 10^-3 10^-2 10^-1 10^0 10])
        xtickangle(60)
        title('Adhesion Frequency vs Spring Constant');
        print('MeanIntegrinfrequency_vs_k_a', '-dpng', '-r400')


        

      
