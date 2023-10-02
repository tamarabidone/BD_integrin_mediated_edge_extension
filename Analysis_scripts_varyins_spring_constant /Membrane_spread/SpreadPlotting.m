clear
close all
clc


FileSaveDirectory = '/Users/remisondaz/Desktop/MATLAB/Varying_spring_constant/Membrane_spread';

load(fullfile(FileSaveDirectory, 'Membranespreadresults.mat'));


k_a = [0.0001 0.001 0.01 0.1 1];
peak = [1 2];
nligands = 400;
runs = 3;

start_time=19001;
end_time=29000;

for i = 1:length(peak)
    for j = 1:length(k_a)
        
        Membrane_spread_per_run=zeros(1,3);
        STD_spread_per_run=zeros(1,3);
        for r = 1:runs

           Membrane_spread_per_run(r)=Spreadresults.(['k_a_', SimFormat(k_a(j)), '_peak_', num2str(peak(i)), '_nligands_', num2str(nligands), 'run_0', num2str(r)]).Membrane_spread;
           STD_spread_per_run(r) = Spreadresults.(['k_a_', SimFormat(k_a(j)), '_peak_', num2str(peak(i)), '_nligands_', num2str(nligands), 'run_0', num2str(r)]).STD;

           
        end

        Membrane_spread(i,j) = mean(Membrane_spread_per_run);
        STD_spread(i,j) = mean(STD_spread_per_run(r));
    end
end

f = figure(2);
        f.Color = 'white';
        h1 = semilogx(k_a, Membrane_spread(1,:), 'o-k', 'linewidth', 4);
        hold on
        errorbar(k_a, Membrane_spread(1,:),STD_spread(1,:), ' -k', 'linewidth', 4)
        hold on
        h2 = semilogx(k_a, Membrane_spread(2,:), '-r', 'linewidth', 4);
        hold on 
        errorbar(k_a, Membrane_spread(2,:),STD_spread(2,:), ' -r', 'linewidth', 4)


        xlabel('Spring Constant (pN/nm)');
        ylabel('Membrane Spread (nm)');
        legend([h1, h2], 'Wild Type', 'Manganese', 'Location','northwest')
        set (gca,'linewidth', 5, 'fontsize', 30)
        xticks([10^-5 10^-4 10^-3 10^-2 10^-1 10^0 10^1])
        % axis([0 10 100 200])
        set (gca,'linewidth', 5, 'fontsize', 30)
        xtickangle(60)
        title('Membrane Spread vs Spring Constant');
        print('Membranespread_vs_spring_constant', '-dpng', '-r400')
