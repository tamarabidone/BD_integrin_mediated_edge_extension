clc
clear
close all

FileSaveDirectory = '/Users/remisondaz/Desktop/MATLAB/Varying_spring_constant/Velocity';

t_max = 29;
load(fullfile(FileSaveDirectory, 'Velocityresults.mat'));

k_a =  [0.0001 0.001 0.01 0.1 1]; 
peak = [1 2];
nligands = 400;

max_peak = length(peak);

for i = 1:length(k_a)
    for j = 1:length(peak)
        Velocity_data(i,j) = Results.(['k_a_', SimFormat(k_a(i)), '_peak_', sprintf('%01d',peak(j)), '_nligands_', sprintf('%03d',nligands)]).Velocity;
        STDVelocity_data(i,j) = Results.(['k_a_', SimFormat(k_a(i)), '_peak_', sprintf('%01d',peak(j)), '_nligands_', sprintf('%03d',nligands)]).STD_Velocity;
        Period_data(i,j) = Results.(['k_a_', SimFormat(k_a(i)), '_peak_', sprintf('%01d',peak(j)), '_nligands_', sprintf('%03d',nligands)]).Periods;
        STDPeriod_data(i,j) = Results.(['k_a_', SimFormat(k_a(i)), '_peak_', sprintf('%01d',peak(j)), '_nligands_', sprintf('%03d',nligands)]).STD_Periods;
        Prominence_data(i,j) = Results.(['k_a_', SimFormat(k_a(i)), '_peak_', sprintf('%01d',peak(j)), '_nligands_', sprintf('%03d',nligands)]).Prominence;
        STDProminence_data(i,j) = Results.(['k_a_', SimFormat(k_a(i)), '_peak_', sprintf('%01d',peak(j)), '_nligands_', sprintf('%03d',nligands)]).STD_Prominence;
        Width_data(i,j) = Results.(['k_a_', SimFormat(k_a(i)), '_peak_', sprintf('%01d',peak(j)), '_nligands_', sprintf('%03d',nligands)]).Width;
        STDWidth_data(i,j) = Results.(['k_a_', SimFormat(k_a(i)), '_peak_', sprintf('%01d',peak(j)), '_nligands_', sprintf('%03d',nligands)]).STD_Width;
        Frequency_data(i,j) = Results.(['k_a_', SimFormat(k_a(i)), '_peak_', sprintf('%01d',peak(j)), '_nligands_', sprintf('%03d',nligands)]).Frequency;
        STDFrequency_data(i, j) = Results.(['k_a_', SimFormat(k_a(i)), '_peak_', sprintf('%01d',peak(j)), '_nligands_', sprintf('%03d',nligands)]).STD_Frequency;
    end

end

f = figure(1);
f.Color = 'white';
h1 = semilogx(k_a, Velocity_data(:,1), 'o-k', 'linewidth', 4);
hold on
errorbar(k_a, Velocity_data(:,1),STDVelocity_data(:,1)./sqrt(3), ' -k', 'linewidth', 4)
hold on
h2 = semilogx(k_a, Velocity_data(:,2), '-r', 'linewidth', 4);
hold on 
errorbar(nligands, Velocity_data(:,2),STDVelocity_data(:,2)./sqrt(3), ' -r', 'linewidth', 4)
xlabel('Spring Constant')
xticks([10^-5 10^-4 10^-3 10^-2 10^-1 10^0 10])
ylabel('V_e_d_g_e (nm/s)')
xtickangle(60)
set (gca,'linewidth', 5, 'fontsize', 30)
axis([10^-5 10 10 19])
f.Position = [655 234 650 530];
legend([h1, h2], 'Wild Type', 'Manganese', 'Location','northwest')
axis square
print('velocity_vs_spring_constant', '-dpng', '-r400')

f = figure(2);
f.Color = 'white';
h3 = semilogx(k_a, Frequency_data(:,1), 'o-k', 'linewidth', 4);
hold on
errorbar(k_a, Frequency_data(:,1),STDFrequency_data(:,1)./sqrt(3), ' -k', 'linewidth', 4)
hold on
h4 = semilogx(k_a, Frequency_data(:,2), '-r', 'linewidth', 4);
hold on 
errorbar(k_a, Frequency_data(:,2),STDFrequency_data(:,2)./sqrt(3), ' -r', 'linewidth', 4)
xlabel('Spring Constant')
xticks([10^-5 10^-4 10^-3 10^-2 10^-1 10^0 10])
ylabel('Frequency (s^-1)')
xtickangle(60)
set (gca,'linewidth', 5, 'fontsize', 30)
axis([10^-5 10 0.4 1.5])
f.Position = [655 234 650 530];
legend([h3, h4], 'Wild Type', 'Manganese', 'Location','northwest')
axis square
print('avg_frequency_vs_spring_constant', '-dpng', '-r400')


f = figure(3);
f.Color = 'white';
h5 = semilogx(k_a, Width_data(:,1), 'o-k', 'linewidth', 4);
hold on
errorbar(k_a, Width_data(:,1),STDWidth_data(:,1)./sqrt(3), ' -k', 'linewidth', 4)
hold on
h6 = semilogx(k_a, Width_data(:,2), '-r', 'linewidth', 4);
hold on 
errorbar(k_a, Width_data(:,2),STDWidth_data(:,2)./sqrt(3), ' -r', 'linewidth', 4)
xlabel('Spring Constant')
xticks([10^-5 10^-4 10^-3 10^-2 10^-1 10^0 10])
ylabel('Peak duration (s)')
xtickangle(60)
set (gca,'linewidth', 5, 'fontsize', 30)
axis([10^-5 10 0.2 1.0])
f.Position = [655 234 650 530];
legend([h5, h6], 'Wild Type', 'Manganese', 'Location','northwest')
axis square
print('width_vs_spring_constant_', '-dpng', '-r400')

f = figure(4);
f.Color = 'white';
h7 = semilogx(k_a, Prominence_data(:,1), 'o-k', 'linewidth', 4);
hold on
errorbar(k_a, Prominence_data(:,1),STDProminence_data(:,1)./sqrt(3), ' -k', 'linewidth', 4)
hold on
h8 = semilogx(k_a, Prominence_data(:,2), '-r', 'linewidth', 4);
hold on 
errorbar(k_a, Prominence_data(:,2),STDProminence_data(:,2)./sqrt(3), ' -r', 'linewidth', 4)
xlabel('Spring Constant')
xticks([10^-5 10^-4 10^-3 10^-2 10^-1 10^0 10])
ylabel('Peak Prominence (nm/s)')
xtickangle(60)
set (gca,'linewidth', 5, 'fontsize', 30)
axis([10^-5 10 5 28])
f.Position = [655 234 650 530];
legend([h7, h8], 'Wild Type', 'Manganese', 'Location','northwest')
axis square
print('peak_prominence_vs_spring_constant', '-dpng', '-r400')
