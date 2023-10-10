clear
close all
clc

%Defining where the plots generated as output will be saved
FileSaveDirectory = '/Users/remisondaz/Desktop/MATLAB/Varying_spring_constant/Actin_retrograde_flow';

load(fullfile(FileSaveDirectory, 'ActinFlowresults.mat'));


k_a = [0.0001 0.001 0.01 0.1 1];
peak = [1 2];
nligands = 400;
runs = 3;

start_time=19001;
end_time=29000;


for i = 1:length(peak)
    for j = 1:length(k_a)

        Actin_force_per_run = NaN(30001, 100);
        Actin_flow_per_run = NaN(30001, 100);
         
        for r = 1:runs

            Actin_f_per_run = Flowresults.(['k_a_', SimFormat(k_a(j)), '_peak_', num2str(peak(i)), '_nligands_', num2str(nligands), 'run_0', num2str(runs)]).Adhesion_force;
            
            Actin_fl_per_run = Flowresults.(['k_a_', SimFormat(k_a(j)), '_peak_', num2str(peak(i)), '_nligands_', num2str(nligands), 'run_0', num2str(runs)]).Retrograde_flow;

            Actin_f_per_run = Actin_f_per_run(start_time:end_time,:);
            Actin_fl_per_run  = Actin_fl_per_run(start_time:end_time,:);
 
            Actin_f_per_run = Actin_f_per_run(~isnan(Actin_f_per_run));
            Actin_fl_per_run = Actin_fl_per_run((Actin_fl_per_run < -1));

            if r==1
      Actin_force_per_run=Actin_f_per_run;
      Actin_flow_per_run=Actin_fl_per_run;
            else
                Actin_force_per_run(end+1:end+length(Actin_f_per_run))=Actin_f_per_run;
                Actin_flow_per_run(end+1:end+length(Actin_fl_per_run))=Actin_fl_per_run;
            end

            clear Actin_f_per_run
            clear Actin_fl_per_run
        end
       

            p = polyfit(Actin_force_per_run, Actin_flow_per_run, 1);
            px = [min(Actin_force_per_run) max(Actin_force_per_run)];
            py = polyval(p, px);
            plot(px, py, 'LineWidth', 2);
            title(['Adhesion force vs Retrograde Flow for k_a = ' num2str(k_a(j)), ' and peak = ' num2str(peak(i))], 'fontsize', 10);
            xlabel('Adhesion force (pN)')
            ylabel('Actin flow (nm/s)')
            % xtickangle(60)
            set (gca,'linewidth', 5, 'fontsize', 15)
            % axis([0 1200 0 200])
            % f.Position = [655 234 650 530];
            % axis square
            print(['Forcevsretrogradeflow_', 'k_a_', SimFormat(k_a(j)), '_peak_', num2str(peak(i))], '-dpng', '-r400')

    end
end


     

