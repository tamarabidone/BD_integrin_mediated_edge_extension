clear
close all
clc


FileSaveDirectory = '/Users/remisondaz/Desktop/MATLAB/Varying_spring_constant/Force_active_integrins';

load(fullfile(FileSaveDirectory, 'Activeintegrinsforce.mat'));


k_a = [0.0001 0.1];
peak = [1 2];
nligands = 400;
runs = 3;

Total_force = NaN(4, 16000);

for i = 1:length(peak)
    for j = 1:length(k_a)

        force_per_run = NaN(30001, 100);
         
        for r = 1:runs

            force_per_run = Activeintegrinsforce.(['k_a_', SimFormat(k_a(j)), '_peak_', num2str(peak(i)), '_nligands_', num2str(nligands), 'run_0', num2str(runs)]).Active_integrins_force;
            
        
           force_per_run = force_per_run(~isnan(force_per_run));
           force_per_run = force_per_run(force_per_run>0);

            if r==1

               Force = force_per_run;
            else
                
                Force(end+1:end+length(force_per_run))=force_per_run;
            end
            
            clear force_per_run

        end

        row_idx = i + length(k_a) * (j-1);
        
        % Store the combined Y_speed values in Total_Y_speed
        Total_force(row_idx, 1:length(Force)) = Force;
        valid_columns = ~all(isnan(Total_force), 1);

       % Remove NaN columns from Total_Y_speed
        Total_force = Total_force(:, valid_columns);

    end
end


bx=boxplot(Total_force','Notch','on', 'symbol', ' ');%,'MarkerStyle','none','Orientation', 'vertical');
axis([0 5 0 5])

% bc.BoxFaceColor=[0.0 0.0 0.8];
% bc.BoxWidth=0.7;
% bc.BoxFaceAlpha=0.1;
colors = [0 0 0; 0.9 0.9 0.9; 0 0 0; 0.9 0.7 0.9];
h = findobj(gca,'Tag','Box');

for j=1:length(h)
     patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.1, 'linewidth',4,'EdgeAlpha',0.3);
end
    set(bx(:,:),'LineWidth',3);
xticklabels({'-', '+', '-', '+'});
ylabel('\it{F_s_u_b} (pN)')
set (gca,'LineWidth',4, 'Fontsize', 45);
  set(gcf,'position',[200,300,500,700])

print('F_sub_3A_k_minus4_minus1', '-dpng', '-r400')