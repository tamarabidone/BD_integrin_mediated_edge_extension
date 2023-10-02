% "prova" is a 4x16000 elements matrix
% 1st row has all values of vpol for WT on a soft substrate
% 2nd row --> Mn soft
% 3rd row --> WT stiff
% 4th row --> Mn stiff

clear
close all
clc


FileSaveDirectory = '/Users/remisondaz/Desktop/MATLAB/Varying_spring_constant/Y_speed_filaments';

load(fullfile(FileSaveDirectory, 'SpeedFilamentresults.mat'));


k_a = 0.1;
peak = [1 2];
nligands = 400;
runs = 3;

Total_Y_speed = NaN(2, 16000);

for j = 1:length(k_a)
    for i = 1:length(peak)
    

        Y_speed_per_run = NaN(30001, 100);
         
        for r = 1:runs

            Y_speed_per_run = SpeedFilamentresults.(['k_a_', SimFormat(k_a(j)), '_peak_', num2str(peak(i)), '_nligands_', num2str(nligands), 'run_0', num2str(runs)]).Y_speed_active_filaments;
            
        
            Y_speed_per_run = Y_speed_per_run(~isnan(Y_speed_per_run));

            if r==1

                Y_speed = Y_speed_per_run;
            else
                
                Y_speed(end+1:end+length(Y_speed_per_run))=Y_speed_per_run;
            end
            
            clear Y_speed_per_run

        end

        row_idx = i + length(k_a) * (j-1);
        
        % Store the combined Y_speed values in Total_Y_speed
        Total_Y_speed(row_idx, 1:length(Y_speed)) = Y_speed;
        valid_columns = ~all(isnan(Total_Y_speed), 1);
        

       % Remove NaN columns from Total_Y_speed
        Total_Y_speed = Total_Y_speed(:, valid_columns);
        

    end
end


bx=boxplot(Total_Y_speed','Notch','on', 'symbol', ' ');%,'MarkerStyle','none','Orientation', 'vertical');
% axis([0 5 0 65])

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
ylabel('\it{v_p_o_l} (nm/s)')
set (gca,'LineWidth',4, 'Fontsize', 45);
  set(gcf,'position',[200,300,500,700])
  title('k_a = 0.1 pN/nm')
  xlabel('Mn^2^+')

print('v_pol_only_large_stiffness', '-dpng', '-r400')