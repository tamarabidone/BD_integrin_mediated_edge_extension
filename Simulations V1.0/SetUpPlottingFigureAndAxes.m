function [FH,AH1,AH2,AH3] = SetUpPlottingFigureAndAxes

    FH = figure(20); 
    clf(FH)
    set(FH,'Units','normalized','Position',[ 0.1626    0.0641    0.8271    0.7648],'Color',[1,1,1])
    % AH1 = axes('Parent',FH,'Position',[0.07,0.30,0.88,0.7]);
    AH1 = axes('Parent',FH,'Position',[0.52      0.0907    0.48      0.8343]);
    AH2 = axes('Parent',FH,'Position',[0.0689    0.5532    0.400     0.3718]);
    AH3 = axes('Parent',FH,'Position',[0.0689    0.0907    0.400     0.3718]);
    %                      Position = [left,     bottom,    width,    height]
end
