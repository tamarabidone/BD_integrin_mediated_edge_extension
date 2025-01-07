function [FILconnections,count] = PlotFilamentsAndMembraneMovie01(nth,count,Filaments,Membrane,Integrins,Ligands,FILconnections,FH,AH1,AH2,AH3,t,nMono,nInteg,MemVel,index,TimeVec,ModelParameters, Data)

    % if (count >= nth) || isequal(t,0)% Plot every 10th, 100th, or nth frame etc.
    %     count = 0;
        figure(FH); clf
        FH.Position = [0.2592    0.2416    0.5265    0.6119];
        AH1 = [];
        AH2 = [];
        AH3 = [];
        TL = tiledlayout(FH,1,1,"TileSpacing","compact","Padding","compact");
        AH1 = nexttile(TL,1);
        idxMF = find(Filaments.Parent == 0); % Find Main Filaments (Filaments that are at the top of the parent/daughter structure)
        nC = length(idxMF);
        C = lines(nC); % Create colors for each of the attached groups of filaments

        %Plot branching region
        % hold(AH1, 'on');
        % Xcoords = Membrane.Nodes.Xcoords;
        % Ycoords_top = Membrane.Nodes.Ycoords;
        % Ycoords_bottom = Membrane.Nodes.Ycoords - 15;
        % X = [Xcoords, fliplr(Xcoords)];  % The X-coordinates are mirrored for top and bottom
        % Y = [Ycoords_top, fliplr(Ycoords_bottom)];  % Top edge and bottom edge
        % fill(X, Y, [0.5059, 0.5059, 0.5059], 'EdgeColor', 'none');

        % For all the attached groups of filaments, plot each group with its own color
        for MF = 1:nC 
            if isequal(MF,1)
                hold(AH1,'on')
            end
            % Find all the filaments that are attached to the main parent
            % filament, and the main parent filament as well.
            idx = find(Filaments.MainIndex == Filaments.Name(idxMF(MF)));
            for f = idx'
                if isequal(Filaments.Parent(f,1),0) %Data.YSpeed(f,1) < 0 
                    % Plotting for Main Filaments
%                     plot(AH1, Filaments.XYCoords{f}(:,1),...
%                               Filaments.XYCoords{f}(:,2),'-','Color',C(MF,:),'LineWidth',2)
                    plot(AH1, Filaments.XYCoords{f}(:,1),...
                              Filaments.XYCoords{f}(:,2),'-','Color',[0, 0.7, 0],'LineWidth',2)
                else
                    %Plotting for daughter filaments
                    idxP  = find(Filaments.Name == Filaments.Parent(f,1));
                    idxXY = find(Filaments.MonomerIndices{idxP} == Filaments.ParentIndex(f,1));
                    XYstart = Filaments.XYCoords{idxP}(idxXY,:);
                    plot(AH1, [XYstart(1,1); Filaments.XYCoords{f}(:,1)],...
                              [XYstart(1,2); Filaments.XYCoords{f}(:,2)],'-','Color',C(MF,:),'LineWidth',2)
                    plot(AH1, [XYstart(1,1); Filaments.XYCoords{f}(:,1)],...
                              [XYstart(1,2); Filaments.XYCoords{f}(:,2)],'-','Color',[0.4902, 0.1529, 0.4902],'LineWidth',2)
                  % plot(AH1, Filaments.XYCoords{f}(:,1),...
                  %             Filaments.XYCoords{f}(:,2),'-','Color',[0.68, 0.85, 0.9],'LineWidth',2)
                end    
            end
        end
        

        % Plot Red dot for capped filaments
        if any(Filaments.IsCapped)
            idxC = find(Filaments.IsCapped);
            for f = 1:length(idxC)
                Xc = Filaments.XYCoords{idxC(f)}(end,1);
                Yc = Filaments.XYCoords{idxC(f)}(end,2);
                plot(AH1,Xc,Yc,'.b','MarkerSize',10)
            end
        end 
        


        % Calculate integrin regions area for density calculation --------------------------
        delX = diff(Membrane.Nodes.Xcoords);
        delY = ModelParameters.ModelDepth;
        AreaOfAllRegions = delX*delY/1000^2; % um^2
        set(AH1,'FontSize',14)
%         Plot inactive integrin indices and plot
         % idx = find(~Integrins.ActiveStatus);
         % indices = randperm(numel(idx), 10);
         % plot(AH1, Integrins.XYpoints(idx,1), Integrins.XYpoints(idx,2), '*k', 'MarkerSize', 3)
        
%        plot active integrins as red
        idx = FILconnections.IntegrinIndex;
        % plot(AH1, Integrins.XYpoints(idx,1), Integrins.XYpoints(idx,2),  '*', 'MarkerSize', 6,'Color',[1,0,0])
        % 
        % if ModelParameters.MolecularClutch_PeakNumber == 1
        %     MCtext = 'WT';
        % else
        %     MCtext = 'Mg^{2+}';
        % end
        % 
        % nBranches = length(Filaments.Parent(Filaments.Parent > 0));
        % 

        idx = max([find(~isnan(MemVel),1,'last'),1]);
        % title(AH1,[{['Time = ', sprintf('%#0.3f',t), ' s      ','FIL connections = ', pad( num2str( numel(FILconnections.IntegrinIndex) ),4,'left'), '     dt = ', num2str(ModelParameters.TimeStep*1000),' ms']};...
        %            {['Filament mass = ', pad( sprintf('%0.0f', nMono (idx,1) ), 6,'left'), ' monomers     ',...
        %             'Membrane speed = ',      pad( sprintf('%0.1f', MemVel(idx,1) ), 7,'left'), ' nm/s']};...
        %            {['Mol. Clutch = ',MCtext,'     k = ',sprintf('%1.4f',ModelParameters.IntegrinSpringConstant), ' pN/nm      nBranches = ',num2str(nBranches)]} ],...
        %              'FontName','monospaced','FontSize', 19, 'FontWeight','bold','HorizontalAlignment','center')
               
       MemYmax = ceil(max(Membrane.Nodes.Ycoords(1)));
       axis(AH1,'equal') 
       %axis(AH1,[-200, 700, (MemYmax-1100), (MemYmax+100)])
       axis(AH1,[-500, 1000, -250, 1000])
       
       box(AH1,'on')
       grid(AH1,'on')
       xlabel(AH1,'X(nm)','FontSize',22,'FontWeight','bold')
       ylabel(AH1,'Y(nm)','FontSize',22,'FontWeight','bold')
       

       % Plot intergin connections --------------------------
       %idx = find(~isnan(FILconnections.IntegrinIndex));
       % for n = 1:length(FILconnections.IntegrinIndex)  %n = 1:nC
       %      Axy = Integrins.XYpoints(FILconnections.IntegrinIndex(n,1),:); 
       %      f = find(Filaments.Name == FILconnections.FilamentName(n,1)); 
       %      Midx = find(Filaments.MonomerIndices{f} == FILconnections.MonomerIndex(n,1));
       %      Mxy  = Filaments.XYCoords{f}(Midx,:);
       %      plot(AH1, [Axy(1);Mxy(1)], [Axy(2);Mxy(2)],':r','LineWidth',2)
       % end
       
       
       % Plot Membrane -------------------------------------------------
             hold(AH1,'on'); 
             plot( AH1, Membrane.Nodes.Xcoords,...
                        Membrane.Nodes.Ycoords, 'k-','Color',[0.55,0.1,0.1],'LineWidth',6)
            

             text(AH1, 0.02, 0.98, ['Time = ', sprintf('%#0.3f', t), ' s'], 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 12, 'Color', 'k', 'FontWeight', 'bold');
      
      % Plot Ligands ---------------------------------------------------
%              plot( AH1, Ligands.XYpoints(:,1),...
%                         Ligands.XYpoints(:,2), 'bs','MarkerSize',8)
      
       set(AH1,'LineWidth',2)
       hold(AH1,'off')
       
     
           
    %end
end



