function [FILconnections,count] = PlotFilamentsAndMembrane(nth,count,Filaments,Membrane,Integrins,Ligands,FILconnections,FH,AH1,AH2,AH3,t,nMono,nAdhes,MemVel,index,TimeVec,ModelParameters)

    if (count >= nth) || isequal(t,0)% Plot every 10th, 100th, or nth frame etc.
        count = 0;
        figure(FH)
        cla(AH1)
        cla(AH2)
        
        idxMF = find(Filaments.Parent == 0); % Find Main Filaments (Filaments that are at the top of the parent/daughter structure)
        nC = length(idxMF);
        C = lines(nC); % Create colors for each of the attached groups of filaments

        % For all the attached groups of filaments, plot each group with its own color
        for MF = 1:nC 
            if isequal(MF,1)
                hold(AH1,'on')
            end
            % Find all the filaments that are attached to the main parent
            % filament, and the main parent filament as well.
            idx = find(Filaments.MainIndex == Filaments.Name(idxMF(MF)));
            for f = idx'
                if isequal(Filaments.Parent(f,1),0)
                    % Plotting for Main Filaments
                    plot(AH1, Filaments.XYCoords{f}(:,1),...
                              Filaments.XYCoords{f}(:,2),'-','Color',C(MF,:),'LineWidth',1)
                else
                    % Plotting for daughter filaments
                    idxP  = find(Filaments.Name == Filaments.Parent(f,1));
                    idxXY = find(Filaments.MonomerIndices{idxP} == Filaments.ParentIndex(f,1));
                    XYstart = Filaments.XYCoords{idxP}(idxXY,:);
                    plot(AH1, [XYstart(1,1); Filaments.XYCoords{f}(:,1)],...
                              [XYstart(1,2); Filaments.XYCoords{f}(:,2)],'-','Color',C(MF,:),'LineWidth',1)
                end    
            end
        end
        

        % Plot Red dot for capped filaments
        if any(Filaments.IsCapped)
            idxC = find(Filaments.IsCapped);
            for f = 1:length(idxC)
                Xc = Filaments.XYCoords{idxC(f)}(end,1);
                Yc = Filaments.XYCoords{idxC(f)}(end,2);
                plot(AH1,Xc,Yc,'.k','MarkerSize',6)
            end
        end 
        


        % Calculate adhesion regions area for density calculation --------------------------
        delX = diff(Membrane.Nodes.Xcoords);
        delY = ModelParameters.ModelDepth;
        AreaOfAllRegions = delX*delY/1000^2; % um^2

        % Plot inactive adhesion indices and plot
        idx = find(~Integrins.ActiveStatus);
        plot(AH1, Integrins.XYpoints(idx,1), Integrins.XYpoints(idx,2), '*k', 'MarkerSize', 6)
       
        % plot active adhesions as red
        idx = FILconnections.AdhesionIndex;
        plot(AH1, Integrins.XYpoints(idx,1), Integrins.XYpoints(idx,2),  '*r', 'MarkerSize', 6)
       
        idx = max([find(~isnan(MemVel),1,'last'),1]);
        title(AH1,[{['Time = ', sprintf('%#0.3f',t), ' s      ',...
                     'FIL connections = ', pad( num2str( numel(FILconnections.AdhesionIndex) ),4,'left'), '     dt = ', num2str(ModelParameters.TimeStep*1000),' ms']};...
                   {['Filament mass = ', pad( sprintf('%0.0f', nMono (idx,1) ), 6,'left'), ' monomers     ',...
                    'Membrane speed = ',      pad( sprintf('%0.1f', MemVel(idx,1) ), 7,'left'), ' nm/s']}],...
                     'FontName','monospaced','FontSize', 14, 'FontWeight','bold','HorizontalAlignment','center')
               
       MemYmax = ceil(max(Membrane.Nodes.Ycoords(1)));
       axis(AH1,'equal') 
       %axis(AH1,[-200, 700, (MemYmax-1100), (MemYmax+100)])
       axis(AH1,[-100, 600, (MemYmax-700), (MemYmax+100)])
       
       box(AH1,'on')
       grid(AH1,'on')
       set(AH1,'FontSize',14)
       xlabel(AH1,'X(nm)','FontSize',22,'FontWeight','bold')
       ylabel(AH1,'Y(nm)','FontSize',22,'FontWeight','bold')
       

       % Plot integrin connections --------------------------
       %idx = find(~isnan(FILconnections.AdhesionIndex));
       for n = 1:length(FILconnections.AdhesionIndex)  %n = 1:nC
            Axy = Integrins.XYpoints(FILconnections.IntegrinIndex(n,1),:); 
            f = find(Filaments.Name == FILconnections.FilamentName(n,1)); 
            Midx = find(Filaments.MonomerIndices{f} == FILconnections.MonomerIndex(n,1));
            Mxy  = Filaments.XYCoords{f}(Midx,:);
            plot(AH1, [Axy(1);Mxy(1)], [Axy(2);Mxy(2)],':k')
       end
       
       
       % Plot Membrane -------------------------------------------------
             hold(AH1,'on'); 
             plot( AH1, Membrane.Nodes.Xcoords,...
                        Membrane.Nodes.Ycoords, 'r-','Color',[0.6,0,0],'LineWidth',6)
      
      % Plot Ligands ---------------------------------------------------
             plot( AH1, Ligands.XYpoints(:,1),...
                        Ligands.XYpoints(:,2), 'bs','MarkerSize',8)
      
       set(AH1,'LineWidth',2)
       hold(AH1,'off')
       
       % Plot Membrane Speed
           idx = find(~isnan(MemVel));
           VelSmooth = NaN(size(MemVel));
           VelSmooth(idx,1) = smoothdata(MemVel(idx),'gaussian',100);
           plot(AH2,TimeVec,MemVel,'b-','LineWidth',0.5) %'Color',[0.8,0.8,0.8]); 
           hold(AH2,'on')
           plot(AH2,TimeVec,VelSmooth,'r-','LineWidth',2);     
           hold(AH2,'off')
           set(AH2,'LineWidth',1,'FontSize',18)
           
           xlim(AH2,[0,max(TimeVec)])
           ylim(AH2,[0,prctile(MemVel(idx),99)])
           xlabel(AH2,'Time (s)','FontSize',18)
           ylabel(AH2,'Membrane Speed (nm/s)','FontSize',22)
           
       
       % Plot Total Monomers
           plot(AH3, TimeVec, nMono,'-','LineWidth',1)
           set(AH3,'LineWidth',1,'FontSize',18)
           xlim(AH3,[0,max(TimeVec)])
           %ylim(AH3,[min([max(nMono),4800]),max(nMono)+100])
           ylim(AH3,[min(nMono), max(nMono)+1])
           xlabel(AH3,'Time (s)','FontSize',18)
           ylabel(AH3,'Filament Mass (monomers)','FontSize',22)
           
    end
end



