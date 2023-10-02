function [Adhesions,Ligands,FALconnections] = ManageAdhesionsAndLigands(Filaments,Adhesions,Ligands,FALconnections,Membrane,ModelParameters)
 


    %% ADHESION SECTION
    
            % Update adhesion regions 
               XY_LT = [ Membrane.Nodes.Xcoords(1), Membrane.Nodes.Ycoords(1) ]; %  Left-Top  xy-coordinates of current adhesion region
               XY_RT = [ Membrane.Nodes.Xcoords(2), Membrane.Nodes.Ycoords(2) ]; %  Right-Top xy-coordinates of current adhesion region
               XY_RB = [ Membrane.Nodes.Xcoords(2), Membrane.Nodes.Ycoords(2) - ModelParameters.ModelDepth]; %  Right-Bottom xy-coordinates of current adhesion region
               XY_LB = [ Membrane.Nodes.Xcoords(1), Membrane.Nodes.Ycoords(1) - ModelParameters.ModelDepth]; %  Left-Bottom  xy-coordinates of current adhesion region
               
            % Record Current Region
               Adhesions.RegionNodes(1,:) = XY_LT;           
               Adhesions.RegionNodes(2,:) = XY_RT;            
               Adhesions.RegionNodes(3,:) = XY_RB; 
               Adhesions.RegionNodes(4,:) = XY_LB; 
            
            % START Periodic boundary for Adhesions below lower adhesion region boundary ---------------------------------------------------------------------------------
                idx = find( Adhesions.XYpoints(:,2) < Adhesions.RegionNodes(4,2) ); 
                if ~isempty(idx)  
                        for idx2 = idx'
                                Adhesions.XYpoints(idx2,2) = Adhesions.XYpoints(idx2,2) + ModelParameters.ModelDepth; 
                                % Find out if adhesion is connected to a filament/ligand
                                [Filaments,Adhesions,Ligands,FALconnections] = DeleteFALconnection('AdhesionIndex',idx2,Filaments,Adhesions,Ligands,FALconnections);
                        end
                end
            
            % START Periodic boundary for Adhesions above upper boundary ---------------------------------------------------------------------------------
                idx = find( Adhesions.XYpoints(:,2) > Adhesions.RegionNodes(2,2) ); 
                if ~isempty(idx)  
                        for idx2 = idx'
                                Adhesions.XYpoints(idx2,2) = Adhesions.XYpoints(idx2,2) - ModelParameters.ModelDepth; 
                                % Find out if adhesion is connected to a filament/ligand
                                [Filaments,Adhesions,Ligands,FALconnections] = DeleteFALconnection('AdhesionIndex',idx2,Filaments,Adhesions,Ligands,FALconnections);
                        end
                end
          
            % START Periodic boundary for Adhesions outside left boundary ---------------------------------------------------------------------------------
                idx = find( Adhesions.XYpoints(:,1) < Adhesions.RegionNodes(1,1) ); 
                if ~isempty(idx)  
                        for idx2 = idx'
                                Adhesions.XYpoints(idx2,1) = Adhesions.XYpoints(idx2,1) + ModelParameters.MembraneWidth; 
                                % Find out if adhesion is connected to a filament/ligand
                                [Filaments,Adhesions,Ligands,FALconnections] = DeleteFALconnection('AdhesionIndex',idx2,Filaments,Adhesions,Ligands,FALconnections);
                        end
                end
            
            
            % START Periodic boundary for Adhesions outside right boundary ---------------------------------------------------------------------------------
                idx = find( Adhesions.XYpoints(:,1) > Adhesions.RegionNodes(2,1) ); 
                if ~isempty(idx)  
                        for idx2 = idx'
                                Adhesions.XYpoints(idx2,1) = Adhesions.XYpoints(idx2,1) - ModelParameters.MembraneWidth; 
                                % Find out if adhesion is connected to a filament/ligand
                                [Filaments,Adhesions,Ligands,FALconnections] = DeleteFALconnection('AdhesionIndex',idx2,Filaments,Adhesions,Ligands,FALconnections);
                        end
                end
            
           
            
    %% LIGAND SECTION
            
            % START Periodic boundary for Ligands below lower region boundary ---------------------------------------------------------------------------------
                idx = find( Ligands.XYpoints(:,2) < Adhesions.RegionNodes(4,2) ); 
                if ~isempty(idx)  
                        for idx2 = idx'
                                Ligands.XYpoints(idx2,2) = Ligands.XYpoints(idx2,2) + ModelParameters.ModelDepth; 
                                % Find out if ligand is connected to a filament/adhesion. If so, delete it
                                [Filaments,Adhesions,Ligands,FALconnections] = DeleteFALconnection('LigandIndex',idx2,Filaments,Adhesions,Ligands,FALconnections);
                        end
                end
            
           
            
            
            
           

end