function [Integrins,Ligands,FILconnections] = ManageIntegrinsAndLigands(Filaments,Integrins,Ligands,FILconnections,Membrane,ModelParameters)
 


    %% INTEGRIN SECTION
    
            % Update adhesion regions 
               XY_LT = [ Membrane.Nodes.Xcoords(1), Membrane.Nodes.Ycoords(1) ]; %  Left-Top  xy-coordinates of current adhesion region
               XY_RT = [ Membrane.Nodes.Xcoords(2), Membrane.Nodes.Ycoords(2) ]; %  Right-Top xy-coordinates of current adhesion region
               XY_RB = [ Membrane.Nodes.Xcoords(2), Membrane.Nodes.Ycoords(2) - ModelParameters.ModelDepth]; %  Right-Bottom xy-coordinates of current adhesion region
               XY_LB = [ Membrane.Nodes.Xcoords(1), Membrane.Nodes.Ycoords(1) - ModelParameters.ModelDepth]; %  Left-Bottom  xy-coordinates of current adhesion region
               
            % Record Current Region
               Integrins.RegionNodes(1,:) = XY_LT;           
               Integrins.RegionNodes(2,:) = XY_RT;            
               Integrins.RegionNodes(3,:) = XY_RB; 
               Integrins.RegionNodes(4,:) = XY_LB; 
            
            % START Periodic boundary for integrins below lower adhesion region boundary ---------------------------------------------------------------------------------
                idx = find( Integrins.XYpoints(:,2) < Integrins.RegionNodes(4,2) ); 
                if ~isempty(idx)  
                        for idx2 = idx'
                                Integrins.XYpoints(idx2,2) = Integrins.XYpoints(idx2,2) + ModelParameters.ModelDepth; 
                                % Find out if integrin is connected to a filament/ligand
                                [Filaments,Integrins,Ligands,FALconnections] = DeleteFILconnection('IntegrinIndex',idx2,Filaments,Integrins,Ligands,FILconnections);
                        end
                end
            
            % START Periodic boundary for Integrins above upper boundary ---------------------------------------------------------------------------------
                idx = find( Integrins.XYpoints(:,2) > Integrins.RegionNodes(2,2) ); 
                if ~isempty(idx)  
                        for idx2 = idx'
                                Integrins.XYpoints(idx2,2) = Integrins.XYpoints(idx2,2) - ModelParameters.ModelDepth; 
                                % Find out if adhesion is connected to a filament/ligand
                                [Filaments,Integrins,Ligands,FILconnections] = DeleteFILconnection('IntegrinIndex',idx2,Filaments,Integrins,Ligands,FILconnections);
                        end
                end
          
            % START Periodic boundary for integrins outside left boundary ---------------------------------------------------------------------------------
                idx = find( Integrins.XYpoints(:,1) < Integrins.RegionNodes(1,1) ); 
                if ~isempty(idx)  
                        for idx2 = idx'
                                Integrins.XYpoints(idx2,1) = Integrins.XYpoints(idx2,1) + ModelParameters.ModelDepth; 
                                % Find out if integrin is connected to a filament/ligand
                                [Filaments,Integrins,Ligands,FILconnections] = DeleteFILconnection('IntegrinIndex',idx2,Filaments,Integrins,Ligands,FILconnections);
                        end
                end
            
            
            % START Periodic boundary for integrins outside right boundary ---------------------------------------------------------------------------------
                idx = find( Integrins.XYpoints(:,1) > Integrins.RegionNodes(2,1) ); 
                if ~isempty(idx)  
                        for idx2 = idx'
                                Integrins.XYpoints(idx2,1) = Integrins.XYpoints(idx2,1) - ModelParameters.MembraneWidth; 
                                % Find out if integrins is connected to a filament/ligand
                                [Filaments,Integrins,Ligands,FILconnections] = DeleteFILconnection('IntegrinIndex',idx2,Filaments,Integrins,Ligands,FILconnections);
                        end
                end
            
           
            
    %% LIGAND SECTION
            
            % START Periodic boundary for Ligands below lower region boundary ---------------------------------------------------------------------------------
                idx = find( Ligands.XYpoints(:,2) < Integrins.RegionNodes(4,2) ); 
                if ~isempty(idx)  
                        for idx2 = idx'
                                Ligands.XYpoints(idx2,2) = Ligands.XYpoints(idx2,2) + ModelParameters.ModelDepth; 
                                % Find out if ligand is connected to a filament/integrin. If so, delete it
                                [Filaments,Integrins,Ligands,FILconnections] = DeleteFILconnection('LigandIndex',idx2,Filaments,Integrins,Ligands,FILconnections);
                        end
                end
            
           
            
            
            
           

end
