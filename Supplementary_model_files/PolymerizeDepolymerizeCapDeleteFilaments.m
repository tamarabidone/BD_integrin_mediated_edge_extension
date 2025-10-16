function [Filaments,Adhesions,Ligands,FALconnections] = PolymerizeDepolymerizeCapDeleteFilaments(nMonomers,Filaments,Adhesions,Ligands,Membrane,FALconnections,ModelParameters)
    %
    % Filaments = PolymerizeDepolymerizeCapDeleteFilaments(Filaments,ModelParamters)
    %
    % KRC 12/07/2020
        
        nMonomers( isnan(nMonomers) ) = 0;
        nF = length(Filaments.XYCoords); % number of filaments
        D = ModelParameters.MonomerLength; %(nm)
        FilamentTips = GetFilamentsTipLocations(Filaments);
        
        % FOR BARBED END OF FILAMENT ------------------------------------------------------------------------------------------------------------------
        for f = 1:nF
                % Cap selected filaments first ----------------------------------------
                if ~Filaments.IsCapped(f,1) 
                    cap_test = rand(1) < ModelParameters.k_cap*ModelParameters.TimeStep;
                    if cap_test
                        Filaments.IsCapped(f,1) = true;
                    end
                end 
                
                % Un-Cap selected filaments -------------------------------------------
                if Filaments.IsCapped(f,1) 
                    cap_test = rand(1) < ModelParameters.k_uncap*ModelParameters.TimeStep;
                    if cap_test
                        Filaments.IsCapped(f,1) = false;
                    end
                end 
                
                % If filaments is not capped, check if it is selected for polymerizing
                if ~Filaments.IsCapped(f,1) 
                        if ModelParameters.BrownianRatchetOn
                            PolymCoeff = CalculatePolymerizationCoefficient(f,Filaments,Membrane,ModelParameters);
                        else
                            PolymCoeff = 1;
                        end
                        
                        y = Membrane.Nodes.Ycoords(1) - Filaments.XYCoords{f}(end,2);

                        % Calculate relative molarity based on an inverse relationship with current total filament mass
                        MolarityCoefficient = 1 - nMonomers / ModelParameters.MaximumFilamentMass;
                        MolarityCoefficient( MolarityCoefficient < 0 ) = 0;
                        r_on = ModelParameters.k_on_barbed * ModelParameters.FreeMonomerMolarity * MolarityCoefficient;
                        
                        % Run random selection test
                        barbed_on_test  =  rand  <  PolymCoeff * r_on * ModelParameters.TimeStep;

                        % Polymerize barbed end for selected filaments ----------------------------------------------------------------------
                        if barbed_on_test 
                            % Add next indices and next XY point on barbed end of filament
                            Filaments.MonomerIndices{f} = [Filaments.MonomerIndices{f}; Filaments.MonomerIndices{f}(end,1)+1]; 
                            Filaments.XYCoords{f} = [Filaments.XYCoords{f}; [Filaments.XYCoords{f}(end,1) + D*Filaments.UnitVector(f,1),...
                                                                             Filaments.XYCoords{f}(end,2) + D*Filaments.UnitVector(f,2)] ];
                        end
                        %--------------------------------------------------------------------------------------------------------------------
                end
        end 
      
        
        
        
        % FOR POINTED END OF FILAMENT------------------------------------------------------------------------------------------------------------------
        
        DeletedFilamentNames = []; % Keep track of all filaments that may have to be deleted
        
        for f = 1:nF
            % De-Polymerize pointed end for selected filaments ----------------------------------------
           
            if Filaments.XYCoords{f}(1,2) < Adhesions.RegionNodes(4,2)
                DF = 100;
            else
                DF = 1;
            end
                pointed_off_test = rand(1) < DF*ModelParameters.k_off_pointed*ModelParameters.TimeStep;
            
           
    
            L = length(Filaments.MonomerIndices{f}); % Get current length of filament
            
            if pointed_off_test && Filaments.Parent(f,1) == 0 % Parent = 0 means filament is a main filament
                        % Record Name and Current Monomer index in case this filament is deleted and we can still find its daughers (BranchIdx) afterwards
                        Name = Filaments.Name(f,1);
                        MonomerIndex = Filaments.MonomerIndices{f}(1); % First, grab index of monomer about to be deleted ( pointed end is position 1 )

                        % If parent filament length is 2 monomers long or less remove it 
                        if L <= 2  
                            dName = Filaments.Name(f,1);
                            DeletedFilamentNames = [DeletedFilamentNames; dName];
                            % Check if there are branches (and branches on branches...)
                            ind = find(  Filaments.MainIndex == dName );
                            if ~isempty(ind)
                                BranchName = Filaments.Name(ind,1);
                                DeletedFilamentNames = [DeletedFilamentNames; BranchName];
                            end
                        end

                        % START: Remove monomer ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                            Filaments.MonomerIndices{f} = Filaments.MonomerIndices{f}(2:end,1);
                            Filaments.XYCoords{f} = Filaments.XYCoords{f}(2:end,:);
                            % Check if adhesion/ligand is attached at this monomer on this filament. If so delete connection.
                            Value = [Filaments.Name(f,1),MonomerIndex];
                            [Filaments,Adhesions,Ligands,FALconnections] = DeleteFALconnection('FilamentNameAndMonomerIndex',Value,Filaments,Adhesions,Ligands,FALconnections);
                            
                            %----------------------------------------------------------------------------------------
                            % Check if there was a branch point there. If so make the branch a main filament (i.e. set Filaments.Parent = 0 );
                            BranchIdx = find( Filaments.Parent == Name );
                            idx1 = [];
                            for k = 1:length(BranchIdx) % If there is more than one branch, find the branch that is attached to the current monomer being depolymerized
                                if Filaments.ParentIndex(BranchIdx(k),1) == MonomerIndex
                                    idx1 = BranchIdx(k); break;
                                end
                            end

                            %----------------------------------------------------------------------------------------
                            % Now that the branch is set as a main filament, find all its daughter filaments and rename their MainIndex (Main Filament) to this branch.
                            if ~isempty(idx1)
                                   % Setting first branch as Main Filament
                                   Filaments.Parent(idx1,1) = 0; 
                                   Filaments.ParentIndex(idx1,1) = 0;
                                   Filaments.MainIndex(idx1,1) = Filaments.Name(idx1,1);
                                   %------------------------------------------------------------------------------
                                   % Find all sub-daughters of the new Main Filament and change their MainIndex to the new Filament name
                                   idx2 = find(Filaments.Parent == Filaments.Name(idx1,1));
                                   if ~isempty(idx2)
                                       idx3 = idx2;
                                       while true
                                            idx4 = find( ismember(Filaments.Parent,Filaments.Name(idx2,1)) );
                                            if isempty(idx4); break; end
                                            idx3 = [idx3; idx4];
                                            idx2 = idx4;
                                       end
                                       Filaments.MainIndex(idx3,1) = Filaments.Name(idx1,1);
                                   end
                                   %------------------------------------------------------------------------------
                            end
                        % END: Remove monomer ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
            end
        end
        
        
        %|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        % Remove deleted filaments before continuing on
        DeletedFilamentNames = unique(DeletedFilamentNames);
        
        for n = 1:length(DeletedFilamentNames)
            
                name = DeletedFilamentNames(n);
                idx = find( Filaments.Name == name );
                
                % Delete filament
                Filaments.Name(idx)           = [];
                Filaments.MonomerIndices(idx) = [];
                Filaments.XYCoords(idx)       = [];
                Filaments.UnitVector(idx,:)   = [];
                Filaments.IsCapped(idx)       = [];
                Filaments.MainIndex(idx)      = [];
                Filaments.Parent(idx)         = [];
                Filaments.ParentIndex(idx)    = [];
                
                % Check for Adhesion/Ligand connections and delete them
                [Filaments,Adhesions,Ligands,FALconnections] = DeleteFALconnection('FilamentName',name,Filaments,Adhesions,Ligands,FALconnections);
        end
        
        
        
end

%======================================================================================================
%======================================================================================================

function PolymCoeff = CalculatePolymerizationCoefficient(f,Filaments,Membrane,ModelParameters)
        
            y = Membrane.Nodes.Ycoords(1) - Filaments.XYCoords{f}(end,2); % Distance of filament tip from membrane
            FilTipX = Filaments.XYCoords{f}(end,1);
            if FilTipX > 0 && FilTipX < 500
                %|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                % Evaluate normalized polymerization velocity eq.3 (Mogliner and Oser 1996) for k_off = 0;
                theta = abs(atand( Filaments.UnitVector(f,1)/Filaments.UnitVector(f,2) ));
                theta(theta < 1 ) = 1; % theta can't be zero otherwise K has 0 in the denominator
                lambda = ModelParameters.PersistenceLength;
                kT = ModelParameters.kT;
                L = length(Filaments.MonomerIndices{f})*ModelParameters.MonomerLength;
                L( L<30  ) = 30;
                L( L>150 ) = 150;
                delta = ModelParameters.MonomerLength*cosd(theta);
                K = 4*lambda*kT/(L.^3*sind(theta).^2); % eq. B.1
                
                if y <= delta && y >= 0 % Three variations of the same calculation but split into three "practical calculations"
                    P_top = sqrt(pi*kT/(2*K))*( erfc(abs((y-delta)*sqrt(K/(2*kT)))) ); 
                    P_bot = sqrt(pi*kT/(2*K))*( erf (y*sqrt(K/(2*kT))) + 1 );
                elseif y < 0 
                    P_top = sqrt(pi*kT/(2*K))*( erfc(abs((y-delta)*sqrt(K/(2*kT)))) ); 
                    P_bot = sqrt(pi*kT/(2*K))*( erfc(abs(y*sqrt(K/(2*kT)))) );
                else % y > delta
                    P_top = sqrt(pi*kT/(2*K))*( erf((y-delta)*sqrt(K/(2*kT))) + 1 ); 
                    P_bot = sqrt(pi*kT/(2*K))*( erf(y*sqrt(K/(2*kT))) + 1 );
                end
                PolymCoeff = P_top./P_bot;
                PolymCoeff( isinf(PolymCoeff) | isnan(PolymCoeff) ) = 0;
            else
                PolymCoeff = 0;
            end

end

