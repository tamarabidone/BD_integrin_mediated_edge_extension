function [FILconnections, Integrins, Ligands] = CreateFILconnections(FILconnections,Filaments,Integrins,Ligands,ModelParameters)

    nF = size(Filaments.XYCoords,1);
    D = ModelParameters.FAL_connection_Distance; % nm
    MinDist = ModelParameters.MonomerLength; % Minimum distance for connection. 
    % Test if Integrins are within a distance D of a filament/Monomer
    if ~isempty(Filaments.XYCoords) % Make sure there are filaments
            %-----------------------------------------------------------------------------------------------------------------------------
            for f = 1:nF
                        % Grab XY points of all monomers making up this filament
                        XY_filament = Filaments.XYCoords{f};
                        idx1 = find( FILconnections.FilamentName  == Filaments.Name(f,1) );
                        if ~isempty(idx1)
                            idx2 = find( ismember( Filaments.MonomerIndices{f},FALconnections.MonomerIndex(idx1)) );
                            XY_filament(idx2,:) = NaN;  % Set monomer positions that are already connected to a ligand/integrin to nan
                        end
                        
                        % Grab XY connections of all unattached (inactive) integrins
                        XY_adhesions = Integrins.XYpoints;
                        XY_adhesions(FILconnections.IntegrinIndex,:) = NaN; % Set adhesions already connected to nan
                        
                        % Grab XY connections of all unattached integrins
                        XY_ligands = Ligands.XYpoints;
                        XY_ligands(FILconnections.LigandIndex,:) = NaN;
                        
                        % Calculate distance from Integrins to filament monomers
                        Distance_IF = sqrt( (XY_integrins(:,1) - XY_filament(:,1)').^2 + (XY_integrins(:,2) - XY_filament(:,2)').^2 );
                        % Calculate distance from Integrins to Ligands
                        Distance_IL = sqrt( (XY_integrins(:,1) - XY_ligands(:,1)' ).^2 + (XY_integrins(:,2) - XY_ligands(:,2)' ).^2 );
                        
                        % Find adesion-filament and adhesion-ligand distances less than D
                        MinInEachRow_IF = false(size(Distance_IF));
                        [~,MinIdx] = min(Distance_IF,[],2,'linear');
                        MinInEachRow_IF(MinIdx) = true;
                        
                        MinInEachRow_IL = false(size(Distance_IL));
                        [~,MinIdx] = min(Distance_IL,[],2,'linear');
                        MinInEachRow_IL(MinIdx) = true;
                        
                        [Iidx1,Midx] = find( Distance_IF <= D & MinInEachRow_IF );
                        [Iidx2,Lidx] = find( Distance_IL <= D & MinInEachRow_IL ) ;
                        
                        % Find integrins that are close to a monomer and ligand
                        [~,idx3,idx4] = intersect(Iidx1,Iidx2);
                        
                        % Create integrin-filament-ligand connection
                        if ~isempty(idx3)
                            for c = 1:length(idx3)
                                p_on = 1 - exp(-ModelParameters.Integrin_ActivationRate * ModelParameters.TimeStep); % Calculate connection probability
                                if rand < p_on % test probability of this interaction happening
                                    a = Iidx1(idx3(c));
                                    l = Lidx (idx4(c)); 
                                    m = Filaments.MonomerIndices{f}( Midx(idx3(c)) );
                                    n = Filaments.Name(f);

                                    FILconnections.FilamentName  = [ FILconnections.FilamentName;  n ];
                                    FILconnections.MonomerIndex  = [ FILconnections.MonomerIndex;  m ];
                                    FILconnections.IntegrinIndex = [ FILconnections.IntegrinIndex; a ];
                                    FILconnections.LigandIndex   = [ FILconnections.LigandIndex;   l ];

                                    Integrins.AttachedFilamentName(a,1) = n;
                                    Integrins.ActiveStatus(a,1) = true;
                                    Integrins.XYpoints(a,:) = Ligands.XYpoints(l,:);

                                    Ligands.AttachedFilamentName(l)  = n;
                                    Ligands.AttachedIntegrinIndex(l) = a;
                                end
                            end
                        end
            end
            %-----------------------------------------------------------------------------------------------------------------------------
    end
    
    
end
