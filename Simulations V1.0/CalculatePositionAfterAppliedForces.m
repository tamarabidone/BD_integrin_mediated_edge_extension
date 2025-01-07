function [Filaments, Membrane, Integrins, Ligands, FILconnections, IntegrinTensions, Data] = ...
                                  CalculatePositionAfterAppliedForces(Filaments,Membrane,Integrins,Ligands,FILconnections,ModelParameters)
                              

    % Find all the main filaments of connected filament structures (group of attached filaments)
    idxMF = find(Filaments.Parent == 0);
    nMF = length(idxMF); % Total number of filament structures
    dt = ModelParameters.TimeStep;
    kT = ModelParameters.kT;
    FilamentTips = GetFilamentsTipLocations(Filaments);
    IntegrinTensions = zeros(size(Integrins.ActiveStatus));
    
    Data.YSpeed    = [];
    Data.FilamentName = [];
    Data.XYPosition = [];
    Data.StructName = [];
    
    
    %% FILAMENT MOTION SECTION 
    if ~isempty(Filaments.Name)
                for MF = 1:nMF % Go through each filament structure and test each of its filaments to see if it is hitting the membrane and calculate the total force acting on the structure.
                        L = 0; % Total length of Filament
                        idx1 = find( Filaments.MainIndex == Filaments.Name(idxMF(MF)) ); % Find all filaments attached to same structure
                        nF = length(idx1); % number of filaments in this structure
                    % Calculate total length of filament structure
                        for f = idx1'  
                            L = L + length(Filaments.MonomerIndices{f}) * ModelParameters.MonomerLength; % L is total length of filament structure in microns
                        end 
                    % CALCULATE FORCES DUE TO FILAMENT-ADHESION CONNECTIONS and break connections where force is greater than threshold
                       [FAx,FAy,Integrins,Ligands,FILconnections,IntegrinTensions] = CalculateForceDueToFILconnections(idx1,Filaments,Integrins,Ligands,IntegrinTensions,FILconnections,ModelParameters);
                        
                    % CALCULATE FORCES ACTING ON FILAMENT FROM MEMBRANE
                        % Calculate sum of the normal forces from filements pushing on membrane (all from the same structure)
                        Fm = TotalForceExertedByMembraneOnFilamentStructure(idx1,Filaments,Membrane,ModelParameters);

                    % CALCULATE NEW POSITION FOR FILAMENT STRUCTURE
                        Nu = ModelParameters.CytoplasmViscosity*(10^-6); % pN s/nm^2
                        r  = 7/2; % nm
                        L( L < r ) = r;
                        gamma1 =  (4*pi*Nu*L) / (0.84+log(L/(2*r) )); % gamma for random fluid motion (pN*s/nm) 
                        SD = sqrt(2*kT*gamma1/dt);                    % sqrt(2*(pN*nm)*(pN*S/um)/S)
                        
                        % Add filament thermal motion if it's turned on ------
                        if ModelParameters.FilamentThermalMotionOn
                            Fx = SD*randn;
                            Fy = SD*randn;
                        else
                            Fx = 0;
                            Fy = 0;
                        end
                        %-----------------------------------------------------
                        
                        tipXpositions = NaN(nF,1);
                        for n = 1:nF % Loop through all the filaments attached to this structure
                                f = idx1(n);
                                yprevious = Filaments.XYCoords{f}(end,2);

                                % START Calculate New Filament Position ------------------------------------------------------------------------------------------
                                Filaments.XYCoords{f}(:,1) = Filaments.XYCoords{f}(:,1)  +  (Fx + FAx     ) * dt/gamma1;   % Compute new position in X direction
                                Filaments.XYCoords{f}(:,2) = Filaments.XYCoords{f}(:,2)  +  (Fy + FAy + Fm) * dt/gamma1;   % Compute new position in Y direction
                                % END Calculate New Filament Position --------------------------------------------------------------------------------------------

                                % Record Filament Speed, Position, Region, etc of each tip -------------------------------------------------
                                YDist = (Filaments.XYCoords{f}(end,2) - yprevious); % Calculate retrogade flow velocity
                                Data.YSpeed       = [ Data.YSpeed;       single(YDist/ModelParameters.TimeStep) ];
                                Data.FilamentName = [ Data.FilamentName; single(Filaments.Name(f,1)) ];
                                Data.XYPosition   = [ Data.XYPosition;   single(Filaments.XYCoords{f}(end,:)) ];
                                Data.StructName   = [ Data.StructName;   single(Filaments.MainIndex(f,1)) ];
                                
                                % Record for mirroring -------------------------------------------------------------------------------------
                                tipXpositions(n,1) = Filaments.XYCoords{f}(end,1); 
                        end
                        

                        % START HORIZONTAL FILAMENT PERIODIC BOUNDARY ---------------------------------------------------------------------------------------------------------------
                        BreakIdx = [];
                        Offset = 0;
                            % Check if any of the filament tips in the structure are crossing the left or right membrane edge
                            if  min(tipXpositions) < Membrane.Nodes.Xcoords(1) 
                                Offset = Membrane.Nodes.Xcoords(2) - max(tipXpositions); 
                            elseif max(tipXpositions) > Membrane.Nodes.Xcoords(2) 
                                Offset = Membrane.Nodes.Xcoords(1) - min(tipXpositions);
                            end
                            % If there are any filaments out of bounds, move the whole structure to the otherside just within the bounds
                            if ~isequal(Offset,0)
                                        % Apply offset to each filament in the structure
                                        for n = 1:nF
                                            f = idx1(n);
                                            Filaments.XYCoords{f}(:,1) = Filaments.XYCoords{f}(:,1) + Offset; % Apply Offset
                                            % Find adhesion connections if they exists, and delete connection
                                            [Filaments,Integrins,Ligands,FALconnections] = DeleteFALconnection('FilamentName',Filaments.Name(f,1),Filaments,Integrins,Ligands,FILconnections);
                                        end
                            end
                        % END HORIZONTAL FILAMENT PERIODIC BOUNDARY ---------------------------------------------------------------------------------------------------------------
                end 
    end
    
    
    
    %% MEMBRANE MOTION SECTION
    
                % Calculate new positions for membrane -------------------------------
                Nu = ModelParameters.CytoplasmViscosity*(10^-6); % pN s/nm^2
                L  = ModelParameters.MembraneWidth;  % membrane width in nm
                r  = ModelParameters.MembraneRadius; % nm 
                ModelParameters.MembraneGamma = (4*pi*Nu*L) / (0.84+log(L/(2*r) )); % gamma for random fluid motion (pN*s/nm) 
                gamma3 = ModelParameters.MembraneGamma; % gamma for random fluid motion (pN*s/nm) 
                
                if ~isempty(Filaments.Name)
                        Ff = CalculateForceFromFilamentsHittingMembraneSegment(Filaments,Membrane,ModelParameters);
                        Ff = sum(Ff,'omitnan');
                        Membrane.Nodes.Ycoords = Membrane.Nodes.Ycoords + Ff*dt/gamma3; % Move y-position of membrane
                end

    
    %% INTEGRIN MOTION SECTION (not connected)
    
                gamma4 = ModelParameters.IntegrinGamma;
                Aidx = find( ~Integrins.ActiveStatus ); % Find all un-activated or unattached adhesion indices
                nA = length(Aidx);
                SD = sqrt(2*kT*gamma4/dt);     

                if ~isempty(Aidx)
                    Fx = SD*randn(nA,1); % Calculate Brownian motion forces
                    Fy = SD*randn(nA,1);
                    Integrins.XYpoints(Aidx,1) = Integrins.XYpoints(Aidx,1) + Fx*dt/gamma4;
                    Integrins.XYpoints(Aidx,2) = Integrins.XYpoints(Aidx,2) + Fy*dt/gamma4;
                end
    
                
end






%% START of Sub-functions ============================================================================================================================================================================================================

function Fm = TotalForceExertedByMembraneOnFilamentStructure(FilInd,Filaments,Membrane,ModelParameters)
            
            Fm = 0;
            for f = FilInd' % Calculate the force of the membrane is exerting on all the filaments in this structure
                    y = Membrane.Nodes.Ycoords(1) - Filaments.XYCoords{f}(end,2); % Distance of filament tip from membrane
                    %|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                    % Evaluate normal force of filament against membrane velocity eq.5 (Mogliner and Oser 1996);
                    theta = abs(atand( Filaments.UnitVector(f,1)/Filaments.UnitVector(f,2) ));
                    theta(theta < 1 ) = 1; % For model stability, keep theta greater than 1 degree
                    lambda = ModelParameters.PersistenceLength;
                    kT = ModelParameters.kT;
                    L = length(Filaments.MonomerIndices{f})*ModelParameters.MonomerLength;
                    L( L<30  ) = 30;
                    L( L>150 ) = 150;
                    delta = ModelParameters.MonomerLength*cosd(theta);
                    K = 4*lambda*kT/(L.^3*sind(theta).^2); % eq. B.1
                    if y >= 0
                        Fn = -sqrt(2*kT*K/pi).*exp(-K*y.^2/(2*kT))./( erf(y*sqrt(K/(2*kT))) + 1 ); % Normal force calculation (force is in -y direction)
                    else
                        Fn = -sqrt(2*kT*K/pi).*exp(-K*y.^2/(2*kT))./( erfc( abs( y*sqrt(K/(2*kT)) ) ) );
                    end
                    %|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                    Fn( isnan(Fn) | isinf(Fn) ) = -100;
                    Fm = Fm + Fn;
            end

end


%======================================================================================================
%======================================================================================================

function Ff = CalculateForceFromFilamentsHittingMembraneSegment(Filaments,Membrane,ModelParameters)
        
        
        % Add up all the filament forces on the membrane
        Ff = []; %0;
        nF = length(Filaments.Name);
        for f = 1:nF
            y = Membrane.Nodes.Ycoords(1) - Filaments.XYCoords{f}(end,2); % Distance of filament tip from membrane
            %|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
            % Evaluate force of filament against membrane  eq.5 (Mogliner and Oser 1996);
            theta = abs( atand( Filaments.UnitVector(f,1) / Filaments.UnitVector(f,2) ) );
            theta(theta < 1 ) = 1; % For model stability, keep theta greater than 1 degree 
            lambda = ModelParameters.PersistenceLength;
            kT = ModelParameters.kT;
            L = length(Filaments.MonomerIndices{f}) * ModelParameters.MonomerLength;
            L( L<30  ) = 30;
            L( L>150 ) = 150;
            delta = ModelParameters.MonomerLength * cosd(theta);
            K = 4*lambda*kT / (L.^3 * sind(theta).^2) ; % eq. B.1
            
            if y >= 0
                Fm = sqrt(2*kT*K/pi).*exp(-K*(y.^2)/(2*kT))./( erf(y*sqrt(K/(2*kT))) + 1 ); % Normal force calculation (force is in +y direction)
            else
                Fm = sqrt(2*kT*K/pi).*exp(-K*(y.^2)/(2*kT))./( erfc( abs( y*sqrt(K/(2*kT)) ) ) ); % Variation of same equation do to bad precision of erf function near 0
            end

            Fm( isnan(Fm) | isinf(Fm) ) = 100;
            %|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
            Ff = [Ff;Fm]; %Ff + Fm;
        end   
       
end

%======================================================================================================
%======================================================================================================
        
function [FAx,FAy,Integrins,Ligands,FILconnections,IntegrinTensions] = CalculateForceDueToFILconnections(idx,Filaments,Integrins,Ligands,IntegrinTensions,FILconnections,ModelParameters)


        [A,B,C,D] = MolecularClutchPeakParameters(ModelParameters.MolecularClutch_PeakNumber);

    % PeakNumber = 1 or 2
    % k_off =  A*exp(B*ConnectionTension) + C*exp(D*ConnectionTension);
        %IntegrinTensions = zeros(size(Integrin.ActiveStatus));
        FAx = 0;
        FAy = 0;
        Aidx = [];
        index = 0;

        if ~isempty(FILconnections.IntegrinIndex)
                    for f = idx'
                         con = find(FILconnections.FilamentName == Filaments.Name(f,1)); % Find adhesions attached to this filament
                         if ~isempty(con)
                                 for c = con'
                                     index = index + 1;
                                     midx = find( Filaments.MonomerIndices{f} == FILconnections.MonomerIndex(c) ); % Find index of attached monomer (used to get XY coords of monomer)
                                     aidx = FILconnections.IntegrinIndex(c,1);                                     % Find index of adhesion in integrins (used to get XY coords of integrins)
                                     xDist = Integrins.XYpoints(aidx,1) - Filaments.XYCoords{f}(midx,1); % X Distance between integrin and attached filament monomer
                                     yDist = Integrins.XYpoints(aidx,2) - Filaments.XYCoords{f}(midx,2); % Y Distance between integrin and attached filament monomer
                                     SeparationDist = sqrt( xDist^2 + yDist^2 );                         % Distance between integrin and attached filament monomer
                                     StretchDist = SeparationDist - ModelParameters.IntegrinSpringEqLength;  % Stretch distance
                                     StretchDist(StretchDist < 0) = 0;                                       % If stretch distance is less than Equilibrium length, StretchDist = 0;
                                     ConnectionTension = StretchDist * ModelParameters.IntegrinSpringConstant; % Fa = ka*x (Calculate spring force between integrin and filament

                                     ForceX = ConnectionTension*(xDist/SeparationDist);   % Fx = F*(x component of separation distance unit vector)
                                     ForceY = ConnectionTension*(yDist/SeparationDist);   % Fy = F*(y component of separation distance unit vector)

                                     
                                     % If Molecular Clutch is turned on, deactivation of adhesion rate is dependent on tension
                                     k_off =  A*exp(B*ConnectionTension) + C*exp(D*ConnectionTension);
                                     if isinf(k_off) 
                                         k_off = 1E20;
                                         % disp('Adhesion tension too high. k_off = Inf'); 
                                     end
                                     
                                     if rand < (1 - exp(-k_off*ModelParameters.TimeStep))    &&   ModelParameters.Adhesion_MolecularClutchOn
                                         Aidx  = [Aidx; aidx]; % Acculumate adhesion indices for connections that were selected to be broken
                                     else
                                         IntegrinTensions(aidx,1) = ConnectionTension; 
                                         FAx = FAx + ForceX;  % Accumulate all the forces from adhesions being exerted on this filament structure
                                         FAy = FAy + ForceY;
                                     end

                                 end
                         end
                    end

                    % Remove connections of deactivated Integrin-Filament connections (only executes if MC is on. Otheriwse BreakIdx is always empty)-----------------
                    if ModelParameters.Adhesion_MolecularClutchOn & ~isempty(Aidx)
                        [Filaments,Integrins,Ligands,FILconnections] = DeleteFILconnection('IntegrinIndex',Aidx,Filaments,Integrins,Ligands,FILconnections);
                    end
           
        end
%         FAx = sum(FAx);  if isempty(FAx); FAx = 0; end  % Sum all the forces from adhesions being exerted on this filament structure
%         FAy = sum(FAy);  if isempty(FAy); FAy = 0; end             

end

