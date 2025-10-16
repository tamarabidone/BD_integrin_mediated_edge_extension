 function [Filaments, Membrane, Adhesions, Ligands, FALconnections, AdhesionTensions, Data] = ...
                                  CalculatePositionAfterAppliedForces(Filaments,Membrane,Adhesions,Ligands,FALconnections,ModelParameters,t)
                              

    % Find all the main filaments of connected filament structures (group of attached filaments)
    idxMF = find(Filaments.Parent == 0);
    nMF = length(idxMF); % Total number of filament structures
    dt = ModelParameters.TimeStep;
    kT = ModelParameters.kT;
    FilamentTips = GetFilamentsTipLocations(Filaments);
    AdhesionTensions = zeros(size(Adhesions.ActiveStatus));
    
    nF_1 = length(Filaments.Name);
    Ysp = NaN(length(Filaments.Name),1);

    Data.YSpeed    = NaN(nF_1,1, 'single');
    Data.FilamentName = NaN(nF_1,1, 'single');
    Data.XYPosition = NaN(nF_1,2, 'single');
    Data.StructName = NaN(nF_1,1, 'single');
    Data.Parent = NaN(nF_1,1, 'single');
    Data.Orientation = NaN(nF_1,2, 'single');
    Data.Main = NaN(nF_1,1, 'single');
    Data.FilamentLength = NaN(nF_1,1, 'single');
    Data.DistanceMembrane = NaN(nF_1,1, 'single');
    Data.ForceonMembrane = NaN(nF_1,1, 'single');
    Data.AdhesionForceComp = []; % [Fax,Fay,Fname]
    idxF = 0;
    %selectedPercentage = ModelParameters.MyosinFraction;
    
    %% FILAMENT MOTION SECTION 
    if ~isempty(Filaments.Name)
            % Apply myosin force to selected filaments *before* looping through structures
            if ModelParameters.ApplyMyosin  %ismember(f, selectedIndices)
                if rem( round(t,10),0.001) == 0
                    for f = 1:length(Filaments.Name)
                        yprevious = Filaments.XYCoords{f}(end,2);
                        %Filaments.XYCoords{f}(:,2) = Filaments.XYCoords{f}(:,2) - (0.002 * rand(1));
                         Nu = ModelParameters.CytoplasmViscosity*(10^-6); % pN s/nm^2
                         r  = 7/2; % nm
                         L = length(Filaments.MonomerIndices{f}) * ModelParameters.MonomerLength;
                         gamma1 =  (4*pi*Nu*L) / (0.84+log(L/(2*r) )); % gamma for random fluid motion (pN*s/nm) 
                         SD = sqrt(2*kT*gamma1/dt);  
                        Filaments.XYCoords{f}(:,2) = Filaments.XYCoords{f}(:,2)  -  (rand * 1) * dt/gamma1;
                        Ysp(f,1) = (Filaments.XYCoords{f}(end,2) - yprevious)/ModelParameters.TimeStep;
                    end
                end
            end

                for MF = 1:nMF % Go through each filament structure and test each of its filaments to see if it is hitting the membrane and calculate the total force acting on the structure.
                        L = 0; % Total length of Filament
                        idx1 = find( Filaments.MainIndex == Filaments.Name(idxMF(MF)) ); % Find all filaments attached to same structure
                        nF = length(idx1); % number of filaments in this structure
                    % Calculate total length of filament structure
                        for f = idx1'  
                            L = L + length(Filaments.MonomerIndices{f}) * ModelParameters.MonomerLength; % L is total length of filament structure in microns
                        end 

                    % CALCULATE FORCES DUE TO FILAMENT-ADHESION CONNECTIONS and break connections where force is greater than threshold
                       [FAx,FAy,FAvalues,Adhesions,Ligands,FALconnections,AdhesionTensions] = CalculateForceDueToFALconnections(idx1,Filaments,Adhesions,Ligands,AdhesionTensions,FALconnections,ModelParameters,Ysp);
                        disp(FAy)
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
                        % if ModelParameters.ApplyMyosin  %ismember(f, selectedIndices)
                        %     MyosForce = 1;
                        % else
                        %     MyosForce = 0;
                        % end

                        tipXpositions = NaN(nF,1);
                        for n = 1:nF % Loop through all the filaments attached to this structure
                                idxF = idxF +1;
                                f = idx1(n);
                                yprevious = Filaments.XYCoords{f}(end,2);

                                % START Calculate New Filament Position ------------------------------------------------------------------------------------------
                                  % Compute new position in X direction
                                  % 
                                    Filaments.XYCoords{f}(:,2) = Filaments.XYCoords{f}(:,2)  +  (Fy + FAy + Fm) * dt/gamma1;   % Compute new position in Y direction
                                    Filaments.XYCoords{f}(:,1) = Filaments.XYCoords{f}(:,1)  +  (Fx + FAx) * dt/gamma1;
                                  % else
                                  %   Filaments.XYCoords{f}(:,2) = Filaments.XYCoords{f}(:,2)  +  (Fy + FAy + Fm) * dt/gamma1;   % Compute new position in Y direction
                                  %   Filaments.XYCoords{f}(:,1) = Filaments.XYCoords{f}(:,1)  +  (Fx + FAx) * dt/gamma1;
                                  % end
 
                                % END Calculate New Filament Position --------------------------------------------------------------------------------------------

                                % Record Filament Speed, Position, Region, etc of each tip -------------------------------------------------
                                YDist = (Filaments.XYCoords{f}(end,2) - yprevious); % Calculate retrogade flow velocity
                                Data.YSpeed(idxF,1)       = single(YDist/ModelParameters.TimeStep);
                                Data.FilamentName(idxF,1) = single(Filaments.Name(f,1)) ;
                                Data.XYPosition(idxF,:)   = single(Filaments.XYCoords{f}(end,:)) ;
                                Data.StructName(idxF,1)   = single(Filaments.MainIndex(f,1)) ;
                                Data.Parent(idxF,1)       = single(Filaments.Parent(f,1)) ;
                                Data.Orientation(idxF,:)  = single(Filaments.UnitVector(f,:)) ;
                                Data.Main(idxF,1)         = single(Filaments.MainIndex(f,:)) ;
                                Data.FilamentLength(idxF,1) = single(length(Filaments.MonomerIndices{f})*ModelParameters.MonomerLength);
                                Data.DistanceMembrane(idxF,1) = single(Membrane.Nodes.Ycoords(1,1)- Filaments.XYCoords{f}(end,2));
                                Data.AdhesionForceComp = [Data.AdhesionForceComp; single(FAvalues)];
                                % Record for mirroring -------------------------------------------------------------------------------------
                                tipXpositions(n,1) = Filaments.XYCoords{f}(end,1); 
                        end
                        

                        % START HORIZONTAL FILAMENT PERIODIC BOUNDARY ---------------------------------------------------------------------------------------------------------------
                        % BreakIdx = [];
                        % Offset = 0;
                        %     % Check if any of the filament tips in the structure are crossing the left or right membrane edge
                        %     if  min(tipXpositions) < Membrane.Nodes.Xcoords(1) 
                        %         Offset = Membrane.Nodes.Xcoords(2) - max(tipXpositions); 
                        %     elseif max(tipXpositions) > Membrane.Nodes.Xcoords(2) 
                        %         Offset = Membrane.Nodes.Xcoords(1) - min(tipXpositions);
                        %     end
                        %     % If there are any filaments out of bounds, move the whole structure to the otherside just within the bounds
                        %     if ~isequal(Offset,0)
                        %                 for n = 1:nF
                        %                     f = idx1(n);
                        %                     Filaments.XYCoords{f}(:,1) = Filaments.XYCoords{f}(:,1) + Offset; % Apply Offset
                        %                     % Find adhesion connections if they exists, and delete connection
                        %                     [Filaments,Adhesions,Ligands,FALconnections] = DeleteFALconnection('FilamentName',Filaments.Name(f,1),Filaments,Adhesions,Ligands,FALconnections);
                        %                 end
                        %     end
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
                %gamma3 = 1;

                if ~isempty(Filaments.Name)
                        Ff = CalculateForceFromFilamentsHittingMembraneSegment(Filaments,Membrane,ModelParameters);
                        Data.ForceonMembrane = Ff;
                        Ff = sum(Ff,'omitnan');
                        Membrane.Nodes.Ycoords = Membrane.Nodes.Ycoords + Ff*dt/gamma3; % Move y-position of membrane
                end

    
    %% ADHESION MOTION SECTION (not connected)
    
                gamma4 = ModelParameters.AdhesionGamma;
                Aidx = find( ~Adhesions.ActiveStatus ); % Find all un-activated or unattached adhesion indices
                nA = length(Aidx);
                SD = sqrt(2*kT*gamma4/dt);     

                if ~isempty(Aidx)
                    Fx = SD*randn(nA,1); % Calculate Brownian motion forces
                    Fy = SD*randn(nA,1);
                    Adhesions.XYpoints(Aidx,1) = Adhesions.XYpoints(Aidx,1) + Fx*dt/gamma4;
                    Adhesions.XYpoints(Aidx,2) = Adhesions.XYpoints(Aidx,2) + Fy*dt/gamma4;
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
                    if ModelParameters.ReducedFlow
                        if y >= 0
                        Fn = (-sqrt(2*kT*K/pi).*exp(-K*(y.^2)/(2*kT))./( erf(y*sqrt(K/(2*kT))) + 1 ))/3; % Normal force calculation (force is in -y direction)
                        else
                        Fn = (-sqrt(2*kT*K/pi).*exp(-K*(y.^2)/(2*kT))./( erfc( abs( y*sqrt(K/(2*kT)) ) ) ))/3; % Variation of same equation do to bad precision of erf function near 0
                        end
                    else
                        if y >= 0
                            Fn = -sqrt(2*kT*K/pi).*exp(-K*y.^2/(2*kT))./( erf(y*sqrt(K/(2*kT))) + 1 ); % Normal force calculation (force is in -y direction)
                        else
                            Fn = -sqrt(2*kT*K/pi).*exp(-K*y.^2/(2*kT))./( erfc( abs( y*sqrt(K/(2*kT)) ) ) );
                        end
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
            
            if ModelParameters.ReducedFlow
                if y >= 0
                Fm = (sqrt(2*kT*K/pi).*exp(-K*(y.^2)/(2*kT))./( erf(y*sqrt(K/(2*kT))) + 1 ))/3; % Normal force calculation (force is in +y direction)
                else
                Fm = (sqrt(2*kT*K/pi).*exp(-K*(y.^2)/(2*kT))./( erfc( abs( y*sqrt(K/(2*kT)) ) ) ))/3; % Variation of same equation do to bad precision of erf function near 0
                end
            else
                if y >= 0
                    Fm = sqrt(2*kT*K/pi).*exp(-K*(y.^2)/(2*kT))./( erf(y*sqrt(K/(2*kT))) + 1 ); % Normal force calculation (force is in +y direction)
                else
                    Fm = sqrt(2*kT*K/pi).*exp(-K*(y.^2)/(2*kT))./( erfc( abs( y*sqrt(K/(2*kT)) ) ) ); % Variation of same equation do to bad precision of erf function near 0
                end
            end

            Fm( isnan(Fm) | isinf(Fm) ) = 100;
            disp(Fm)
            %|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
            Ff = [Ff;Fm]; %Ff + Fm;
        end   
       
end

%======================================================================================================
%======================================================================================================

% function [] = CalculateTensionPerAdhesion(Filaments,Adhesions,Ligands,AdhesionTensions,FALconnections,ModelParameters)
% 
% [a,b,c,d] = MolecularClutchPeakParameters(ModelParameters.MolecularClutch_PeakNumber);
% 
% 
% 
% 
% end

%======================================================================================================
%======================================================================================================
        
function [FAx,FAy,FAvalues,Adhesions,Ligands,FALconnections,AdhesionTensions] = CalculateForceDueToFALconnections(idx,Filaments,Adhesions,Ligands,AdhesionTensions,FALconnections,ModelParameters,Ysp)


        [A,B,C,D] = MolecularClutchPeakParameters(ModelParameters.MolecularClutch_PeakNumber);
        % [a,b,c,d] = MyosinPeakParameters();

    % PeakNumber = 1 or 2
    % k_off =  A*exp(B*ConnectionTension) + C*exp(D*ConnectionTension);
        %AdhesionTensions = zeros(size(Adhesions.ActiveStatus));
        FAx = 0;
        FAy = 0;
        Aidx = [];
        index = 0;
        FAvalues = [];
        nAdhesions = length(find(Adhesions.ActiveStatus));
        v_opt = 8;
        sigma = 0.2;  
        k = 0.05;     

        if ~isempty(FALconnections.AdhesionIndex)
                    for f = idx'
                         con = find(FALconnections.FilamentName == Filaments.Name(f,1)); % Find adhesions attached to this filament
                         %v = Ysp(f,1);
                         %MyosScalingFactor = (d + a .* (1 ./ ((abs(v) + 0.01) .* c .* sqrt(2 .* pi))) .* exp(-((log(abs(v) + 0.01) - b).^2) ./ (2 .* c.^2)))/48.9227;
                         if ~isempty(con)
                                 for c = con'
                                     index = index + 1;
                                     midx = find( Filaments.MonomerIndices{f} == FALconnections.MonomerIndex(c) ); % Find index of attached monomer (used to get XY coords of monomer)
                                     aidx = FALconnections.AdhesionIndex(c,1);                                     % Find index of adhesion in Adhesions (used to get XY coords of adhesion)
                                     xDist = Adhesions.XYpoints(aidx,1) - Filaments.XYCoords{f}(midx,1); % X Distance between adhesion and attached filament monomer
                                     yDist = Adhesions.XYpoints(aidx,2) - Filaments.XYCoords{f}(midx,2); % Y Distance between adhesion and attached filament monomer
                                     SeparationDist = sqrt( xDist^2 + yDist^2 );                         % Distance between adhesion and attached filament monomer
                                     StretchDist = SeparationDist - ModelParameters.AdhesionSpringEqLength;  % Stretch distance
                                     StretchDist(StretchDist < 0) = 0;                                       % If stretch distance is less than Equilibrium length, StretchDist = 0;
                                     %ConnectionTension = StretchDist * ModelParameters.AdhesionSpringConstant;
                                     %ConnectionTension = StretchDist * ModelParameters.AdhesionSpringConstant .*(1 + exp(abs(StretchDist) / 27.5))* (v * exp(-v / v_opt)); % Fa = ka*x (Calculate spring force between adhesion and filament
                                     ConnectionTension = StretchDist * ModelParameters.AdhesionSpringConstant; %* MyosScalingFactor;
                                     ForceX = ConnectionTension*(xDist/SeparationDist);   % Fx = F*(x component of separation distance unit vector)
                                     ForceY = ConnectionTension*(yDist/SeparationDist);   % Fy = F*(y component of separation distance unit vector)

                                     FAvalues = [FAvalues; [ForceX,ForceY,Filaments.Name(f,1)] ];
                                     % If Molecular Clutch is turned on, deactivation of adhesion rate is dependent on tension
                                     %k_off =  A*exp(B*ConnectionTension) + C*exp(D*ConnectionTension);
                                     if ModelParameters.ExponentialLifetime == 2
                                        L0 = 6; 
                                        k_off = 1./(L0 * exp(-k * ConnectionTension));
                                     elseif ModelParameters.ExponentialLifetime == 3
                                        L0 = 0.2; 
                                        k_off = 1./(L0 * exp(-k * ConnectionTension));
                                     else
                                        k_off =  A*exp(B*ConnectionTension) + C*exp(D*ConnectionTension);
                                     end
                                        if isinf(k_off) 
                                         k_off = 1E20;
                                         % disp('Adhesion tension too high. k_off = Inf'); 
                                     end
                                     
                                     if rand < (1 - exp(-k_off*ModelParameters.TimeStep))    &&   ModelParameters.Adhesion_MolecularClutchOn
                                         Aidx  = [Aidx; aidx]; % Acculumate adhesion indices for connections that were selected to be broken
                                     else
                                         AdhesionTensions(aidx,1) = ConnectionTension; 
                                         FAx = FAx + ForceX;  % Accumulate all the forces from adhesions being exerted on this filament structure
                                         FAy = FAy + ForceY;
                                     end

                                 end
                         end
                    end

                    % Remove connections of deactivated Adhesion-Filament connections (only executes if MC is on. Otheriwse BreakIdx is always empty)-----------------
                    if ModelParameters.Adhesion_MolecularClutchOn & ~isempty(Aidx)
                        [Filaments,Adhesions,Ligands,FALconnections] = DeleteFALconnection('AdhesionIndex',Aidx,Filaments,Adhesions,Ligands,FALconnections);
                    end
           
        end
%         FAx = sum(FAx);  if isempty(FAx); FAx = 0; end  % Sum all the forces from adhesions being exerted on this filament structure
%         FAy = sum(FAy);  if isempty(FAy); FAy = 0; end             

end

