function Filaments = InitializeActinFilaments(ModelParameters,Membrane)
    % This function creates all the initial actin filaments
    N = ModelParameters.StartingNumberOfFilaments;
    theta = 45+90*rand(N,1); %180*rand(N,1);

    %----------------------------------------------------------
    SD  = 30;
    InitialLengths = abs( round(ModelParameters.FilamentsInitialLength   +   SD*randn(N,1)) ); % In units of monomers
    VerticalOffSets = ModelParameters.VerticalOffSet * rand(N,1);
    Filaments.Name = (1:N)';
    Filaments.MonomerIndices = cell(N,1);
    Filaments.XYCoords = cell(N,1);
    Filaments.UnitVector = zeros(N,2);
    D = ModelParameters.MonomerLength; %(nm)
    
    for f = 1:N
        % Create initial two points of filament as a horizontal line (origin [0,0] is implied to be the second point).
            x = ModelParameters.MonomerLength;
            y = 0;
        % Now randomly rotate between 0 and 180 degrees
            R = [ [cosd(theta(f,1)), -sind(theta(f,1))];... % create rotation matrix
                  [sind(theta(f,1)),  cosd(theta(f,1))] ];
             
            M = R*[x;y]; 
            xr = M(1); % Where second point is the end closest to the membrane
            yr = M(2); % Where second point is the end closest to the membrane
            
            Filaments.MonomerIndices{f,1}  = 1; %% Index values for the first point of the filament
            Filaments.XYCoords{f,1} = [0,0]; % coordinate of first point
            Filaments.UnitVector(f,:) = [xr,yr]./vecnorm([xr;yr]);
            % plot(x,y,'-b',xr,yr,'-r',x(2),y(2),'.b',xr(2),yr(2),'.r'); axis equal; axis([-3,3,-3,3]); pause
            
            % Make Filaments the desired length
            for k = 1:InitialLengths(f) % Each loops adds a monomer to filament f
                Filaments.MonomerIndices{f} = [Filaments.MonomerIndices{f}; Filaments.MonomerIndices{f}(end,1)+1;]; 
                Filaments.XYCoords{f} = [Filaments.XYCoords{f}; [Filaments.XYCoords{f}(end,1) + D*Filaments.UnitVector(f,1),...
                                                                 Filaments.XYCoords{f}(end,2) + D*Filaments.UnitVector(f,2)] ];
            end
            
            % Make the filament's tip at [0,0] 
            Filaments.XYCoords{f}(:,1) = Filaments.XYCoords{f}(:,1) - Filaments.XYCoords{f}(end,1);
            Filaments.XYCoords{f}(:,2) = Filaments.XYCoords{f}(:,2) - Filaments.XYCoords{f}(end,2);
            
            % Spread out Filaments in the x direction and apply y-offset --------------------------------------
            Xrand = ModelParameters.MembraneWidth*rand;
            Filaments.XYCoords{f}(:,1) = Filaments.XYCoords{f}(:,1) + Xrand;
            Filaments.XYCoords{f}(:,2) = Filaments.XYCoords{f}(:,2) + VerticalOffSets(f,1); 
    end
    
    % Setup starting/default parameters
    Filaments.IsCapped = false(N,1);          % Has the filament been capped?
    Filaments.MainIndex = (1:N)';             % Main Filament group that this filament is a part of. (based on parent filament name: Filament structure)
    Filaments.Parent = zeros(N,1);            % The filament parent name that each filament branched off of (zero if it's a main filament)
    Filaments.ParentIndex  = zeros(N,1);      % The Monomer Index on the parent filament where this filament branched from.
  

% FULL FILAMENT DESCRIPTION:    
% Combine all the parameters for this filament ----------------------------------------------------
%         Filaments.Name           = [ Filaments.Name; NewFilamentName];                 % This is a unique number assigned to each filament regardless if it's a main filament or branch
%         Filaments.MonomerIndices = [ Filaments.MonomerIndices; {MonomerIndices} ];     % A unique number assigned to each monomer
%         Filaments.XYCoords       = [ Filaments.XYCoords; {XYCoords}];                  % The XY coordinates of each monomer in this filament (the barded-end (tip) is at XY(end,:)
%         Filaments.UnitVector     = [ Filaments.UnitVector; UnitVector ];               % Unit vector giving the direction the tip is pointing
%         Filaments.IsCapped       = [ Filaments.IsCapped; false ];                      % Logical array: true if filament is capped
%         Filaments.MainIndex      = [ Filaments.MainIndex; NewFilamentName ];           % All filaments belonging to the same structure have the MainIndex value (which is the name of the main filament all the branches and branches of branches branched from)
%         Filaments.Parent         = [ Filaments.Parent; 0 ];                            % Name of the filament each branch branched from ( equals 0 if it's main filament)
%         Filaments.ParentIndex    = [ Filaments.ParentIndex; 0 ];                       % Monomer index number of the parent filament the branch branched from (0 if it's a main filament)
    
end 


