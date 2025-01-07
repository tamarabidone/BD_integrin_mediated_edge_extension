function Integrins = InitializeIntegrins(ModelParameters,Membrane)
    
%     ModelParameters.MemWidth = 500; % nm
%     ModelParameters.ModelDepth = 1000; % nm
%     ModelParameters.IntegrinTotal = 200; 
    
    nA = ModelParameters.IntegrinTotal;  
    
    Integrins.XYpoints      = [ModelParameters.MembraneWidth*rand(nA,1), -ModelParameters.ModelDepth*rand(nA,1)];
    Integrins.ActiveStatus  = false(nA,1);
    Integrins.AttachedFilamentName = NaN(nA,1);
    Integrins.AttachedLigandIndex  = NaN(nA,1);
    Integrins.Orientation = 360*rand(nA,1);
    Integrins.nBonds = zeros(nA,1);
    
    Integrins.RegionNodes = zeros(4,2); % Since there is one membrane segment, create one adhesion region
    

       XY_LT = [ Membrane.Nodes.Xcoords(1), Membrane.Nodes.Ycoords(1) ]; %  Left-Top  xy-coordinates of current adhesion region
       XY_RT = [ Membrane.Nodes.Xcoords(2), Membrane.Nodes.Ycoords(2) ]; %  Right-Top xy-coordinates of current adhesion region
       XY_RB = [ Membrane.Nodes.Xcoords(2), Membrane.Nodes.Ycoords(2) - ModelParameters.ModelDepth]; %  Right-Bottom xy-coordinates of current adhesion region
       XY_LB = [ Membrane.Nodes.Xcoords(1), Membrane.Nodes.Ycoords(1) - ModelParameters.ModelDepth]; %  Left-Bottom  xy-coordinates of current adhesion region


       % Record Current Region
       Integrins.RegionNodes(1,:) = XY_LT;           
       Integrins.RegionNodes(2,:) = XY_RT;            
       Integrins.RegionNodes(3,:) = XY_RB; 
       Integrins.RegionNodes(4,:) = XY_LB; 

            
end

