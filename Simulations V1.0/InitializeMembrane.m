function Membrane = InitializeMembrane(ModelParameters)
   
    
    % ModelParameters.MembraneWidth = 500; % nm
    % ModelParameters.ModelDepth = 1000; % nm
    
    Membrane.Nodes.Xcoords = [0, ModelParameters.MembraneWidth];
    Membrane.Nodes.Ycoords = [0, 0];
    
end
