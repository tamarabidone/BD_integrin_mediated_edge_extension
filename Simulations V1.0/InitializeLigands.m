function Ligands = InitializeLigands(ModelParameters)
    
%    
    
    nL = ModelParameters.LigandTotal;  
            
    Ligands.XYpoints               = [ ModelParameters.MembraneWidth * rand(nL,1),... % These are the current XY positions of the ligands. It is updated each model iteration.
                                       -ModelParameters.ModelDepth * rand(nL,1) ];
                                                               
    Ligands.AttachedFilamentName   = NaN(nL,1);
    Ligands.AttachedIntegrinIndex  = NaN(nL,1);
            
end
