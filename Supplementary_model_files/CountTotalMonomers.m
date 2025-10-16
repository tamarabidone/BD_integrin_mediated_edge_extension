function [nMonomers, DeletedFilamentNames]  = CountTotalMonomers(Filaments, Membrane, ModelParameters)

    %nF = length(MVAR.Filaments.MonomerIndices);
    DeletedFilamentNames = [];
    nMonomers = 0;
    idxMF = find(Filaments.Parent == 0);
    nMF = length(idxMF); % Total number of filament structures
    BoundaryEdge = Membrane.Nodes.Ycoords(end) - ModelParameters.ModelDepth;
    for MF = 1:nMF
        Structmonomers = 0;
        Outsidemonomers = 0;
        idx1 = find( Filaments.MainIndex == Filaments.Name(idxMF(MF)) ); 
        nF = length(idx1);
    for f = 1:nF
        Structmonomers = Structmonomers + length(Filaments.XYCoords{f}(:,2));
        outRange = sum( (Filaments.XYCoords{f}(:,2)) < BoundaryEdge );
        Outsidemonomers = Outsidemonomers + outRange;
        inRange =  sum( (Filaments.XYCoords{f}(:,2)) > BoundaryEdge );
        nMonomers = nMonomers + inRange;
    end
        if Outsidemonomers == Structmonomers
            names = Filaments.Name(idx1,1);
            DeletedFilamentNames = [DeletedFilamentNames; names];
        end

    end