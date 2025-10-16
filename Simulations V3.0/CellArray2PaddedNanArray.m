function [ArrayOut,nCells,Ns] = CellArray2PaddedNanArray(ArrayIn)

        nCells = size(ArrayIn,1);
        Ns = cellfun(@length ,ArrayIn);
        MaxDataLength = max(Ns); % Find the length of the longest vector
        ArrayOut = NaN(MaxDataLength,nCells);
        
        for n = 1:nCells
            L = length(ArrayIn{n});
            ArrayOut(1:L,n) = ArrayIn{n};
        end

end