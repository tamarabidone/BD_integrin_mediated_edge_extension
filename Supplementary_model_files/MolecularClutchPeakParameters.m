function [A,B,C,D] = MolecularClutchPeakParameters(PeakNumber)

    % PeakNumber = 1 or 2
    % k_off =  A*exp(B*ConnectionTension) + C*exp(D*ConnectionTension);

    
        switch PeakNumber
            case 1 
                A =  2;
                B = -0.064;
                C =  0.0005;
                D =  0.288;
            case 2 
                A =  0.5;
                B = -0.0488;
                C =  0.00005;
                D =  0.261;
            case 3 
                A =  1.25;
                B = -0.0357;
                C =  0.0005;
                D =  0.2657;
            case 4
               A =  0.9;
               B = -0.0155;
               C = 0.0005;
               D = 0.2317;

        end
   
end
