function [A,B,C,D] = MolecularClutchPeakParameters(PeakNumber)

    % PeakNumber = 1 or 2
    % k_off =  A*exp(B*ConnectionTension) + C*exp(D*ConnectionTension);
    
        switch PeakNumber
              case 1 % Peak 1 (Tau_max = 3 2)
                A =  2;
                B = -0.064;
                C =  0.0005;
                D =  0.288;
            case 2 % Peak 2 (Tau_max = 7.5 s)
                A =  0.5;
                B = -0.0488;
                C =  0.00005;
                D =  0.261;
        end
   
end
