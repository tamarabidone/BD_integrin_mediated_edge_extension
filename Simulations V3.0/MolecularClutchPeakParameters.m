function [A,Alpha,B,Beta] = MolecularClutchPeakParameters(PeakNumber)

    % PeakNumber = 1, 2 or 3
    % k_off =  A*exp(Alpha*ConnectionTension) + B*exp(Beta*ConnectionTension);
    
        switch PeakNumber
            case 1 % Wild Type 
                A     =  2;
                Alpha = -0.064;
                B     =  0.00005;
                Beta  =  0.26;
            case 2 % koff- 
                A     =  2;
                Alpha = -0.46;
                B     =  0.00005;
                Beta  =  0.78;
            case 3 % koff+
                A     =  4;
                Alpha = -0.28;
                B     =  0.00005;
                Beta  =  0.95;
        end
   
end
