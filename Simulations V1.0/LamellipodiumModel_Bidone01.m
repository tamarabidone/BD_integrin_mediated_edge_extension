function LamellipodiumModel_Bidone01(ModelParameters,SaveDir,SaveName)

            % Initialize model variables for Filaments, Adhesions, Ligands, and Filament/Adhesion/Ligand connections
            Membrane       = InitializeMembrane(ModelParameters);
            Filaments      = InitializeActinFilaments(ModelParameters,Membrane); 
            Adhesions      = InitializeAdhesions(ModelParameters,Membrane);
            Ligands        = InitializeLigands(ModelParameters);
            FALconnections = InitializeFALconnections;
%             ShowPlot       = true;
            
%             if ShowPlot  
%                 [FH,AH1,AH2,AH3] = SetUpPlottingFigureAndAxes; 
%             end
            
            % Pre-allocate space and create parameters to be saved after model completes (variable used for recording data is: SimData)
            TimeVec = 0:ModelParameters.TimeStep:ModelParameters.TotalSimulationTime;
            TV2     = 0:0.001:ModelParameters.TotalSimulationTime; % This time vector is for sampling only at 1ms intervals
            nMono   = NaN(length(TV2),1);
            nAdhes  = NaN(length(TV2),1);
            MemVel  = NaN(length(TV2),1);
            
            SimData.ModelParameters = ModelParameters;
            SimData.TimeVector = TV2;
            SimData.MembranePosition = NaN(length(TV2),1);
            SimData.AdhesionData     = NaN(ModelParameters.AdhesionTotal,5,length(TV2));
            
            DATA = cell(length(TV2),1);
            MembranePrevious = Membrane;
            
            % Initialize parameters to control the "nth" frame to plot
            index = 0;  
            count = 0;
            %nth = 100; % plot every nth frame
            
            %disp('Starting model....')
            
    %% START Model 
            tic
            for t = TimeVec
                % Main model calculations ----------------------------------------------------------------------------------
                [Filaments, Adhesions, Ligands, FALconnections] = PolymerizeDepolymerizeCapDeleteFilaments(CountTotalMonomers(Filaments),Filaments,Adhesions,Ligands,Membrane,FALconnections,ModelParameters);
                 Filaments = BranchFilamentsInBranchWindowIfSelected(Filaments,ModelParameters,Membrane);
                [FALconnections, Adhesions, Ligands] = CreateFALconnections(FALconnections,Filaments,Adhesions,Ligands,ModelParameters);
                [Filaments, Membrane, Adhesions, Ligands, FALconnections, AdhesionTensions, Data] = CalculatePositionAfterAppliedForces(Filaments,Membrane,Adhesions,Ligands,FALconnections,ModelParameters);
                [Adhesions, Ligands, FALconnections] = ManageAdhesionsAndLigands(Filaments,Adhesions,Ligands,FALconnections,Membrane,ModelParameters);    
                
                % Calculate speed, mass, and add random filament if necessarry ---------------------------------------------
                if rem( round(t,10), 0.1) == 0 % Only record in 1 ms intervals
                    index = index + 1;
                    Data.Timepoint  = t;
                    DATA{index,1}   = Data; % Data contains retrograde flow values for all filaments at each timepoint
                    nMono(index,1)  = CountTotalMonomers(Filaments);
                    nAdhes(index,1) = length(find(Adhesions.ActiveStatus));
                    MemVel(index,1) = (Membrane.Nodes.Ycoords(1) - MembranePrevious.Nodes.Ycoords(1)) / ModelParameters.TimeStep; 
                    % Record membrane segment position, and Adhesion data --------------------------------------------------
                    SimData.MembranePosition(index,1) = Membrane.Nodes.Ycoords(1); 
                    SimData.AdhesionData(:,:,index)   = [Adhesions.XYpoints, Adhesions.ActiveStatus, AdhesionTensions, Adhesions.AttachedFilamentName];
                    SimData.nMonomers   = nMono;
                    SimData.nAdhesions  = nAdhes;
                    SimData.MemVelocity = MemVel;
                    SimData.Data        = DATA; 
                end
                %-----------------------------------------------------------------------------------------------------------
                count = count + 1;
                MembranePrevious = Membrane;
            end
            toc
            
    %% END Model 


    try
        save(fullfile(SaveDir,SaveName),'SimData','-mat','-v7.3')
    catch
        disp(['File not saved: ',SaveName])
    end
    
end

