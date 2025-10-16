function LamellipodiumModel_Bidone01(ModelParameters,SaveDir,SaveName)

            % Initialize model variables for Filaments, Adhesions, Ligands, and Filament/Adhesion/Ligand connections
            Membrane       = InitializeMembrane(ModelParameters);
            Filaments      = InitializeActinFilaments(ModelParameters,Membrane); 
            Adhesions      = InitializeAdhesions(ModelParameters,Membrane);
            Ligands        = InitializeLigands(ModelParameters);
            FALconnections = InitializeFALconnections;

            ModelParameters.SampleTimeStep = 1E-2;
            ModelParameters.SaveTimeStep = 1;
            
            % Pre-allocate space and create parameters to be saved after model completes (variable used for recording data is: SimData)
            TimeVec = 0:ModelParameters.TimeStep:ModelParameters.TotalSimulationTime;
            TV2     = 0:0.001:ModelParameters.TotalSimulationTime; % This time vector is for sampling only at 1ms intervals
            
            %N = length(TV2);
            N = ModelParameters.SaveTimeStep/ModelParameters.SampleTimeStep;
            nMono   = NaN(N,1);
            nAdhes  = NaN(N,1);
            MemVel  = NaN(N,1);

            SimData.ModelParameters = ModelParameters;
            SimData.TimeVector = TV2;
            SimData.MembranePosition = NaN(N,1);
            SimData.AdhesionData     = NaN(ModelParameters.AdhesionTotal,5,N);
            
            DATA = cell(N,1);
            MembranePrevious = Membrane;
            
            % Initialize parameters to control the "nth" frame to plot
            index = 0;  
            count = 0;
            SaveIndex = 0;
            WriteIndex = 0;
            %nth = 100; % plot every nth frame
            
            %disp('Starting model....')
            
    %% START Model 
            
            for t = TimeVec
                %tic
                % Main model calculations ----------------------------------------------------------------------------------
                [nMonomers, DeletedFilamentNames] = CountTotalMonomers(Filaments, Membrane, ModelParameters);
                [Filaments, Adhesions, Ligands, FALconnections] = PolymerizeDepolymerizeCapDeleteFilaments(nMonomers,Filaments,Adhesions,Ligands,Membrane,FALconnections,ModelParameters);
                 Filaments = BranchFilamentsInBranchWindowIfSelected(Filaments,ModelParameters,Membrane);
                [FALconnections, Adhesions, Ligands] = CreateFALconnections(FALconnections,Filaments,Adhesions,Ligands,ModelParameters);
                [Filaments, Membrane, Adhesions, Ligands, FALconnections, AdhesionTensions, Data] = CalculatePositionAfterAppliedForces(Filaments,Membrane,Adhesions,Ligands,FALconnections,ModelParameters, t);
                [Adhesions, Ligands, FALconnections] = ManageAdhesionsAndLigands(Filaments,Adhesions,Ligands,FALconnections,Membrane,ModelParameters);    
                
                % Calculate speed, mass, and add random filament if necessarry ---------------------------------------------
                if rem( round(t,10), ModelParameters.SampleTimeStep) == 0 % Only record in 1 ms intervals
                    index = index + 1;
                    SaveIndex = SaveIndex + 1;
                    Data.Timepoint  = t;
                    DATA{SaveIndex,1}   = Data; % Data contains retrograde flow values for all filaments at each timepoint
                    nMono(index,1)  = CountTotalMonomers(Filaments, Membrane, ModelParameters);
                    nAdhes(index,1) = length(find(Adhesions.ActiveStatus));
                    MemVel(index,1) = (Membrane.Nodes.Ycoords(1) - MembranePrevious.Nodes.Ycoords(1)) / ModelParameters.TimeStep; 
                    % Record membrane segment position, and Adhesion data --------------------------------------------------
                    SimData.MembranePosition(index,1) = Membrane.Nodes.Ycoords(1); 
                    SimData.AdhesionData(:,:,index)   = [Adhesions.XYpoints, Adhesions.ActiveStatus, AdhesionTensions, Adhesions.AttachedFilamentName];
                    SimData.nMonomers   = nMono;
                    SimData.nAdhesions  = nAdhes;
                    SimData.MemVelocity = MemVel;
                    %SimData.Data        = DATA; 
                    
                    if SaveIndex >= N
                        WriteIndex = WriteIndex + 1;
                        disp([num2str(SaveIndex),' of ',num2str(10000)])
                        TempSaveDir = fullfile(SaveDir, ['TempSaveFolder_', SaveName]);
                        TempSaveName = [SaveName(1,1:end-4), sprintf('_n%04d',WriteIndex),'.mat'];
                        if ~isfolder(TempSaveDir)
                            mkdir(TempSaveDir)
                        end
                        save(fullfile(TempSaveDir,TempSaveName),'DATA','-mat','-v7.3')
                        % Reset sampling variables --------------
                        DATA = cell(N,1);
                        SaveIndex = 0;
                        % ---------------------------------------
                    end

                end
                %-----------------------------------------------------------------------------------------------------------
                count = count + 1;
                MembranePrevious = Membrane;
                %toc
            end

             %% Concatenate temporary data files
            tempFiles = dir(fullfile(TempSaveDir, '*.mat'));
            AllData = {};
            for k = 1:length(tempFiles)
                tempFilePath = fullfile(tempFiles(k).folder, tempFiles(k).name);
                tempData = load(tempFilePath, 'DATA');
                AllData = [AllData; tempData.DATA];
            end
            SimData.Data = AllData;
    
    %% Delete temporary save folder
            rmdir(TempSaveDir, 's');

            
            
    %% END Model 



    try
        save(fullfile(SaveDir,SaveName),'SimData','-mat','-v7.3')
    catch
        disp(['File not saved: ',SaveName])
    end
    
end