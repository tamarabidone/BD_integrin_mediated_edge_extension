clc
ModelParameters = InitializeModelParameters;  % Initialize Model paramaeters

% Special Conditions (or edits to ModelParameters) ==============================================================
    ModelParameters.TimeStep = 1e-4;
    ModelParameters.TotalSimulationTime = 10; % seconds
    ModelParameters.FilamentThermalMotionOn = false; 
    ModelParameters.CytoplasmViscosity = 1e5 * 0.0001; % Pascal * seconds
    ModelParameters.VerticalOffSet = -200;
    ModelParameters.StartingNumberOfFilaments = 32; % This is a ballpark value based on MaximumFilamentMass
    ModelParameters.AdhesionSpringConstant = 0.01; % Change (substrate rigidity)
    ModelParameters.k_off_pointed = 7; % s^-1
    ModelParameters.k_branch = 2.2; % s^-1
    ModelParameters.FAL_connection_Distance = 10*2.75; % nm
    ModelParameters.MaximumFilamentMass = 4000; % monomers
    ModelParameters.MolecularClutch_PeakNumber = 1; % Catch bond (1 = WildType or 2 = Manganese)
    ModelParameters.MyosinFraction = 1;
    ModelParameters.ApplyMyosin = false;
    ModelParameters.ReducedFlow = false;
    ModelParameters.ExponentialLifetime = false;
% ===============================================================================================================
            %save ('SimData_0001pN_wt_1.mat', '-struct', 'SimData');
            Membrane       = InitializeMembrane(ModelParameters);
            Filaments      = InitializeActinFilaments(ModelParameters,Membrane); 
            Adhesions      = InitializeAdhesions(ModelParameters,Membrane);
            Ligands        = InitializeLigands(ModelParameters);
            FALconnections = InitializeFALconnections;
            ShowPlot       = true;
            
            if ShowPlot  
                [FH,AH1,AH2,AH3] = SetUpPlottingFigureAndAxes; 
            end
            
            % Pre-allocate space and create parameters to be saved after model completes (variable used for recording data is: SimData)
            TimeVec = 0:ModelParameters.TimeStep:ModelParameters.TotalSimulationTime;
            TV2 = 0:0.001:ModelParameters.TotalSimulationTime; % This time vector is for sampling only at 1ms intervals
            nMono   = NaN(length(TimeVec),1);
            nAdhes  = NaN(length(TimeVec),1);
            MemVel  = NaN(length(TimeVec),1);
            
            SimData.ModelParameters = ModelParameters;
            SimData.TimeVector = TimeVec;
            SimData.MembranePosition = NaN(length(TimeVec),1);
            SimData.AdhesionData     = NaN(ModelParameters.AdhesionTotal,5,length(TV2));
            
            DATA = cell(length(TV2),1);
            MembranePrevious = Membrane;
            
            % Initialize parameters to control the "nth" frame to plot
            index = 0; 
            index2 = 0;
            count = 0;
            nth = 1000; % plot every nth frame
            MeanTensions = [];
            
            disp('Starting model....')

%% START Movie

% MovieDirectory = '/Users/remisondaz/Desktop/MATLAB/Simulations_V1_DifferentClutch';
% v = VideoWriter(fullfile(MovieDirectory, 'Movie01.avi'), 'Motion JPEG AVI');
% v.Quality = 95;
% v.FrameRate = 10; %Frames per second
% open(v)
% MovieIdx = 0;   


%% START Model 

            for t = TimeVec
                index = index + 1;
                % Main model calculations ----------------------------------------------------------------------------------
                [nMonomers, DeletedFilamentNames]  = CountTotalMonomers(Filaments, Membrane, ModelParameters);
                [Filaments, Adhesions, Ligands, FALconnections] = PolymerizeDepolymerizeCapDeleteFilaments(nMonomers, Filaments,Adhesions,Ligands,Membrane,FALconnections,ModelParameters);
                 Filaments = BranchFilamentsInBranchWindowIfSelected(Filaments,ModelParameters,Membrane);
                [FALconnections, Adhesions, Ligands] = CreateFALconnections(FALconnections,Filaments,Adhesions,Ligands,ModelParameters);
                [Filaments, Membrane, Adhesions, Ligands, FALconnections, AdhesionTensions, Data] = CalculatePositionAfterAppliedForces(Filaments,Membrane,Adhesions,Ligands,FALconnections,ModelParameters, t);
                [Adhesions, Ligands, FALconnections] = ManageAdhesionsAndLigands(Filaments,Adhesions,Ligands,FALconnections,Membrane,ModelParameters); 

                nMono(index,1)  = nMonomers;
                nAdhes(index,1) = length(find(Adhesions.ActiveStatus));
                MemVel(index,1) = (Membrane.Nodes.Ycoords(1) - MembranePrevious.Nodes.Ycoords(1)) / ModelParameters.TimeStep; 
                disp(MemVel(index,1))
                
                % Calculate speed, mass, and add random filament if necessarry ---------------------------------------------
                if rem( round(t,10),0.1) == 0 % Only record in 1ms intervals
                    disp(t)
                    index2 = index2 + 1;
                    Data.Timepoint  = t;
                    DATA{index2,1}   = Data; % Data contains retrograde flow values for all filaments at each timepoint
                     
                    % Record membrane segment position, and Adhesion data ------------------
                    SimData.MembranePosition(index2,1) = Membrane.Nodes.Ycoords(1); 
                    SimData.AdhesionData(:,:,index2)   = [Adhesions.XYpoints, Adhesions.ActiveStatus, AdhesionTensions, Adhesions.AttachedFilamentName];
                    SimData.nMonomers   = nMono;
                    SimData.nAdhesions  = nAdhes;
                    SimData.MemVelocity = MemVel;
                    SimData.Data        = DATA; 
                end
                
                %[Filaments] = AddRandomFilaments(Filaments,Membrane,ModelParameters,nMono(index,1)); % If neccesary add a new filament
                
                % Create plot ----------------------------------------------------------------------------------------------
                if ShowPlot
                    %MovieIdx = MovieIdx + 1;
                    idx = find(AdhesionTensions > 0.001);
                    MeanTensions = [MeanTensions; mean(AdhesionTensions(idx))];
                    [FALconnections,count] = PlotFilamentsAndMembrane(TimeVec,nth,count,Filaments,Membrane,Adhesions,Ligands,FALconnections,FH,AH1,AH2,AH3,t,nMono,nAdhes,MemVel,index,TV2,ModelParameters, AdhesionTensions,MeanTensions);
                    %[FALconnections,count] = PlotFilamentsAndMembraneMovie01(nth,count,Filaments,Membrane,Adhesions,Ligands,FALconnections,FH,AH1,AH2,AH3,t,nMono,nAdhes,MemVel,index,TV2,ModelParameters);
                    %Fimage = getframe(FH);
                    %writeVideo(v,Fimage.cdata)
                    %imwrite( Fimage.cdata, fullfile(MovieDirectory,['Frame_',sprintf('%06d',MovieIdx),'.tif']) )
                end
                
                if round(t,10) == 0.5
                    disp('pause')
                elseif round(t,10) == 10
                    disp('pause')
                elseif round(t,10) == 20
                    disp('pause')
                end

                count = count + 1;
                MembranePrevious = Membrane;
                drawnow
            end
            
%% END Model 


% close(v)



         %%save(fullfile(SaveDir,SaveName),'SimData','-mat','-v7.3')
    %% catch
       %%  disp(['File not saved: ',SaveName])
    %% end


