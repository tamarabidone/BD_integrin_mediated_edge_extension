clc
ModelParameters = InitializeModelParameters;  % Initialize Model paramaeters

% Special Conditions (or edits to ModelParameters) ==============================================================
    ModelParameters.TimeStep = 1e-4;
    ModelParameters.TotalSimulationTime = 50; % seconds
    ModelParameters.FilamentThermalMotionOn = false; 
    ModelParameters.CytoplasmViscosity = 1e5 * 0.0001; % Pascal * seconds
    ModelParameters.VerticalOffSet = -200;
    ModelParameters.StartingNumberOfFilaments = 32; % This is a ballpark value based on MaximumFilamentMass
    ModelParameters.IntegrinSpringConstant = 0.0001; % Change (substrate rigidity)
    ModelParameters.k_off_pointed = 7; % s^-1
    ModelParameters.k_branch = 2.2; % s^-1
    ModelParameters.FAL_connection_Distance = 10*2.75; % nm
    ModelParameters.MaximumFilamentMass = 4000; % monomers
    ModelParameters.MolecularClutch_PeakNumber = 1; % Catch bond (1 = WildType or 2 = Manganese)
% ===============================================================================================================
            %save ('SimData_0001pN_wt_1.mat', '-struct', 'SimData');
            Membrane       = InitializeMembrane(ModelParameters);
            Filaments      = InitializeActinFilaments(ModelParameters,Membrane); 
            Integrins      = InitializeIntegrins(ModelParameters,Membrane);
            Ligands        = InitializeLigands(ModelParameters);
            FILconnections = InitializeFILconnections;
            ShowPlot       = true;
            
            if ShowPlot  
                [FH,AH1,AH2,AH3] = SetUpPlottingFigureAndAxes; 
            end
            
            % Pre-allocate space and create parameters to be saved after model completes (variable used for recording data is: SimData)
            TimeVec = 0:ModelParameters.TimeStep:ModelParameters.TotalSimulationTime;
            TV2 = 0:0.001:ModelParameters.TotalSimulationTime; % This time vector is for sampling only at 1ms intervals
            nMono   = NaN(length(TV2),1);
            nInteg  = NaN(length(TV2),1);
            MemVel  = NaN(length(TV2),1);
            
            SimData.ModelParameters = ModelParameters;
            SimData.TimeVector = TV2;
            SimData.MembranePosition = NaN(length(TV2),1);
            SimData.IntegrinData     = NaN(ModelParameters.IntegrinTotal,5,length(TV2));
            
            DATA = cell(length(TV2),1);
            MembranePrevious = Membrane;
            
            % Initialize parameters to control the "nth" frame to plot
            index = 0;  
            count = 0;
            nth = 1000; % plot every nth frame
            
            disp('Starting model....')

%% START Movie

MovieDirectory = '/uufs/chpc.utah.edu/common/home/bidone-group3/Remi/Soft_WT_fig5';
v = VideoWriter(fullfile(MovieDirectory, 'Movie01.avi'), 'Motion JPEG AVI');
v.Quality = 95;
v.FrameRate = 10; %Frames per second
open(v)
MovieIdx = 0;   


%% START Model 

            for t = TimeVec
                % Main model calculations ----------------------------------------------------------------------------------
                [Filaments, Integrins, Ligands, FILconnections] = PolymerizeDepolymerizeCapDeleteFilaments(CountTotalMonomers(Filaments),Filaments,Integrins,Ligands,Membrane,FILconnections,ModelParameters);
                 Filaments = BranchFilamentsInBranchWindowIfSelected(Filaments,ModelParameters,Membrane);
                [FILconnections, Integrins, Ligands] = CreateFILconnections(FILconnections,Filaments,Integrins,Ligands,ModelParameters);
                [Filaments, Membrane, Integrins, Ligands, FILconnections, IntegrinTensions, Data] = CalculatePositionAfterAppliedForces(Filaments,Membrane,Integrins,Ligands,FILconnections,ModelParameters);
                [Integrins, Ligands, FALconnections] = ManageIntegrinsAndLigands(Filaments,Integrins,Ligands,FILconnections,Membrane,ModelParameters);    
                
                % Calculate speed, mass, and add random filament if necessarry ---------------------------------------------
                if rem( round(t,10),0.1) == 0 % Only record in 1ms intervals
                    disp(t)
                    index = index + 1;
                    Data.Timepoint  = t;
                    DATA{index,1}   = Data; % Data contains retrograde flow values for all filaments at each timepoint
                    nMono(index,1)  = CountTotalMonomers(Filaments);
                    nAdhes(index,1) = length(find(Integrins.ActiveStatus));
                    MemVel(index,1) = (Membrane.Nodes.Ycoords(1) - MembranePrevious.Nodes.Ycoords(1)) / ModelParameters.TimeStep; 
                    
                    % Record membrane segment position, and Adhesion data ------------------
                    SimData.MembranePosition(index,1) = Membrane.Nodes.Ycoords(1); 
                    SimData.IntegrinData(:,:,index)   = [Integrins.XYpoints, Integrins.ActiveStatus, IntegrinTensions, Integrins.AttachedFilamentName];
                    SimData.nMonomers   = nMono;
                    SimData.nAdhesions  = nAdhes;
                    SimData.MemVelocity = MemVel;
                    SimData.Data        = DATA; 
                end
                
                %[Filaments] = AddRandomFilaments(Filaments,Membrane,ModelParameters,nMono(index,1)); % If neccesary add a new filament
                
                % Create plot ----------------------------------------------------------------------------------------------
                if ShowPlot
                    MovieIdx = MovieIdx + 1;
                    %[FILconnections,count] = PlotFilamentsAndMembrane(nth,count,Filaments,Membrane,Integrins,Ligands,FILconnections,FH,AH1,AH2,AH3,t,nMono,nInteg,MemVel,index,TV2,ModelParameters);
                    [FILconnections,count] = PlotFilamentsAndMembraneMovie01(nth,count,Filaments,Membrane,Integrins,Ligands,FILconnections,FH,AH1,AH2,AH3,t,nMono,nInteg,MemVel,index,TV2,ModelParameters);
                    Fimage = getframe(FH);
                    writeVideo(v,Fimage.cdata)
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
