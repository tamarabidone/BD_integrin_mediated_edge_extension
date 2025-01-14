# BD_integrin_mediated_edge_extension
BD model of integrin-mediated edge extension


The LamellipodiumModel simulates the dynamics of actin filament polymerization and membrane deformation in response to force interactions. It models the behavior of filaments, adhesions, and ligands on a membrane and tracks the membrane position, monomer counts, adhesion status, and filament dynamics over time.

The simulation includes several components:

Filament polymerization/depolymerization and branching.
Adhesion formation and management, including periodic boundary conditions.
Ligand interaction with actin filaments and adhesion sites.
Membrane displacement and velocity calculations.


Subfunctions Used in the Simulation
1. InitializeMembrane(ModelParameters)
Purpose: Initializes the membrane structure with the specified model parameters. Sets up membrane nodes and geometry.

2. InitializeActinFilaments(ModelParameters, Membrane)
Purpose: Initializes the actin filaments according to the model parameters and the membrane geometry. Sets up filament properties like length, position, and connectivity.

3. InitializeAdhesions(ModelParameters, Membrane)
Purpose: Initializes the adhesion sites on the membrane. Sets up the initial adhesion points where actin filaments can bind.

4. InitializeLigands(ModelParameters)
Purpose: Initializes the ligands involved in adhesion interactions.

Sets up ligand positions and properties.
5. InitializeFALconnections()
Purpose: Initializes the connections between filaments, adhesions, and ligands.

This structure keeps track of which filaments are connected to which adhesions and ligands.
6. PolymerizeDepolymerizeCapDeleteFilaments(CountTotalMonomers, Filaments, Adhesions, Ligands, Membrane, FALconnections, ModelParameters)
Purpose: Handles the polymerization, depolymerization, capping, and deletion of actin filaments.

Updates filament lengths and positions based on polymerization/depolymerization dynamics.
7. BranchFilamentsInBranchWindowIfSelected(Filaments, ModelParameters, Membrane)
Purpose: Manages the branching of actin filaments if the branching process is enabled in the model.

Updates the filament structure to reflect new branches.
8. CreateFALconnections(FALconnections, Filaments, Adhesions, Ligands, ModelParameters)
Purpose: Establishes new connections between actin filaments, adhesions, and ligands.

Updates the connections structure to reflect new interactions.
9. CalculatePositionAfterAppliedForces(Filaments, Membrane, Adhesions, Ligands, FALconnections, ModelParameters)
Purpose: Calculates the new positions of filaments, adhesions, and ligands after applying forces to the system.

Updates the positions of all model components based on force interactions.
10. ManageAdhesionsAndLigands(Filaments, Adhesions, Ligands, FALconnections, Membrane, ModelParameters)
Purpose: Handles the periodic boundary conditions for adhesions and ligands.

Ensures that the adhesions and ligands wrap around the boundaries of the model domain correctly, updating their positions accordingly.
11. CountTotalMonomers(Filaments)
Purpose: Counts the total number of monomers in the system at a given time step.

Iterates over all filaments and sums the number of monomers.
Data Structure: SimData
The SimData structure contains all the data recorded during the simulation. Key fields include:

ModelParameters: The input parameters for the simulation.
TimeVector: A time vector of simulation steps.
MembranePosition: The position of the membrane at each time step.
AdhesionData: Data for each adhesion, including position, status, tension, and the filaments attached.
nMonomers: The total number of monomers in the system at each time step.
nAdhesions: The number of active adhesions at each time step.
MemVelocity: The velocity of the membrane at each time step.
Data: A cell array containing detailed data for each time point.
Saving Results
The simulation results are saved as a .mat file, including all recorded data in SimData. This can be loaded into MATLAB for further analysis.

To save the results:

save(fullfile(SaveDir, SaveName), 'SimData', '-mat', '-v7.3')

How to Run the Model
To run the model, simply call the LamellipodiumModel_Bidone01 function with the appropriate inputs:

ModelParameters = struct(); % Define your model parameters
SaveDir = 'path/to/save/results';
SaveName = 'simulation_results.mat';
LamellipodiumModel_Bidone01(ModelParameters, SaveDir, SaveName);
Make sure the input parameters (ModelParameters, SaveDir, SaveName) are set properly before running the simulation.

MODEL OUTPUT ANALYSIS:
  Assuming model code is run and raw output files are created (i.e. assuming you have run the script: SimulationTaskList_001.m),
  in Analysis_scripts directory run analyis files in each folder followed by plotting code.
  Adjust raw data file read directory variables when necessary.
  
  
