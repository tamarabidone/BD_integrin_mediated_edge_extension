# BD_integrin_mediated_edge_extension
BD model of integrin-mediated edge extension

RUNNING MODEL:
  To run a visual version of the model, run the script: LamellipodiumModel_Bidone01_NonFunction.m.
  To run multiple runs while varying parameters, run or create a script like SimulationTaskList_001.m. 
  Be sure to update set SaveDirectory on line 4.
  Also, engage MATLAB parallel processing before running, or comment out line 25 and relace with line 26.
  Upon launching the code (visual or not), the code runs through the LamellipodiumModel_Bidone01.m file, which contains all of the sub 
  functions located within the Simulation V1.0 folder.
  At first, the code cycles through a number of functions to establish the properties of the elements of the model before running the 
  simulation:
  InitializeMembrane.m : This code initializes a membrane structure by defining its node coordinates along the x-axis (from 0 to 500 nm) and 
  setting all y-coordinates to zero, based on the provided model parameters.
  InitializeActinFilaments.m : This code initializes a set of actin filaments by creating their starting positions, orientations, lengths, and 
  other properties based on model parameters, while also randomly distributing them in space and applying offsets to simulate biological 
  variability.
  InitializeAdhesions.m : This code initializes adhesion points in a model by randomly generating their positions within a defined model 
  domain (500 x 500 nm), assigning default properties like active status, orientation, and attachment information, and defining the boundaries 
  of the adhesion region based on the membrane and model depth.
  InitializeLigands.m : This code initializes ligand molecules by randomly generating their positions within a defined model 
  domain (500 x 500 nm) and setting their default attachment properties to indicate that they are not yet connected to any filaments or 
  adhesions.
  InitializeFALconnections.m : This code initializes the FALconnections structure, which represents connections between filaments, integrins, 
  and ligands, by creating empty arrays for indices of integrins, ligands, filaments, and monomers, starting with no connections 
  established.

  Then, the code runs through a set of functions that are repeated after every timestep:

  PolymerizeDepolymerizeCapDeleteFilaments.m : manages the dynamic processes of actin filament polymerization, depolymerization, capping, and 
    deletion within a simulation framework. Here's a breakdown of its functionality:
  
    The function modifies the filament structure over time, considering:
    Polymerization: Growth of filaments at their barbed ends. Determines if a filament grows based on a PolymCoeff, Brownian ratchet effects,      and monomer availability.
    The CalculatePolymerizationCoefficient function computes the polymerization probability using mechanical models of membrane force 
    interactions and filament flexibility.
    Depolymerization: Shrinkage of filaments at their pointed ends. Filaments shrink at their pointed ends based on a probability defined by 
    k_off_pointed.
    Capping: Temporary halting of polymerization. Filaments are capped or uncapped based on random probability and rates (k_cap and k_uncap).
    Deletion: Removal of filaments that become too short or are otherwise flagged for removal. Filaments are deleted if their length falls     
    below 
    a threshold (e.g., 2 monomers). Deleted filaments are removed from all related structures, including: Filament-specific arrays (e.g., 
    MonomerIndices, XYCoords). Adhesion and ligand connections using the DeleteFALconnection function.
    
    Helper Functions
    CalculatePolymerizationCoefficient.m = This function calculates the coefficient affecting filament polymerization rates using: Distance 
    from 
    the membrane (y), Filament orientation (theta), Persistence length, thermal energy, and geometric considerations.


  BranchFilamentsInBranchWindowIfSelected.m : Provides functionality for simulating the polymerization and depolymerization of actin filaments     and for branching filaments under specific conditions. Below is a breakdown of the different key components and their roles:
    
    Polymerization (Barbed End): For each filament, a random process is simulated to cap or uncaps the filament. If uncapped, the filament may     undergo polymerization at its barbed end based on the available monomer concentration and the specific polymerization coefficient.
    The PolymerizationCoefficient is calculated based on various factors, including membrane distance and filament geometry.
    
    Depolymerization (Pointed End):The filament may also undergo depolymerization at its pointed end, where the oldest monomers are removed.
    If the filament is small (2 or fewer monomers), it is removed entirely. The DeleteFALconnection function is used to delete any filament-       adhesion-ligand (FAL) connections associated with the depolymerized monomer.
    
    Filament Deletion: After polymerization or depolymerization, filaments that meet specific criteria (e.g., reaching a certain length or         undergoing certain events) are removed from the system.
    
    
    Branching Window: The code checks if the filament's tip is within a defined "branching window" on the membrane.
    If the filament is within this window and a branch is not already present at the tip, a branching event can occur based on a probability       determined by the branching rate (k_branch).
    
    Branching Criteria: A minimum separation between branches is enforced, meaning branches cannot be created too close to each other.
    The angle of the branch is randomly determined based on a Gaussian distribution with a standard deviation (branchAngleSTD) and a mean          branching angle (branchAngle).
    
    Branch Creation: If a branch is selected, a new filament is created starting from the filament tip in a direction determined by the parent     filament’s orientation and the random angle. New properties (e.g., name, monomer indices, coordinates, and parent-child relationships) are     set for the new branch.
    
    Utility Functions: IsFilamentWithinTheBranchingWindowOfAMembraneSegment: Checks if the filament's tip is within the branching window of        any membrane segment, based on the filament's coordinates and the geometry of the membrane.
    
    CalculatePolymerizationCoefficient: Computes the polymerization coefficient considering various factors, including the filament's distance     from the membrane, its length, and its orientation.
    
    Key Concepts and Parameters:
    Filament Geometry: Filaments are described by their coordinates, monomer indices, and orientation vectors.
    Model Parameters: Includes key parameters such as monomer length, membrane properties, branching rates, and polymerization coefficients.
    Membrane Geometry: The membrane is described by a set of nodes and segments, and the branching window is a region near the membrane where      branching is allowed.
    This structure provides a simulation of filament dynamics, allowing the study of polymerization, depolymerization, and branching processes     in a simulated environment with adjustable parameters for the system.



MODEL OUTPUT ANALYSIS:
  Assuming model code is run and raw output files are created (i.e. assuming you have run the script: SimulationTaskList_001.m),
  in Analysis_scripts directory run analyis files in each folder followed by plotting code.
  Adjust raw data file read directory variables when necessary.
  
  
