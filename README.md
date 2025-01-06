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

    CreateFALconnections.m : establishes connections between filaments, integrins, and ligands based on spatial proximity and activation          probabilities. Here's a detailed breakdown of its key components:

      These connections, called FAL connections (Filament-Adhesion-Ligand connections), represent interactions that may influence the     
      dynamics of filament growth and cell adhesion.
      
      Input Arguments:
      FALconnections: Contains existing connections between filaments, adhesions, and ligands.
      Filaments: Contains information about the filaments, including their coordinates (XYCoords), names, and monomer indices.
      Adhesions: Contains information about the adhesions, including their coordinates (XYpoints).
      Ligands: Contains information about the ligands, including their coordinates (XYpoints).
      ModelParameters: Contains model parameters such as the connection distance (FAL_connection_Distance), monomer length, and adhesion 
      activation rate.
      
      Loop Over Filaments: For each filament, the code calculates the distances between the filament's monomers and the adhesions, as well as 
      between the adhesions and ligands.
      XY_filament: Coordinates of the filament's monomers.
      XY_adhesions: Coordinates of adhesions that are not yet connected to a filament or ligand.
      XY_ligands: Coordinates of ligands that are not yet connected to an adhesion.
      
      Distance Calculation: The function calculates the Euclidean distance from each adhesion to the filament monomers (Distance_AF) and from 
      each adhesion to the ligands (Distance_AL).
      
      Find Minimum Distances: For each adhesion, it finds the nearest filament monomer and the nearest ligand. This is done by identifying 
      the minimum distance between the adhesion and all monomers, as well as the minimum distance between the adhesion and all ligands.
      
      Filter Connections Based on Distance: Only connections where the distance is less than or equal to the predefined connection distance 
      (D) are considered. This ensures that only the adhesions within range of a filament and ligand can potentially connect.
      
      Probability of Connection: For each valid connection (adhesion-filament-ligand triplet), the probability of the connection forming is 
      calculated using an exponential decay based on the adhesion activation rate and the simulation time step.
      p_on: The probability that a connection will occur at the current time step.
      
      Create FAL Connections: If the connection is activated (based on the calculated probability), a new FAL connection is established by 
      adding the connection details to the FALconnections structure.
      Additionally, the adhesion and ligand structures are updated to reflect the new connection:
      The adhesion is marked as attached to the filament, and its position is updated to the position of the ligand.
      The ligand is marked as attached to the filament and the specific adhesion.
      Output: The updated FALconnections, Adhesions, and Ligands structures are returned, reflecting the newly established connections.

    CalculatePositionAfterAppliedForces.m : Filament Motion Section:

      Filament structures are calculated by finding the tips of each filament and calculating the forces acting on them (due to adhesion, 
      membrane, and Brownian motion).
      Forces due to connections between the filaments and adhesions are computed using CalculateForceDueToFALconnections.
      The motion of the filaments is updated according to these forces, using a model of random fluid motion (Brownian motion).
      Membrane Motion Section:
      
      This section calculates the force exerted by the filaments on the membrane. The membrane's movement is then computed based on the total 
      force.
      
      Adhesion Motion Section: Adhesions that are not yet activated undergo Brownian motion, where their positions are updated based on 
      random forces.
      
      Sub-functions: 
      TotalForceExertedByMembraneOnFilamentStructure: Calculates the normal force exerted by the membrane on the filament structure.
      CalculateForceFromFilamentsHittingMembraneSegment: Calculates the total force from all filaments interacting with the membrane.
      CalculateForceDueToFALconnections: Calculates the forces between filaments and adhesions (via Focal Adhesion-Ligand (FAL) connections), 
      including molecular clutch deactivation based on tension.
      Parameters Involved:
      ModelParameters: Contains constants such as CytoplasmViscosity, kT (Boltzmann constant), PersistenceLength, and adhesion-related 
      parameters.
      Filaments: Contains the information about filament names, positions, unit vectors, and lengths.
      Adhesions: Stores information about adhesion sites and their active status.
      Ligands: Stores ligand-related data for adhesion interactions.
      FALconnections: Stores information about the connections between filaments and adhesions.
      Membrane: Contains membrane properties, including its position and boundary conditions.

      Force Calculations:
      Forces acting on filaments from adhesions (FAx, FAy).
      Forces exerted by the membrane on filaments (Fm).
      Forces due to Brownian motion are included in the filament's movement.
      
      Thermal Motion:
      If enabled, thermal motion is added to the filament motion using random Gaussian noise scaled by the diffusion coefficient.
      
      Molecular Clutch: 
      The code implements a molecular clutch mechanism, where adhesion connections can break based on a tension-dependent rate (k_off), which 
      is modeled using exponential decay.
      
      Boundary Conditions: The code checks if any filament tips go out of bounds (left or right of the membrane) and applies periodic 
      boundary conditions by shifting the filaments to the opposite side.

  `ManageAdhesionsAndLigands.m : Performs periodic boundary conditions for adhesion and ligand points. The periodic boundary conditions 
     ensure that if the coordinates of adhesions or ligands fall outside the defined region of the membrane (based on RegionNodes), they are 
     wrapped around to the opposite side, effectively creating a continuous system.

    Key Operations:
      Adhesion Handling: For each boundary (top, bottom, left, right), if adhesion points move outside the defined region, their positions 
      are adjusted by adding or subtracting ModelDepth or MembraneWidth (depending on the boundary). Then, the function checks whether the 
      adhesion point is connected to any filament or ligand. If so, the connection is removed by calling the DeleteFALconnection function.

    Ligand Handling: Similar to the adhesion section, ligands are handled with periodic boundary conditions, and their connections to 
      filaments or adhesions are removed if they go out of bounds.

    Helper Function (DeleteFALconnection): This function is used to delete connections between filaments, adhesions, and ligands whenever an 
      adhesion or ligand moves across the boundary.

  CountTotalMonomers.m : Counts the total number of monomers across all filaments. Here's a breakdown of the function:


MODEL OUTPUT ANALYSIS:
  Assuming model code is run and raw output files are created (i.e. assuming you have run the script: SimulationTaskList_001.m),
  in Analysis_scripts directory run analyis files in each folder followed by plotting code.
  Adjust raw data file read directory variables when necessary.
  
  
