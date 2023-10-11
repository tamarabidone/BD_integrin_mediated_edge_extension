The "Adhesionfrequencyanalysis.m" file is used to extract the frequency on the adhesions from the raw files obtained from running the model. 
Thus, the model output obtained needs to be accessible through a specific path specified in the Adhesionfrequencyanalysis file. 
Upon running the Adhesionfrequencyanalysis file, a .mat file will be generated containing the necessary information needed for the Frequencyplotting file. 
Running the Frequencyplotting file will generate figures that will be saved in the directory of your choice specified in the beginning of the code.
The Frequencyplotting files contains a segment of code that concatenates the data from 3 separate runs with the same conditions of substrate stiffness, and will then compute the average frequency active of integrins through the 30 seconds, and generate plots of those average frequencies for the different substrate stiffnesses.

