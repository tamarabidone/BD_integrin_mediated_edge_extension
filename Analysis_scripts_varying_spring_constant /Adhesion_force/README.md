The "Adhesionforceanalysis.m" file is used to extract the force on the adhesions from the raw files obtained from running the model. 
Thus, the model output obtained needs to be accessible through a specific path in the Adhesionforceanalysis file.
Upon running the Adhesionforceanalysis file, a .mat file will be generated containing the necessary information needed for the Forceplotting file. 
This .mat file needs to be in the same folder as the Forceplotting file. 
Running the Forceplotting file will generate figures that will be saved in the same directory as the Adhesionforceanalysis and Forceplotting files. The Forceplotting files contains two segments of code. The first segment will generate a histogram plot of the distribution of forces on adhesions for 3 different spring constants. The other segment of code will generate a figure plotting the mean and standard deviation of force on the adhesions for all 5 different spring constants. 
