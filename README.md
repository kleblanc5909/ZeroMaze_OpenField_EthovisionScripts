# ZeroMaze_OpenField_EthovisionScripts
Python scripts for analyzing zero maze and open field data from Noldus Ethovision raw data files.

These scripts use pandas dataframes to import raw data files exported from Ethovision. 

Input Files:These are excel documents with a header of approximately 40 rows containing the independent experimental variables. Following the header is the beginning of the data, starting with a row of the dependent variables. The data block is approximately 25 columns wide and ranges from 600 rows to over 200,000 rows depending on the length of the video recording. 

Output Files: Excel files consisting of a row of column titles with the selected output variable names, both dependent and independent variables from the header, followed by separate rows of the data for each subject. In general, these scripts analyze the number of movements, the average duration of movements, the speed while moving, average velocity, the percent of time in the open or center, and the number of movements into the open or center.

The different scripts depend on the experimental task (zero maze or open field), and the variables that are analyzed. Optogenetic experiments take into account the LED status, DREADD experiments include variables for drug. 
