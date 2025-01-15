# extracellularEIS

This github repo includes:
- file structure and code for loading in lookup tables and raw EIS data, and outputting TER, TEC, membrane ratio, and resnorm in a summary table for each raw EIS sweep.
- Matlab code for calculating mean absolute error (MAE) on autolab data (calcMAE.m)
- Python code for mac and PC automatically generate lookup tables (20241212_makeLookupTableGeneralizedMac.py and 20241212_makeLookupTableGeneralizedPC.py)

Using calcMAE.m to calculate the mean absolute error of an EIS sweep an RCRC fit:
calcMAE.m - change line 2 and 5 to reference the folder with the raw data and the summary table, correspondingly, and run the matlab code to output a file with the resnorm and RMSE as 2 columns, comma delimited.
dataLocation = '/Users/Athena/GaTech Dropbox/Athena Chien/WPI EIS Project/Biological Sample Data/9-23-24 Testing/20240923_WPIbioExp2/';
SummaryTable = readtable('/Users/Athena/GaTech Dropbox/Athena Chien/PBL Hanna-Athena Files/dataforBPSabstract/newSummary/20240923_WPIBIoExp2_newSettings_summary_table.csv');

20241212_makeLookupTableGeneralizedMac.py and 20241212_makeLookupTableGeneralizedPC.py are python scripts to automatically generate the lookup tables the matlab code needs to run fitting, on either a mac or PC, accordingly.
These files require editing the experiment name and changing the filepaths in lines (18, 20, 26, 27).
The provided makeLookupTable code requires all data files have the text _freq in the filename (this prevents attempts to process dotfiles in the folder).
- Replace "workingFolder" with the file and path to direct to the desired folder with the raw data files. Replace "summaryTable" name with the desired name of the summary table. Replace "raw directory" with the path to the "raw data" folder.
- This code automatically detects any files with "modelCell" to set cross-sectional area to 1, and all remaining to the cross-sectional area of a Corning 3460 Transwell, 1.12 $cm^{2}$.
- This code copies the raw data files and the final summary table to the folder with the fit code. Additional columns can be filled out manually for additional notes, and then should be manually copied to the fit folder under "lookup table".
