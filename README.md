# extracellularEIS

This github repo includes:
Matlab code for calculating RMSE (and other error metrics) on autolab data.

Change working folder and other filepaths for the lookup table code.

calcRMSE.m - change line 2 and 5 to reference the folder with the raw data and the summary table, correspondingly, and run the matlab code to output a file with the resnorm and RMSE as 2 columns, comma delimited.
dataLocation = '/Users/Athena/GaTech Dropbox/Athena Chien/WPI EIS Project/Biological Sample Data/9-23-24 Testing/20240923_WPIbioExp2/';
SummaryTable = readtable('/Users/Athena/GaTech Dropbox/Athena Chien/PBL Hanna-Athena Files/dataforBPSabstract/newSummary/20240923_WPIBIoExp2_newSettings_summary_table.csv');


20241212_makeLookupTableGeneralizedMac.py and 20241212_makeLookupTableGeneralizedPC.py are python scripts to automatically generate the lookup tables the matlab code needs to run fitting, on either a mac or PC, accordingly.
These files require editing the experiment name and changing the filepaths in lines (18, 20, 26, 27).
The provided makeLookupTable code requires all data files have the text _freq in the filename (this prevents attempts to process dotfiles in the folder).
