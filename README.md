# extracellularEIS Repo

This github repo includes:
- file structure and code for fitting raw EIS data with the RCRC model (as specified in Protocols paper: XX)
  - lookup table must be loaded into ``lookup table`` folder, raw EIS datafiles into ``raw data`` folder
  - then when the matlab code ``NOVA_batch_20240125.m`` is run, TER, TEC, membrane ratio, and resnorm are output as a summary table for each raw EIS sweep in ``FIT`` folder
- Matlab code for calculating mean absolute error (MAE) on EIS data data (calcMAE.m)
- Python code for mac and PC to automatically generate lookup tables (20241212_makeLookupTableGeneralizedMac.py and 20241212_makeLookupTableGeneralizedPC.py)

## Fitting raw EIS data with RCRC model:
1. Open the the downloaded ``extracellularEIS-main`` folder from Github. Copy and paste the template lookup table in the ``lookup table`` folder titled ``Template Lookup Table.xlsx`` and rename the copy with your experiment name (e.g. ``20240706Exp1 Lookup Table.xlsx``).
2. In the renamed lookup table, copy the name(s) of the raw data text file(s) (e.g. data1\_freq.txt) you would like to fit into Column B, under ``plateID.``
3. Type the cross-sectional area of the sample for each raw file in Column K, under ``measArea`` (1.12 $\mathrm{cm^2}$ for Corning 3460 cell culture inserts).
4. Number the files in column A with integer values starting at 1.
5. Save and close the lookup table.
6. Copy all the raw impedance data files output by the NOVA software into the ``raw data" folder in the ``extracellularEIS-main`` folder.
7. Open the ``NOVA\_batch\_20240125.m`` file in Matlab and change the file name and path in line 51 to point to the desired lookup table (e.g. ``20240706Exp1 Lookup Table.xlsx``).
8. In the Matlab software, click the ``Editor`` menu in the top bar, and then the ``Run`` button to begin fitting. The software runs multiple fits in parallel, and should show ``parallel processing`` in the lower left corner. A figure will pop up and update with the raw data Nyquist plots and the overlaid fit line as the fit code is running through each sweep.
9. When completed, a file with the format ``\_summary\_table.csv`` will be generated in the ``FIT`` folder. 

## Using calcMAE.m to calculate the mean absolute error (MAE) of an EIS sweep an RCRC fit:
- change line 2 and 5 to reference the folder with the raw data and the summary table, correspondingly, and run the matlab code to output a file with the resnorm and MAE as 2 columns, comma delimited:
dataLocation = '/Users/Athena/GaTech Dropbox/Athena Chien/WPI EIS Project/Biological Sample Data/9-23-24 Testing/20240923_WPIbioExp2/';
SummaryTable = readtable('/Users/Athena/GaTech Dropbox/Athena Chien/PBL Hanna-Athena Files/dataforBPSabstract/newSummary/20240923_WPIBIoExp2_newSettings_summary_table.csv');

## Automatically generating lookup tables:
20241212_makeLookupTableGeneralizedMac.py and 20241212_makeLookupTableGeneralizedPC.py are python scripts to automatically generate the lookup tables the matlab code needs to run fitting, on either a mac or PC, accordingly.
These files require editing the experiment name and changing the filepaths in lines (18, 20, 26, 27).
The provided makeLookupTable code requires all data files have the text _freq in the filename (this prevents attempts to process dotfiles in the folder).
- Replace "workingFolder" with the file and path to direct to the desired folder with the raw data files. Replace "summaryTable" name with the desired name of the summary table. Replace "raw directory" with the path to the "raw data" folder.
- This code automatically detects any files with "modelCell" to set cross-sectional area to 1, and all remaining to the cross-sectional area of a Corning 3460 Transwell, 1.12 $cm^{2}$.
- This code copies the raw data files and the final summary table to the folder with the fit code. Additional columns can be filled out manually for additional notes, and then should be manually copied to the fit folder under "lookup table".
