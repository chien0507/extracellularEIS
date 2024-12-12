# extracellularEIS

This github repo includes:
Matlab code for calculating RMSE (and other error metrics) on autolab data.

Change working folder and other filepaths for the lookup table code.

20241212_makeLookupTableGeneralizedMac.py and 20241212_makeLookupTableGeneralizedPC.py are python scripts to automatically generate the lookup tables the matlab code needs to run fitting, on either a mac or PC, accordingly.
These files require editing the experiment name and changing the filepaths in lines (18, 20, 26, 27).
The provided makeLookupTable code requires all data files have the text _freq in the filename (this prevents attempts to process dotfiles in the folder).
