
# run on python 3.9.7
# Code written by Athena Chien to automatically make lookup table from raw EIS files in a folder on a PC 
# (see 20241212_makeLookupTableGeneralizedMac.py if you are using a mac)
# reads any file with _freq in it as a raw data file to put in the lookup table
# sets cross-sectional area to 1 if filename has modelCell in it
# otherwise sets to 1.12 cm^2 (area for corning 3460 cell culture insert)

import pandas as pd
import os
import shutil
import csv

# initialize dataframe
fileList = pd.DataFrame(columns = ["plateID", "plateWellIDs", "treatment", "treatDate",	"treatTime", "cultMediaID", "cultMediaDate", "cultApVol", "cultBaVol", "measArea", "deviceID", "devApVol", "devBaVol", "devTempC", "devMediaID", "devMediaDate", "carrier", "carrierREF", "cellOrigin",	"cellType", "cellLine", "cellVariant", "cellClone", "cellDisease", "cellPassage", "diffDayFrozen", "diffFrozenID", "diffBy", "seedDensity", "seedDate", "notes"])

# change line BELOW for the experiment name
expName = '20240706_eisExp1'
# change line BELOW for the correct path to your raw data files
directoryToSplit = "C:\\Users\\myuser\\rawData\\" + expName + '\\'

# make table name
tableName = expName + ' lookup table'

# change lines BELOW for the correct path to the 'raw data' and 'lookup table' folders downloaded from the Github.
rawDestinationFold = 'C:\\Users\\myuser\\extracellularEIS\\raw data\\' # where the fit is occurring
lookupDestinationFold = 'C:\\Users\\myuser\\extracellularEIS\\lookup table\\'

# make dataframe for lookup table
for file in os.listdir(directoryToSplit):
    if "_freq" in file:
        print(file)
        currentIdx = len(fileList.index)
        # add check here later to look for duplicate files?
        fileList.at[currentIdx, 'plateID'] = file
        # if contains modelCell use 1, otherwise use 1.12
        # fileList.loc[len(fileList.index)]["plateID"] = file
        if "modelCell" in file:
            fileList.at[currentIdx, "measArea"] = 1
        else: # not model cell, using 12well TW
            fileList.at[currentIdx,"measArea"] = 1.12

# make excel file of lookup table and save to directoryToSplit
print(fileList)
# write excel file
fileList.to_excel(directoryToSplit + tableName+".xlsx")  

# copies all raw files to the rawDestinationFold
for file in os.listdir(directoryToSplit):
    if "_freq" in file:

        # open the source  file in binary mode with 'rb'
        source_file =  open(directoryToSplit+file, 'rb')
        # open the destination file in binary mode with 'wb'
        destination_file = open(rawDestinationFold+file, 'wb')

        # use the shutil.copyobj() method to copy the contents of source_file to destination_file
        shutil.copyfileobj(source_file, destination_file)
        print("copied from "+directoryToSplit+file+" to " + rawDestinationFold + file)

# copy lookup table from directoryToSplit to the folder

# open the source file in binary mode with 'rb'
source_file =  open(directoryToSplit + tableName+ ".xlsx", 'rb')
# open the destination file in binary mode with 'wb'
destination_file = open(lookupDestinationFold + tableName+".xlsx", 'wb')

# use the shutil.copyobj() method to copy the contents of source_file to destination_file
shutil.copyfileobj(source_file, destination_file)
print("copied "+directoryToSplit+ tableName+".xlsx to " + lookupDestinationFold + tableName+".xlsx")

# then run fits manually in the matlab script