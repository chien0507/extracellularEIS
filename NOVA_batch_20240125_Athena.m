%% SELECT LOOKUP TABLES
clear all; close all; clc
warning('off')

%% define default fit settings
settings = table();

settings.initial_guess_method = 'grid'; % Options: 'random', 'grid'
settings.grid_size = 3; % Applicable if initial_guess_method is 'grid', normally set to 3
settings.total_random_guesses = 15; % [15 seems like good compromise] total number of initial guesses to make for each circuit using the random method
settings.show_raw_data_and_fit = 1; % show each meas and resulting fit
settings.show_all_fits = 0; % toggle to show each surface of best fits

% settings.fit_circuit = 'RCRC_penalty' ; % select the circuit to fit 'RCRC','RCRCalpha',RCRCpRC,'RCRCdistributed'
settings.fit_circuit = 'RCRC' ; % select the circuit to fit 'RCRC','RCRCalpha',RCRCpRC,'RCRCdistributed'
settings.normalize_data = 1; % boolean to select if you want to normalize the data between -1 and 1. 
settings.target_tau_ratio = 1; % for RCRC_penalty fit_circuit, you can specify a custom target for the tau ratio. Must be double (e.g., 0.1, 1, 10)
settings.lambda_penalty = 1; % for RCRC_penalty fit_curcuit, you can scale the magnitude of penalty by lambda (default 1e-1)

settings.remove_positive_imag_impedance = 0; % optionally remove the positive imaginary impedance values
settings.R_max = 5000; % [Ohms.cm2] maximum expected tissue resistance
settings.C_max = 1e-4; % [F/cm2] maximum expected tissue capacitance
settings.C_min = 1e-8; % [F/cm2] maximum expected tissue capacitance
settings.R_solution_max = 1000; % [Ohms.cm2] maximum possible expected solution resistance - used only if non RC circuit

settings.make_tau_ratio_greater_than_1 = 1; % either 1 or 0. maes tau ratios always greater than or less than 1. 

settings.remove_excess_control_meas = 1; % optionally remove extra control measurments (performed on days that the cells were not measured)
settings.control_IDs = ['blank_transwell', 'model_cell', 'empty_scaffold']; % string names for control measurments
settings.signal_threshold = 3*(-3.2647); % (-3.2647*3) is the median minimum value of y for all blank transwells as of 20220615
settings.phase_angle_threshold = 89; % (89) [degrees] the steepest phase angle expected for the data

settings.max_allow_freq = 1e6; % [Hz] maximum allowable measurment frequency
settings.min_allow_freq = 0.1; % [Hz] minimum allowable measurment frequency

settings.smooth_freq_data = 0; % use the matlab smooth function on the real and imaginary data
settings.smooth_percentage = 0.001; % what percentage of the range of data recorded should be "smoothed"
if settings.smooth_freq_data == 0; settings.smooth_percentage = 0; end

% set the evom measurment frequency, the resistance that the evom would
% have measured is included in the summary table for comparison's sake
settings.f_evom = 12.5;  % [Hz]

% default settings
settings.LT_date_format = 'dd-MMM-yyyy';
settings.data_folder_name = 'raw data'; % relative to pwd
settings.lookup_folder_name = 'lookup table'; % relative to pwd
settings.fit_folder_name = 'FIT'; % relative to pwd

%% run program
lookup_table = cell(1,3);
i = 0;

% i = i+1; lookup_table(i,:) = {'Athena 20221130 lookup table 2',0,settings}; 
% i = i+1; lookup_table(i,:) = {['Athena 20220202_Glu7_Tox3 lookup table.xlsx'],1,settings}; 
% i = i+1; lookup_table(i,:) = {['Athena20230228 lookup table.xlsx'],1,settings}; 
% i = i+1; lookup_table(i,:) = {'clinical scaffold comparison 20230413 lookup table.xlsx',1,settings}; 
%  i = i+1; lookup_table(i,:) = {'Jiwon_Devika_Treatment lookup table.xlsx',1,settings}; 
% i = i+1; lookup_table(i,:) = {'Connie_Collection_20Jul23 lookup table.xlsx',1,settings};
% i = i+1; lookup_table(i,:) = {'Jiwon_Devika_Treatment_D2C lookup table.xlsx',1,settings};
% i = i+1; lookup_table(i,:) = {'Jair_Misc_3 lookup table.xlsx',1,settings};
% i = i+1; lookup_table(i,:) = {'CellBox_11Oct23 lookup table.xlsx',1,settings};
% i = i+1; lookup_table(i,:) = {'TSP 20231214 lookup table.xlsx',1,settings};
% i = i+1; lookup_table(i,:) = {'0124Jair misc lookup table.xlsx',1,settings};
% i = i+1; lookup_table(i,:) = {'Jiwon Colby Jair Starvation lookup table 2.xlsx',1,settings};
% i = i+1; lookup_table(i,:) = {'helena hts 20240125 lookup table.xlsx',1,settings};
% i = i+1; lookup_table(i,:) = {'ciliopathy lookup table 7.xlsx',1,settings};
% i = i+1; lookup_table(i,:) = {'MigrExp1New-fullSweep lookup table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'AthenaLowCaExp3 Lookup Table',0,settings};
% i = i+1; lookup_table(i,:) = {'AthenaLowCaExp3Troubleshooting Lookup Table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'30sMeasExp4 Lookup Table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'AthenaMediaChangeExp5-cut Lookup Table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'AthenaMediaChangeExp5-abridged Lookup Table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'AthenaMigrExp3 Lookup Table',0,settings};
% i = i+1; lookup_table(i,:) = {'AthenaElecMat10 Lookup Table',0,settings};
% i = i+1; lookup_table(i,:) = {'ARPEExp1 Lookup Table',0,settings};
% i = i+1; lookup_table(i,:) = {'ARPEExp2 Lookup Table',0,settings};
% i = i+1; lookup_table(i,:) = {'ElecMat9 Lookup Table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'ElecMat11 Lookup Table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'ElecMat12 Lookup Table',0,settings};
% i = i+1; lookup_table(i,:) = {'ElecMat13 Lookup Table',0,settings};
% i = i+1; lookup_table(i,:) = {'ElecMat14 Lookup Table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'LowCaExp4-2 Lookup Table.xlsx',0,settings}; 
% i = i+1; lookup_table(i,:) = {'20240730_MigrExp4 Lookup Table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'20240806_WPIbioExp1 Lookup Table',0,settings};
% i = i+1; lookup_table(i,:) = {'20240819_ApicalExp1 Lookup Table',0,settings};
% i = i+1; lookup_table(i,:) = {'20240827_LowCaExp5 Lookup Table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'20240827_ElecMat16 Lookup Table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'20240904_LowCaExp6 Lookup Table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'20240910_ApicalExp2 Lookup Table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'20240911_BasalExp1 Lookup Table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'20240917_LowCaExp7 Lookup Table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'20240924_LPSExp2 Lookup Table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'20240923_WPIbioExp2 Lookup Table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'20240228_dualMeasExp2 Lookup Table',0,settings};
% i = i+1; lookup_table(i,:) = {'20241001_ApicalExp3 Lookup Table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'20241001_BasalExp2 Lookup Table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'20241008_ApicalExp4 Lookup Table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'20241008_BasalExp3 Lookup Table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'20241105_WPIBioExp3 Lookup Table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'20241107_WPIBioExp4 Lookup Table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'20241112_BasalExp4 lookup table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'20241114_WPIBioExp5 lookup table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'20241120_WPIBioExp6_16HBE lookup table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'20241120_chirpBioExp1 lookup table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'ChirpBioExp1json lookup table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'ChirpBioExp1jsonGal lookup table',0,settings}; 

% i = i+1; lookup_table(i,:) = {'20240923_WPIbioExp22 Lookup Table',0,settings}; 

% i = i+1; lookup_table(i,:) = {'20241203_BasalExp5 Lookup Table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'20241204_BasalExp6 Lookup Table',0,settings};
% i = i+1; lookup_table(i,:) = {'ChirpBioExp2jsonGal lookup table',0,settings};
% i = i+1; lookup_table(i,:) = {'ChirpBioExp2jsonPot lookup table',0,settings};

% i = i+1; lookup_table(i,:) = {'ChirpBioExp2Cellsjson Lookup Table',0,settings}; 
% i = i+1; lookup_table(i,:) = {'LowCaExp4 Lookup Table',0,settings}; 

% i = i+1; lookup_table(i,:) = {'20241210_bubbleData Lookup Table',0,settings}; 
i = i+1; lookup_table(i,:) = {'20250122_LPSExp3 Lookup Table',0,settings}; 

tic
for idx = 1:i
    [T_master,T_final,T_summary] = NOVA_function_20240125_Athena(lookup_table{idx,:});
end
toc

% check if seed date exists in teh plate id before appending seed date to
% unique ID!