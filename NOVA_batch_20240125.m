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
i = i+1; lookup_table(i,:) = {'20240923_Exp2 Lookup Table',0,settings}; 


tic
for idx = 1:i
    [T_master,T_final,T_summary] = NOVA_function_20240125(lookup_table{idx,:});
end
toc
