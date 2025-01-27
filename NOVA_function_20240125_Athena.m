%% NOVA SUMMARY TABLE
% NOTE: %$ means Athena has edited this line

function [T_master,T_final,T_summary] = NOVA_function_20240125_Athena(lookup_table,remove_all_but_last_meas_of_day,settings)
    % this function imports lookup tables (of a particular format/naming 
    % convention) and performs fits on all associated raw data. After
    % fitting has completed, this script stores the data in a summary
    % table.

    % INPUT PROPERTIES:
    only_keep_last = remove_all_but_last_meas_of_day;
    
    %% SETTINGS: define and unpack

    total_random_guesses = settings.total_random_guesses; % [30 seems like good compromise] total number of initial guesses to make for each circuit
    show_raw_data_and_fit = settings.show_raw_data_and_fit; % show each meas and resulting fit
    show_all_fits = settings.show_all_fits; % toggle to show each surface of best fits
    
    remove_positive_imag_impedance = settings.remove_positive_imag_impedance; % optionally remove the positive imaginary impedance values
    R_max = settings.R_max; % [Ohms.cm2] maximum expected tissue resistance
    C_max = settings.C_max; % [F/cm2] maximum expected tissue capacitance
    R_solution_max = settings.R_solution_max; % [Ohms.cm2] maximum possible expected solution resistance - used only if non RC circuit
    
    make_tau_ratio_greater_than_1 = settings.make_tau_ratio_greater_than_1; % either 1 or 0. maes tau ratios always greater than or less than 1. 
    
    remove_excess_control_meas = settings.remove_excess_control_meas; % optionally remove extra control measurments (performed on days that the cells were not measured)
    control_IDs = settings.control_IDs; % string names for control measurments
    signal_threshold = settings.signal_threshold; % -3.2647 is the median minimum value of y for all blank transwells as of 20220615
    phase_angle_threshold = settings.phase_angle_threshold; % [degrees] the steepest phase angle expected for the data
    
    max_allow_freq = settings.max_allow_freq; % [Hz] maximum allowable measurment frequency
    min_allow_freq = settings.min_allow_freq; % [Hz] minimum allowable measurment frequency
    
    smooth_freq_data = settings.smooth_freq_data; % use the matlab smooth function on the real and imaginary data
    smooth_percentage = settings.smooth_percentage; % what percentage of the range of data recorded should be "smoothed"
    
    % set the evom measurment frequency, the resistance that the evom would
    % have measured is included in the summary table for comparison's sake
    f_evom = settings.f_evom;  % [Hz]
    
    % default settings
    LT_date_format = settings.LT_date_format;
    data_folder_name = settings.data_folder_name;
    lookup_folder_name = settings.lookup_folder_name;
    fit_folder_name = settings.fit_folder_name;

    % collect time stamp for this fitDate
    fitDate = datetime("now");
    fitDate.Format =  strcat(LT_date_format,' HH:mm:ss');
    
    %% LOOKUP TABLE: load
    
    % define the lookup table name and path
    lookup_table_name = lookup_table;
    lookup_table_name_and_path = fullfile(pwd,lookup_folder_name,lookup_table_name);

    % define the output summary table name based on the main lookup table
    % name
    T_summary_name = extractBefore(lookup_table_name," Lookup");
    T_summary_name = replace(T_summary_name," ","_");

    % define import options for lookup table to work correctly
    LTopts = setImportOptions(lookup_table_name_and_path);

    % read in the lookup table
    LT = readtable(lookup_table_name_and_path,LTopts); % LT == lookup table

    % move the controls to the end of the document. This is useful for
    % quickly checking that the actual data fits look good. Do the controls
    % last is just more convenient for debugging purposes.
    [n_LT_rows,~] = size(LT);
    for i = 1:n_LT_rows
        this_plateID = LT.plateID{1};
        if contains(control_IDs,this_plateID)
            this_LT_row = LT(1,:);
            LT(1,:) = [];
            LT(end+1,:) = this_LT_row;
        end
    end

    % change all date time arrays to cell strings for reading and writing
    % data
    LT = convertvars(LT, @isdatetime, @(t) cellstr(t));
    LT = convertvars(LT, @isduration, @(t) cellstr(t,'dd:hh:mm:ss'));

    % count the number of IDs in the lookup table
    [nRows_LT,~] = size(LT); 

    % the column called "notes" can have multiple tags. We need to load the
    % notes for each row in the LT, find the maximum number of unique
    % entries, and rebuild the notes column to be a cell array with any
    % missing values set to "none."
    no_tag_str = 'none';
    tag_delimiter = ';';
    tag = LT.notes;

    % replace empty cells with a no note tag
    idx_no_tag = cellfun(@isempty,tag);
    tag(idx_no_tag) = {no_tag_str};
    
    % if the number of tags is not equal to the maximum entry, then append
    % the no note tag, proportional to the number of missing entries
    delimiterN = count(tag,tag_delimiter); % get the number of delimiters for all entries
    maxN_delimiter = max(delimiterN); % get the maximum possible number of entries
    if maxN_delimiter > 0 % if there are multiple tags possible
        for idx_this_tag = 1:length(tag)
            this_tag = tag{idx_this_tag}; % get the current tag
            this_delimiterN = count(this_tag,tag_delimiter); % count the number of delimiters
            if this_delimiterN < maxN_delimiter % if this entry has fewer delimieters than the max possible...
                this_tag_appendix = ''; % init a tag to be appended
                for n_missing_tag = 1:(maxN_delimiter-this_delimiterN) % recursively build the suffix until the number of entries match the max possible
                    this_tag_appendix = strcat(this_tag_appendix,tag_delimiter,32,no_tag_str); % add delimiter followed by space followed by empty array
                end
                this_tag = strcat(this_tag,this_tag_appendix);
                tag{idx_this_tag} = this_tag; % update the tag value
            end
        end
    end
    % now that all tags have the same number of entries, split thepx_ cells at
    % the delimiter (if more than 0)
    if maxN_delimiter > 0
        tag = strtrim(split(tag,tag_delimiter));
    end

    % save the new values as tags in the LT
    LT.tag = tag;

    % in historical experiments the seeded date property for the cells
    % was called "seedDate" and now it is called "cellSeedDate." To
    % handle old lookup tables, we have added this method to add the
    % missing column to the lookup table by duplicating the seedDate
    % column to the cellSeedDate column
    if any(contains(LT.Properties.VariableNames,"seedDate")) && ~any(contains(LT.Properties.VariableNames,"cellSeedDate"))
        LT.cellSeedDate = LT.seedDate;
    end

    % if there are no specific wells set in the lookup table, the well array
    % returns an array of 'double' NaNs which has brace indexing errors. To
    % work around this, check if the type is double and if it is, change it to
    % an empty array of characters.
    if strcmp(class(LT.plateWellIDs),'double'); LT.plateWellIDs = cell(nRows_LT,1); end
     
    %% RAW DATA: read and format dir

    % detect all raw data files
    % specify folder with raw data and build path to folder
    data_path = fullfile(pwd,data_folder_name);
    
    % get the directory information of all raw data file in raw data folder
    dir_struct = struct2cell(dir(data_path));
    data_file_names = dir_struct(1,:);
    
    %% RAW DATA: match each LT plateID+wellID with raw data file
    % find the indicies of the raw data files and fix LT entries, when
    % needed

    matching_LTID = cell(nRows_LT,2);
    for i = 1:nRows_LT
        % get the ith unqiue cell line ID
        this_plateID = LT.plateID{i};

        % if the treatment date is empty or NaT, set it to the cell seeding 
        % date. if the entire column is empty, then the variable assignment
        % is not consistent with tables that have dates already
        if cellfun(@isempty,LT.treatDate(i)) || contains(LT.treatDate{i},'NaT')
            LT.treatDate(i) = LT.cellSeedDate(i);
        end

        % get the wells in this treatment. Must be a character array of
        % IDs. For example, 'A1 A2 A3'. These must match the suffix of the
        % raw data file name for it to match
        these_wells = LT.plateWellIDs{i};
    
        % loop through all raw files and find those with matching IDs
        for j = 1:length(data_file_names)
            this_file_name = data_file_names{j};
            this_file_name_prefix = extractBefore(this_file_name,' ');
            this_file_name_suffix = extractAfter(this_file_name,' ');
    
            % if the file name has the .txt suffix, remove it here to match
            % with the lookup table indicies
            if contains(this_file_name_suffix,'.')
                this_file_name_suffix = extractBefore(this_file_name_suffix,'.');
            end
            % if the string is empty after the loop above, re-write the entire
            % filename to the file prefix variable
            if isempty(this_file_name_prefix)
                this_file_name_prefix = this_file_name;
            end
            % not sure that this line is necessary
            if isempty(this_file_name_suffix)
                this_file_name_suffix = 'A1';
            end
            % if no specific wells are returned, then assume that the entire
            % transwell contains the treatment of interest.
            if isempty(these_wells)
                these_wells = ['A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 ' ...
                    'B1 B2 B3 B4 B5 B6 B7 B8 B9 B10 B11 B12 ' ...
                    'C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 C11 C12 ' ...
                    'D1 D2 D3 D4 D5 D6 D7 D8 D9 D10 D11 D12 ' ...
                    'E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11 E12 ' ...
                    'F1 F2 F3 F4 F5 F6 F7 F8 F9 F10 F11 F12 ' ...
                    'G1 G2 G3 G4 G5 G6 G7 G8 G9 G10 G11 G12 ' ...
                    'H1 H2 H3 H4 H5 H6 H7 H8 H9 H10 H11 H12'];
                LT.plateWellIDs{i} = these_wells;
            end
            
            % -----------------------------------------------------------------
            % critical logic for matching LT entries with ID and well. 
            if strcmp(this_plateID,this_file_name_prefix)
                if any(matches(split(these_wells),this_file_name_suffix))
                    matching_LTID{i,1} = [matching_LTID{i,1}; j];
                    matching_LTID{i,2} = [matching_LTID{i,2}; {this_file_name}];
                end
            end
            % -----------------------------------------------------------------
    
        end
    end
    
    % check for redundant files inside the idx_all_files and matching_LTID
    % scripts and keep each unique entry
    varnames = {'i','fileName'};
    matching_LTID = unique(table(vertcat(matching_LTID{:,1}),vertcat(matching_LTID{:,2}),'VariableNames',varnames),'rows');
%     % CHANGED^ to allow rerun
%     matching_LTID = table(vertcat(matching_LTID{:,1}),vertcat(matching_LTID{:,2}),'VariableNames',varnames);

    % count the total number of files to process
    [n_rawFiles,~] = size(matching_LTID);
    
    %% PROCESS FILES 
    % for each file, import the data, process it, and store it in a summary 
    % table

    % make a reduced form of the lookup table to identify biological
    % repeats within a lookup table
    LT_reduced = removevars(LT,["plateID","plateWellIDs","cellSeedDate","treatDate","treatTime","notes","tag"]); % removes subjective/relative columns from lookup table such as plate name. diffFrozenID and techRepeat are the same thing so ignore it here to work correctly
    [~,~,repeat_idx] = unique(LT_reduced);
    unique_grid = zeros(length(repeat_idx),length(repeat_idx));
    for this_repeat_idx = 1:max(repeat_idx)
        unique_grid(:,this_repeat_idx) = repeat_idx == this_repeat_idx;
    end
    % the technical repeat array assigns a number to each technical
    % replicate. This makes assigning properties in plotting functions
    % easier since it is simply a number, independent of cell line and
    % frozen date, etc.
    techRepeat_array = sum( (cumsum(unique_grid).*unique_grid)' )'; 

    % initialize ouput large tables for storing in after this loop
    T_frequency = [];
    T_summary = []; 
    fit_table = []; T_summary_abp = []; % used for displaying results in command window only
    lookup_table_name = {lookup_table_name}; % needed to be a cell for saving to T_summary

    for idx_rawFile = 1:n_rawFiles
        % select a file to load into memory 
        this_idx = matching_LTID.i(idx_rawFile);
        this_file_name = data_file_names{this_idx};
        this_file = fullfile(pwd,data_folder_name,this_file_name);
        
        % use readcell to import
        fprintf('loading EIS data: %s ...',this_file)
        N = readcell(this_file,'DatetimeType','datetime');
        fprintf(' done!\n')
    
        % measure the data length
        [data_length,~] = size(N);
    
        % remove any data rows that are listed as infinite or NaN. This is
        % a strange behavior of some NOVA outputs and needs to be removed.
        remove_this_row = [];
        for i = 1:data_length
            this_row = N(i,:);
            for ii = 1:length(this_row)
                this_entry = this_row{ii};
                if ischar(this_entry)
                    if contains(this_entry,'∞')
                        % fprintf('contains a infinite entry!\n')
                        remove_this_row = [remove_this_row; i];
                    end
                end
            end
        end
        if length(remove_this_row)>1
            remove_this_row = unique(remove_this_row);
            N(remove_this_row,:) = [];
        end
        % RE-measure the data length!!!!!!!!!!!!!!
        [data_length,~] = size(N);
    
        % find all cells that contain the remark that starts each 
        % measurment. For NOVA, this remark is a date array.
        idx = [];
        for i = 1:data_length
            if isdatetime(N{i,1}); idx = [idx i]; end
        end
        % transpose idx array for legibility
        idx = idx';    
    
        % get the column numbers for each entry that we expect from the NOVA
        % file
        column_names = N(2,:);
        t_col = []; f_col = []; x_col = []; y_col = []; ocp_col = []; A_col = []; meas_start_col = [];
        for i = 1:length(column_names)
           if strcmp(column_names{i},'Time (s)') || strcmp(column_names{i},'Time')%$ changed bc potentiostat code not outputting files with right col name
               t_col = i;
               if strcmp(column_names{i},'Time')%$
                   disp("Time column is labelled 'Time' not 'Time (s)'");%$
               end%$
           elseif strcmp(column_names{i},'Frequency (Hz)')
               f_col = i;
           elseif strcmp(column_names{i},'Z'' (Ω)') 
               x_col = i;
           elseif strcmp(column_names{i},'-Z'''' (Ω)')
               y_col = i;
           elseif strcmp(column_names{i},{'OCP value (V)'})
               ocp_col = i;
           elseif strcmp(column_names{i},{'A ratio'})
               A_col = i;
           elseif strcmp(column_names{i},{'meas start time (s)'})
               meas_start_col = i;
           elseif strcmp(column_names{i},{'Z (Ω)'})
               magnitude_col = i;
           elseif strcmp(column_names{i},{'-Phase (°)'})
               phase_col = i;
           end 
        end
        
        % build a 3D array of the data, where the 3rd dimension is equal to the
        % number of measurements recorded
        num_meas = length(idx);
        each_meas_length = diff([idx; data_length+1]); % +1 becuase this algorithm is designed to have a date below the last measurment. All lengths are, thus, 1 greater than they should be, this +1 makes the final measurment match this error
        max_meas_length = max(each_meas_length);
        [~,n_meas_cols] = size(N);
        D = cell(max_meas_length,n_meas_cols,num_meas);
        for i = 1:num_meas
            if i<num_meas
                D(1:each_meas_length(i),:,i) = N(idx(i):idx(i+1)-1,:);
            else
                D(1:each_meas_length(i),:,i) = N(idx(i):end,:);
            end
        end
    
        % OPTIONAL 
        % remove all but the last measurment of the day
        % find and remove redundant measurments
        all_exp_dates = N(idx,1);
        % keep only the last unique measurment on each day (asuming that
        % bad measurments came before the final "good" measurment
        if only_keep_last
            idx_remove = [];
            for i = 1:length(idx)  
                all_exp_dates_string{i,1} = datestr(all_exp_dates{i});
                if i>1
                    if strcmp(all_exp_dates_string{i,1},all_exp_dates_string{i-1,1})
                        idx_remove = [idx_remove i-1]; % 
                    end
                end
            end
            if length(all_exp_dates)>1
                idx(idx_remove) = [];
                all_exp_dates(idx_remove) = [];
                D(:,:,idx_remove) = [];
            end
        end
    
        % if this measurment is a control measurment, it is possible that the
        % control was measured on many days - other than the days that this
        % particular experiment was performed. As such, it is useful to remove
        % the unecessary control measurments
        filename_before_space = extractBefore(data_file_names{this_idx},' ');
        if isempty(filename_before_space); filename_before_space = data_file_names{this_idx}; end
        is_control_ID = contains(control_IDs,filename_before_space);
        if and(is_control_ID,remove_excess_control_meas)
            % find all the unique dates that measurements were performed
            fprintf('filtering control measurement...\n')
            if exist('T_summary')
                idx_keep = [];
                all_T_summary_measDate_NOVA = unique(datetime(T_summary.measDate_NOVA));
                for exp_idx = 1:length(all_T_summary_measDate_NOVA)
                    this_exp_date = all_T_summary_measDate_NOVA(exp_idx);
                    % check to see if this_exp_date matches any dates in
                    % all_exp_dates array (the array of all exp dates for this
                    % specific file). If it does not, then delete all of these
                    % dates from the control file 
    
                    % this line creates a time table array which can be indexed
                    % to see if a date exists in the array
                    tmp = timetable(vertcat(all_exp_dates{:}),linspace(1,length(all_exp_dates),length(all_exp_dates))');
                    idx_keep = [idx_keep; tmp.Var1(this_exp_date)];
                end
    
                % reduce the control arrray to only the indicies to keep
                if length(idx_keep)>1
                    idx = idx(idx_keep);
                    all_exp_dates = all_exp_dates(idx_keep);
                    D = D(:,:,idx_keep);
                else
                    idx = idx(end);
                    all_exp_dates = all_exp_dates(end);
                    D = D(:,:,end);
                end
            end
        end

        % for each measurment in the D array, get the date, and process the
        % data
        num_meas = length(idx);
        sdSignal = []; % array of all signal standard deviations
        ySignalMin = []; % temporary array to find maximum y in the signal
        for n = 1:num_meas
            % check if T_summary entry exists already
    
            measDate_NOVA = all_exp_dates{n};
            t = cell2mat(D(3:each_meas_length(n),t_col,n)); %$% most common error, t(1), because t = []; t is the start times for each freq, usually means measurement didn't complete.
            fprintf('meas start time: %1.3f seconds\n',t(1))
            f = cell2mat(D(3:each_meas_length(n),f_col,n));
            x = cell2mat(D(3:each_meas_length(n),x_col,n));
            y = -cell2mat(D(3:each_meas_length(n),y_col,n));
            m = cell2mat(D(3:each_meas_length(n),magnitude_col,n)); % signal magnitude
            fprintf('filename: %s\n',data_file_names{this_idx})
            p = -cell2mat(D(3:each_meas_length(n),phase_col,n)); % signal phase in degrees
            ocp = cell2mat(D(3,ocp_col,n));
    
            % track the absolute minimum y in the signal for establishing a
            % crossover threshold
            ySignalMin = [ySignalMin; min(y)];
    
            % A col and meas start time columns in the txt file will only exist if the
            % labview code that synchronizes EIS data was run in parallel with the
            % nova. To make this code compatible with measurments performed without
            % intracellular recordings, check to see if A_col and meas_start_col were
            % detected. If they were not, set them to empty values
            if ~isempty(A_col)
                A = cell2mat(D(3,A_col,i));
            else
                A = 0;
            end
            if ~isempty(meas_start_col)
                measStart = cell2mat(D(3,meas_start_col,i));
            else
                measStart = t(1);
            end
            if isempty(ocp_col)
                ocp = NaN;
            end
    
            % remove positive outliers
            if remove_positive_imag_impedance
                idx_remove = find(y>0);
                t(idx_remove) = [];
                f(idx_remove) = [];
                x(idx_remove) = [];
                y(idx_remove) = [];
            end
        
            % remove frequences above max
            idx_remove = find(f>max_allow_freq);
            max_freq = max(f);
            t(idx_remove) = [];
            f(idx_remove) = [];
            x(idx_remove) = [];
            y(idx_remove) = [];
        
            % remove frequences below max
            idx_remove = find(f<min_allow_freq);
            min_freq = min(f);
            t(idx_remove) = [];
            f(idx_remove) = [];
            x(idx_remove) = [];
            y(idx_remove) = [];
    
            % build omega and zdata arrays for data analysis and fitting
            w = 2.*pi.*f;
            zdata = [x y];
            dt_freq = t-t(1);
    
            % smooth data if desired
            if smooth_freq_data
                x = smooth(log10(f),x,smooth_percentage,'lowess');
                y = smooth(log10(f),y,smooth_percentage,'lowess');
            end
        
            % % make the guesses for this measurement
            % clear p0 % ensure that the any previous initial guesses are cleared from memory
            % p0 = cell(1,1);
            % ter = round(max(abs(zdata(:,1)))); if ter<=0; ter=1; end
            % tec = C_max; % maximum expected tec
            % inverse_tec = round(1/tec); % taking the inverse allows us to use the randi function (which only generates random integers. if tec == e-6, thats too small w/out inverting)
            % for num_guess = 1:total_random_guesses
            %     p0{num_guess} = [randi(ceil(min(abs(zdata(:,1)))))*1.1 ...
            %         randi(ter)*10 2/randi(inverse_tec)...
            %         randi(ter)*10 2/randi(inverse_tec) ...
            %         1]; % 1 is for the alpha term
            % end
            % Make the guesses for this measurement
            clear p0 % Ensure that any previous initial guesses are cleared from memory
            p0 = cell(1,1);
            ter = round(max(abs(zdata(:,1)))); if ter <= 0; ter = 1; end
            tec = settings.C_max; % Maximum expected tec
            inverse_tec = round(1/tec); % Taking the inverse allows us to use the randi function
            
            % Define initial guesses based on the selected method
            switch settings.initial_guess_method
                case 'random'
                    % Existing random guessing strategy
                    for num_guess = 1:settings.total_init_guesses
                        p0{num_guess} = [randi(ceil(min(abs(zdata(:,1)))))*1.1 ...
                                         randi(ter)*10 2/randi(inverse_tec)...
                                         randi(ter)*10 2/randi(inverse_tec) ...
                                         1]; % 1 is for the alpha term
                    end
            
                case 'grid'
                    % Define ranges for each parameter based on settings
                    Rsol_range = linspace(1, settings.R_solution_max, settings.grid_size);
                    R_range = linspace(1, settings.R_max, settings.grid_size);
                    C_range = linspace(settings.C_min, settings.C_max, settings.grid_size);
                    
                    % Generate a grid of initial guesses
                    [Rsol, R1, C1, R2, C2] = ndgrid(Rsol_range, R_range, C_range, R_range, C_range);
                    Rsol = Rsol(:); R1 = R1(:); C1 = C1(:); R2 = R2(:); C2 = C2(:); % Flatten the arrays
                    
                    % Limit the number of initial guesses to the user-defined total
                    total_combinations = length(Rsol);
                    %p0 = cell(total_combinations,1);  % Pre-allocate the cell array for initial guesses
                    
                    for idx = 1:total_combinations
                        % Assign the grid values to p0, ensuring C values are in the correct format
                        p0{idx} = [Rsol(idx), R1(idx), C1(idx), R2(idx), C2(idx), 1];  % Note the inversion of C values, 1 is for alpha
                    end
                otherwise
                    error('Invalid initial_guess_method. Choose "random" or "grid".');
            end

    
            % find the metadata for this cell line in the look up table and
            % find the index in the lookup table that matches this ID
            fileName = {this_file_name};
            plateID = {extractBefore(data_file_names{this_idx},' ')};
            wellID = {extractAfter(data_file_names{this_idx},' ')};
    
            % check for file extensions and ignore if it exists
            if contains(wellID,'.')
                wellID = extractBefore(wellID,'.');
            end
    
            if isempty(plateID{:})
                plateID = fileName;
                wellID = {'A1'};
            end
    
            idx_LT_ID = strcmp(vertcat(LT.plateID(:)),plateID);
            idx_LT_well = [];
            for idx_LT_row = 1:length(idx_LT_ID)
                % check each row to see if the well name was tested
                idx_LT_well(idx_LT_row) = any(strcmp(split(LT.plateWellIDs{idx_LT_row}),wellID));
            end
            idx_LT = find(idx_LT_ID.*idx_LT_well');
            if length(idx_LT)>1
                % if more than one LT entry with the same ID exsists, compare
                % the measurment date in the NOVA file with the measurment date
                % in the lookup table. Optionally group this data by the
                % closest datetime?????????????????
                tmp_LT = LT.treatDate(idx_LT);
                for idx_each_LT = 1:length(tmp_LT)
                    if contains(tmp_LT(idx_each_LT),'NaT')
                        tmp_LT(idx_each_LT) = cellstr(measDate_NOVA);
                    end
                end
                LT.treatDate(idx_LT) = tmp_LT;
                unique_LT_treatmentDates = dateshift(datetime(unique(LT.treatDate(idx_LT)),'InputFormat',LT_date_format),'start','day');
                duration_between_LT_dates = between(measDate_NOVA,unique_LT_treatmentDates,'Days');
                duration_between_LT_dates = caldays(duration_between_LT_dates); % convert the duration array to type double
                % remove all LT datetimes that are before this meas date
                idx_LT_dates_remove = duration_between_LT_dates>0;
                % get the idx of the date in the unique_LT_treatmentDates
                % folder
                dates_to_remove = unique_LT_treatmentDates(idx_LT_dates_remove);
                % correct the idx of dates to remove
                idx_LT_dates_remove = contains(LT.treatDate(idx_LT),cellstr(dates_to_remove));
                % remove the dates idx from the lookup table options
                idx_LT(idx_LT_dates_remove) = [];
                % determine the durations of the viable rows in the lookup
                % table
                duration_between_LT_dates = caldays(between(measDate_NOVA,datetime(LT.treatDate(idx_LT)),'Days'));
                [~,idx_LT_keep] = min(abs(duration_between_LT_dates));
                idx_LT = idx_LT(idx_LT_keep);  
                if length(idx_LT) > 1
                    warning('Duplicate entries for this ID and well have been found!\nUsing the first entry only!\n')
                    idx_LT = idx_LT(1);
                end
            end
            this_LT = LT(idx_LT,:); %% CHECK THAT THIS MATCHES UP THE INDEX AND ROW CORRECTLY
            measArea = this_LT.measArea;
    
    
            % ----------------------------------------------------------------
            % fit data
            % ----------------------------------------------------------------
            [pfit,resnorm,idx_best_fit,info,all_results,fun] = fitZ12_20240125(p0,w,zdata,settings,R_max,C_max,R_solution_max,measArea);
            firstorderopt = all_results.output{idx_best_fit,1}.firstorderopt;
            residual =  all_results.residual{idx_best_fit};
            exitflag = all_results.exitflag{idx_best_fit};
            % ----------------------------------------------------------------
            
    
            % assign outputs
            switch settings.fit_circuit
                case 'RCRCdistributed'
                    Rblank = pfit(1);
                    % two RCRC circuits in parallel
                    R1 = pfit(2)+pfit(4);
                    R2 = pfit(6)+pfit(8);
                    C1 = 1/(1/pfit(3)+1/pfit(5));
                    C2 = 1/(1/pfit(7)+1/pfit(9));
                    alpha = pfit(10);
                case 'RCRCpRC'
                    Rblank = pfit(1);
                    R1 = pfit(2);
                    C1 = pfit(3);
                    R2 = pfit(4);
                    C2 = pfit(5);
                    Rp = pfit(6);
                    Cp = pfit(7);
                    alpha = pfit(8);
                otherwise
                    Rblank = pfit(1);
                    R1 = pfit(2);
                    C1 = pfit(3);
                    R2 = pfit(4);
                    C2 = pfit(5);
                    alpha = pfit(6);
            end
    
            % check the SNR of the signal here to decide if it makes sense to
            % include the sep R1 C1 R2 C2 portion of the circuit. 
            distance_btwn_points = sqrt(diff(x).^2+diff(y).^2);
            sdSignal = [sdSignal std(distance_btwn_points)]; % not stored or used anywhere at the moment
    
            % run fit in model
            %z_fit = fun(log(pfit),w); % remember, funRCRC takes exp of input variables to allow solvers to provide equal contribution for R and C despite magnitude differences
            wfit = 2.*pi().*logspace(log10(min_allow_freq),log10(max_allow_freq),1000);
            z_fit = fun(log10(pfit),wfit); % remember, funRCRC takes 10^x of input variables to allow solvers to provide equal contribution for R and C despite magnitude differences
            x_fit = z_fit(:,1);
            y_fit = z_fit(:,2);
    
            % the evom 3 takes measurments at 12.5 hertz, use this frequency to
            % see what resistance you would have measured with the EVOM
            z_evom = fun(log10(pfit),f_evom.*pi.*2); % remember, funRCRC takes 10^x of input variables to allow solvers to provide equal contribution for R and C despite magnitude differences
            x_evom = z_evom(:,1);
            y_evom = z_evom(:,2);
    
            if show_raw_data_and_fit
                figure(1)
                cla
                plot(x,y,'ob','MarkerSize',10)
                hold on
                plot(x_fit,y_fit,'-k','Linewidth',1)
                plot(x_evom,y_evom,'sr','MarkerSize',15)
                hold off
                if max(x) > 0
                    xlim([0 max(x)])
                end
                axis equal
                xlabel('Real')
                ylabel('Imag')
                title(strcat('Z12 fit: ', makeTissueID(this_LT,wellID),'-rep',num2str(techRepeat_array(idx_LT)),'-n',num2str(n)))   
            end
    
            if show_all_fits
                tmp_all_fits = vertcat(all_results.pfit{:});
                all_r1 = tmp_all_fits(:,2);
                all_r2 = tmp_all_fits(:,4);
                all_c1 = tmp_all_fits(:,3);
                all_c2 = tmp_all_fits(:,5);
                all_resnorm = vertcat(all_results.resnorm{:});
    
                figure(2)
                cla
                %scatter3(all_r1,all_r2,all_resnorm,'filled','AlphaData',0.5)
                scatter(all_r1+all_r2,all_resnorm,'filled','AlphaData',0.2)
                xlabel('Rt')
                ylabel('resnorm')
                %zlabel('resnorm')
                set(gca,'FontSize',12)
                
                figure(3)
                cla
                %scatter3(all_c1,all_c2,all_resnorm,'filled','AlphaData',0.5)
                scatter(1./(1./all_c1+1./all_c2),all_resnorm,'filled','AlphaData',0.2)
                xlabel('Ct')
                ylabel('resnorm')
                %zlabel('resnorm')
                set(gca,'FontSize',12)
    
                %drawnow()
            end
    
            % calculate the tau ratio of the fit
            t1 = R1*C1;
            t2 = R2*C2;
            if make_tau_ratio_greater_than_1
                if t1>t2
                    tau_r = t1/t2;
                else
                    tau_r = t2/t1;
                end
            else
                if t2>t1
                    tau_r = t1/t2;
                else
                    tau_r = t2/t1;
                end
            end
    
            % CELL AGE DURATIONS
            cellSeedDate = this_LT.cellSeedDate;
            cellDaysInVitro = getDaysInVitro(this_LT);
            if ~(cellfun(@isempty,cellSeedDate) || contains(cellSeedDate,'NaT'))
                seedAgeDays = between(cellSeedDate,measDate_NOVA,'Days');
            else
                seedAgeDays = calendarDuration([0 0 0]);
            end
            seedAgeDays = caldays(seedAgeDays);
            cellAgeDays = seedAgeDays + cellDaysInVitro;
            seedAgeWeeks = manualAgeWeeks(seedAgeDays);
            cellAgeWeeks = manualAgeWeeks(cellAgeDays);

            % CULTURE MEDIA DURATIONS
            cultMediaDate = this_LT.cultMediaDate;
            if ~(cellfun(@isempty,cultMediaDate) || contains(cultMediaDate,'NaT'))
                cultMediaAgeDays = between(cultMediaDate,measDate_NOVA,'Days');
                cultMediaAgeDays = caldays(cultMediaAgeDays);
            else
                cultMediaAgeDays = NaN;
            end

            % MEASURMENT DEVICE MEDIA DURATIONS
            devMediaDate = this_LT.devMediaDate;
            if ~(cellfun(@isempty,devMediaDate) || contains(devMediaDate,'NaT'))
                devMediaAgeDays = between(devMediaDate,measDate_NOVA,'Days');
                devMediaAgeDays = caldays(devMediaAgeDays);
            else
                devMediaAgeDays = NaN;
            end

            % define variable names for table metadata
            t0 = t(1);
            fitCircuit = {settings.fit_circuit};
            target_tau_ratio = settings.target_tau_ratio;
            lambda_penalty = settings.lambda_penalty;
            initialGuessMethod = {settings.initial_guess_method};
            gridSize = settings.grid_size;
            totalRandomGuesses = settings.total_random_guesses;
            normalizeData = {settings.normalize_data};
            techRepeatNum = techRepeat_array(idx_LT);
            medianSdSignal = median(sdSignal);
            Rt = (R1+R2);
            Ct = 1/(1/C1+1/C2);
            evomResistance = x_evom;
            fMin = min(f);
            fMax = max(f);
            nFreq = length(f);
            dy = max(y)-min(y);
            dx = max(x)-min(x);
            dydx = dy/dx;
            phaseAngleMax = max(abs(p));
            minY = min(y);
            R1pR2 = R1+R2;
            measuredResistanceOnly = ~((minY<signal_threshold) & (R1pR2>abs(signal_threshold)));
            % measuredResistanceOnly = 1;
            
            if measuredResistanceOnly
                % try a fit to a single point resistor
                normalization_array_R = max(abs(x));
                funR = @(p,w) p(1)./normalization_array_R.*ones(length(w),1);
        
                %lsqcurvefit(funR,Rblank,w,x./normalization_array,0,R_solution_max);
                pfit_rsol = lsqcurvefit(funR,Rblank,w,x./normalization_array_R,0,R_solution_max,info.opts);
    
                % replace RCRC outputs with single point resistor solution
                Rt = NaN;
                %evomResistance = Rblank;
                R1 = NaN;
                R2 = NaN;
                C1 = NaN;
                C2 = NaN;
                Rblank = pfit_rsol;
                Ct = NaN;
                tau_r = NaN;
            end
    
            % add cross sectional area to all electrical values
            R1Area = R1*measArea;
            R2Area = R2*measArea;
            C1Area = C1/measArea;
            C2Area = C2/measArea;
            RblankArea = Rblank*measArea;
            evomResistanceArea = evomResistance*measArea;
            RtArea = Rt*measArea;
            CtArea = Ct/measArea;
            tau1 = R1Area*C1Area;
            tau2 = R2Area*C2Area;
            tauRtCt = RtArea*CtArea;
            tauMEM = tau1*tau2/(RtArea*CtArea);
            
            tissueID = makeTissueID(this_LT,wellID);
            measID = strcat(tissueID,'-rep',num2str(techRepeatNum),'-n',num2str(n));
            treatment = this_LT.treatment;
            treatmentDateTime = datetime(this_LT.treatDate)+duration(this_LT.treatTime);
            cellDensity = calculateCellDensity(this_LT);
            %measID = strcat(tissueID,'-',datestr(treatmentDateTime,'mmm-dd-yyyy-HH:MM:SS'),'-',num2str(n));
            TEP = ocp; % for nova ocp is Apical side positive TEP which is standard convention

            % show results in command window
            % fit_table = [fit_table; table(n,RblankArea,R1Area,C1Area,R2Area,C2Area,exitflag,RtArea,CtArea,tauRtCt,tauMEM)];
            fit_table = [fit_table; table(plateID,wellID,cellAgeDays,RblankArea,RtArea,CtArea,tau_r,R1Area,C1Area,R2Area,C2Area,exitflag)];
            fprintf('\n')
            %uit = uitable(figTABLE,'Data',fit_table);
            %drawnow limitrate
            %pause(0.01);
            disp(fit_table)
            resnorm_info = strcat('fit:',32,num2str(n),32,'out of',32,num2str(num_meas),...
                32,'[resnorm =',32,num2str(resnorm),', firstorderopt=',32,num2str(firstorderopt),']');
            fprintf('most recent fit: %s, %1.0f out of %1.0f files (%1.1f%%) ...\n', ...
                data_file_names{this_idx}, ...
                idx_rawFile,...
                n_rawFiles,...
                idx_rawFile/n_rawFiles*100)
            disp(resnorm_info)
            fprintf('-------------------------------------------------------\n')    

            tmp_lhs = table(plateID,wellID,treatment,cellAgeDays,cellAgeWeeks,seedAgeDays,seedAgeWeeks,techRepeatNum,...
                RtArea,RblankArea,CtArea,TEP,...
                evomResistanceArea,...
                measuredResistanceOnly,...
                R1Area,R2Area,C1Area,C2Area,tau1,tau2,tau_r,tauRtCt,tauMEM);

            tmp_rhs = table(treatmentDateTime,fileName,measDate_NOVA,...
                cellDensity,cultMediaAgeDays,devMediaAgeDays,...
                t0,measStart,...
                n,nFreq,fMin,fMax,...
                fitCircuit,normalizeData,target_tau_ratio,lambda_penalty,initialGuessMethod,gridSize,totalRandomGuesses,alpha,exitflag,resnorm,firstorderopt,...
                idx_best_fit,total_random_guesses,...
                smooth_freq_data,smooth_percentage,...
                signal_threshold,minY,R1pR2,medianSdSignal,phaseAngleMax,fitDate,lookup_table_name,...
                measID,tissueID);

            this_LT_reduced = removevars(this_LT,{'plateID','treatment'}); % if want to reorder lookup table, and include in tmp tables, need to remove here
            tmp = [tmp_lhs, tmp_rhs, this_LT_reduced];
            T_summary = [T_summary; tmp];
    
        end 
    
    end
    
    %% Set relative time, measurement IDX, normalized ephys values
    % find all unique datetime entries in the summary table, for each unique
    % entry, subtract the minimum time value from the array
    unique_meas_in_T_summary = unique(table(datetime(T_summary.treatDate),T_summary.fileName,'VariableNames',{'treatDate','fileName'}),'rows');
%     unique_meas_in_T_summary = table(datetime(T_summary.treatDate),T_summary.fileName,'VariableNames',{'treatDate','fileName'}); % CHANGED because want to rerun
    [n_unique_meas,~] = size(unique_meas_in_T_summary);
    [n_entries,~] = size(T_summary);
    meas_number = zeros(n_entries,1);
    for i = 1:n_unique_meas
        this_treatmentDate = unique_meas_in_T_summary.treatDate(i);
        this_fileName = unique_meas_in_T_summary.fileName(i);
    
        idx_this_treatmentDate = datetime(T_summary.treatDate) == this_treatmentDate;
        idx_this_fileName = strcmp(T_summary.fileName,this_fileName);
        idx_both_match = idx_this_treatmentDate.*idx_this_fileName;
    
        % find all indicies that match requirments for unique, continuous
        % recording
        idxs_unique_recording = find(idx_both_match);
    
        % get value and index of start time
        [mint0,idxt0] = min(T_summary.t0(idxs_unique_recording));
        idxt0 = min(idxs_unique_recording)+idxt0-1;
    
        % normalized values
        T_summary.nRtArea(idxs_unique_recording) = (T_summary.RtArea(idxs_unique_recording) - T_summary.RtArea(idxt0))./T_summary.RtArea(idxt0);
        T_summary.nCtArea(idxs_unique_recording) = (T_summary.CtArea(idxs_unique_recording) - T_summary.CtArea(idxt0))./T_summary.CtArea(idxt0);
        T_summary.nRblankArea(idxs_unique_recording) = (T_summary.RblankArea(idxs_unique_recording) - T_summary.RblankArea(idxt0))./T_summary.RblankArea(idxt0);
        T_summary.nR1Area(idxs_unique_recording) = (T_summary.R1Area(idxs_unique_recording) - T_summary.R1Area(idxt0))./T_summary.R1Area(idxt0);
        T_summary.nR2Area(idxs_unique_recording) = (T_summary.R2Area(idxs_unique_recording) - T_summary.R2Area(idxt0))./T_summary.R2Area(idxt0);
        T_summary.nC1Area(idxs_unique_recording) = (T_summary.C1Area(idxs_unique_recording) - T_summary.C1Area(idxt0))./T_summary.C1Area(idxt0);
        T_summary.nC2Area(idxs_unique_recording) = (T_summary.C2Area(idxs_unique_recording) - T_summary.C2Area(idxt0))./T_summary.C2Area(idxt0);
    
        % relative to start values
        T_summary.t0(idxs_unique_recording) = T_summary.t0(idxs_unique_recording) - mint0;
        T_summary.rTEP(idxs_unique_recording) = T_summary.TEP(idxs_unique_recording) - T_summary.TEP(idxt0);
    
        meas_number = meas_number + idx_both_match.*i;
    end
    T_summary.meas_number = meas_number;
    
    
    %% save the results to a FIT folder

    if ~isfolder(fit_folder_name); mkdir(fit_folder_name); end

    T_summary = label_final_meas_of_day(T_summary);

    fprintf('writing summary table to %s ... ',fit_folder_name)
    summary_table_name = strcat(T_summary_name,'_summary_table.csv');
    writetable(T_summary,fullfile(pwd,fit_folder_name,summary_table_name));
    fprintf('done!\n')

    fprintf('median signal standard deviation (used for controls): %1.4f\n',medianSdSignal) 
    

    T_master = []; T_final = [];
    % update the master file
    % [T_master,master_folder] = updateMaster(T_summary);
    
    % update the final meas file
    % T_final = updateFinal(T_master,master_folder);
    
    disp('all done!')
    fprintf('\n')
    
    %% Label final meas of each day
    function T_summary = label_final_meas_of_day(T_summary)
        % get the unique tissue IDs
        fprintf('labeling final meas of each day in summary table...\n')
        [uniqueTissues,~,uniqueTissues_ic] = unique(T_summary.tissueID);
        T_summary.is_final_meas_of_day = zeros(length(T_summary{:,1}),1);

        for this_unique = 1:length(uniqueTissues)
            all_entries = T_summary(uniqueTissues_ic==this_unique,:);
            % get all unique measurment days for this tissue
            all_seedAgeDays = unique(all_entries.seedAgeDays);
            for idx_thisDay = 1:length(all_seedAgeDays)
                all_idx_for_this_day = find(and(uniqueTissues_ic==this_unique,all_seedAgeDays(idx_thisDay)==T_summary.seedAgeDays));
                T_summary.is_final_meas_of_day(max(all_idx_for_this_day))=1;
            end
            
        end

    end

    %% IMPORT OPTIONS FOR LTs
    function importTableOpts = setImportOptions(table_path)
        % Detect initial import options from the file
        importTableOpts = detectImportOptions(table_path);
        
        % Retrieve the names of all columns in the table
        varNames = importTableOpts.VariableNames;
        
        % Define settings for each column in a structured array
        columnSettings = struct(...
            'Name', {'treatTime', 'cultMediaDate', 'devMediaDate', 'carrierREF', 'cellOrigin', 'cellDisease', 'cellLine', 'cellVariant', 'cellClone', 'cellType', 'diffFrozenID', 'diffDayFrozen', 'cellDaysInVitro', 'diffBy','nCellsSeeded'}, ...
            'Type', {'duration', 'datetime', 'datetime', 'char', 'char', 'char', 'char', 'char', 'char', 'char', 'char', 'char', 'char', 'char', 'char'}, ...
            'FillValue', {hours(0), NaT, NaT, '', '', '', '', '', '', '', '', '', '', '', ''} ...
        );
        
        % Iterate over the column settings
        for i = 1:length(columnSettings)
            % Check if the current column exists in the table
            if ismember(columnSettings(i).Name, varNames)
                % Set the variable type for the existing column
                importTableOpts = setvartype(importTableOpts, {columnSettings(i).Name}, columnSettings(i).Type);
                
                % Check if a fill value is defined and needs to be set
                if ~isempty(columnSettings(i).FillValue)
                    % Set the fill value for the existing column
                    importTableOpts = setvaropts(importTableOpts, {columnSettings(i).Name}, 'FillValue', columnSettings(i).FillValue);
                end
            end
        end
    end


    %% small function to make custom tissueID
    function tid = makeTissueID(LT_entry,wi)
        % format the cell origin in case of multiple user inputs
        co = LT_entry.cellOrigin{:};
        co = replace(co,' ','');

        patterns = {'human','Human'}; % use full words first
        abbrev = 'h';
        co = replace(co,patterns,abbrev);

        patterns = {'pig','Pig'}; % use full words first
        abbrev = 'p';
        co = replace(co,patterns,abbrev);

        patterns = {'mouse','Mouse'}; % use full words first
        abbrev = 'm';
        co = replace(co,patterns,abbrev);

        patterns = {'rat','Rat'}; % use full words first
        abbrev = 'r';
        co = replace(co,patterns,abbrev);

        patterns = {'induced pluripotent stem cell','induced pluripotent stem cells','pluripotent stem cells','pluripotent stem cell','iPSC','iPSCs','ipscs','ipsc','iPSs','iPS','ipss'};
        abbrev = 'ips';
        co = replace(co,patterns,abbrev);

        patterns = {'fetal stem cell','fetal stem cells','FSCs','FSC','fSCs','fSC','fscs','fsc','Fetal','fetal'};
        abbrev = 'f';
        co = replace(co,patterns,abbrev);

        patterns = {'Immortal Stem Cells','Immortal Stem cells','Immortal stem cells','immortal stem cells',...
            'Immortal Stem Cell','Immortal Stem cell','Immortal stem cell','immortal stem cell',...
            'Immortal Cells','Immortal cells','immortal cells',...
            'Immortal','immortal'};
        abbrev = 'imm';
        co = replace(co,patterns,abbrev);

        patterns = {'adult stem cell','adult stem cells','ASCs','ASC','aSCs','aSC','ascs','asc','Adult','adult'};
        abbrev = 'a';
        co = replace(co,patterns,abbrev);

        patterns = {'brochial','Bronchial','Retinal','retinal','pigment','Pigment','pigmented','Pigmented','epithelia','epithelium','Epithelia','epithelia','epithelial','Epithelial'};
        abbrev = '';
        co = replace(co,patterns,abbrev);

        cl = LT_entry.cellLine{:};

        cv = LT_entry.cellVariant{:};

        cc = LT_entry.cellClone{:};

        ct = LT_entry.cellType{:};

        sd = LT_entry.cellSeedDate;
        % if the table is missing a seeded on date, use an array of 0s for
        % the seeding date. Else, format the seeded date to a year month
        % day format
        if cellfun(@isempty,sd) || contains(sd,'NaT')
            sd = '00000000';
        else
            sd = datestr(sd,'yyyymmdd');
        end

        % get the plate id, if it has the same seed date in the name,
        % remove it here
        pid = LT_entry.plateID;
        pid = replace(pid,sd,'');
        pid = replace(pid,'_','');
        pid = replace(pid,'-','');

        % build the complete plateID

        % build the unique tissueID
        tid = strcat(sd,'-',pid,'-',co,ct,'-',cl,cv,cc,'-',wi);

    end

    %% CONVERT MASTER to CREATE FINAL TABLE
    % take a master table and remove all but the final value for each unique
    % tissueID
    function T_final = updateFinal(T_master,master_folder)
        % get the unique tissue IDs
        fprintf('updating final summary table with new results...\n')
        [uniqueTissues,~,uniqueTissues_ic] = unique(T_master.tissueID);
        T_final = [];
        for this_unique = 1:length(uniqueTissues)
            all_entries = T_master(uniqueTissues_ic==this_unique,:);
            
            % find the maximum n for this tissueID and append to the master
            % table
            final_entry = all_entries(all_entries.n==max(all_entries.n),:);
            %final_entry = all_entries(cell2mat(all_entries.n)==max(cell2mat(all_entries.n)),:);
            % version above is needed if n is stored as cells
    
            % append
            T_final = [T_final; final_entry];
            
        end
    
        writetable(T_final,fullfile(pwd,master_folder,"final_meas_summary_table.csv"))
    end
    
    %% MASTER SUMMARY TABLE
    % read and write to a master summary table for cross reference between many
    % experiments of a variety of types
    function [T_master,master_folder] = updateMaster(T_summary)
        % make a directory for master summary data
        fprintf('-------------------------------------------------------\n')
        master_folder = 'Master Fit';
        if ~isfolder(master_folder); mkdir(master_folder); end
    
        % set the master summary table filename
        master_table_name = strcat('master_summary_table.csv');
        master_path_full = fullfile(pwd,master_folder,master_table_name);
        master_path_exists = isfile(master_path_full);
    
        % check if the file exists already
        master_path_exists = 0; % default overwrite master for now, in the future change to update only
        if master_path_exists
            fprintf('master file detected... loading file...\n')
            % define importing options
            masterOpts = setImportOptions(master_path_full);
            T_master = readtable(master_path_full,masterOpts);
    
            % matlab struggles with updating datetime arrays, store the
            % datetime array into a cell to save computation effort
            T1 = convertvars(T_master, @isdatetime, @(t) cellstr(t));
            T2 = convertvars(T_summary, @isdatetime, @(t) cellstr(t));

            % matlab struggles with updating duration arrays, store the
            % duration array into a cell to save computation effort
            T1 = convertvars(T1, @isduration, @(t) cellstr(t));
            T2 = convertvars(T2, @isduration, @(t) cellstr(t));
    
            % store all numbers as cells in case labels are primarly numeric
            % but in one lookup table they are strings. This avoids the
            % conversion between double and cell not possible error in matlab
            % table append
            T1 = convertvars(T1, @isnumeric, @(t) num2cell(t));
            T2 = convertvars(T2, @isnumeric, @(t) num2cell(t));
    
            % for each entry in the lookup table, find the matching uniqueID
            % and date time entries in the master table, if they exist
            all_summary_measIDs = T_summary.measID;
            nEntries = length(all_summary_measIDs);
            updated_indicies = 0;
            appended_indicies = 0;
            for idx_this_entry = 1:nEntries
                idx_of_match = find(strcmp(T_master.measID,all_summary_measIDs{idx_this_entry}));
                if ~isempty(idx_of_match)
                    % if there is already a matching index, update it in the
                    % table
                    updated_indicies = updated_indicies+1;
                else
                    % if no matching measID is in the master table, append it
                    % to the end of the master table.
                    appended_indicies = appended_indicies+1;
                    idx_of_match = height(T1)+1;
                end
                
                % update/append the matching data
                T1(idx_of_match,T2.Properties.VariableNames)=T2(idx_this_entry,:);
            end
              
            % update the master table
            T_master = T1;
            fprintf('updated %1.0f indicies and appended %1.0f new indicies in %s\n',updated_indicies,appended_indicies,master_table_name);
    
        else
            fprintf('no master table detected... starting new file!\n')
            T_master = convertvars(T_summary, @isdatetime, @(t) cellstr(t));
        end
    
        T_master = sortrows(T_master,'measID');
        % update/write the new master table
        writetable(T_master,master_path_full);
        pause(1);
        % return the table in a format that is consistent for the final
        % meas summary file. To enforce this, re-read the written file just
        % like how the final meas summary table will.
        masterOpts = setImportOptions(master_path_full);
        T_master = readtable(master_path_full,masterOpts);
    end
    
    %% custom function for grouping by cell age
    function seedAgeWeeks = manualAgeWeeks(seedAgeDays)
        days_offset = 3; % days
        week_duration = 7; % days
        if seedAgeDays<=days_offset
            seedAgeWeeks = 0;
        else
            seedAgeWeeks = ceil((seedAgeDays-days_offset)/week_duration);
        end
    end
    
    %% Funciton to determine the days in vitro
    function cellDaysInVitro = getDaysInVitro(this_LT_row)
        % Check if the 'cellDaysInVitro' column exists in the row
        if ~any(contains(this_LT_row.Properties.VariableNames, "cellDaysInVitro"))
            cellDaysInVitro = 0;  % Default value if the column doesn't exist
        else
            % Extract the data from the 'cellDaysInVitro' column
            columnData = this_LT_row.cellDaysInVitro;
    
            % Check if the data is already numeric
            if isnumeric(columnData)
                cellDaysInVitro = columnData;
            elseif iscell(columnData) && ~isempty(columnData{1})
                % Extract numeric part from the string using regular expression
                numericStr = regexp(columnData{1}, '\d+', 'match');
                
                % Check if any numeric part was found
                if ~isempty(numericStr)
                    % Convert the first numeric string found to a number
                    cellDaysInVitro = str2double(numericStr{1});
                else
                    % Default to 0 if no numeric part is found
                    cellDaysInVitro = 0;
                end
            else
                % Default to 0 if the cell is empty or contains non-numeric, non-string data
                cellDaysInVitro = 0;
            end
        end
    end

    %% SEED DENSITY FUNCTION
    function cellDensity = calculateCellDensity(this_LT)
        % Initialize cellDensity to NaN to handle cases where necessary data is missing
        cellDensity = NaN;
    
        % Check for 'nCellsSeeded' or 'seedDensity' columns and select the appropriate one
        if any(contains(this_LT.Properties.VariableNames, "nCellsSeeded"))
            cellCountStr = this_LT.nCellsSeeded;
        elseif any(contains(this_LT.Properties.VariableNames, "seedDensity"))
            cellCountStr = this_LT.seedDensity;
        else
            return;  % Exit the function if neither column exists
        end
    
        % Ensure cellCountStr is a character array
        if ~ischar(cellCountStr) && ~iscell(cellCountStr)
            cellCountStr = num2str(cellCountStr);
        elseif iscell(cellCountStr)
            cellCountStr = char(cellCountStr{1});  % Convert cell array to char
        end
    
        % Check for 'measArea' column existence
        if ~any(contains(this_LT.Properties.VariableNames, "measArea"))
            return;  % Exit the function if 'measArea' column does not exist
        end
    
        % Extract the numeric value and scale from cellCountStr using regular expression
        tokens = regexp(cellCountStr, '^(\d+\.?\d*)([kMmc]?)$', 'tokens');
        if isempty(tokens)
            return;  % Exit the function if parsing fails
        end
    
        numbers = str2double(tokens{1}{1});
        scale = tokens{1}{2};
    
        % Apply the scale if present
        switch scale
            case 'k'
                numbers = numbers * 1e3;
            case 'M'
                numbers = numbers * 1e6;
            case 'c'
                numbers = numbers * 1e-2;
            % Handle other cases as necessary
        end
    
        % Get the measured area, ensuring it's numeric
        measArea = this_LT.measArea;
        if iscell(measArea)
            measArea = str2double(measArea{1});
        elseif ischar(measArea)
            measArea = str2double(measArea);
        end
    
        % Calculate the cell density
        cellDensity = numbers / measArea;
    end



end