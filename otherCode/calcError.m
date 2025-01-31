%% gaussian noise
% Read in Summary Table
% SummaryTable = readtable('/Users/hannakhor/GaTech Dropbox/Hanna Khor/WPI EIS Project/Biological Sample Data/9-23-24 Testing/20240923_WPIbioExp2/20240923_WPIbioExp2_summary_table.csv');
% [row,col] = size(SummaryTable);
% oldFilename = {'Placeholder'};
% dataLocation = '/Users/hannakhor/GaTech Dropbox/Hanna Khor/WPI EIS Project/Biological Sample Data/9-23-24 Testing/20240923_WPIbioExp2/';

% SummaryTable = readtable('/Users/Athena/GaTech Dropbox/Athena Chien/PBL Hanna-Athena Files/dataforBPSabstract/newSummary/20240923_WPIBioExp2_newSettings_summary_table.csv');
SummaryTable = readtable('/Users/Athena/GaTech Dropbox/Athena Chien/PBL Hanna-Athena Files/dataforBPSabstract/newSummary/20240923_WPIBioExp2test_subsetnewSettings_summary_table.csv');
[row,col] = size(SummaryTable);
oldFilename = {'Placeholder'};
dataLocation = '/Users/Athena/GaTech Dropbox/Athena Chien/WPI EIS Project/Biological Sample Data/9-23-24 Testing/20240923_WPIbioExp2/';

mae_vec = []; rmse_vec = []; mape_vec = []; chi_sq_vec = []; 
chi_sq_norm_vec = []; resnorm_vec = [];

for index = 1:row % read through every row of summary table
    num_freqs = SummaryTable.nFreq(index);

    % Calculate the impedance points from the model
    p = [SummaryTable.RblankArea(index)/SummaryTable.measArea(index) SummaryTable.R1Area(index)/SummaryTable.measArea(index) SummaryTable.C1Area(index)*SummaryTable.measArea(index) SummaryTable.R2Area(index)/SummaryTable.measArea(index) SummaryTable.C2Area(index)*SummaryTable.measArea(index)];
    % if p has NANs - skip to next iteration of the for loop - this applies
    % to blank TW measurements
    if any(isnan(p))
        continue;
    end
    filename = SummaryTable.plateID(index);
    DataTable = readtable([dataLocation char(filename)]);
    w = [DataTable.Frequency_Hz_(1:num_freqs)]*2*pi;
    zdata = funRCRC(p,w); % ZDATA HERE IS THE ESTIMATED VALUES
    % plot(zdata(:,1),-zdata(:,2),'Color', [.7 .7 .7],'LineWidth',3)

    % Determine if the next data sample is in the same data file
    if isequal(index,row)
        nextFilename = {'Placeholder'};
    elseif (index+1) > row
        nextFilename = {'Placeholder'};
    else
        nextFilename = SummaryTable.plateID(index+1);
    end

    if strcmp(oldFilename, filename) || strcmp(filename, nextFilename)
        
        if strcmp(oldFilename, filename)
            occurance = occurance + 1;
            lower = 1 + ((occurance-1)*(num_freqs+2));
            upper = num_freqs + ((occurance-1)*(num_freqs+2));
        else
            occurance = 1;
            lower = 1;
            upper = num_freqs;
        end
  
        % Read in the raw data to compare with model
        DataTable = readtable([dataLocation char(filename)]);
        reactance_data = (DataTable.x_Z_____(lower:upper));
        resistance_data = (DataTable.Z____(lower:upper));
             
    else % next sweep in a new datafile
        DataTable = readtable([dataLocation char(filename)]);
        lower = 1; upper = num_freqs;
        reactance_data = (DataTable.x_Z_____(lower:upper));
        resistance_data = (DataTable.Z____(lower:upper));
    end

    % Checking goodness of fit metrics
    reactance_data = -1*reactance_data; % don't use absolute values here, use negative reactance for both
    estimate_X = zdata(:,2);
    estimate_R = zdata(:,1);
    nn = length(zdata);
    resid_R = resistance_data(:)-estimate_R(:);
    resid_X = reactance_data(:)-estimate_X(:);
    rmse_x=sqrt(sum(resid_X.^2)/nn);
    rmse_r=sqrt(sum(resid_R.^2)/nn);
    rmse = rmse_r + rmse_x;
    mae_x=(1/nn)*(abs(sum(resid_X)));
    mae_r=(1/nn)*(abs(sum(resid_R)));
    mae = mae_r + mae_x;
    mape_x=(100/nn)*(sum((resid_X./reactance_data(:))));
    mape_r=(100/nn)*(sum((resid_R./resistance_data(:))));
    mape = mape_r + mape_x;
    chi_sq_r=sum((resid_R.^2)./estimate_R(:));
    chi_sq_x=sum((resid_X.^2)./estimate_X(:));
    chi_sq=chi_sq_r+chi_sq_x;

    chi_sq_r_norm=sum( (resid_R./max(resistance_data(:))).^2 ./( estimate_R(:)/max(abs(resistance_data(:))) ) );
    chi_sq_x_norm=sum( (resid_X./max(reactance_data(:))).^2 ./( estimate_X(:)/max(abs(reactance_data(:))) ) ); 
    chi_sq_norm=chi_sq_r_norm+chi_sq_x_norm;
    %resnorm = SummaryTable.resnorm(index);
    normalizationArray = [max(abs(resistance_data)), max(abs(reactance_data))];
    resnorm_zdata = [resistance_data, reactance_data];
    resnorm_zdata_fitting = resnorm_zdata ./ normalizationArray;
    calcZ = originalfunRCRC(log10(p),w);
    calcZ_norm(:, 1) = calcZ(:, 1) ./ normalizationArray(1); % normalization array(1) is the max of the RAW real data
    calcZ_norm(:, 2) = calcZ(:, 2) ./ normalizationArray(2); % normalization array(1) is the max of the RAW imag data
    normZres = calcZ_norm - resnorm_zdata_fitting;
    % Residuals
    resnorm = sum(normZres(:, 1).^2) + sum(normZres(:, 2).^2);
    oldFilename = filename; 

    % Prepping Histogram Data
    mae_vec = [mae_vec mae];
    rmse_vec = [rmse_vec rmse];
    mape_vec = [mape_vec mape];
    chi_sq_vec = [chi_sq_vec chi_sq];
    resnorm_vec = [resnorm_vec resnorm];
    chi_sq_norm_vec = [chi_sq_norm_vec chi_sq_norm];
end

% output the error values as text files
% error = [{'mae_vec'} {'rmse_vec'} {'mape_vec'} {'chi_sq_vec'} {'resnorm_vec'} {'chi_sq_norm_vec'}];
data = [];
error = [{'mae_vec'} {'resnorm_vec'} ];
for j = 1:length(error)
    data = [data, eval([error{j}])'];
end
dataWfilenames = [SummaryTable.plateID, num2cell(data)];
writecell(dataWfilenames, 'LowCaExp4_RCRC_errors.txt');
