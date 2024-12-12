% Read in Summary Table
% SummaryTable = readtable('/Users/hannakhor/GaTech Dropbox/Hanna Khor/WPI EIS Project/Biological Sample Data/9-23-24 Testing/20240923_WPIbioExp2/20240923_WPIbioExp2_summary_table.csv');
% [row,col] = size(SummaryTable);
% oldFilename = {'Placeholder'};
% dataLocation = '/Users/hannakhor/GaTech Dropbox/Hanna Khor/WPI EIS Project/Biological Sample Data/9-23-24 Testing/20240923_WPIbioExp2/';

SummaryTable = readtable('/Users/Athena/GaTech Dropbox/Athena Chien/PBL Hanna-Athena Files/dataforBPSabstract/newSummary/20240923_WPIBIoExp2_newSettings_summary_table.csv');
[row,col] = size(SummaryTable);
oldFilename = {'Placeholder'};
dataLocation = '/Users/Athena/GaTech Dropbox/Athena Chien/WPI EIS Project/Biological Sample Data/9-23-24 Testing/20240923_WPIbioExp2/';

rmse_vec = [];
resnorm_vec = [];


for i = 1:row
    p = [SummaryTable.RblankArea(index)/SummaryTable.measArea(index) SummaryTable.R1Area(index)/SummaryTable.measArea(index) SummaryTable.C1Area(index)*SummaryTable.measArea(index) SummaryTable.R2Area(index)/SummaryTable.measArea(index) SummaryTable.C2Area(index)*SummaryTable.measArea(index)];
    filename = SummaryTable.plateID(index);
    DataTable = readtable([dataLocation char(filename)]);
    w = [DataTable.Frequency_Hz_(1:66)]*2*pi;
    zdata = funRCRC(p,w);
    
    hold on
    
    if isequal(n,row)
        nextFilename = {'Placeholder'};
    elseif (index+1) > row
        nextFilename = {'Placeholder'};
    else
        nextFilename = SummaryTable.plateID(index+1);
    end
    
    if strcmp(oldFilename, filename) | strcmp(filename, nextFilename)
        
        if strcmp(oldFilename, filename)
            occurance = occurance + 1;
            lower = 1 + ((occurance-1)*68);
            upper = 66 + ((occurance-1)*68);
        else
            occurance = 1;
            lower = 1;
            upper = 66;
        end
    
        DataTable = readtable([dataLocation char(filename)]);
        reactance_data = (DataTable.x_Z_____(lower:upper))+noise_settings(:,m);
        resistance_data = (DataTable.Z____(lower:upper))+noise_settings(:,m);
    else 
        DataTable = readtable([dataLocation char(filename)]);
        lower = 1; upper = 66;
        reactance_data = (DataTable.x_Z_____(lower:upper))+noise_settings(:,m);
        resistance_data = (DataTable.Z____(lower:upper))+noise_settings(:,m);
    end
    
    % Checking goodness of fit metrics
    reactance_data = reactance_data; 
    estimate_y = abs(zdata(:,2));
    resistance_data = resistance_data; 
    estimate_x = abs(zdata(:,1));
    rmse_y=sqrt(sum((reactance_data(:)-estimate_y(:)).^2)/(numel(reactance_data)));
    rmse_x=sqrt(sum((resistance_data(:)-estimate_x(:)).^2)/(numel(resistance_data)));
    rmse = rmse_y + rmse_x;
    % calculate resnorm
    calcZ = funRCRC(log10(p),w);
    normalization_array = max(abs(zdata), [], 1);
    calcZ_norm(:, 1) = calcZ(:, 1) / normalization_array(1); % normalization array(1) is the max of the RAW real data
    calcZ_norm(:, 2) = calcZ(:, 2) / normalization_array(2); % normalization array(1) is the max of the RAW imag data
    normZres = calcZ_norm - zdata_fitting;
    % Residuals
    resnorm = sum(normZres(:, 1).^2) + sum(normZres(:, 2).^2);
    oldFilename = filename; 
    
    % Prepping Histogram Data
    rmse_vec = [rmse_vec rmse];
    resnorm_vec = [resnorm_vec resnorm];
end

assignin('base',strcat('rmse',num2str(m)),rmse_vec);
assignin('base',strcat('resnorm',num2str(m)),resnorm_vec);

error = [{'rmse'} {'resnorm'} ];
for j = 1:length(error)
    data = [eval([error{j} '1']); eval([error{j} '2']); eval([error{j} '3'])]';
    writematrix(data,error{j});
end
