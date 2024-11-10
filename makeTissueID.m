function [tid, lookup_tid] = makeTissueID(LT_entry, wi)
    % Validate LT_entry structure
    validateLTEntry(LT_entry);
    if nargin < 2 || isempty(wi)
        error('The wellID (wi) is missing or empty, which is required for generating a unique tissue ID.');
    end

    
    sd = formatDate(LT_entry.cellSeedDate);
    pid = cleanPlateID(LT_entry.plateID, sd);
    treat = cleanTreatment(LT_entry.treatment);

    % Reverse compatibility for cell name property. cellName used to be
    % cellLine
    if any(ismember(LT_entry.Properties.VariableNames,'cellName'))
        cn = cleanCellName(LT_entry.cellName);
    elseif any(ismember(LT_entry.Properties.VariableNames,'cellLine'))
        cn = cleanCellName(LT_entry.cellLine);
    else
        cn = 'NA';
    end



    % Assemble unique tissueID
    tid = buildTissueID(cn, sd, wi, treat);

    % make a complete tissue id for making unique properties
    lookup_tid = strcat(pid,'-',tid);
end

function validateLTEntry(LT_entry)
    requiredFields = {'cellSeedDate', 'plateID'};
    missingFields = setdiff(requiredFields, fieldnames(LT_entry));
    if ~isempty(missingFields)
        error('LT_entry is missing one or more required fields: %s', strjoin(missingFields, ', '));
    end
end

function sd = formatDate(sd)
    if isempty(sd) || contains(sd, 'NaT')
        sd = '00000000';
    else
        sd = datestr(sd, 'yyyymmdd');
    end
end

function pid = cleanPlateID(pid, sd)
    pidLower = pid;
    %pidLower = lower(pid); % Convert to lowercase for consistency
    pidCleaned = replace(pidLower, {lower(sd), '_', '-'}, '');
    pid = pidCleaned;
end

function cnUpper = cleanCellName(cn)
    cnUpper = upper(cn);
end

function treatmentClean = cleanTreatment(treatment)
    treatmentClean = lower(treatment);
end

function tid = buildTissueID(cn, sd, wi, treat)
    tid = strcat(cn, '-', sd, '-', wi);
    if ~isempty(treat{:})
        tid = strcat(tid, '-', treat);
    end
end