function [p_best,resnorm_best,idx_best,info,all_results,fun] = fitZ12_20240125(p0,w,zdata,settings,Rt_max,Ct_max,R_solution_max,cross_section)
% this function is used to fit EIS data to an RCRC circuit. Specifically,
% for epithelia.
%
% Inputs: 
% fun = reference to a function file that contains the model equation for
%   fitting. Generally, this should be a function with 4 inputs to match
%   the 4 parameters in the RCRC circuit
% p0 = mx4 array of the initial guesses for each parameter in the impedance
%   model. The data should be in the following form [R1 C1 R2 C2]. m is the
%   number of initial guesses to try. 
% w = the frequency at each impedance measurment (rad/s).
% zdata = the actual impedance data at each measured frequency (w).
% settings = struct containing properties about the fit process. Including
%   fit_circuit
% Output:
% pfit = the model's best fit for the 4 input parameters (p). The data is
%   in the following form [R1 C1 R2 C2];
% pmin = the minimum bound for each parameter during the fit
% pmax = the maximum bound for each parameter during the fit

%% fit_circuit function declaration
% Define the function to fit based on the selected circuit model
switch settings.fit_circuit
    case 'RCRC'
        fun = @funRCRC;

    case 'RCRC_penalty'
        % Check if settings contain target_tau_ratio and lambda_penalty
        if ismember('target_tau_ratio',settings.Properties.VariableNames) && ismember('lambda_penalty',settings.Properties.VariableNames)
            % If both parameters are available, pass them to the function
            fun = @(p, w) funRCRC_penalty(p, w, settings.target_tau_ratio, settings.lambda_penalty);
        else
            % Default values or a warning could be issued if parameters are missing
            default_lambda = 1e-1; % Example default value
            default_tau_ratio = 1; % Default tau ratio (1:1)
            fun = @(p, w) funRCRC_penalty(p, w, default_tau_ratio, default_lambda);
            warning('Default values for lambda_penalty and target_tau_ratio are used as they are not specified in settings.');
        end
        % Initialize the third column for penalty
        zdata(:,3) = zeros(size(zdata(:,1)));

    case 'RCRCalpha'
        fun = @funRCRCalpha;

    case 'RCRCdistributed'
        fun = @funRCRCdistributed;

    case 'RCRCpRC'
        fun = @funRCRCpRC;

    otherwise
        warning('Invalid fit_circuit selected in function fitZ12!');
end

%% Data normalization
% Normalize the function so that resistance and reactance contribute equally
if settings.normalize_data
    % Calculate normalization factors based on the maximum absolute value in each column of zdata
    normalization_array = max(abs(zdata), [], 1);
    % Replace 0 with 1 to avoid division by zero
    normalization_array(normalization_array == 0) = 1;
    
    % Adjust the function handle to apply normalization
    original_fun = fun; % Store the original function handle
    fun = @(p, w) original_fun(p, w) ./ normalization_array; % Apply normalization within the function call
    zdata_fitting = zdata ./ normalization_array; % Normalize zdata for fitting
else
    normalization_array = ones(1, size(zdata, 2));
    zdata_fitting = zdata; % Use unnormalized data if normalization is not enabled
end

%% Upper and lower bounds
% define the maximum and minimum possible values for the function
% parameters
rsolmax = R_solution_max/cross_section;
rtmax = Rt_max/cross_section; % most likely
ctmax = Ct_max*cross_section;
switch settings.fit_circuit
    case 'RCRCpRC'
        info.pmin = [0 0 0 0 0 0 0 0];
        info.pmax = [rsolmax rtmax ctmax rtmax ctmax rtmax ctmax 1];
    otherwise
        info.pmin = [0 0 0 0 0 0];
        info.pmax = [rsolmax rtmax ctmax rtmax ctmax 1];
end
info.pmin = log10(info.pmin);
info.pmax = log10(info.pmax);

%% FIT options
info.opts = optimoptions('lsqcurvefit',...
    'FunctionTolerance',1e-9,...
    'StepTolerance',1e-9,...
    'MaxFunctionEvaluations',3e3,...
    'MaxIterations', 3e3,...
    'OptimalityTolerance', 5e-7, ...
    'Display','off');

%% INIT FIT
% get the number of intitial gueses to try to fit to the data.
[~,m] = size(p0);

% g for guess
pfit = cell(m,1);
resnorm = cell(m,1);
residual = cell(m,1);
exitflag = cell(m,1);
output = cell(m,1);
lambda = cell(m,1);
jacobian = cell(m,1);

%% FIT
parfor g = 1:m
% for g = 1:m
    this_p0 = log10(p0{g}); % use the log to equalize the contribution of the capacitance and resistances to the overall fit.
    [pfit{g},resnorm{g},residual{g},exitflag{g},output{g},lambda{g},jacobian{g}] = lsqcurvefit(fun,this_p0,w,zdata_fitting,info.pmin,info.pmax,info.opts);
    pfit{g} = 10.^(pfit{g}); % return the exponentials back to a standard form 
end

%% SELECT BEST RESULT
% get the optimal result and store it in the output
p_best = pfit{1}; 
resnorm_best = resnorm{1}; 
idx_best = 1;
    
% use the resnorm to evaluate the quality of the fit. Therefore, select
% the smallest resnorm as the best fit for the data.
% optional add: && exitflag{g}>0 to if statment inside loop
for g = 1:m
    if resnorm{g} <= resnorm_best 
        p_best = pfit{g};
        resnorm_best = resnorm{g};
        idx_best = g;
    end
end

%% RETURN ALL RESULTS
all_results = table(pfit,resnorm,residual,exitflag,output,lambda,jacobian);

% un-normalize the function file
if settings.normalize_data
    fun = @(p,w)fun(p,w).*normalization_array;
end

% exit flag values
% 1. Function converged to a solution x.
% 
% 2. Change in x was less than the specified tolerance.
% 
% 3. Change in the residual was less than the specified tolerance.
% 
% 4. Relative magnitude of search direction was smaller than the step tolerance.
% 
% 0. Number of iterations exceeded options.MaxIterations or number of function evaluations exceeded options.MaxFunctionEvaluations.
% 
% -1. A plot function or output function stopped the solver.
% 
% -2. Problem is infeasible: the bounds lb and ub are inconsistent.
