% =====================
% Import and Initialize
% =====================
format long g;
arr = readmatrix("PITT1.txt"); %text file was created using ECLab Software custom txt export function
data = arr;

% Define the threshold for step changes (1 mV = 0.001 V)
threshold = 0.001;
EPSILON = 1e-12; %safeguard variable to prevent dividing by zero
MILLIAMPS_TO_AMPS = 1000;

% Initialize variables
steps = {}; % Cell array to store each step
current_step = 1;
steps{current_step} = data(1, :); % Start with the first row of data
R_constant = 6e-5; % assume 500 nm radius particles for some calculations

% Initial guesses for the parameters [Lambda, tau, Delta_Q]
Lambda_tau_DeltaQ_InitialGuess =  [0.5, 1, 0.003];


% =====================
% Process data row-by-row into voltage steps
% =====================

for i = 2:size(data, 1)
    % Check if the difference in potential (Column 2) exceeds the threshold
    if abs(data(i, 2) - data(i - 1, 2)) > threshold
        % Start a new step
        current_step = current_step + 1;
        steps{current_step} = []; % Initialize new step
        first_time = data(i, 1); % Set the new step's reference time
    end

    % Append the current row to the current step
    steps{current_step} = [steps{current_step}; data(i, :)];

    % Adjust time for the current row in the step
    steps{current_step}(end, 1) = data(i, 1) - first_time;

    % Adjust the current from mA to Amps:
    steps{current_step}(end, 5) = (data(i, 5))/MILLIAMPS_TO_AMPS;
end


% =====================
% Define the model
% =====================
% as described in Montella, Claude. "Apparent diffusion coefficient of intercalated species measured with PITT: A simple formulation." Electrochimica acta 51.15 (2006): 3102-3111.
% and J. Phys. Chem. Lett. 2022, 13, 3165−3172
I_model_Montella = @(params, t) compute_I_Montella(t, params); 

% =====================
% Fit the model for all potential steps
% =====================

Lambda_tau_DeltaQ_array = zeros(110, 3); % Preallocate for B, D, Q values
resnorm_array = zeros(110, 1); % Preallocate for resnorm values
totalChargeArray = zeros(110, 1); 

for i = 1:110 % Iterate over all steps
    [stepOfInterest, stepOfInterest_time_current] = extractStep(steps, i);
    stepSizeInfo = size(stepOfInterest_time_current);
    stepDataSize = stepSizeInfo(1);
    Lambda_tau_DeltaQ_InitialGuess_modified = Lambda_tau_DeltaQ_InitialGuess;

    if (stepDataSize > 5) %only normalize if there are enough points to do so
        %calculate total charge by integrating using trapezoidal method
        totalCharge = trapz(stepOfInterest_time_current(:,1), stepOfInterest_time_current(:,2)); 
        % fprintf('totalCharge: %g\n', totalCharge);
        totalChargeArray(i) = totalCharge;
        % Normalize the third parameter by total charge
        Lambda_tau_DeltaQ_InitialGuess_modified(3) = totalCharge;
    end
    %[Lambda_tau_DeltaQ, resnorm] = PittFit(stepOfInterest_time_current(:,1), stepOfInterest_time_current(:,2), ...
     %                          R_constant, Lambda_tau_DeltaQ_InitialGuess_modified, I_t_Li);
    [Lambda_tau_DeltaQ, resnorm] = PittFitMontella(stepOfInterest_time_current(:,1), stepOfInterest_time_current(:,2), ...
                               Lambda_tau_DeltaQ_InitialGuess_modified, I_model_Montella);

    Lambda_tau_DeltaQ_array(i, :) = Lambda_tau_DeltaQ; % Store Lambda_tau_DeltaQ values
    resnorm_array(i) = resnorm; % Store resnorm value
end

plotResults(Lambda_tau_DeltaQ_array, resnorm_array, totalChargeArray, steps, I_model_Montella, indexesToPlot)
    


%take ratio of fit charge to expected charge from integrating current
%during step (should eb equal to 1)
sanityCheck = Lambda_tau_DeltaQ_array(:,3)./totalChargeArray;
figure;
plot(sanityCheck)
xlabel('Step #'); % Label x-axis
ylabel('Sanity ratio (nominally 1))'); % Label y-axis
title('Sanity check of measured charge)');

% =====================
% Helper functions
% =====================

function [Lambda_tau_DeltaQ, resnorm] = PittFitMontella(t_data, I_data, Lambda_tau_DeltaQ_InitialGuess, I_model_Montella)
% PittFit fits the transient current I(t) to the Li equation using lsqcurvefit
% with bounds to prevent singularities (such as B=1) or non-physical parameters.
%
% Usage:
%   [Lambda_tau_DeltaQ, resnorm] = PittFitMontella(t_data, I_data, ...
%           Lambda_tau_DeltaQ_InitialGuess, I_model_Montella)
%
% Inputs:
%   t_data  : time vector (column) for the step
%   I_data  : current vector (column) for the step
%   Lambda_tau_DeltaQ_InitialGuess : 1x3 initial guess for [Lambda, tau, DeltaQ]
%   I_model_Montella  : function handle for I(t) = I_model_Montella(params, t)
%
% Outputs:
%   Lambda_tau_DeltaQ   : best-fit parameters [B, D, Q]
%   resnorm : residual norm from lsqcurvefit

    % -----------------------------
    % 1) Quick checks and data cleanup
    % -----------------------------
    % Ensure time data is strictly > 0 so we don’t get 1/sqrt(0). 
    % Shift tiny or negative times if needed.
    idxZeroOrNeg = (t_data <= 0);
    if any(idxZeroOrNeg)
        t_data(idxZeroOrNeg) = 1e-9;  % Shift them up slightly
    end

    % Remove any NaN/Inf from inputs
    idxBad = isnan(t_data) | isnan(I_data) | isinf(t_data) | isinf(I_data);
    if any(idxBad)
        t_data(idxBad) = [];
        I_data(idxBad) = [];
    end

    % If there's not enough data left to fit, return defaults
    if length(t_data) < 3
        Lambda_tau_DeltaQ   = Lambda_tau_DeltaQ_InitialGuess; 
        resnorm = NaN;
        return;
    end

    % If needed, also ensure initial guess won't blow up (e.g., B ~ 1).
    % This is just a tiny shift if B=1 exactly.
    if abs(Lambda_tau_DeltaQ_InitialGuess(1) - 1) < 1e-5
        Lambda_tau_DeltaQ_InitialGuess(1) = 1.001; % or 0.999, whichever side is physically valid
    end

    % Also ensure positivity on D, Q
    if Lambda_tau_DeltaQ_InitialGuess(2) <= 0
        Lambda_tau_DeltaQ_InitialGuess(2) = 1e-14;  % small positive
    end
    if Lambda_tau_DeltaQ_InitialGuess(3) <= 0
        Lambda_tau_DeltaQ_InitialGuess(3) = 1e-9;   % small positive
    end

    % -----------------------------
    % 2) Define bounds for [Lambda, tau, Q]
    % -----------------------------
    % The bounds below are EXAMPLES. Adjust them for your system:
    %   - B != 1, but can be near 0.5 to 2.0 (example)
    %   - D in [1e-14, 1e6] (example)
    %   - Q in [1e-9, 1e-2]  (example)
    lb = [0.01,   1e-14, 1e-9 ];
    ub = [1e3,   1e6 , 1e6 ];

    % now let's fit only for Lambda, tau, taking deltaQ as fixed
    %lb = [1e-4,   1e-6];
    %ub = [1e4,   1e6 ];
    % Update initial guess to only include Lambda and tau
    %Lambda_tau_InitialGuess = Lambda_tau_DeltaQ_InitialGuess(1:2); % Extract initial guess for Lambda and tau
    %fixedDeltaQ = Lambda_tau_DeltaQ_InitialGuess(3); % Example value, replace with the desired fixed value

    % Update the objective function to include the fixed DeltaQ
    %I_model_Montella_FixedDeltaQ = @(params, t) compute_I_Montella_FixedDeltaQ(params, t, fixedDeltaQ);


    % In particular, B=1 leads to a singularity in your I_t_Li expression
    % so we do not allow it to land exactly at 1. If your B can
    % physically be >1 or <1, define a range that excludes 1 or includes it
    % in a way that your code can handle. (You might re-factor the equation
    % if B=1 is physically possible.)

    % -----------------------------
    % 3) Define lsqcurvefit Options
    % -----------------------------
    options = optimoptions('lsqcurvefit', ...
        'Algorithm',          'trust-region-reflective', ...
        'Display',            'off',           ...  % or 'iter', 'iter-detailed'
        'MaxIterations',      300,            ...
        'MaxFunctionEvaluations', 2000,       ...
        'FunctionTolerance',  1e-9,          ...
        'StepTolerance',      1e-9);

    % -----------------------------
    % 4) Call lsqcurvefit
    % -----------------------------
    % If your function is well-behaved, the solver should not produce NaNs.
    % However, if it does, the bounding constraints and checks above
    % help keep it in a valid region.
    %
    % Because we wrote:  I_t_Li = @(params, t) ...
    % we simply pass it in. Notice we specify lb, ub for bounding.
    %
    try
        % [Lambda_tau, resnorm] = lsqcurvefit( ...
        %     @(params, t) I_model_Montella_FixedDeltaQ(params, t), ...
        %     Lambda_tau_InitialGuess, ...
        %     t_data, ...
        %     I_data, ...
        %     lb, ...
        %     ub, ...
        %     options);
        % % Add the fixed DeltaQ back to the output for consistency
        % Lambda_tau_DeltaQ = [Lambda_tau, fixedDeltaQ];
        % 

        [Lambda_tau_DeltaQ, resnorm] = lsqcurvefit( ...
            @(params, t) I_model_Montella(params, t), ...
            Lambda_tau_DeltaQ_InitialGuess, ...
            t_data, ...
            I_data, ...
            lb, ...
            ub, ...
            options);
        

    catch ME
        % If any unexpected numerical error occurs, handle it gracefully
        warning('PittFit:SolverError','lsqcurvefit encountered an error: %s', ME.message);
        Lambda_tau_DeltaQ   = [NaN, NaN, NaN];
        resnorm = NaN;
    end
    
    % -----------------------------
    % 5) (Optional) Final checks
    % -----------------------------
    % If the solver solution ended near a singularity, you could warn:
    if abs(Lambda_tau_DeltaQ(1) - 1) < 1e-5
        warning('PittFit:NearSingularity','B is near 1 => Potential singularity in I_model_Montella');
    end
    
end


% =====================
% Helper functions
% =====================

function I = compute_I_Montella_FixedDeltaQ(params, t, fixedDeltaQ)
    Lambda = params(1);
    tau = params(2);
    Delta_Q = fixedDeltaQ; % Use the fixed value
    I = compute_I_Montella(t, [Lambda, tau, Delta_Q]);
end


function I = compute_I_Montella(t, params)  % from J. Phys. Chem. Lett. 2022, 13, 3165−3172
    Lambda = params(1);
    tau = params(2);
    Delta_Q = params(3);
    MAXIMUM_N = 100; % could be 100... but 10 is faster!
    I = zeros(size(t)); % Initialize output
    for n = 1:MAXIMUM_N % Truncate summation at n = MAXIMUM_N for numerical stability
        b_n = solve_transcendental(n, Lambda); % Solve transcendental equation
        I = I + (6 * Delta_Q / tau) * (Lambda^2 / (Lambda^2 - Lambda + b_n^2)) .* exp(-b_n^2 .* t / tau);
    end
end

function b_n = solve_transcendental(n, Lambda)
    % Solve b*cot(b) + Lambda - 1 = 0 for b using fzero
    if Lambda <= 0
        error('Lambda must be positive for valid solutions.');
    end
    
    f = @(b) real(b .* cot(b) + Lambda - 1); % Ensure real output


    b_guess = n * pi; % Initial guess for b

    % Safeguard the interval to avoid singularities
    lower_bound = b_guess - pi + 1e-9;
    upper_bound = b_guess - 1e-9;

    % Try solving the equation, handle failure gracefully
    try
        b_n = fzero(f, [lower_bound, upper_bound]); % Find root
    catch ME
        warning('fzero failed for n = %d, Lambda = %.4f: %s', n, Lambda, ME.message);
        b_n = NaN; % Assign NaN to indicate failure
    end
end
