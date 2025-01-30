% =====================
% Import and Initialize
% =====================
format long g;
data = readmatrix("PITT1.txt"); %text file was created using ECLab Software custom txt export function

% Define the threshold for step changes (1 mV = 0.001 V)
threshold = 0.001;
EPSILON = 1e-12; %safeguard variable to prevent dividing by zero
MILLIAMPS_TO_AMPS = 1000;
FUDGE_FACTOR = 0; %-1.7e-7;
% Initialize variables
steps = {}; % Cell array to store each step
current_step = 1;
steps{current_step} = data(1, :); % Start with the first row of data
R_constant = 6e-5; % assume 500 nm radius particles for some calculations

% Initial guesses for the parameters [Lambda, tau, Delta_Q]
Lambda_tau_DeltaQ_InitialGuess =  [1, 1, 0.003];


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
    steps{current_step}(end, 5) = (data(i, 5))/MILLIAMPS_TO_AMPS + FUDGE_FACTOR;
end



% =====================
% Fit the model for all potential steps
% =====================

Lambda_tau_DeltaQ_array = zeros(110, 3); % Preallocate for B, D, Q values
resnorm_array = zeros(110, 1); % Preallocate for resnorm values
totalChargeArray = zeros(110, 1);




totalSteps = 110;

numSegments = 10; % Number of segments for the progress bar

% Initialize progress bar
progress('_start')

for i = 1:totalSteps
    % Your processing logic
    [stepOfInterest, stepOfInterest_time_current] = extractStep(steps, i);
    stepSizeInfo = size(stepOfInterest_time_current);
    stepDataSize = stepSizeInfo(1);
    Lambda_tau_DeltaQ_InitialGuess_modified = Lambda_tau_DeltaQ_InitialGuess;

    if (stepDataSize > 5)
        % Calculate total charge
        totalCharge = trapz(stepOfInterest_time_current(:, 1), stepOfInterest_time_current(:, 2)); 
        totalChargeArray(i) = totalCharge;

        % Normalize the third parameter by total charge
        Lambda_tau_DeltaQ_InitialGuess_modified(3) = totalCharge;
    end

    % Perform fitting
    [Lambda_tau_DeltaQ, resnorm] = PittFitMontella(stepOfInterest_time_current(:, 1), stepOfInterest_time_current(:, 2), Lambda_tau_DeltaQ_InitialGuess_modified);

    Lambda_tau_DeltaQ_array(i, :) = Lambda_tau_DeltaQ; % Store Lambda_tau_DeltaQ values
    resnorm_array(i) = resnorm; % Store resnorm value

    % Update the progress bar
    progress(100*i / totalSteps) % updates the progress bar. i is a percentage.

end

% Complete the progress bar
progress('_end') % ends the progress bar.



plotResults(Lambda_tau_DeltaQ_array, resnorm_array, totalChargeArray, steps, @I_model_Montella)
 