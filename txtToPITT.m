% =====================
% Import and Initialize
% =====================

format long g;
arr = readmatrix("PITT1.txt"); %text file was created using ECLab Software custom txt export function
data = arr;

% Define the threshold for step changes (1 mV = 0.001 V)
threshold = 0.001;

% Initialize variables
steps = {}; % Cell array to store each step
current_step = 1;
steps{current_step} = data(1, :); % Start with the first row of data

R_constant = 5e-5; % assume 500 nm radius particles
B_D_Q_InitialGuess =  [1.08, 4e-10, 0.00025858];
epsilon = 1e-12;


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
end


% do a test!
elementToObtain = 25;
[stepOfInterest, lsqcurvefit] = extractStep(steps, elementToObtain);


% Define the equation for I(t) as an anonymous function with safeguards
%% defined according to Li as: 
I_t_Li = @(params, t) -3 * params(2) * params(3) / max(R_constant^2, epsilon) * ( ...
    (-params(1)/(params(1) - 1)) * (1 - erfc(R_constant ./ sqrt(max(params(2) * t, epsilon)))) + ...
    (params(1)^2 / (params(1) - 1)) * exp(max(params(2) * t / (R_constant^2) * (params(1) - 1)^2, -700)) .* ...
    erfc((params(1) - 1) * sqrt(max(params(2) * t / (R_constant^2), epsilon)) + ...
    (R_constant ./ sqrt(max(params(2) * t, epsilon))) + (params(1) - 1) * sqrt(max(params(2) * t / (R_constant^2), epsilon))) ...
    );



% =====================
% Iterate over all data
% =====================

for i=1:1:60 % over all data (just 60 for now)
    [stepOfInterest, stepOfInterest_time_current] = extractStep(steps, i);
    [B_D_Q, resnorm] = PittFit(stepOfInterest_time_current(:,1), stepOfInterest_time_current(:,2), R_constant, B_D_Q_InitialGuess, I_t_Li);
    B_D_Q_array(i,:) = B_D_Q;
end

figure;
plot(B_D_Q_array(:,2))
hold on;
ylabel('Diffusion coefficient (cm^{2}/s)');
xlabel('step number');
title('Diffusion coefficient as a function of step number');
grid on;


%[B_D_Q, resnorm] = PittFit(stepOfInterest_time_current(:,1), stepOfInterest_time_current(:,2), R_constant, B_D_Q_InitialGuess, I_t_Li);

indexToPlot = 50

% Plot the results from indexToPlot
I_fitted = I_t_Li(B_D_Q, stepOfInterest_time_current(:,1));
t_inv_sqrt = 1 ./ sqrt(stepOfInterest_time_current(:,1)); % Compute t^{-1/2} for plotting

figure;
plot(t_inv_sqrt, stepOfInterest_time_current(:,2), 'ro', 'MarkerSize', 8, 'DisplayName', 'Data'); % Original data
hold on;
plot(t_inv_sqrt, I_fitted, 'b-', 'LineWidth', 2, 'DisplayName', 'Fitted Curve'); % Fitted curve
xlabel('t^{-1/2} (s^{-1/2})');
ylabel('Current (A)');
title('Transient Current vs. t^{-1/2}');
xlim([0, 3]); % Limit x-axis range to 0 to 1.4
grid on;
legend('show');




% =====================
% Helper functions
% =====================


function [B_D_Q, resnorm] = PittFit(t_data, I_data, R_constant, initial_guess, I_t_Li)
    % PittFit: Fits the given data using the Pitt and a model for I(t)
    % 
    % Inputs:
    %   t_data        - Time data (column vector)
    %   I_data        - Current data (column vector)
    %   R_constant    - Constant value for R
    %   initial_guess - Initial guesses for [B, D, Q]
    %
    % Outputs:
    %   B_D_Q      - an array with Fitted value of parameters B, D, Q
    %   resnorm  - Residual norm from the fitting process

    % Use lsqcurvefit to optimize parameters
    % options = optimset('Display', 'iter'); % Show optimization progress
    options= optimset('Display', 'none'); % don't show progress

    [B_D_Q, resnorm] = lsqcurvefit(@(params, t) I_t_Li(params, t), initial_guess, t_data, I_data, [], [], options);
    
    % Extract fitted parameters
    B = B_D_Q(1);
    D = B_D_Q(2);
    Q = B_D_Q(3);
    
    disp(['Fitted Parameters: B = ', num2str(B), ', D = ', num2str(D), ', Q = ', num2str(Q)]);
    disp(['Constant Parameter: R = ', num2str(R_constant)]);

end



function [stepOfInterest, stepOfInterest_time_current] = extractStep(steps, elementToObtain)
% Access specific steps

    if elementToObtain <= length(steps) && ~isempty(steps{elementToObtain})
        stepOfInterest = steps{elementToObtain}; % Access the specific step
        stepOfInterest_time_current = [stepOfInterest(:,1), stepOfInterest(:,5)];
        stepOfInterest_time = stepOfInterest(:,1);
        t_inv_sqrt = stepOfInterest_time.^(-0.5);
    
        %figure;
        %plot(t_inv_sqrt, stepOfInterest_time_current(:,2), '.', 'MarkerSize', 12); % Plot points only
    
        %xlabel('t^{-1/2} (s^{-1/2})');
        %ylabel('Current (A)');
        %title('Transient Current vs. t^{-1/2}');
        %grid on;
    
        % Set graph range
        %xlim([0, 1.4]); % Limit x-axis range to 0 to 1.4
    
        % Plotting semilogy for step_50_time_current
        %figure;
        %semilogy(stepOfInterest_time_current(:,1), stepOfInterest_time_current(:,2));
    
        %disp('Step of interest  stepOfInterest_time_current:');
        %disp(stepOfInterest_time_current);
    else
        fprintf('Step %d does not exist.\n', elementToObtain);
    end
end
