% =====================
% This function takes a txt file of a rotating disk electrode and obtains
% diffusion coefficient of the analyte in solvent
% kinetics of electrode coming soon!
% =====================


% =====================
% Import and Initialize
% =====================
format long g;
data = readmatrix("data/rde1.txt"); %text file was created using ECLab Software custom txt export function

% Define the threshold for step changes (1 mV = 0.001 V)
threshold = 0.005
MILLIAMPS_TO_AMPS = 1000;

% Initialize variables
steps = {}; % Cell array to store each step
current_step = 1;
steps{current_step} = data(1, :); % Start with the first row of data
R_constant = 6e-5; % assume 500 nm radius particles for some calculations
first_time = data(1, 1); % Set the new step's reference time

% =====================
% Process data row-by-row into voltage steps
% =====================
previousVoltageDifference = 0;

for i = 2:size(data, 1)
    % Check if the rotator potential (Column 7) changes
    % (indicating the beginning of a cycle)
    if abs(data(i, 7) - data(i - 1, 7)) > threshold
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
    %steps{current_step}(end, 5) = (data(i, 5))/MILLIAMPS_TO_AMPS + FUDGE_FACTOR;
end

totalSteps = length(steps)
figure; 
hold on;

for i = 1:totalSteps
    % Extract the data for this step
    step = extractStep(steps, i);

    % Compute average RPM
    rpm = step(:,7);
    rpm_avg = mean(rpm);

    % Extract voltage and current
    voltage = step(:,2);
    current = step(:,5)/MILLIAMPS_TO_AMPS;

    % Plot, and set 'DisplayName' for the legend
    plot(voltage, current, ...
         'DisplayName', sprintf('RPM = %.0f', rpm_avg));
end

% Automatically create legend from 'DisplayName' properties
legend show;

hold off;
xlabel('Voltage (V)');
ylabel('Current (A)');
title('CV Curves in Rotating Disk Electrode for various RPMs');
grid on;




% Define the representative voltages (in mV) at which to sample
voltagesToAnalyze = [50, 0, -50, -100, -150, -200, -250, -300];

figure; 
hold on; 
colors = lines(length(voltagesToAnalyze)); % Distinct colors for each voltage

for vIndex = 1:length(voltagesToAnalyze)
    targetVoltage = voltagesToAnalyze(vIndex);

    % Arrays to store data from each step
    I_values = zeros(totalSteps, 1);
    rpm_values = zeros(totalSteps, 1);

    for stepIndex = 1:totalSteps
        % Extract the data for this step
        step = extractStep(steps, stepIndex);

        % Column 2: Voltage, Column 5: Current, Column 7: RPM
        voltData = step(:,2);
        currData = step(:,5)/MILLIAMPS_TO_AMPS;
        rpmData  = step(:,7);

        % Find index of the point closest to targetVoltage
        [~, idxClosest] = min(abs(1000*voltData - targetVoltage));

        % Store the corresponding current & rpm
        I_values(stepIndex)   = currData(idxClosest);
        rpm_values(stepIndex) = rpmData(idxClosest);
    end

    % K–L plot: typically 1/I vs. 1/sqrt(RPM)
    plot(1./sqrt(rpm_values), 1./I_values, ...
         'o-', 'Color', colors(vIndex,:), ...
         'DisplayName', sprintf('V = %d mV', targetVoltage));
end

% Create legend from 'DisplayName' values
legend('Location','best');
xlabel('1 / \sqrt{RPM}', 'Interpreter','none');
ylabel('1 / I (A^{-1})', 'Interpreter','none');
title('Koutechy–Levich Plot');
grid on;
hold off;


figure; 
hold on; 
colors = lines(length(voltagesToAnalyze)); % Distinct colors for each voltage

% Create arrays to store slope and intercept values for each voltage
slopes = zeros(length(voltagesToAnalyze), 1);
intercepts = zeros(length(voltagesToAnalyze), 1);

for vIndex = 1:length(voltagesToAnalyze)
    targetVoltage = voltagesToAnalyze(vIndex);

    % Arrays to store data from each step
    I_values = zeros(totalSteps, 1);
    %rpm_values = zeros(totalSteps, 1);
    angularVelocity_values = zeros(totalSteps, 1);
    for stepIndex = 1:totalSteps
        % Extract the data for this step
        step = extractStep(steps, stepIndex);

        % Column 2: Voltage, Column 5: Current, Column 7: RPM
        voltData = step(:,2);
        currData = step(:,5)/MILLIAMPS_TO_AMPS;
        %rpmData  = step(:,7);
        angularVelocityData  = step(:,7)*2*pi/60; % convert from RPM to rad/s

        % Find index of the point closest to targetVoltage
        [~, idxClosest] = min(abs(1000*voltData - targetVoltage));

        % Store current and rpm
        I_values(stepIndex)   = currData(idxClosest);
        %rpm_values(stepIndex) = rpmData(idxClosest);
        angularVelocity_values(stepIndex) = angularVelocityData(idxClosest);
    end

    % The x-values for K-L plotting
    xData = 1./sqrt(angularVelocity_values);  % e.g., 1/sqrt(rad/s)
    yData = 1./I_values;          % 1/I

    % === Plot original data points ===
    hData = plot(xData, yData, 'o', ...
        'Color', colors(vIndex,:), ...
        'MarkerFaceColor', colors(vIndex,:), ...
        'DisplayName', sprintf('Data: V = %d mV', targetVoltage));

    % === Fit a linear curve y = m*x + b ===
    p = polyfit(xData, yData, 1);
    m = p(1);      % slope
    b = p(2);      % intercept

    % Store for later (if needed)
    slopes(vIndex)     = m;
    intercepts(vIndex) = b;

    % === Evaluate the fit line for plotting ===
    xFit = linspace(min(xData), max(xData), 50);  % More points for a smooth line
    yFit = polyval(p, xFit);

    % === Plot the fit line ===
    plot(xFit, yFit, '--', ...
        'Color', colors(vIndex,:), ...
        'HandleVisibility','off');  % <-- This hides the line from the legend
end

legend('Location','best');
xlabel('1 / \sqrt{rad/s}', 'Interpreter','none');
ylabel('1 / I (A^{-1})', 'Interpreter','none');
title('Koutechy–Levich Plot with Linear Fits');
grid on;
hold off;

% -- Optional: Display the slope/intercept arrays in the Command Window --
disp(table(voltagesToAnalyze(:), slopes, intercepts, ...
    'VariableNames', {'Voltage_mV','Slope','Intercept'}));




% --- Define constants ---
F   = 96485;        % Faraday's constant [C/mol]
A = pi * (0.5/2)^2;   % Electrode area in cm^2 (for a 5 mm diameter electrode)
C   = 2.5e-7;       % Bulk concentration [mol/cm^3] (dissolved oxygen)
nu  = 0.01;         % Kinematic viscosity [cm^2/s] (example)
n   = 4;            % Number of electrons transferred (oxygen reduction - could maybe be 2 also)

% slopes is the array you already computed via polyfit 
%slopes   = ...   % the array from your KL fits, slope = m
nVoltages = length(slopes);

% We'll evaluate for n=2 and n=4
nValues = [2, 4];

% Pre-allocate a matrix to store D for each n
diffusionCoeffsMatrix = zeros(nVoltages, length(nValues));

% Loop over n = 2 and n = 4
for j = 1:length(nValues)
    nCurrent = nValues(j);

    for i = 1:nVoltages
        m = slopes(i);        % slope from K–L plot
        m = abs(m);           % ensure positivity (ORR currents can be negative)
        Bprime = 1 / m;       % B' = 0.62 * n * F * A * D^(2/3) * nu^(-1/6) * C

        % Rearrange to solve for D:
        %   B' = 0.62 * n * F * A * D^(2/3) * nu^(-1/6) * C
        % => D^(2/3) = B' / (0.62 * n * F * A * C) * nu^(1/6)
        D_23 = ( Bprime / (0.62 * nCurrent * F * A * C ) ) * ( nu^(1/6) );
        D_val = (D_23)^(3/2);

        diffusionCoeffsMatrix(i, j) = D_val;
    end
end

% Create separate tables for clarity (optional)
resultsTable_n2 = table(voltagesToAnalyze(:), slopes(:), diffusionCoeffsMatrix(:,1), ...
    'VariableNames', {'Voltage_mV','Slope','D_cm2_s_n2'});
resultsTable_n4 = table(voltagesToAnalyze(:), slopes(:), diffusionCoeffsMatrix(:,2), ...
    'VariableNames', {'Voltage_mV','Slope','D_cm2_s_n4'});

disp('--- Results assuming n = 2 ---');
disp(resultsTable_n2);

disp('--- Results assuming n = 4 ---');
disp(resultsTable_n4);

% Plot both sets of results on the same figure
figure; hold on; grid on;

plot(voltagesToAnalyze, diffusionCoeffsMatrix(:,1), 'o-', ...
    'DisplayName','n = 2');
plot(voltagesToAnalyze, diffusionCoeffsMatrix(:,2), 's--', ...
    'DisplayName','n = 4');

% Add a horizontal line at 1.4e-5 cm^2/s
hLine = yline(1.4e-5, '--k', 'DisplayName','Literature Value');

xlabel('Voltage (mV)');
ylabel('Diffusion Coefficient (cm^2/s)');
title('Comparison of Diffusion Coefficients (n = 2 vs. n = 4)');
legend('Location','best');




% --- 1) Kinetic current from intercept ---
I_k_values = 1 ./ intercepts;   % [A]

% --- 2) Convert mV -> V (adjust if you have a reference offset) ---
E = voltagesToAnalyze / 1000;   % [V]

% --- 3) Tafel plot: E vs. log10(I_k) ---
figure; hold on; grid on;
logIk = log10(abs(I_k_values));
plot(logIk, E, 'o','DisplayName','Kinetic Data');

% Fit a subset or all the points
idxFit = 1:length(logIk);
pTafel = polyfit(logIk(idxFit), E(idxFit), 1);
E_fit  = polyval(pTafel, logIk(idxFit));

plot(logIk(idxFit), E_fit, '--','DisplayName','Linear Fit');
xlabel('log_{10}(I_k) [A]');
ylabel('E (V vs Ag/AgCl;)');  % specify your reference
title('Tafel Plot');
legend('Location','best');

% --- 4) Slope in V/dec ---
slope_Vdec = pTafel(1);
% Convert slope to mV/dec
slope_mVdec = slope_Vdec * 1000;

% Display in Command Window
disp('================ Tafel Analysis ================');
disp(['Tafel slope: ', num2str(slope_mVdec, '%.2f'), ' mV/dec']);
disp('================================================');
