% test_Zfit.m
close all;
clc;

% Set up test file path (relative path to /tests folder)
testFile = fullfile(pwd, 'testData/EIS DC3.mpr');

% Expected nominal parameter values
nominal_p = [499, 1000, 10e-9, 3.57e3, 2.2e-6];

% Run Zfit on the test .mpr file
disp("Running Zfit test...");

% Parse MPR file
ans = parseBiologicMPR(testFile);

% Extract data
data = ans.Modules{2}.Data.DataPoints;
columns = ans.Modules{2}.Data.Columns; % Extract column names

% Identify relevant column indices
freqIdx = find(strcmp(columns, 'freq'), 1);
ReIdx   = find(strcmp(columns, 'Re(Z)'), 1);
ImIdx   = find(strcmp(columns, '-Im(Z)'), 1);

% Check for required columns
if isempty(freqIdx) || isempty(ReIdx) || isempty(ImIdx)
    error('Missing required columns: freq, ReIdx, or ImIdx not found.');
end

% Extract data
freqData = data(:, freqIdx);  
ReData = data(:, ReIdx);      
ImData = data(:, ImIdx);      

% Prepare input for Zfit
varargin = [freqData, ReData, -ImData];

% Circuit model for fitting
circuit='s(R1,s(p(R1,C1),p(R1,C1)))';  
indexes=[];
param=[100, 1000, 1e-9, 1e3, 1e-6];  % Initial guess
LB=[0,0,0,0,0];

% Run Zfit
[p, ~, ~] = Zfit(varargin, '', circuit, param, indexes, 'fitNP', LB);

% Display results
disp("Zfit Results:");
disp(p);

% Compute percentage error for each parameter
error_percent = abs((p - nominal_p) ./ nominal_p) * 100;

% Print errors
disp("Percentage Error for Each Parameter:");
disp(error_percent);

% Check if parameters are within 5% of nominal values
tolerance = 10; % percent
if all(error_percent < tolerance)
    disp("✅ Test PASSED: Parameters are within 10% of nominal values.");
else
    disp("❌ Test FAILED: Some parameters deviate more than 10%.");
end
