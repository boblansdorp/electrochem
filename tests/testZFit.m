classdef testZFit < matlab.unittest.TestCase
    methods (Test)
        function testZfitAccuracy(testCase)
            % Set up test file path
            testFile = fullfile(pwd, 'testData', 'EIS DC3.mpr')

            % Expected nominal parameter values
            nominal_p = [499, 1000, 10e-9, 3.57e3, 2.2e-6];

            % Parse MPR file
            ans = parseBiologicMPR(testFile);
            data = ans.Modules{2}.Data.DataPoints;
            columns = ans.Modules{2}.Data.Columns; 

            % Identify relevant column indices
            freqIdx = find(strcmp(columns, 'freq'), 1);
            ReIdx   = find(strcmp(columns, 'Re(Z)'), 1);
            ImIdx   = find(strcmp(columns, '-Im(Z)'), 1);

            % Verify required columns exist
            testCase.verifyNotEmpty(freqIdx, "Missing 'freq' column");
            testCase.verifyNotEmpty(ReIdx, "Missing 'Re(Z)' column");
            testCase.verifyNotEmpty(ImIdx, "Missing '-Im(Z)' column");

            % Extract data
            freqData = data(:, freqIdx);  
            ReData = data(:, ReIdx);      
            ImData = data(:, ImIdx);      

            % Prepare input for Zfit
            varargin = [freqData, ReData, -ImData];

            % Circuit model for fitting
            circuit='s(R1,s(p(R1,C1),p(R1,C1)))';  
            indexes=[];
            param=[100, 1000, 1e-9,1e3, 1e-6];  
            LB=[0,0,0,0,0];

            % Run Zfit WITHOUT PLOTS
            [p, ~, ~] = Zfit(varargin, '', circuit, param, indexes, 'fitNP', LB);
            
            % Close hidden figures if Zfit still generates any
            close all hidden;

            % Compute percentage error
            error_percent = abs((p - nominal_p) ./ nominal_p) * 100;

            % Verify each parameter is within 5% tolerance
            tolerance = 10; % percent
            testCase.verifyLessThan(error_percent, tolerance, "Zfit parameters exceed tolerance");
        end
    end
end
