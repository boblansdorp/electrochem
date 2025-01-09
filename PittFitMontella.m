function [bestParams, bestResnorm] = PittFitMontella(t_data, I_data, Lambda_tau_DeltaQ_InitialGuess)
% Fits the transient current I(t) to the Montella model using lsqcurvefit.
% First checks residues for multiple initial guesses and uses the best one.

    % -----------------------------
    % 1) Quick checks and data cleanup
    % -----------------------------
    idxZeroOrNeg = (t_data <= 0);
    if any(idxZeroOrNeg)
        t_data(idxZeroOrNeg) = 1e-9; % Shift tiny or negative times
    end

    idxBad = isnan(t_data) | isnan(I_data) | isinf(t_data) | isinf(I_data);
    if any(idxBad)
        t_data(idxBad) = [];
        I_data(idxBad) = [];
    end

    if length(t_data) < 3
        bestParams = Lambda_tau_DeltaQ_InitialGuess;
        bestResnorm = NaN;
        return;
    end

    Lambda_tau_DeltaQ_InitialGuess(1:3) = max(Lambda_tau_DeltaQ_InitialGuess(1:3), [1e-9, 1e-9, 1e-9]);

    % -----------------------------
    % 2) Generate Initial Guess Array
    % -----------------------------
    scalingFactors = [1/4, 1/2, 1, 2, 4]; % Define scaling factors
    InitialGuess_Array = []; % Initialize array to store combinations

    % Generate all combinations using nested for loops
    for scale1 = scalingFactors
        for scale2 = scalingFactors
            for scale3 = scalingFactors
                scaledGuess = Lambda_tau_DeltaQ_InitialGuess .* [scale1, scale2, scale3];
                InitialGuess_Array = [InitialGuess_Array; scaledGuess];
            end
        end
    end

    % -----------------------------
    % 3) Identify Best Initial Guess by Residue
    % -----------------------------
    numGuesses = size(InitialGuess_Array, 1); % Number of initial guesses
    bestResnorm = Inf; % Initialize best residual norm
    bestInitialGuess = NaN(1, 3); % Initialize best initial guess

    for i = 1:numGuesses
        initialGuess = InitialGuess_Array(i, :);

        % Compute model values using current initial guess
        modelValues = I_model_Montella(initialGuess, t_data);

        % Calculate residue as the sum of squared differences
        residue = sum((I_data - modelValues).^2);

        % Update best guess if this residue is smaller
        if residue < bestResnorm
            bestResnorm = residue;
            bestInitialGuess = initialGuess;
        end
    end

    % -----------------------------
    % 4) Define Bounds for [Lambda, tau, Delta_Q]
    % -----------------------------
    lb = [1e-4, 1e-4, 1e-9];
    ub = [1e4, 1e4, 1e6];

    % -----------------------------
    % 5) Perform Optimization Using Best Initial Guess
    % -----------------------------
    options = optimoptions('lsqcurvefit', ...
        'Algorithm', 'trust-region-reflective', ...
        'Display', 'off', ...
        'MaxIterations', 300, ...
        'MaxFunctionEvaluations', 2000, ...
        'FunctionTolerance', 1e-10, ...
        'StepTolerance', 1e-10);

    try
        [bestParams, bestResnorm] = lsqcurvefit( ...
            @(params, t) I_model_Montella(params, t), ...
            bestInitialGuess, ...
            t_data, ...
            I_data, ...
            lb, ...
            ub, ...
            options);
    catch ME
        warning('PittFitMontella:SolverError', 'Error during lsqcurvefit: %s', ME.message);
        bestParams = [NaN, NaN, NaN];
        bestResnorm = NaN;
    end
end