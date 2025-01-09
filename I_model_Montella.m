function I = I_model_Montella(params, t)
% Computes the current response based on the Montella model.
    I = compute_I_Montella(params, t); 
end

function I = compute_I_Montella(params, t)
% Computes the transient current for the Montella model.
    Lambda = params(1);
    tau = params(2);
    Delta_Q = params(3);
    MAXIMUM_N = 100;

    I = zeros(size(t));
    for n = 1:MAXIMUM_N
        b_n = solve_transcendental(n, Lambda);
        I = I + (6 * Delta_Q / tau) * (Lambda^2 / (Lambda^2 - Lambda + b_n^2)) .* exp(-b_n^2 .* t / tau);
    end
end

function b_n = solve_transcendental(n, Lambda)
% Solves b*cot(b) + Lambda - 1 = 0 for b using fzero.
    if Lambda <= 0
        error('Lambda must be positive.');
    end

    f = @(b) b .* cot(b) + Lambda - 1;
    b_guess = n * pi;
    lb = b_guess - pi + 1e-9;
    ub = b_guess - 1e-9;
    
    % Use correct fzero interval syntax
    try
        b_n = fzero(f, [lb, ub]);
    catch
        warning('solve_transcendental:SolverError', 'Failed for n = %d, Lambda = %.4f', n, Lambda);
        b_n = NaN;
    end
end
