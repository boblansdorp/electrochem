
function plotResults(Lambda_tau_DeltaQ_array, resnorm_array, totalChargeArray, steps, I_model_Komayko)
    indexesToPlot = [20, 40, 60, 80, 100, 120, 140];

    figure;
    
    % Plot diffusion parameter Lambda
    plot(Lambda_tau_DeltaQ_array(:, 1), 'b-o', 'LineWidth', 1.5, 'DisplayName', 'Diffusion Coefficient');
    ylabel('\Lambda (unitless)');
    xlabel('Step Number');
    grid on;
    
    % Set the left y-axis range
    ylim([10^-2 100]); % Set the range for the left y-axis
    set(gca, 'YScale', 'log'); % Set left y-axis to logarithmic scale
    
    % Customize plot
    title('Resistance ratio Lambda vs. Step Number');
    legend('show');
        
    
    figure;
    
    % Plot diffusion parameter tau
    yyaxis left; % Use left y-axis
    plot(Lambda_tau_DeltaQ_array(:, 2), 'b-o', 'LineWidth', 1.5, 'DisplayName', 'Diffusion Coefficient');
    ylabel('\tau (s)'); % Use LaTeX-like syntax for special characters
    xlabel('Step Number');
    grid on;
    
    % Set the left y-axis range
    ylim([0 400]); % Set the range for the left y-axis
    set(gca, 'YScale', 'log'); % Set left y-axis to logarithmic scale
    
    
    % Plot resnorm
    % yyaxis right; % Use right y-axis
    % plot(resnorm_array, 'r-s', 'LineWidth', 1.5, 'DisplayName', 'Residual Norm');
    % ylabel('Residual Norm');
    
    % Customize plot
    title('Diffusion Coefficient Parameter τ and Residual Norm vs. Step Number');
    legend('show');
    


    figure;
    
    % Plot diffusion coefficient
    yyaxis left; % Use left y-axis
    %10^-12 cm^2 per particle
    plot(10^-12./Lambda_tau_DeltaQ_array(:, 2), 'b-o', 'LineWidth', 1.5, 'DisplayName', 'Diffusion Coefficient');
    ylabel('D (cm^2/s)'); % Use LaTeX-like syntax for special characters
    xlabel('Step Number');
    grid on;
    
    % Set the left y-axis range
    ylim([8*10^(-14) 4*10^(-12)]); % Set the range for the left y-axis
    set(gca, 'YScale', 'log'); % Set left y-axis to logarithmic scale
    
    
    % Plot resnorm
    % yyaxis right; % Use right y-axis
    % plot(resnorm_array, 'r-s', 'LineWidth', 1.5, 'DisplayName', 'Residual Norm');
    % ylabel('Residual Norm');
    
    % Customize plot
    title('Diffusion Coefficient D vs. Step Number');
    legend('show');
    

    figure;
    
    
    % Plot DeltaQ
    yyaxis left; % Use left y-axis
    plot(Lambda_tau_DeltaQ_array(:, 3), 'b-o', 'LineWidth', 1.5, 'DisplayName', 'DeltaQ (Fitted)');
    hold on; % Keep the left y-axis active for the next plot
    
    % Plot TotalChargeArray
    plot(totalChargeArray, 'g--s', 'LineWidth', 1.5, 'DisplayName', 'Total Charge (Experimental)');
    
    % Left y-axis labels and logarithmic scale
    ylabel('DeltaQ / Total Charge (Coulombs?)');
    xlabel('Step Number');
    grid on;
    set(gca, 'YScale', 'log'); % Set left y-axis to logarithmic scale
    
    % Adjust the left y-axis range (ensure it doesn't include zero or negative values)
    ylim([1e-6 1e-2]); % Adjust based on your data
    
    % Plot resnorm on the right y-axis
    yyaxis right; % Switch to right y-axis
    plot(resnorm_array, 'r-s', 'LineWidth', 1.5, 'DisplayName', 'Residual Norm');
    ylabel('Residual Norm');
    
    % Customize plot
    title('DeltaQ, Total Charge, and Residual Norm vs. Step Number');
    legend('show'); % Display the legend
    hold off;
    
    
    %[B_D_Q, resnorm] = PittFit(stepOfInterest_time_current(:,1), stepOfInterest_time_current(:,2), R_constant, B_D_Q_InitialGuess, I_t_Li);
    colors = lines(length(indexesToPlot)); % Generate distinct colors for each index
    
    figure; % Create a single figure for all plots
    hold on; % Enable overlaying multiple plots
    
    for i = 1:length(indexesToPlot)
        indexToPlot = indexesToPlot(i);
        
        % Extract step data
        [stepOfInterest, stepOfInterest_time_current] = extractStep(steps, indexToPlot);
        
        % Compute theoretical and experimental data
        I_fitted = I_model_Montella(Lambda_tau_DeltaQ_array(indexToPlot,:), stepOfInterest_time_current(:, 1));
        t_inv_sqrt = 1 ./ sqrt(stepOfInterest_time_current(:, 1)); % Compute t^{-1/2}
        
        % Plot experimental data
        plot(t_inv_sqrt, stepOfInterest_time_current(:, 2), 'o', 'Color', colors(i, :), ...
            'MarkerSize', 8, 'DisplayName', sprintf('Experimental (Index %d)', indexToPlot));
        
        % Plot theoretical data
        plot(t_inv_sqrt, I_fitted, '-', 'Color', colors(i, :), ...
            'LineWidth', 2, 'DisplayName', sprintf('Theoretical (Index %d)', indexToPlot));
    end
    
    % Customize the graph
    xlabel('t^{-1/2} (s^{-1/2})');
    ylabel('Current (A)');
    title('Transient Current vs. t^{-1/2}');
    xlim([0, 3]); % Adjust x-axis limits
    grid on;
    legend('show'); % Display the legend
    hold off; % Release the hold on the current figure
    
    
    
    
    figure; % Create a single figure with time linear for all plots
    hold on; % Enable overlaying multiple plots
    
    for i = 1:length(indexesToPlot)
        indexToPlot = indexesToPlot(i);
        
        % Extract step data
        [stepOfInterest, stepOfInterest_time_current] = extractStep(steps, indexToPlot);
        
        % Compute theoretical and experimental data
        I_fitted = I_model_Komayko(Lambda_tau_DeltaQ_array(indexToPlot,:), stepOfInterest_time_current(:, 1));
        %t_inv_sqrt = 1 ./ sqrt(stepOfInterest_time_current(:, 1)); % Compute t^{-1/2}
        
        % Plot experimental data
        semilogx(stepOfInterest_time_current(:, 1), stepOfInterest_time_current(:, 2), 'o', 'Color', colors(i, :), ...
            'MarkerSize', 8, 'DisplayName', sprintf('Experimental (Index %d)', indexToPlot));
        
        % Plot theoretical data
        semilogx(stepOfInterest_time_current(:, 1), I_fitted, '-', 'Color', colors(i, :), ...
            'LineWidth', 2, 'DisplayName', sprintf('Theoretical (Index %d)', indexToPlot));
    end
    
    % Customize the graph
    xlabel('t (s)');
    set(gca, 'XScale', 'log'); % Explicitly set the x-axis to logarithmic
    
    ylabel('Current (A)');
    title('Transient Current vs. time');
    %xlim([0, 3]); % Adjust x-axis limits
    ylim([0, 2E-6]); % Adjust x-axis limits

    grid on;
    legend('show'); % Display the legend
    hold off; % Release the hold on the current figure
    


    
    figure; % Create a single figure with time linear for all plots
    hold on; % Enable overlaying multiple plots
    
    for i = 1:length(indexesToPlot)
        indexToPlot = indexesToPlot(i);
        
        % Extract step data
        [stepOfInterest, stepOfInterest_time_current] = extractStep(steps, indexToPlot);
        
        % Compute theoretical and experimental data
        I_fitted = I_model_Komayko(Lambda_tau_DeltaQ_array(indexToPlot,:), stepOfInterest_time_current(:, 1));
        %t_inv_sqrt = 1 ./ sqrt(stepOfInterest_time_current(:, 1)); % Compute t^{-1/2}
        
        % Plot experimental data
        semilogx(stepOfInterest_time_current(:, 1), stepOfInterest_time_current(:, 2), 'o', 'Color', colors(i, :), ...
            'MarkerSize', 8, 'DisplayName', sprintf('Experimental (Index %d)', indexToPlot));
        
        % Plot theoretical data
        semilogx(stepOfInterest_time_current(:, 1), I_fitted, '-', 'Color', colors(i, :), ...
            'LineWidth', 2, 'DisplayName', sprintf('Theoretical (Index %d)', indexToPlot));
    end
    
    % Customize the graph
    xlabel('t (s)');
    set(gca, 'XScale', 'log'); % Explicitly set the x-axis to logarithmic
    
    ylabel('Current (A)');
    title('Transient Current vs. time');
    %xlim([0, 3]); % Adjust x-axis limits
    ylim([0, 2E-6]); % Adjust x-axis limits

    grid on;
    legend('show'); % Display the legend
    hold off; % Release the hold on the current figure
    

    
    %check that I*sqrt(t) vs log(t) has a curve similar to predictions
    figure; % Create a new figure with a logarithmic x-axis
    hold on; % Enable overlaying multiple plots
    
    for i = 1:length(indexesToPlot)
        indexToPlot = indexesToPlot(i);
        
        % Extract step data
        [stepOfInterest, stepOfInterest_time_current] = extractStep(steps, indexToPlot);
        
        % Ensure strictly positive x-axis values
        validIdx = stepOfInterest_time_current(:, 1) > 0; % Filter out non-positive time values
        t_valid = stepOfInterest_time_current(validIdx, 1); % Filtered time values
        current_valid = stepOfInterest_time_current(validIdx, 2); % Corresponding current values
        
        % Compute theoretical and experimental data
        I_fitted_valid = I_model_Komayko(Lambda_tau_DeltaQ_array(indexToPlot, :), t_valid);
        t_sqrt_valid = sqrt(t_valid); % Compute sqrt(t)
        
        % Plot experimental data
        semilogx(t_valid, current_valid .* t_sqrt_valid, 'o', ...
            'Color', colors(i, :), 'MarkerSize', 8, ...
            'DisplayName', sprintf('Experimental (Index %d)', indexToPlot));
        
        % Plot theoretical data
        semilogx(t_valid, I_fitted_valid .* t_sqrt_valid, '-', ...
            'Color', colors(i, :), 'LineWidth', 2, ...
            'DisplayName', sprintf('Theoretical (Index %d)', indexToPlot));
    end
    
    % Customize the graph
    xlabel('t (s)'); % Label x-axis
    ylabel('Current * sqrt(t) (A sqrt(s))'); % Label y-axis
    title('Transient Current vs. Time (Logarithmic Scale)');
    set(gca, 'XScale', 'log'); % Explicitly set the x-axis to logarithmic
     ylim([0, 2E-6]); % Adjust x-axis limits

    grid on; % Enable grid
    legend('show'); % Show legend
    hold off; % Release hold on current figure

    
    %take ratio of fit charge to expected charge from integrating current
    %during step (should eb equal to 1)
    sanityCheck = Lambda_tau_DeltaQ_array(:,3)./totalChargeArray;
    figure;
    hold on;
    plot(sanityCheck)
    ylim([0 3]); % Adjust based on your data

    xlabel('Step #'); % Label x-axis
    ylabel('Sanity ratio (nominally 1))'); % Label y-axis
    title('Sanity check of measured charge)');
    


end