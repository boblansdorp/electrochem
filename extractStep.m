function [stepOfInterest, stepOfInterest_time_current] = extractStep(steps, elementToObtain)
% Access specific steps and retain rows based on logarithmic spacing
    MINIMUM_STEP_RATIO = 1.1;
    if elementToObtain <= length(steps) && ~isempty(steps{elementToObtain})
        % Access the specific step
        stepOfInterest = steps{elementToObtain};
        stepOfInterest_time_current = [];

        % Initialize variables
        lastTime = stepOfInterest(1, 1); % Start with the first time point
        stepOfInterest_time_current = [lastTime, stepOfInterest(1, 5)]; % Include the first row

        % Loop through the rest of the rows
        for n = 2:size(stepOfInterest, 1)
            currentTime = stepOfInterest(n, 1);

            % Include the row if the time difference is > 30%
            if currentTime > (MINIMUM_STEP_RATIO * lastTime)
                stepOfInterest_time_current = [stepOfInterest_time_current; currentTime, stepOfInterest(n, 5)];
                lastTime = currentTime; % Update the last time
            end
        end
    else
        fprintf('Step %d does not exist or is empty.\n', elementToObtain);
        stepOfInterest = [];
        stepOfInterest_time_current = [];
    end
end
