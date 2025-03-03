function time_autocorrelation = get_time_autocorrelation(ensemble)
    % This function computes the time autocorrelation function for an ensemble of signals.

    [number_of_realizations, number_of_samples] = size(ensemble);
    
    max_tau = number_of_samples-1;  % Maximum lag
    time_autocorrelation = zeros(1, max_tau + 1); % Stores final autocorrelation values

    % Compute time autocorrelation function
    for tau = 0:max_tau
        accumulate = 0;
        
        % Loop over all realizations
        for realization = 1:number_of_realizations
            for count = 1:(number_of_samples - tau)
                accumulate = accumulate + (ensemble(realization, count) * ensemble(realization, count + tau));
            end
        end

        % Normalize by the number of valid samples
        time_autocorrelation(tau+1) = accumulate / ((number_of_samples - tau) * number_of_realizations);
    end
    time_autocorrelation_function = [fliplr(time_autocorrelation(2:end)) time_autocorrelation];

end
