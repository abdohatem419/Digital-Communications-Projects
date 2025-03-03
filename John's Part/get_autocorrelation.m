function [autocorrelation,autocorrelation_average] = get_autocorrelation(ensemble)
%   This function computes the autocorrelation function Rx(t1, tau).
%   The output is a 2D array where each row represents Rx(t1, tau)
%   for different values of t1 (from 1 to number_of_samples).

if ~ismatrix(ensemble)
    error('Invalid Input to Function, please provide a 2D array that represents an ensemble of realizations of a random process');
end    

number_of_samples = size(ensemble, 2);
number_of_realizations = size(ensemble, 1);

max_tau = number_of_samples - 1;  % Maximum time shift = 699
autocorrelation = zeros(number_of_samples, 2 * max_tau + 1);  

% Compute autocorrelation for different values of t1 and tau
for t1 = 1:number_of_samples
    for tau = -max_tau:max_tau
        t2 = t1 + tau;
        
        if t2 >= 1 && t2 <= number_of_samples
            for realization_index = 1:number_of_realizations
                autocorrelation(t1, tau + max_tau + 1) = autocorrelation(t1, tau + max_tau + 1) + (ensemble(realization_index, t1) * ensemble(realization_index, t2)) / number_of_realizations;     
            end    
        end
    end    
end

autocorrelation_average = get_statistical_mean(autocorrelation);

end
