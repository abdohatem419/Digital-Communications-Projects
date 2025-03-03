function statistical_mean = get_statistical_mean(ensemble)
%   This function is to compute the statistical mean of a given random
%   process that is represented as an ensemble of multiple realizations.

if ~ismatrix(ensemble)
    error('Invalid Input to Function, the input must be a 2D array representing an ensemble');
end    

number_of_samples = size(ensemble,2);
number_of_realizations = size(ensemble,1);

statistical_mean = zeros(1,number_of_samples);

for sample_index = 1:number_of_samples
    for realization_index = 1:number_of_realizations
        statistical_mean(sample_index) = statistical_mean(sample_index) + (ensemble(realization_index,sample_index)/number_of_realizations);
    end    
end    

end