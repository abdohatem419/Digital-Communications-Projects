function time_mean = get_time_mean(ensemble)
%   This function computes the time mean for realizations.
%   The output is a 2D array where each row represents time mean of each
%   waveform 

number_of_samples = size(ensemble, 2);
number_of_realizations = size(ensemble, 1);
time_mean=zeros(number_of_realizations,1);

for realization_index =1:number_of_realizations
    sum=0;
    for sample_index=1:number_of_samples
        sum=ensemble(realization_index,sample_index)+sum;
    end
    time_mean(realization_index)=sum/number_of_samples;
end


end