function ensemble = generate_ensemble(choice)
%   This is an ensemble generator function.
%   The input to this function will be a user's choice (1,2,3)
%   1. Unipolar     2. Polar NRZ    3. Polar RZ

number_of_realizations = 500;
bits_per_realization = 103;
required_bits_per_realization = 100;
samples_per_bit = 7;

ensemble_bits = randi([0,1],number_of_realizations,bits_per_realization);
delay_time = randi([0,6],number_of_realizations,1);
A = 4;

ensemble = zeros(number_of_realizations, required_bits_per_realization* samples_per_bit);

switch (choice)

    case (1)
        for i = 1:number_of_realizations
            realization_bits = ensemble_bits (i,:);
            realization_values = realization_bits*A;
            realization_samples = repelem(realization_values,samples_per_bit);
            realization_samples_delayed = realization_samples(delay_time(i)+1:delay_time(i)+(samples_per_bit*required_bits_per_realization));
            ensemble(i,:) = realization_samples_delayed;
        end
    case (2)
        for i = 1:number_of_realizations
            realization_bits = ensemble_bits(i,:);
            realization_values = A*(2*realization_bits - 1);
            realization_samples = repelem(realization_values,samples_per_bit);
            realization_samples_delayed = realization_samples(delay_time(i)+1:delay_time(i)+(samples_per_bit*required_bits_per_realization));
            ensemble(i,:) = realization_samples_delayed;
        end    
    case (3)
        for i = 1:number_of_realizations
            realization_bits = ensemble_bits(i,:);
            realization_values = A*(2*realization_bits - 1);

            realization_samples = zeros(1,bits_per_realization*samples_per_bit);
            for j = 1:bits_per_realization
                realization_samples((j-1)*samples_per_bit+1:(j-1)*samples_per_bit+4) = realization_values(j);
            end

            realization_samples_delayed = realization_samples(delay_time(i)+1:delay_time(i)+(samples_per_bit*required_bits_per_realization));
            ensemble(i,:) = realization_samples_delayed;
        end   
    otherwise
            error('Incorrect Value for "Choice". Please Choose Either 1, 2 or 3');

end