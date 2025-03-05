%variables definition
number_of_realizations = 500;
bits_per_realization = 103;
required_bits_per_realization = 100;
samples_per_bit = 7;
realization_index = 1;
%ensemble information
ensemble_bits = randi([0,1],number_of_realizations,bits_per_realization);
delay_time = randi([0,6],number_of_realizations,1);
A = 4;
ensemble = zeros(number_of_realizations, required_bits_per_realization* samples_per_bit);
%starting point
choice = input("Enter The Required Line Code:\n1_ Unipolar NRZ\n2_ Polar NRZ\n3_ Polar RZ\n");
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
% Generate the time vector
total_samples = required_bits_per_realization * samples_per_bit;
time = (0:total_samples - 1); % Sample indices as time
% Extract the realization
realization_signal = ensemble(realization_index, :);
% Plot
figure;
plot(time, realization_signal, 'LineWidth', 2);
xlabel('Time (samples)');
ylabel('Amplitude');
title(['Realization ', num2str(realization_index), ' Line Code']);
grid on;

%statistical mean calculation
number_of_samples = size(ensemble,2);
number_of_realizations = size(ensemble,1);
statistical_mean = zeros(1,number_of_samples);
for sample_index = 1:number_of_samples
    for realization_index = 1:number_of_realizations
        statistical_mean(sample_index) = statistical_mean(sample_index) + (ensemble(realization_index,sample_index)/number_of_realizations);
    end    
end

%statistical mean plotting
figure;
plot(statistical_mean, 'LineWidth', 1.5);
xlabel('Time (samples)');
ylabel('Mean Value');
title('Statistical Mean of the Ensemble');
grid on;
xlim([1 length(statistical_mean)]);

%statistical autocorrelation calculation
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
number_of_samples = size(autocorrelation,2);
number_of_realizations = size(autocorrelation,1);
autocorrelation_average = zeros(1,number_of_samples);
for sample_index = 1:number_of_samples
    for realization_index = 1:number_of_realizations
        autocorrelation_average(sample_index) = autocorrelation_average(sample_index) + (autocorrelation(realization_index,sample_index)/number_of_realizations);
    end    
end    

%statistical autocorrelation plotting 
% Create the tau axis ranging from -max_tau to +max_tau
tau_values = -max_tau:max_tau;
t1_values  = [1,100,120];
figure;
hold on;
% Plot for each selected t1
 for i = 1:length(t1_values)
     t1_index = t1_values(i);
     plot(tau_values, autocorrelation(t1_index, :), 'LineWidth', 1.5, 'DisplayName', ['t_1 = ', num2str(t1_index)]);
 end
% Plot the averaged autocorrelation in bold
plot(tau_values, autocorrelation_average, 'k', 'LineWidth', 2.5, 'DisplayName', 'Average R_x(\tau)');
xlabel('\tau (Time Shift)');
ylabel('Autocorrelation R_x(t_1, \tau)');
title('Autocorrelation Function for Different t_1 Values and Average');
legend show;
grid on;
hold off;

%time mean calculation
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

%time mean plotting
figure;
plot(time_mean, 'LineWidth', 1.5);
xlabel('Time (samples)');
ylabel('Mean Value');
title('time Mean of the realizations');
grid on;
xlim([1 length(time_mean)]);

%time autocorrelation calculation
max_tau = number_of_samples-1;  % Maximum lag
time_autocorrelation = zeros(1, max_tau + 1); % Stores final autocorrelation values
% Compute time autocorrelation function
for tau = 0:max_tau
    sum = 0;
    % Loop over all realizations
    for realization = 1:number_of_realizations
        for count = 1:(number_of_samples - tau)
            sum = sum + (ensemble(realization, count) * ensemble(realization, count + tau));
        end
    end
    % Normalize by the number of valid samples
    time_autocorrelation(tau+1) = sum / ((number_of_samples - tau) * number_of_realizations);
end
time_autocorrelation_function = [fliplr(time_autocorrelation(2:end)) time_autocorrelation];

%time autocorrelation plotting
tau_values = -699:699;  % Range for symmetric shifts
% Plot the autocorrelation function
figure;
plot(tau_values, time_autocorrelation_function, 'k', 'LineWidth', 2.5, 'DisplayName', '<R_x(\tau)>');
xlabel('Time Shift (tau)');
ylabel('Autocorrelation');
title('Time Autocorrelation Function');
legend show;
grid on;
hold off;

%Bandwidth calculation
k = -number_of_samples + 1: number_of_samples - 1;
fs = 100;  % Sampling frequency
psd = abs(fftshift(fft(autocorrelation_average)));
figure;
plot(k * fs / (2 * number_of_samples), psd, 'b', 'LineWidth', 1.5);
title("PSD");
xlabel("Frequency (Hz)");
ylabel("PSD");
% Zoom into relevant range
% Highlight zero-crossing point
[~, zero_idx] = min(abs(psd)); % Find nearest zero-crossing index
hold on;
plot(k(zero_idx) * fs / (2 * number_of_samples), psd(zero_idx), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
hold off;