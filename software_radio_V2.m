%%variables definition%%
number_of_realizations = 500;
bits_per_realization = 103;
required_bits_per_realization = 100;
samples_per_bit = 7;
realization_index_first  = 1;
realization_index_second = 2;
%%User Interface%%
choice = input("Enter The Required Line Code:\n1_ Unipolar NRZ\n2_ Polar NRZ\n3_ Polar RZ\n");
%%Function calls%%
ensemble_array = generate_ensemble(choice,number_of_realizations,bits_per_realization,required_bits_per_realization,samples_per_bit);
plot_realization(ensemble_array, realization_index_first, realization_index_second, samples_per_bit, required_bits_per_realization);
sm = get_statistical_mean(ensemble_array);
plot_statistical_mean(sm);
[autocorrelation_R,autocorrelation_average_Ravg] = get_autocorrelation(ensemble_array);
t1_values  = [1,100,120];
plot_autocorrelation(autocorrelation_R, autocorrelation_average_Ravg, size(ensemble_array,2), t1_values);
tm = get_time_mean(ensemble_array);
plot_time_mean(tm);
[tm_autocorrelation,tm_autocorrelation_function] = get_time_autocorrelation(ensemble_array);
plot_time_autocorrelation(tm_autocorrelation_function,ensemble_array);
calculate_BW(ensemble_array,autocorrelation_average_Ravg,size(ensemble_array,2),choice);
is_process_stationary(sm,choice,autocorrelation_average_Ravg);
is_process_ergoidic(sm,tm,tm_autocorrelation_function,autocorrelation_average_Ravg);
%%Functions%%
function ensemble = generate_ensemble(choice,number_of_realizations,bits_per_realization,required_bits_per_realization,samples_per_bit)
    %   This is an ensemble generator function.
    %   The input to this function will be a user's choice (1,2,3)
    %   1. Unipolar     2. Polar NRZ    3. Polar RZ

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
end

function plot_realization(ensemble, realization_index_first,realization_index_second, samples_per_bit, required_bits_per_realization)
% Generate the time vector
    total_samples = required_bits_per_realization * samples_per_bit;
    time = (0:total_samples - 1); % Sample indices as time

    % Extract the realization
    realization_signal_first  = ensemble(realization_index_first, :);
    realization_signal_second = ensemble(realization_index_second, :);

    % Plot
    figure;
    subplot(2,1,1);
    plot(time, realization_signal_first, 'LineWidth', 2);
    xlabel('Time (samples)');
    ylabel('Amplitude');
    title(['Realization ', num2str(realization_index_first), ' Line Code']);
    grid on;
    subplot(2,1,2);
    plot(time, realization_signal_second, 'LineWidth', 2);
    xlabel('Time (samples)');
    ylabel('Amplitude');
    title(['Realization ', num2str(realization_index_second), ' Line Code']);
    grid on;
end

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

function plot_statistical_mean(statistical_mean)
% This function plots the statistical mean of a given random process.
    figure;
    plot(statistical_mean, 'LineWidth', 1.5);
    xlabel('Time (samples)');
    ylabel('Mean Value');
    title('Statistical Mean of the Ensemble');
    grid on;
    xlim([1 length(statistical_mean)]);
end

function [autocorrelation,autocorrelation_average] = get_autocorrelation(ensemble)
%   This function computes the autocorrelation function Rx(t1, tau).
%   The output is a 2D array where each row represents Rx(t1, tau)
%   for different values of t1 (from 1 to number_of_samples).

if ~ismatrix(ensemble)
    error('Invalid Input to Function, please provide a 2D array that represents an ensemble of realizations of a random process');
end    

number_of_samples = size(ensemble, 2);
number_of_realizations = size(ensemble, 1);

max_tau = number_of_samples;  % Maximum time shift = 700
autocorrelation = zeros(number_of_samples, 2 * max_tau + 1);  

% Compute autocorrelation for different values of t1 and tau

for t1 = 1:number_of_samples
    for tau = -max_tau:max_tau 

        t2 = t1 + tau;
        %if (t2 <= 0)
            %t2 = t2 + max_tau;
        if ((t2 > 700) || (t2 <= 0))
            t2 = randi([1,700],1,1); % To avoid any discontinuities due to going out of bounds
        else
            t2 = t2;
        end

        for realization_index = 1:number_of_realizations 

                autocorrelation(t1,tau + max_tau + 1) = autocorrelation(t1,tau + max_tau + 1) + ((ensemble(realization_index,t1) * ensemble(realization_index,t2))/number_of_realizations);
        end
    end
end

autocorrelation_average = get_statistical_mean(autocorrelation);

end

function plot_autocorrelation(autocorrelation, autocorrelation_average, max_tau, t1_values)
    %   This function plots the autocorrelation function R_x(?) for selected values of t1,
    %   along with the averaged autocorrelation across all t1 values.
    %
    %   Inputs:
    %       autocorrelation          - The 2D autocorrelation matrix (rows = different t1 values)
    %       autocorrelation_average  - The mean autocorrelation function across all t1 values
    %       max_tau                  - The maximum time shift (?) considered
    %       t1_values                - A vector containing the values of t1 to plot (e.g., [1, 100, 300])
    %
    %   Output:
    %       A plot of R_x(?) against real values of ?, including both selected t1 values and the average.

    % Create the tau axis ranging from -max_tau to +max_tau
    tau_values = -max_tau:max_tau;

    % Create figure
    figure;
    hold on;

    % Plot for each selected t1
    for i = 1:length(t1_values)
        t1_index = t1_values(i);
        plot(tau_values, autocorrelation(t1_index, :), 'LineWidth', 1.5, 'DisplayName', ['t_1 = ', num2str(t1_index)]);
    end

    % Plot the averaged autocorrelation in bold
    plot(tau_values, autocorrelation_average, 'k', 'LineWidth', 2.5, 'DisplayName', 'Average R_x(\tau)');

    % Labels and title
    xlabel('\tau (Time Shift)');
    ylabel('Autocorrelation R_x(t_1, \tau)');
    title('Autocorrelation Function for Different t_1 Values and Average');
    legend show;
    grid on;
    hold off;
end

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

function plot_time_mean(time_mean)
    % This function plots the time mean of a given ensemble.
    figure;
    plot(time_mean, 'LineWidth', 1.5);
    xlabel('Time (samples)');
    ylabel('Mean Value');
    title('time Mean of the realizations');
    grid on;
    xlim([1 length(time_mean)]);
end

function [time_autocorrelation,time_autocorrelation_function] = get_time_autocorrelation(ensemble)
    % This function computes the time autocorrelation function for an ensemble of signals.
    number_of_samples = size(ensemble,2);
    number_of_realizations = size(ensemble,1);

    max_tau = number_of_samples;  % Maximum lag
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

end

function plot_time_autocorrelation(time_autocorrelation_function,ensemble)
    %time autocorrelation plotting
    tau_values = -700:700;  % Range for symmetric shifts
    % Plot the autocorrelation function
    figure;
    plot(tau_values, time_autocorrelation_function, 'k', 'LineWidth', 2.5, 'DisplayName', '<R_x(\tau)>');
    xlabel('Time Shift (tau)');
    ylabel('Autocorrelation');
    title('Time Autocorrelation Function');
    legend show;
    grid on;
    hold off;
end

function calculate_BW(ensemble, autocorrelation_average, number_of_samples, choice)
    % Bandwidth calculation
    k = -number_of_samples : number_of_samples;
    fs = 100;  % Sampling frequency
    
    % Compute the Power Spectral Density (PSD)
    % Note that we divide by fs to normalize the effect of fourier
    % transform on the amplitude
    psd = abs(fftshift(fft(autocorrelation_average)))/fs;
    
    % Define parameters
    A = 4;
    T = 7;
    
    % Define the triangular function in time domain
    Tri = @(x, a) max(1 - abs(x) / (a / 2), 0);
    tau_values = -number_of_samples:number_of_samples;
    
    % Frequency axis
    freq_axis = k * fs / (2 * number_of_samples);
    
    % Plot PSD
    figure;
    hold on;
    plot(freq_axis, psd, 'k', 'LineWidth', 2, 'DisplayName', 'Average PSD');
    
    % Compute and plot based on choice
    switch choice
        case 1 % Unipolar
            Rx_unipolar = (A^2 / 4) * Tri(tau_values, 2*T) + (A^2 / 4);
            psd_unipolar = abs(fftshift(fft(Rx_unipolar)))/fs;
            plot(freq_axis, psd_unipolar, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Unipolar PSD');
        case 2 % Polar NRZ
            Rx_polar_nrz = (A^2) * Tri(tau_values, 2*T);
            psd_polar_nrz = abs(fftshift(fft(Rx_polar_nrz)))/fs;
            plot(freq_axis, psd_polar_nrz, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Polar NRZ PSD');
        case 3 % Polar RZ
            Rx_polar_rz = (A^2 / 2) * Tri(tau_values, T);
            psd_polar_rz = abs(fftshift(fft(Rx_polar_rz)))/fs;
            plot(freq_axis, psd_polar_rz, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Polar RZ PSD');
    end
    title("Power Spectral Density (PSD)");
    xlabel("Frequency (Hz)");
    ylabel("PSD");
    legend show;
    grid on;
    hold off;
end

function is_process_stationary(statistical_mean, choice, autocorrelation)
    % Define parameters
    A = 4;  % Amplitude
    T = 7;  % Bit Period in Samples
    tau_values = -700:700;
    
    % Define the triangular function
    Tri = @(x, a) max(1 - abs(x) / (a / 2), 0); % Tri function with width a
    
    % Compute the average of statistical_mean
    avg_value = mean(statistical_mean);
    figure;
    plot(statistical_mean, 'LineWidth', 1.5, 'DisplayName', ['Average mean = ', num2str(avg_value)]);
    hold on;
    if choice == 1
        yline(2, '--r', 'LineWidth', 1.5, 'DisplayName', ['Theoretical mean = ', num2str(2)]);
    elseif choice == 2 || choice == 3
        yline(0, '--r', 'LineWidth', 1.5, 'DisplayName', ['Theoretical mean = ', num2str(0)]);
    end
    hold off;
    xlabel('Time (samples)');
    ylabel('Mean Value');
    title('Statistical Mean of the Ensemble');
    legend show;
    grid on;
    xlim([1 length(statistical_mean)]);
    
    % Create figure for autocorrelation function
    figure;
    hold on;
    
    % Plot the averaged autocorrelation in bold
    plot(tau_values, autocorrelation, 'k', 'LineWidth', 2.5, 'DisplayName', ...
        sprintf('Average R_x(\tau), R(0) = %.2f , DC = %.2f', autocorrelation(701), autocorrelation(1401)));
    
    % Switch case for different choices to sketch the theoretical
    % autocorrelation function:
    switch choice
        case 1
            Rx_tau = ((A^2)/4) * Tri(tau_values, 2*T) + (A^2)/4;
            plot(tau_values, Rx_tau, 'b--', 'LineWidth', 2, 'DisplayName', 'R_x(\tau) theoretical for Unipolar NRZ');
        case 2
            Rx_tau = (A^2) * Tri(tau_values, 2*T);
            plot(tau_values, Rx_tau, 'r--', 'LineWidth', 2, 'DisplayName', 'R_x(\tau) theoretical for Polar NRZ');
        case 3
            Rx_tau = ((A^2)/2) * Tri(tau_values, T);
            plot(tau_values, Rx_tau, 'g--', 'LineWidth', 2, 'DisplayName', 'R_x(\tau) theoretical for Polar RZ');
    end
    
    % Labels and title
    xlabel('\tau (Time Shift)');
    ylabel('Autocorrelation R_x(t_1, \tau)');
    title('Autocorrelation Function Average');
    legend show;
    grid on;
    hold off;
end

function is_process_ergoidic(statistical_mean,time_mean,time_autocorrelation,autocorrelation)
    tau_values = -700:700;
    avg_value_sm = mean(statistical_mean);
    avg_value_tm = mean(time_mean);
    figure;
    subplot(1,2,1);
    plot(statistical_mean, 'LineWidth', 1.5, 'DisplayName', ['Average mean = ', num2str(avg_value_sm)]);
    xlabel('Time (samples)');
    ylabel('Mean Value');
    title('Statistical Mean of the Ensemble');
    legend show;
    grid on;
    xlim([1 length(statistical_mean)]);
    subplot(1,2,2);
    plot(time_mean, 'LineWidth', 1.5, 'DisplayName', ['Average mean = ', num2str(avg_value_tm)]);
    xlabel('Realizations');
    ylabel('Mean Value');
    title('time Mean of the Realizations');
    legend show;
    grid on;
    xlim([1 length(time_mean)]);
    
    figure;
    hold on;
    subplot(1,2,1);
    plot(tau_values, autocorrelation, 'k', 'LineWidth', 2.5, 'DisplayName', ...
        sprintf('Average R_x(\tau), R(0) = %.2f , DC = %.2f', autocorrelation(701), autocorrelation(1401)));
    xlabel('\tau (Time Shift)');
    ylabel('Autocorrelation R_x(t_1, \tau)');
    title('Autocorrelation Function Average');
    legend show;
    grid on;
    hold off;
    
    hold on;
    subplot(1,2,2);
    plot(tau_values, time_autocorrelation, 'k', 'LineWidth', 2.5, 'DisplayName', ...
        sprintf('Average R_t(\tau), R(0) = %.2f , DC = %.2f', time_autocorrelation(701), time_autocorrelation(1399)));
    xlabel('\tau (Time Shift)');
    ylabel('Autocorrelation R_x(t_1, \tau)');
    title('Time Autocorrelation Function');
    legend show;
    grid on;
    hold off;
    
end