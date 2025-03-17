%%variables definition%%
number_of_realizations = 500;
bits_per_realization = 103;
required_bits_per_realization = 100;
samples_per_bit = 7;
realization_index = 1;
%%User Interface%%
choice = input("Enter The Required Line Code:\n1_ Unipolar NRZ\n2_ Polar NRZ\n3_ Polar RZ\n");
%%Function calls%%
ensemble_array = generate_ensemble(choice,number_of_realizations,bits_per_realization,required_bits_per_realization,samples_per_bit);
plot_realization(ensemble_array, realization_index, samples_per_bit, required_bits_per_realization);
sm = get_statistical_mean(ensemble_array);
plot_statistical_mean(sm);
[autocorrelation_R,autocorrelation_average_Ravg] = get_autocorrelation(ensemble_array);
t1_values  = [1,100,120];
plot_autocorrelation(autocorrelation_R, autocorrelation_average_Ravg, (size(ensemble_array,2)-1), t1_values);
tm = get_time_mean(ensemble_array);
plot_time_mean(tm);
[tm_autocorrelation,tm_autocorrelation_function] = get_time_autocorrelation(ensemble_array);
plot_time_autocorrelation(tm_autocorrelation_function,ensemble_array);
calculate_BW(ensemble_array,autocorrelation_average_Ravg,size(ensemble_array,2));
is_process_stationary(sm,choice,autocorrelation_average_Ravg);
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

function plot_realization(ensemble, realization_index, samples_per_bit, required_bits_per_realization)
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

end

function plot_time_autocorrelation(time_autocorrelation_function,ensemble)
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
end

function calculate_BW(ensemble,autocorrelation_average,number_of_samples)
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
end

function is_process_stationary(statistical_mean,choice,autocorrelation)
    % Compute the average of statistical_mean
    avg_value = mean(statistical_mean);
    figure;
    %plot(statistical_mean, 'LineWidth', 1.5);
    plot(statistical_mean, 'LineWidth', 1.5, 'DisplayName', ['Average mean = ', num2str(avg_value)]);
        if choice == 1
            yline(2, '--r', 'LineWidth', 1.5, 'DisplayName', ['Theoritical mean = ', num2str(2)]);
        elseif choice == 2
            yline(0, '--r', 'LineWidth', 1.5, 'DisplayName', ['Theoritical mean = ', num2str(0)]);
        elseif choice == 3
           yline(0, '--r', 'LineWidth', 1.5, 'DisplayName', ['Theoritical mean = ', num2str(0)]);
        end
    xlabel('Time (samples)');
    ylabel('Mean Value');
    title('Statistical Mean of the Ensemble');
    legend show;
    grid on;
    xlim([1 length(statistical_mean)]);
    tau_values = -699:699;
    % Create figure
    figure;
    hold on;
    % Plot the averaged autocorrelation in bold
    plot(tau_values, autocorrelation, 'k', 'LineWidth', 2.5, 'DisplayName', ...
    sprintf('Average R_x(\\tau), R(0) = %.2f , DC = %.2f', autocorrelation(1,1),mean(autocorrelation)));
    % Labels and title
    xlabel('\tau (Time Shift)');
    ylabel('Autocorrelation R_x(t_1, \tau)');
    title('Autocorrelation Function Average');
    legend show;
    grid on;
    hold off;
end