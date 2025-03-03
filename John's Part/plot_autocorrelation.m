function plot_autocorrelation(autocorrelation, autocorrelation_average, max_tau, t1_values)
%   This function plots the autocorrelation function R_x(τ) for selected values of t1,
%   along with the averaged autocorrelation across all t1 values.
%
%   Inputs:
%       autocorrelation          - The 2D autocorrelation matrix (rows = different t1 values)
%       autocorrelation_average  - The mean autocorrelation function across all t1 values
%       max_tau                  - The maximum time shift (τ) considered
%       t1_values                - A vector containing the values of t1 to plot (e.g., [1, 100, 300])
%
%   Output:
%       A plot of R_x(τ) against real values of τ, including both selected t1 values and the average.

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

