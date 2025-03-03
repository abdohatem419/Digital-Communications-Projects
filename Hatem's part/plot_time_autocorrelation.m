function plot_time_autocorrelation(time_autocorrelation_function)
% This function plots the time autocorrelation of a given ensemble.

    % Define time lag axis
    tau_values = -699:699;  % Range for symmetric shifts

    % Plot the autocorrelation function
    figure;
    plot(tau_values, time_autocorrelation_function, 'b', 'LineWidth', 2);
    xlabel('Time Shift (tau)');
    ylabel('Autocorrelation');
    title('Time Autocorrelation Function');
end
