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

