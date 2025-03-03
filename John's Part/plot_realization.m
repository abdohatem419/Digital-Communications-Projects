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

