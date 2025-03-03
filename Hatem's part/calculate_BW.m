function calculate_BW(ensemble,autocorrelation)

    number_of_samples = size(ensemble, 2);
    
    k = -number_of_samples + 1: number_of_samples - 1;
    fs = 100;  %Sampling frequency
    figure;
    plot(k*fs/(2*number_of_samples),abs(fftshift(fft(autocorrelation))));
    title("PSD");
    xlabel("Frequency(Hz)");
    ylabel("PSD");

end