function plot_time_mean(time_mean)
% This function plots the time mean of a given ensemble.
    figure;
    time=1:500; 
    figure;
    plot(time,time_mean(1:500));
    axis([1 500 -5 5]); 
    title('Time Mean '); 
    xlabel('time'); 
    ylabel('Time Mean');
end
