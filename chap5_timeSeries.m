%% Chapter 5 Time series Analysis


% if the time series doesn't exist, then run the LOADfiles script.
if exist('t_series', 'var') == 0
    run LOADfiles.m;
end


% nfft = 2^nextpow2(length(t_series.t));  % next power of 2 to use for fft
dt = t_series.t(2) - t_series.t(1);     % sampling period [day]
fs = 1/dt;                              % sampling frequency [1/day]

% --- PERIODOGRAM VALUES ------
% quiet time
stop = find(t_series.t == 80);
denq = t_series.Den(1:stop);
denq_dt = detrend(denq);                % detrend it
nfft = 2^nextpow2(length(denq_dt));     % next power of 2 to use for fft points
[Pxxq, fq] = periodogram(denq_dt, [], nfft, fs);

% storm time
start = find(t_series.t == 83);
denst = t_series.Den(start:end);
denst_dt = detrend(denst);                  % detrend it
nfft = 2^nextpow2(length(denst_dt));        % next power of 2 to use for fft points
[Pxxst, fst] = periodogram(denst_dt, [], nfft, fs);





%%

figure('units', 'normalized', 'outerposition',[0 1 .8 1]);

subplot(221)
plot(t_series.t(1:stop), t_series.Den(1:stop), 'b'); hold on;
plot(t_series.t(start:end), t_series.Den(start:end), 'r'); 
plot(t_series.t(stop:start), t_series.Den(stop:start), 'k');
legend('Quiet Time', 'Storm time')
ylim([.99*min(t_series.Den), 1.01*max(t_series.Den)]);
title(['Average Global Neutral Denisty, Zp = ', num2str(t_series.Zlvl)])
vline([t_series.t(2), t_series.t(stop), t_series.t(start), t_series.t(end-1)], {'b', 'b', 'r', 'r'});
xlabel('Model Time [Days]')
ylabel('Global Averaged Neutral Density [kg/m3]');
grid on;
hold off;

subplot(222)
x1_norm = (Pxxq -  min(Pxxq))/range(Pxxq);
plot(fq, x1_norm, 'b');
t1 = 'Autospectrum for Quiet Time Density';
t2 = ['(', num2str(t_series.t(1)), ' < t < ', num2str(t_series.t(stop)), ')'];
xlim([0 5]);
title({t1; t2});
xlabel('Frequency [1/days]')
ylabel('Normalized Power');
grid on;

subplot(223)
x2_norm = Pxxst/max(Pxxq);
plot(fst, x2_norm, 'r')
t1 = 'Autospectrum for Storm Time Density';
t2 = ['(', num2str(t_series.t(start)), ' < t < ', num2str(t_series.t(end)), ')'];
title({t1; t2});
xlim([0 5])
xlabel('Frequency [1/days]')
ylabel('Normalized Power');
grid on;



subplot(224)
test = t_series.Den;
s = spectrogram(test);
spectrogram(test_norm, 'yaxis')



