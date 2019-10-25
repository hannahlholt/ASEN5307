% this program tries to understand how wavelets work.

% % %  ---- REAL DATA ----
% T = readtable('t_series_data.csv');
% 
% x_tot = T{:,2};
% t = T{:,1};
% dt = T{2,1} - T{1,1};    % dt is in hours
% fs = 1/dt;
% 
% % quiet time
% stop = find(t == 79);
% xq = x_tot(1:stop);
% xq_dt = detrend(xq);                % detrend it
% 
% % storm time
% start = find(t == 84);
% xst = x_tot(start:end);
% xst_dt = detrend(xst);                  % detrend it
% 
% % now splice together the storm and quiet time
% x = [xq_dt; xst_dt];
% 
% name1 = './Figures/Powerspectrum_obs.png';
% name2 = './Figures/Morlet_obs.png';
% xlimits = [0, 1.3];



% ---- SYNTHETIC DATA----
dt = 0.5;
fs = 1/dt;
t = [0:dt:10000];
ind = find(t == 3000);
A = 40;
B = 100;
P1 = 50;
P2 = 100;

x = A * sin(2*pi*t/P1);
x(ind:end) = x(ind:end) +  B * sin(2*pi*t(ind:end)/P2);
n = std(x) *  randn(size(t));
x = x + n;
name1 = './Figures/Powerspectrum_syn.png';
name2 = './Figures/Morlet_syn.png';
xlimits = [0, 200];


nfft = 2^nextpow2(length(x));     % next power of 2 to use for fft points
[Px, fq] = periodogram(x, [], nfft, fs);
xnorm = Px/sum(Px(:));

h = figure();
subplot(211);
plot(t, x);
title('Signal');

subplot(212);
plot(1./fq, xnorm);
set(gca, 'YScale', 'log');
title('Power spectrum')
xlabel('Time');
xlim(xlimits)
grid on;
saveas(gcf, name1);


% waitfor(h)

cwt(x, 'amor', seconds(dt))
saveas(gcf, name2);

