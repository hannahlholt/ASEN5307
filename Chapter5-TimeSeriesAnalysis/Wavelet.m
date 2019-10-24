% this program tries to understand how sucky wavelets work.
% cool cool


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
nfft = 2^nextpow2(length(x));     % next power of 2 to use for fft points

[Px, fq] = periodogram(x, [], nfft, fs);
xnorm = Px/sum(Px(:));

figure();
subplot(211);
plot(t, x);
title('Signal');

subplot(212);
plot(1./fq, xnorm);
set(gca, 'YScale', 'log');
title('Power spectrum')
xlabel('Time');
xlim([0 200])
grid on;

saveas(gcf, './Powerspectrum.png');


cwt(x, 'amor', seconds(dt))
saveas(gcf, './wavelet.png');

% tms = (0:numel(x)-1)/fs;
% 
% figure
% subplot(2,1,1)
% plot(t,x)
% axis tight
% title('Signal')
% xlabel('Time (s)')
% ylabel('Amplitude')
% 
% 
% subplot(2,1,2)
% surface(t,1./frq,abs(cfs))
% axis tight
% shading flat
% xlabel('Time (s)')
% ylabel('Period (sec)')
% % set(gca,'yscale','log')
% 