% this program tries to understand how sucky wavelets work.
% cool cool


dt = 0.5;
fs = 1/dt;
t = [0:dt:5000];

A = 23;
B = 10;
f1 = 1/50;
f2 = 1/100;

x = A * sin(2*pi*f1*t)  + B * sin(2*pi*f2*t) + randn(size(t));

% x = x + x_noise;

% [cfs,frq]  = ;
cwt(x, 'amor', seconds(dt))


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