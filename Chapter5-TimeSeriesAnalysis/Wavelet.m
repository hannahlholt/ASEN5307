% this program tries to understand how wavelets work.
close all;
clear all;

%%  ---- REAL DATA ----
T = readtable('t_series_data.csv');

x_tot = T{:,2};
t = T{:,1};
UT = T{:,3};            % universal time
dt = T{2,1} - T{1,1};    % dt is in hours
fs = 1/dt;
xname = '(days)';

% quiet time
stop = find(t == 80);
xq = x_tot(1:stop);
tq = t(1:stop);
xq_dt = detrend(xq, 2);    % detrend it

% storm time
start = find(t == 82);
xst = x_tot(start:end);
tst = t(start:end);
xst_dt = detrend(xst, 3);                  % detrend it

% now splice together the storm and quiet time
x = [xq_dt; xst_dt];
tsplice = [tq; tst];
UTsplice = [UT(1:stop); UT(start:end)];

name1 = './Figures/Powerspectrum_obs.png';
name2 = './Figures/Morlet_obs.png';
name3 = './Figures/MorletPhase_obs.png';
name4 = './Figures/RecreatedSignal.png';
name5 = './Figures/RecreatedSignalZOOM.png';
xlimits = [0, 1.3];

%% ---- SYNTHETIC DATA----
% dt = 0.5;
% fs = 1/dt;
% t = [0:dt:10000];
% xname = '(sec)';
% ind = find(t == 3000);
% A = 40;
% B = 100;
% P1 = 50;
% P2 = 100;
% 
% x = A * sin(2*pi*t/P1);
% x(ind:end) = x(ind:end) +  B * sin(2*pi*t(ind:end)/P2);
% n = std(x) *  randn(size(t));
% x = x + n;
% name1 = './Figures/Powerspectrum_syn.png';
% name2 = './Figures/Morlet_syn.png';
% xlimits = [0, 200];

%% ----- COMPUTATIONS ----------
% ----- COMPUTE POWER SPECTRUM ---------------
nfft = 2^nextpow2(length(x));     % next power of 2 to use for fft points
[Px, fq] = periodogram(x, [], nfft, fs);
xnorm = Px/sum(Px(:));
 
% ----- COMPUTE WAVELET TRANSFORM -------------
[wt, period, coi] = cwt(x, 'amor', days(dt));

Z_og = abs(wt);
Z = Z_og;
[X,Y] = meshgrid(tsplice, days(period));

% for plotting cone of influence (coi)
Z1 = repmat(days(coi'), length(Z(:,1)), 1);
Z(Y > Z1) = NaN;

% ---- COMPUTE PHASE -------
% 1) one day oscillation
n = days(1);

% index in time space (rows) of day oscillation
[~, idx1] = min(abs(period-n));


% 2) half day oscillation
n = days(0.5);
% index in time space (rows) of half day oscillation
[~, idx2] = min(abs(period-n));

phi = 180/pi*(angle(wt));            % phase angle [-180, 180]

% --- RECREATE SIGNAL -----
x_day = Z_og(idx1,:) .* sin( 2*pi./days(period(idx1)) + pi/180.*phi(idx1,:));
x_halfday = Z_og(idx2,:) .* sin( 2*pi./days(period(idx2)) + pi/180.*phi(idx2,:));


% --- create Phase of the North pole w.r.t noon
lon_Npole = 190;                %  degrees E  , N pole at 2018??
dt = lon_Npole/15 - 12;
phi0 = mod(2*pi*dt/24, 2*pi);

% sine wave of pole w.r.t noon. (MAX when pole is at noon!!)
x_pole_noon = cos(2*pi/24*UTsplice + phi0); 

%% ------ PLOTTING -------

%% PLOT 1D POWER SPECTRUM
h0 = figure();
subplot(211);
plot(tsplice, x);
title('Signal');

subplot(212);
plot(1./fq, xnorm);
set(gca, 'YScale', 'log');
title('Power spectrum')
xlabel(['Time ', xname]);
xlim(xlimits)
grid on;
saveas(h0, name1);
%

%% PLOT WAVELET ANALYSIS
h1 = figure();
contourf(X, Y, Z, 100, 'linecolor', 'none')
hold on;
colormap jet;
cbar = colorbar();
set(gca, 'YScale', 'log')
ylabel(['Period ', xname]);
xlabel(['Time ', xname]);
cbar.Label.String = 'Magnitude';
grid on;
saveas(h1, name2);


%% PLOT PHASE
h3 = figure('units', 'normalized', 'position', [0 .5 1 1], 'visible', 'on');
daystop = 72;
B = 12;

% -----------------------------------------------
subplot(211)
plot(tsplice, phi(idx1,:)); hold on;
plot(tsplice, zeros(length(tsplice),1), 'k')
axis([70 daystop -180 180])

ax = gca;
ax.XAxis.TickValues = [70:1/B:daystop];
oldtick = ax.XAxis.TickValues;
[~, indx] = min(abs(tsplice-oldtick));
ax.XTickLabel = num2str([0; UTsplice(indx(2:end))]);
title('One Day Periodicity');
xlabel('UT')
ylabel('Phase Angle (deg)')
grid on;
% -----------------------------------------------
subplot(212)
plot(tsplice, phi(idx2,:)); hold on;
plot(tsplice, zeros(length(tsplice),1), 'k')
axis([70 daystop -180 180])

ax = gca;
ax.XAxis.TickValues = [70:1/B:daystop];
oldtick = ax.XAxis.TickValues;
[~, indx] = min(abs(tsplice-oldtick));
ax.XTickLabel = num2str([0; UTsplice(indx(2:end))]);
title('Half Day Periodicity');
xlabel('UT')
ylabel('Phase Angle (deg)')
grid on;
saveas(h3, name3);


%% PLOT Recreated day and half-day signal
h4 = figure('units', 'normalized', 'position', [0 .5 1 1], 'visible', 'on');
B = 12;

% -----------------------------------------------
plot(tsplice, x_day, tsplice, x_halfday); hold on;
plot(tsplice, zeros(length(tsplice),1), 'k')
legend('x_{day}', 'x_{half day}');

title('Recreated TIEGCM Density Signal');
xlabel('Model Time [days]')
ylabel('x(t)')
grid on;
 saveas(h4, name4);
% -----------------------------------------------'

% zoomed in plots 
h5 = figure('units', 'normalized', 'position', [0 .5 1 1], 'visible', 'on');

% BEFORE STORM
daystart = 72;
daystop = 75;

subplot(211)
hold on;
yyaxis left
p1 = plot(tsplice, x_day, tsplice, x_halfday);
plot(tsplice, zeros(length(tsplice),1), 'k')
xlim([daystart daystop])

ax = gca;
ax.XAxis.TickValues = [70:1/B:daystop];
oldtick = ax.XAxis.TickValues;
[~, indx] = min(abs(tsplice-oldtick));
ax.XTickLabel = num2str([0; UTsplice(indx(2:end))]);
xlabel('UT')
ylabel('x(t)')
grid on;


yyaxis right 
ax2 = gca;
ax2.YColor = 'k';
ax2.YLim = [-1, 1];
ax2.YAxis(2).TickValues = linspace(-1,1,9);
ticks1 = linspace(12, 21, 4);
ticks2 = linspace(0, 12, 5);
ax2.YTickLabel = [ticks1, ticks2];
p2 = plot(tsplice, x_pole_noon, 'k');
ylabel('Position of N Mag Pole (LT)')
hold off;

legend([p1; p2], 'x_{day}', 'x_{half day}', 'Phase of N Mag Pole');
title({['Recreated TIEGCM Signal ZOOMED']; ['Model Day ', num2str(daystart), ' to ', num2str(daystop)]});


% AFTER STORM
daystart = 87;
daystop = 90;

subplot(212)
hold on;
yyaxis left
p1 = plot(tsplice, x_day, tsplice, x_halfday);
plot(tsplice, zeros(length(tsplice),1), 'k')
xlim([daystart daystop])

ax = gca;
ax.XAxis.TickValues = [70:1/B:daystop];
oldtick = ax.XAxis.TickValues;
[~, indx] = min(abs(tsplice-oldtick));
ax.XTickLabel = num2str([0; UTsplice(indx(2:end))]);
xlabel('UT')
ylabel('x(t)')
grid on;


yyaxis right 
ax2 = gca;
ax2.YColor = 'k';
ax2.YLim = [-1, 1];
ax2.YAxis(2).TickValues = linspace(-1,1,9);
ticks1 = linspace(12, 21, 4);
ticks2 = linspace(0, 12, 5);
ax2.YTickLabel = [ticks1, ticks2];
p2 = plot(tsplice, x_pole_noon, 'k');
ylabel('Position of N Mag Pole (LT)')
hold off;

legend([p1; p2], 'x_{day}', 'x_{half day}', 'Phase of N Mag Pole');
title({['Recreated TIEGCM Signal ZOOMED']; ['Model Day ', num2str(daystart), ' to ', num2str(daystop)]});

 saveas(h5, name5);






