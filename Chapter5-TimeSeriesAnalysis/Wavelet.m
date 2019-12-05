% this program tries to understand how wavelets work.
% close all;
% clear all;

% if the time series doesn't exist, then run the LOADfiles script.
load t_series.mat

%% CALL IGRF TO GET MAG. FIELD.
% create datenum for the model times 
DateNum = datenum(t_series.date(1,:), 1, t_series.date(2,:), t_series.date(3,:), t_series.date(4,:), 0);

Nlon = length(t_series.lon);
Nlat = length(t_series.lat);
Ntime = length(DateNum);
Bx = zeros(Nlat, Nlon);
By = Bx; Bz = Bx;

% Find Magnetic Field Coords for every lat/lon (just do for first time)
for t = 1:1
    for i = 1:Nlat
        for j = 1:Nlon
            [Bx(i,j), By(i,j), Bz(i,j)] = igrf(DateNum(t), t_series.lat(i), t_series.lon(j), t_series.Zp_latlon(j,i)./1000);
        end
    end   
end

%%
I = atand(Bz./hypot(Bx, By));  % inclination of field [deg]

[latSP_ind, lonSP_ind] = find( I == min(I(:)) );    % lat and lon indx of South mag pole (-90 inc)
[latNP_ind, lonNP_ind] = find( I == max(I(:)) );    % lat and lon indx of North mag pole (+90 inc)

t_series.lat(latSP_ind)
t_series.lat(latNP_ind)




%%  ---- REAL DATA ----
% TIEGCM is from year 2003.
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

name1 = './Figures/Powerspectrum_obs.png';
name2 = './Figures/Morlet_obs.png';
name3 = './Figures/MorletPhase_obs.png';
name4 = './Figures/RecreatedSignal.png';
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
h3 = figure('units', 'normalized', 'position', [0 .5 1 1], 'visible', 'off');
daystop = 72;
B = 12;
UTplt = [UT(1:stop); UT(start:end)];

% -----------------------------------------------
subplot(211)
plot(tsplice, phi(idx1,:)); hold on;
plot(tsplice, zeros(length(tsplice),1), 'k')
axis([70 daystop -180 180])

ax = gca;
ax.XAxis.TickValues = [70:1/B:daystop];
oldtick = ax.XAxis.TickValues;
[~, indx] = min(abs(tsplice-oldtick));
ax.XTickLabel = num2str([0; UTplt(indx(2:end))]);
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
ax.XTickLabel = num2str([0; UTplt(indx(2:end))]);
title('Half Day Periodicity');
xlabel('UT')
ylabel('Phase Angle (deg)')
grid on;
saveas(h3, name3);


%% PLOT Recreated day and half-day signal
h4 = figure('units', 'normalized', 'position', [0 .5 1 1], 'visible', 'on');
B = 12;
UTplt = [UT(1:stop); UT(start:end)];

% -----------------------------------------------
subplot(211)
plot(tsplice, x_day, tsplice, x_halfday); hold on;
plot(tsplice, zeros(length(tsplice),1), 'k')
legend('x_{day}', 'x_{half day}');

title('Recreated TIEGCM Density Signal');
xlabel('Model Time [days]')
ylabel('x(t)')
grid on;
% -----------------------------------------------
% zoomed in
daystart = 72;
daystop = 75;

subplot(212)
plot(tsplice, x_day, tsplice, x_halfday); hold on;
plot(tsplice, zeros(length(tsplice),1), 'k')
xlim([daystart daystop])
legend('x_{day}', 'x_{half day}');
ax = gca;
ax.XAxis.TickValues = [70:1/B:daystop];
oldtick = ax.XAxis.TickValues;
[~, indx] = min(abs(tsplice-oldtick));
ax.XTickLabel = num2str([0; UTplt(indx(2:end))]);
title({['Recreated TIEGCM Signal ZOOMED']; ['Model Day ', num2str(daystart), ' to ', num2str(daystop)]});
xlabel('UT')
ylabel('x(t)')
grid on;


 saveas(h4, name4);





