% this program tries to understand how wavelets work.
close all;
% clear all;

addpath(genpath('~/Documents/MATLAB/ASEN5307/Utilities'))
var = 'Density';

try 
    load t_series.mat
catch
    LOADfiles
end

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

%% FIND LAT/LON OF N AND S MAG FIELD
I = atand(Bz./hypot(Bx, By));  % inclination of field [deg]

[latSP_ind, lonSP_ind] = find( I == min(I(:)) );    % lat and lon indx of South mag pole (-90 inc)
[latNP_ind, lonNP_ind] = find( I == max(I(:)) );    % lat and lon indx of North mag pole (+90 inc)

% lat and lon of North and South geomag pole in 2003. 
% lat = [90, -90]
% lon = [0, 360] deg East
Bnorth = [t_series.lat(latNP_ind), t_series.lon(lonNP_ind)];
Bsouth = [t_series.lat(latSP_ind), t_series.lon(lonSP_ind)];


%%  ---- REAL DATA ----
% TIEGCM is from year 2003.
DAT = t_series;

switch var
    case 'Density' 
        x_tot = DAT.Den;       % total density (global average)
    case 'Temperature'
        x_tot = DAT.T;         % temperature
    case 'Joule Heating'
        x_tot = DAT.QJoule;    % temperature
    otherwise
        error('Bad variable name.');
end

t = DAT.t;
UT = DAT.UT;            % universal time
dtN = t(2) - t(1);   % dt is in hours
fs = 1/dtN;
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

name0 = ['./Figures/RealVsSyn_',var,'.png'];
name1 = ['./Figures/Powerspectrum_obs_',var,'.png'];
name2 = ['./Figures/Morlet_obs_',var,'.png'];
name3 = ['./Figures/MorletPhase_obs_',var,'.png'];
name4 = ['./Figures/SynSignal_',var,'.png'];
name5 = ['./Figures/SynSignalZOOM_',var,'.png'];
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
[wt, period, coi] = cwt(x, 'amor', days(dtN));

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
dT = 3;
[~, idx1] = min(abs(period-n));
% [~, trueMAXind] = max(Z_og(idx1-dT:idx1+dT,1));
% idx1 = idx1 + (-dT + trueMAXind);
idx1 = idx1 - 3;

% 2) half day oscillation
n = days(0.5);
% index in time space (rows) of half day oscillation
[~, idx2] = min(abs(period-n));
idx2 = idx2 - 2;

phi = 180/pi*(angle(wt));            % phase angle [-180, 180]

% --- RECREATE SIGNAL -----

x_day = Z_og(idx1,:) .* sin( 2*pi./days(period(idx1)) + pi/180.*phi(idx1,:));
x_halfday = Z_og(idx2,:) .* sin( 2*pi./days(period(idx2)) + pi/180.*phi(idx2,:));

% --- create Phase of the North pole w.r.t noon
lon_Npole = mod(Bnorth(2), 360);                %  degrees E, N pole at 2003
dtN = lon_Npole/15 - 12;
phiN = mod(2*pi*dtN/24, 2*pi);

% sine wave of N pole w.r.t noon. (MAX when pole is at noon!!)
N_pole_noon = cos(2*pi/24*UTsplice + phiN); 

% --- create Phase of the South pole w.r.t noon
lon_Spole = mod(Bsouth(2), 360);                %  degrees E, S pole at 2003
dtS = lon_Spole/15 - 12;
phiS = mod(2*pi*dtS/24, 2*pi);

% sine wave of S pole w.r.t noon. (MAX when pole is at noon!!)
S_pole_noon = cos(2*pi/24*UTsplice + phiS); 


%% ------ PLOTTING -------
h0 = figure('units', 'normalized', 'position', [0 .5 1 1], 'visible', 'off');
subplot(211)
plot(tsplice, x)
hold on
plot(tsplice, x_day + x_halfday)
xlim([70, 80])
legend('TIEGCM', 'Synthetic')
title([var, ' Oscilations Before Storm'])
grid on;

subplot(212)
plot(tsplice, x)
hold on
plot(tsplice, x_day + x_halfday)
xlim([82, 100]);
legend('TIEGCM', 'Synthetic')
title([var, ' Oscilations After Storm'])
grid on;
saveas(h0, name0);


%% PLOT 1D POWER SPECTRUM
h1 = figure('visible', 'off');
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
saveas(h1, name1);
%

%% PLOT WAVELET ANALYSIS
h2 = figure('visible', 'off');
contourf(X, Y, Z, 100, 'linecolor', 'none')
hold on;
colormap jet;
cbar = colorbar();
set(gca, 'YScale', 'log')
ylabel(['Period ', xname]);
xlabel(['Time ', xname]);
cbar.Label.String = 'Magnitude';
grid on;
saveas(h2, name2);


%% PLOT PHASE
h3 = figure('units', 'normalized', 'position', [0 .5 1 1], 'visible', 'off');
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
h4 = figure('units', 'normalized', 'position', [0 .5 1 1], 'visible', 'off');
B = 12;

% -----------------------------------------------
plot(tsplice, x_day, tsplice, x_halfday); hold on;
plot(tsplice, zeros(length(tsplice),1), 'k')
legend('x_{day}', 'x_{half day}');

title('Recreated TIEGCM Density Signal');
xlabel('Model Time [days]')
ylabel(var)
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
ylabel(var)
grid on;


yyaxis right 
ax2 = gca;
ax2.YColor = 'k';
ax2.YLim = [-1, 1];
ax2.YAxis(2).TickValues = linspace(-1,1,9);
ticks1 = linspace(12, 21, 4);
ticks2 = linspace(0, 12, 5);
ax2.YTickLabel = [ticks1, ticks2];
p2 = plot(tsplice, N_pole_noon, 'k', tsplice, S_pole_noon, 'k');
ylabel('Position of N/S Mag Pole (LT)')
hold off;

legend([p1; p2], 'x_{day}', 'x_{half day}', 'Phase of N Mag Pole', 'Phase of S Mag Pole');
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
ylabel(var)
grid on;


yyaxis right 
ax2 = gca;
ax2.YColor = 'k';
ax2.YLim = [-1, 1];
ax2.YAxis(2).TickValues = linspace(-1,1,9);
ticks1 = linspace(12, 21, 4);
ticks2 = linspace(0, 12, 5);
ax2.YTickLabel = [ticks1, ticks2];
p2 = plot(tsplice, N_pole_noon, 'k', tsplice, S_pole_noon, 'k');
ylabel('Position of N/S Mag Pole (LT)')
hold off;

legend([p1; p2], 'x_{day}', 'x_{half day}', 'Phase of N Mag Pole', 'Phase of S Mag Pole');
title({['Recreated TIEGCM Signal ZOOMED']; ['Model Day ', num2str(daystart), ' to ', num2str(daystop)]});

 saveas(h5, name5);






