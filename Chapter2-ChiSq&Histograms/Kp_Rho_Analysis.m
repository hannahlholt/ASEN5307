% Using a version of ASEN5307_test.m to check 

close all;

output = '~/MATLAB/TIEGCM_files/';
runtype = 'lowF107.lowtohighKp/';
files = dir([output, runtype, '*.nc']);
FigFolder = './Figures/';
addpath(FigFolder);

day_want = 80;
%%
modeltime = [];
Kp = [];
he = [];
den = [];
Z = [];
UT = [];
wn = [];
Tn = [];

% since this is the 5 deg res, lon pts = 72, lat pts = 36, and Z = 29
% den(lon, lat, Z, time) = (72, 36, 29, 72)

for i=1:length(files)
    name = [output, runtype, files(i).name];

    modeltime = [modeltime, ncread(name,'mtime')];  % (day, hr, min)
    Kp = [Kp; ncread(name, 'Kp')];                  % Kp index
    
    den = cat(4, den, ncread(name,'DEN')/1e3);          % total density [kg/m^3], ILEV
%     he = cat(4, he, ncread(name, 'HE'));              % helium mmmr
    Tn = cat(4, Tn, ncread(name, 'TN'));

    % These dont change for the runs
    if i == 1
        Z = ncread(name, 'ilev');              % interface pressure level
        Zp = ncread(name, 'Z')/100;            % geopotential height [m]
        UT = ncread(name, 'ut');               % UT time from model time hr and minute [hrs]  
        lat = ncread(name, 'lat');
        lon = ncread(name, 'lon');
    end
end

%% condense density to specific pressure level
z_want = 22;                % ~400 km
day_want = 81;
hr_want = 12;

t_ind = find(modeltime(1,:) == day_want & modeltime(2,:) == hr_want);

den_z = squeeze(den(:,:,z_want,:));
den_z_ut = squeeze(den(:,:,z_want,t_ind(1)));
tn_z_ut = squeeze(Tn(:,:,z_want,t_ind(1)));


%% create a synthetic normal gaussian and normalize
den_data = den_z_ut(:);
num_bins = round(sqrt(length(den_data)));
xmin = min(den_data);
xmax = max(den_data);
binspace = linspace(xmin, xmax, num_bins);

figure();
h = histogram(den_data, binspace);
grid on;
line1 = ['Global density distribution @ Z = ', num2str(Z(z_want)), ', (~400 km)'];
line2 = ['Day ', num2str(day_want), ', ',  num2str(hr_want), ':00'];
title({line1, line2});
xlabel('Total Mass Denisty, \rho [kg/m^3]');
n_meas = h.Values;
saveas(gcf, [FigFolder, 'GlobalDensityHistogram.png'])


midpoints = binspace(1:end-1)+h.BinWidth;

% evaluate Gaussian at the midpoints of each bin
n_syn = normpdf(midpoints, mean(den_data), std(den_data));
n_syn = n_syn / sum(n_syn);                % normalize to a total of one
n_syn = sum(n_meas) * n_syn;               % scale to data meas

figure()
subplot(1,2,1), bar(midpoints, n_syn, 'r');
title('Synthetic Gaussian for \chi^2')
subplot(1,2,2), bar(midpoints, n_meas, 'b');
title('Actual Data for \chi^2')

%% Chi^2 test
chi2calc = sum(((n_meas - n_syn).^2) ./ n_syn);
DegofFr = h.NumBins - (2 + 1);
ConfLvl = 0.95;  %(2 sigma confidence lvl)
chi2crit = chi2inv(ConfLvl, DegofFr)
if chi2calc > chi2crit
   Conclusion = 'Reject null hypothesis!'
else
    Conclusion = 'Do not reject null hypothesis without another cause!'
end

[h, p, stats] = chi2gof(den_data);


%% Compare global temp and density at specific UT and pressure lvl (=
% 22) to see if they are associated
temp_data = tn_z_ut(:);

figure('Visible', 'on');
subplot(121);
num_bins = round(sqrt(length(den_data)));
xmin = min(den_data);
xmax = max(den_data);
binspace = linspace(xmin, xmax, num_bins);
h1 = histogram(den_data, binspace);
grid on;
title(['Day ', num2str(day_want), ', ',  num2str(hr_want), ':00']);
xlabel('Total Mass Denisty, \rho [kg/m^3]');
n_meas1 = h1.Values;
bin_mid1 = h1.BinEdges(1:end-1) + h1.BinWidth;

subplot(122);
xmin = min(temp_data);
xmax = max(temp_data);
binspace = linspace(xmin, xmax, num_bins);
h2 = histogram(temp_data, binspace);
grid on;
title(['Day ', num2str(day_want), ', ',  num2str(hr_want), ':00']);
xlabel('Neutral Temperature, T [K]');
n_meas2 = h2.Values;
bin_mid2 = h2.BinEdges(1:end-1) + h2.BinWidth;

sgtitle(['Global distributions @ Z = ', num2str(Z(z_want)), ', (~400 km)'])
saveas(gcf, [FigFolder, 'Dens&TempHistogram.png'])


%% 2D Chi Squared Test - not working...
% now we want to make a table containing the number of counts for a
% specific temp/den point
matrix_den = repmat(n_meas1', 1, h2.NumBins);
matrix_temp = repmat(n_meas2, h1.NumBins, 1);

den_temp_counts = matrix_den + matrix_temp; % number of counts in each den/temp bin;
row_tot = sum(den_temp_counts, 2);
col_tot = sum(den_temp_counts, 1);
total = sum(row_tot);

% the expected counts in each den/temp point if they were independent of each other
den_temp_exp = ones(size(den_temp_counts)) .* round(row_tot/total .* col_tot/total .* total);     

chi2stat = sum( (den_temp_exp - den_temp_counts).^2 ./ den_temp_exp, 'all')
DOF = (h1.NumBins - 1)*(h2.NumBins - 1);        % degrees of freedom!
alpha = 0.05;
ConfLvl = 1-alpha;
chi2crit = chi2inv(ConfLvl, DOF)

if chi2stat > chi2crit
   Conclusion = 'Reject null hypothesis! The variables are associated.'
else
    Conclusion = 'Do not reject null hypothesis without another cause!. Variables are independent.'
end


% cross tab tests that x,y are independent of each other
% - Null Hypoth - they are independent.  If P-value < 0.05, then reject the
% null!!
% they are inversely correlated as shown by r12 and r21 values.
[r, p] =  corrcoef(den_data, temp_data);


%% How density and Kp at z = 3.5 changes as the storm approaches

day_want = [70:1:99];

samps = length(modeltime(1,:));
TOD = double(modeltime(1,:))+double(modeltime(2,:))/24+double(modeltime(3,:))/1440;       % time of model day

% find global mean of density at the z = 3.5 level and compare to Kp index
rho_mean = squeeze(den(:,:,z_want,:));
rho_mean = squeeze(mean(rho_mean, [1 2 ]));

figure()
subplot(2,1,1)
plot(TOD, rho_mean)
xlim([day_want(1), day_want(end)]);
title(['Global Mean Density @ Z = ', num2str(Z(z_want)), ', (~400 km)']);
ylabel('Total Mass Density , \rho [kg/m^3]')
grid on;

subplot(2,1,2)
plot(TOD, Kp)
axis([day_want(1), day_want(end), 0, 5]);
title('Kp Index');
xlabel('Model Day');
grid on;

saveas(gcf, [FigFolder, 'Kp_DenVsTime.png'])









