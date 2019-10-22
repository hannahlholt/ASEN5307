%% Chap 3 fitting and analyzing residuals of the ideal gas equation

close all;
% clear all;
% run startup.m;

day_want = 79;

output = '~/MATLAB/TIEGCM_files/';
runtype = 'lowF107.lowtohighKp/';
UtilitiesFolder = '/home/haho3703/MATLAB/ASEN5307/Utilities';
FigFolder = './Figures/';
filenames = dir([output, runtype, '*.nc']);
addpath(UtilitiesFolder, FigFolder);
%%
% since this is the 5 deg res, lon pts = 72, lat pts = 36, and Z = 29
% den(lon, lat, Z, time) = (72, 36, 29, 72)
modeltime = [];
Kp = [];
den = [];
Tn = [];
n2 = [];
o2 = [];
o1 = [];
he = [];

for i=1:length(filenames)
    name = [output, runtype, filenames(i).name]

    modeltime = [modeltime, ncread(name,'mtime')];  % (day, hr, min)
    Kp = [Kp; ncread(name, 'Kp')];                  % Kp index
    
    den = cat(4, den, ncread(name,'DEN')/1e3);       % total density [kg/m^3], ILEV
    Tn = cat(4, Tn, ncread(name, 'TN'));             % Temperature [K], LEV
    n2 = cat(4, n2, ncread(name, 'N2'));             % N2 mmr, LEV
    o2 = cat(4, o2, ncread(name, 'O2'));              
    o1 = cat(4, o1, ncread(name, 'O1'));              
    he = cat(4, he, ncread(name, 'HE'));             
    
    % These dont change for the runs
    if i == 1
        Z = ncread(name, 'ilev');           % interface pressure level
        UT = ncread(name, 'ut');                % UT time from model time hr and minute [hrs]  
        lat = ncread(name, 'lat');
        lon = ncread(name, 'lon');
        p0 = ncread(name, 'p0_model') / 10; % p0 reference pressure used by model [Pascals]
        g0 = ncread(name, 'grav') / 100;    % const. gravitational acceleration [m/s]
    end
end
 
P = p0 .* exp(-Z);                      % pressure array on ilevs [Pa]

%% condense density to specific pressure level
z_want = 22;
hr_want = 12;

Z_mat = repmat(Z, [1, 72, 36]);   
Z_mat = permute(Z_mat, [2 3 1]);

% sample every 2 hrs
skip = 5;
% [Temp, den, mbar] in each column for every time step
globe_avg = zeros(length(modeltime)/skip, 3);
ind = 1;
for t_ind = [1:skip:length(modeltime)]
    t_ind
    denx = squeeze(den(:, :, :, t_ind));        
    n2x = squeeze(n2(:, : , :, t_ind));
    o2x = squeeze(o2(:, :, :, t_ind));
    o1x = squeeze(o1(:, :, :, t_ind));
    hex = squeeze(he(:, :, :, t_ind));
    Tnx = squeeze(Tn(:, :, :, t_ind));            

    % put Temperature on ilevs
    lonPts = length(Tnx(:,1,1));
    latPts = length(Tnx(1,:,1));
    altPts = length(Tnx(1,1,:));
    slicePts = [lonPts, latPts];
    Tnx = CONVERT2ILEV(Tnx, slicePts, altPts);

    % convert everything to ilevs
    %----- Initialize your OBJECTS ---------
    N2 = TIEGCMspecies('N2', 28.01, n2x, denx, Z_mat);
    O2 = TIEGCMspecies('O2', 32, o2x, denx, Z_mat);
    O1 = TIEGCMspecies('O1', 16, o1x, denx, Z_mat);
    He = TIEGCMspecies('He', 4, hex, denx, Z_mat);
    mbarx = (N2.mmr/N2.weight + O2.mmr/O2.weight + O1.mmr/O1.weight + He.mmr/He.weight).^-1;     

    % now find global mean of Temperature, mbar and den at ~ 400 km
    
    globe_avg(ind,1) = mean(Tnx(:,:,z_want), 'all');      % take mean over lon and lat (i.e. rows and cols)
    globe_avg(ind,2) = mean(denx(:,:,z_want), 'all');
    globe_avg(ind,3) = mean(mbarx(:,:,z_want), 'all');
    ind = ind+1;
   
end

%%
% Calculate fits and residuals of Ideal Gas equation
x = 1./globe_avg(:,1);
y1 = globe_avg(:,2);

y2 = globe_avg(:,2)./globe_avg(:,3);

[p1, s1] = polyfit(x, y1, 1);
[p2, s2] = polyfit(x, y2, 1);
[y1_fit, std1] = polyconf(p1, x, s1, 'alpha', 0.05);
[y2_fit, std2] = polyconf(p2, x, s2, 'alpha', 0.05);
res1 = y1 - y1_fit;
res2 = y2 - y2_fit; 

thisfig = figure('units', 'normalized', 'outerposition',[0 1 .6 1]);
subplot(221)
plot(x, y1, 'r.'); hold on;
plot(x, y1_fit, 'k-');
plot(x, y1_fit+std1, 'k--', x, y1_fit-std1, 'k--')
xlabel('Global 1/T average [1/K]');
ylabel('\rho [kg/m^3]');
legend('Data', ['y = ', sprintf('%.3e', p1(1)), 'x + ',  sprintf('%.2f', p1(2))], '+/- 1\sigma', 'location', 'best' )
title('Global \rho Average');
grid on;

subplot(222);
plot(x, y2, 'b.'); hold on;
plot(x, y2_fit, 'k-');
plot(x, y2_fit+std2, 'k--', x, y2_fit-std2, 'k--');
legend('Data', ['y = ', sprintf('%.3e', p2(1)), 'x + ',  sprintf('%.2f', p2(2))], '+/- 1\sigma', 'location', 'best' )
xlabel('Global Average of 1/T [1/K]');
ylabel('\rho / mbar');
title('Global Average of \rho / mbar')
grid on;

subplot(223);
num_bins = round(sqrt(length(res1)));
xmin = min(res1);
xmax = max(res1);
binspace = linspace(xmin, xmax, num_bins);
h1 = histogram(res1, binspace);
title('Normalized Postfit Residuals');
% xlim([xmin xmax])
grid on;

subplot(224);
xmin = min(res2);
xmax = max(res2);
binspace = linspace(xmin, xmax, num_bins);
h2 = histogram(res2, binspace);
title('Normalized Postfit Residuals');
% xlim([xmin xmax]);
grid on;

saveas(thisfig, [FigFolder, 'IdealGasFits.png'])

% perform chi2 test on residuals
bin_mid1 = h1.BinEdges(1:end-1) + h1.BinWidth;
bin_mid2 = h2.BinEdges(1:end-1) + h2.BinWidth;

n_syn1 = normpdf(bin_mid1, mean(res1), std(res1));
n_syn1 = (n_syn1 ./ sum(n_syn1)) .* sum(h1.Values);

n_syn2 = normpdf(bin_mid2, mean(res2), std(res2));
n_syn2 = (n_syn2 ./ sum(n_syn2)) .* sum(h2.Values);

%%
chicalc2 = sum((h2.Values - n_syn2).^2 ./n_syn2 )
chi2crit = chi2inv(0.95, num_bins -(1 + 2))

%%


