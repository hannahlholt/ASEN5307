
close all;


TIEGCMoutput = '/home/haho3703/MATLAB/TIEGCM_files/';
UtilitiesFolder = '/home/haho3703/MATLAB/ASEN5307/Utilities';
FigFolder = './Figures/';
addpath(TIEGCMoutput, FigFolder, UtilitiesFolder);
run startup.m

%----------------
ut_want = 1;        % what time segment desired from simulation
feat = 7;           % Select latitude and longitude desired
pdrag = 1;          % 1 if pdrag file used, 0 if not, 2 = special case
day_want = 81;      % storm hits on day 80 at 00:00 hrs - use on Erics Runs
res = 2.5;          % simulation resolution

% ----- Global Features -------
% for paper
if feat == 7
    lon_want = [90, -72.5];     % two points for paper -[MAX, MIN] latitude is in between south max and min
    lat_want = -53.75; 
    savename = 'S_He_Extrema';
    plotname = 'HELIUM EXTREMA AT 400 KM';
end

if feat == 1
lon_want = 95;          % North nighttime maximum feature
lat_want = 61.25;
savename = 'N_He_max';
plotname = 'HELIUM ENHANCEMENT AT 400 KM';
end
if feat == 2 
lon_want = 100;          % South nighttime maximum feature
lat_want = -58.75;
savename = 'S_He_max';
plotname = 'HELIUM ENHANCEMENT AT 400 KM';
end
if feat == 3 
lon_want = -77.5;       % North Daytime minimum feature
lat_want = 18.75;
savename = 'N_He_min';
plotname = 'HELIUM DEPLETION AT 400 KM';
end
if feat == 4
lon_want = -77.5;       % South Daytime minimum feature
lat_want = -46.25;      % true minimum is at -13.75 in the Southern hemi!!
savename = 'S_He_min';
plotname = 'HELIUM DEPLETION AT 400 KM';
end
if feat == 5        % N helium max in the 5deg resolution
lon_want = 100;
lat_want = 47.5;
savename = 'N_He_max';
plotname = 'HELIUM ENHANCEMENT AT 400 KM';   
end

%-----Loading Viki's tiegcm simulation-----
% format is (lon,lat,ilev,UT) 
if pdrag == 1
    filename = [TIEGCMoutput, 'HSUVW.tiegcm2.0_dres.pdrag_f107_180_001.nc'];
    id = 'pdrag';
end
if pdrag == 0
    filename = [TIEGCMoutput, 'HSUVW.tiegcm2.0_dres.nodragtest_ctrSS_f107_180_001.nc'];
    id = 'no Ion Drag';
end
if pdrag == 2
    filename = [TIEGCMoutput, 'lowF107.lowtohighKp/s0', num2str(day_want),'.nc'];
end
% ----------- Parse Data from file -------------- 
den = ncread(filename,'DEN') / 1e3;     % total density [kg/m^3] ILEV
zp = ncread(filename, 'Z') / 100;       % geopotential height [m] on ILEV
p0 = ncread(filename, 'p0_model') / 10; % p0 reference pressure used by model [Pascals]
g0 = ncread(filename, 'grav') / 100;    % const. gravitational acceleration [m/s]
z_ilev = ncread(filename, 'ilev');      % interface pressure level
Tn = ncread(filename,'TN');             % temperature on levs [K] 
lat = ncread(filename,'lat');           % latitude points
lon = ncread(filename,'lon');           % longitude points [-180, 180]
wn = ncread(filename,'WN')/100;         % neutral vertical winds on ilev [m/s]

% make a 144x57 matrix with each whole column being one of the pressure ilevs. 
% (i.e. rows are all the same for a given column)
Z = repmat(z_ilev', size(den,1), 1 );   

P = p0 .* exp(-Z);                      % pressure array on ilevs [Pa]
altPts = length(z_ilev);
lonPts = length(lon);

he = ncread(filename,'HE');             % Units of mass mixing ratio
n2 = ncread(filename,'N2');
o2 = ncread(filename,'O2');
o1 = ncread(filename,'O1');

he2 = squeeze(he(:, :, 22, ut_want+1));     % to plot helium at 400 km

% Condense to UT time and latitude  
i = find(lat == lat_want);

den = squeeze(den(:, i, :, ut_want+1));        
n2 = squeeze(n2(:, i , :, ut_want+1));
o2 = squeeze(o2(:, i, :, ut_want+1));
o1 = squeeze(o1(:, i, :, ut_want+1));
he = squeeze(he(:, i, :, ut_want+1));
zp = squeeze(zp(:, i, :, ut_want+1));
Tn = squeeze(Tn(:, i, :, ut_want+1));            
wn = squeeze(wn(:, i, :, ut_want+1)); 

% put Temperature on ilevs
Tn = CONVERT2ILEV(Tn, lonPts, altPts);

%----- Initialize your OBJECTS ---------
N2 = TIEGCMspecies('N2', 28.01, n2, den, Z);
O2 = TIEGCMspecies('O2', 32, o2, den, Z);
O1 = TIEGCMspecies('O1', 16, o1, den, Z);
He = TIEGCMspecies('He', 4, he, den, Z);

% Calculate Diffusion Coefficients for each species
DiffCoeff(Tn, P, N2, O2, O1, He)

% Thermal diffusion coefficient
D_Ti = -0.38 * He.Di; 

% ---------- Calculations --------------

% mean molecular mass on ilevs [kg/mol]
mbar = (N2.mmr/N2.weight + O2.mmr/O2.weight + O1.mmr/O1.weight + He.mmr/He.weight).^-1;        

% Scale Heights for each species and for general atm. 
[H_P, H_T, H_T_He, H_rho_diff, H_rho_star] = ScaleHeight(N2, O2, O1, He, Tn, mbar, zp, den, g0, Av, kb);

% omega "winds" and needed gradients
omega = wn./H_P;                   % TIEGCM omega for every lon/alt [1/s]
omegaGrad = zeros(size(omega));    % gradient of omega w.r.t Z
omegaExpGrad = omegaGrad;          % gradient of (omega * exp(-Z)) w.r.t Z
for l = 1:lonPts
   omegaGrad(l,:) = ThreePtGrad( Z(l,:), omega(l,:) );
   omegaExpGrad(l,:) = ThreePtGrad( Z(l,:), exp(-Z(l,:)) .* omega(l,:) );
end

% Vertical Advection Terms and assign to Objects
VertAdvection(omega, Z, p0, g0, N2, O2, O1, He)

% Horizontal Mass Flux Divergence Terms for species and total gas 
Tot_Mdiv = HorMassFluxDivergence(omegaExpGrad, Z, p0, g0, N2, O2, O1, He);       

% -------------------------------------

%% ------------- PLOTTING -------------------------------------------------
% -------------------------------------------------------------------------

%% 1 - Plot Total Gas Features
x_label = 'Longitude';
% saveFig = '0';
% PLOT_TotalGas(res, x_label, zp, z_ilev, lon, lon_want, lat_want, omega, omegaGrad, Tot_Mdiv, plotname, saveFig)

%% 2 - Plot Specific Species Behavior
close all;
saveFig = savename;
factor = P/g0; 
PLOT_Species(res, x_label, zp, z_ilev, lon, lon_want, lat_want, factor, omega, omegaGrad, plotname, saveFig, He)

%% 3 - Plot MMRs of N2, O2, O, and He
i = find(lon == lon_want(1));

thisfig = figure();
plot(N2.mmr(i,:), zp(i,:)/1000, 'Linewidth', 1.5);
hold on
plot(O2.mmr(i,:), zp(i,:)/1000, 'Linewidth', 1.5);
hold on
plot(O1.mmr(i,:), zp(i,:)/1000, 'Linewidth', 1.5);
hold on
plot(He.mmr(i,:), zp(i,:)/1000, 'Linewidth', 1.5);

legend('N_2','O_2','O','He', 'location', 'Northwest');
title(['TIEGCM Mass Mixing Ratios for Lat = ', num2str(lat_want) , ', Lon = ', num2str(lon_want(1))]);
ylim([120 550]);
set(gca, 'XScale', 'log');
xlabel('Mass Mixing Ratio');
ylabel('Geopotential Altitude [km]');
grid on;
saveas(thisfig, ['./Figures/Helium_Paper/MMRs_', savename, '.svg'])

%% 4 - Plot global Helium Density at ~400 km on Pressure lvl
y = lat;        % row vector with average of each column.
x = lon;
[X, Y] = meshgrid(x, y);        % mesh grid used for every subplot
Z = he2';

mag_eq = importdata('Magnetic_equator_lat_lon.txt'); % first column is lon, second column is lat

figure();
contourf(X, Y, Z, 25)
colormap(jet(300))
title(['Helium Denisty at Zp = ', num2str(z_ilev(22)) ,', UT = 0, Storm Day ', num2str(day_want)])
cbar = colorbar();
cbar.Label.String = '\rho_{He} [kg/m^3]';
    
hold on;   
plot(mag_eq(:,1), mag_eq(:,2), 'r', 'Linewidth', 2.5)

xlabel('LT [Hr]')
ylabel('Lat [Deg]')
xticks(linspace(-180, 180, 9))
xticklabels({'12', '15', '18', '21', '0', '3', '6', '9', '12'})

%% Extras 
% make sure the lat and lon max is correct!
[lon_max, lat_max] = find(he2 == max(he2(:,15), [], 'all'))
lon(lon_max)
lat(lat_max)


