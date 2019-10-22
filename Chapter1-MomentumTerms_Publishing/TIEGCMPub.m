%% TIEGCM Continuity and Momentum Equation Analysis Program
% The programs and functions used to evaluate the different terms in the momentum equation for the thermosphere. Written by Hannah Holt. Last Updated 9/4/2019

%% -------------- MAIN ----------------------
%% Define Program Characteristics
close all;

output = '~/Documents/MATLAB/TIEGCM/TIEGCM_output/';
FigFolder = './Figures/';
addpath(output, FigFolder);

%----------------
ut_want = 1;        % what time segment desired from simulation
feat = 1;           % Select latitude and longitude desired
pdrag = 1;          % 1 if pdrag file used, 0 if not

% ----- Global Features -------
if feat == 1
lon_want = 95;          % North nighttime maximum feature
lat_want = 61.25;
savename = 'N_He_max';
plotname = 'HELIUM ENHANCEMENT AT 400 KM';
end
if feat == 2 
lon_want = 85;          % South nighttime maximum feature
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

%% Loading Viki's TIEGCM Simulation
if pdrag == 1
    filename = [output, 'HSUVW.tiegcm2.0_dres.pdrag_f107_180_001.nc'];
    id = 'pdrag';
end
if pdrag == 0
    filename = [output, 'HSUVW.tiegcm2.0_dres.nodragtest_ctrSS_f107_180_001.nc'];
    id = 'no Ion Drag';
end

% Parse Data from file
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
altPts = size(z_ilev);
lonPts = size(lon);

he = ncread(filename,'HE');             % Units of mass mixing ratio
n2 = ncread(filename,'N2');
o2 = ncread(filename,'O2');
o1 = ncread(filename,'O1');

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

%% Initialize Objects
N2 = TIEGCMspecies('N2', 28.01, n2, den, Z);
O2 = TIEGCMspecies('O2', 32, o2, den, Z);
O1 = TIEGCMspecies('O1', 16, o1, den, Z);
He = TIEGCMspecies('He', 4, he, den, Z);


%% Main Function Calls
% Diffusion Coefficients for each species
DiffCoeff(Tn, P, N2, O2, O1, He)

% Thermal diffusion coefficient
D_Ti = -0.36 * He.Di; 

% mean molecular mass on ilevs [kg/mol]
mbar = (N2.mmr/N2.weight + O2.mmr/O2.weight + O1.mmr/O1.weight + He.mmr/He.weight).^-1;        

% Scale Heights for each species and for general atm. 
[H_P, H_T, H_T_He, H_rho_diff, H_rho_star] = ScaleHeight(N2, O2, O1, He, Tn, mbar, zp, den, g0, Av, kb);

% Omega "winds" and needed gradients
omega = wn./H_P;                   % TIEGCM omega for every lon/alt [1/s]
omegaGrad = zeros(size(omega));    % gradient of omega w.r.t Z
omegaExpGrad = omegaGrad;          % gradient of omega * exp(-Z) w.r.t Z
for l = 1:lonPts
   omegaGrad(l,:) = ThreePtGrad( Z(l,:), omega(l,:) );
   omegaExpGrad(l,:) = ThreePtGrad( Z(l,:), exp(-Z(l,:)) .* omega(l,:) );
end

% Vertical Advection Terms and assign to Objects
VertAdvection(omega, Z, p0, g0, N2, O2, O1, He)

% Horizontal Mass Flux Divergence Terms for objecta and total gas 
Tot_Mdiv = HorMassFluxDivergence(omegaExpGrad, Z, p0, g0, N2, O2, O1, He);       


%% Plot Total Gas Features
saveFig = '0';
PLOT_TotalGas(zp, z_ilev, lon, omega, omegaGrad, Tot_Mdiv, lon_want, lat_want, plotname, saveFig)

%% Plot Specific Species Behavior
PLOT_Species(zp, z_ilev, lon, lon_want, lat_want, omega, plotname, savename, He)

%% -------------- CLASSES ---------------------
%% A 'TIEGCMspecies' object
% <include>TIEGCMspecies.m</include>
%

%% -------------- FUNCTIONS ------------------
%% -- Diffusion Coefficient Calculation
%
% <include>DiffCoeff.m</include>
%
% <include>Di_TIEGCM.m</include>
%


%% -- Horizontal Total Mass Flux Divergence Calculation
%
% <include>HorMassFluxDivergence.m</include>
%


%% -- Vertical Advection Calculation
%
% <include>VertAdvection.m</include>
%

%% -- Scale Height Calculation
%
% <include>ScaleHeight.m</include>
%

%% -- Converting from lev to ilevs
%
% <include>CONVERT2ILEV.m</include>
%


%%
% Hannah Holt 
% CU BOULDER Aerospace Engineering Sciences.
% Thank you.