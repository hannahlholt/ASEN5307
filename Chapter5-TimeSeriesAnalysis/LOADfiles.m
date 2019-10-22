%% Program to load all the TIEGCM files 

close all;
% clear all;
% run startup.m;
UtilitiesFolder = '/home/haho3703/MATLAB/ASEN5307/Utilities';

addpath('/home/haho3703/MATLAB/ASEN5307', UtilitiesFolder);
output = '~/MATLAB/TIEGCM_files/';
runtype = 'lowF107.lowtohighKp/';
filenames = dir([output, runtype, '*.nc']);
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
    name = [output, runtype, filenames(i).name];
    
    fprintf(['Loading File: ', name, '\n\n']);

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
step = 1;
% [Temp, den, mbar] in each column for every time step
globe_avg = zeros(length(modeltime)/step, 3);
mtime = zeros(length(modeltime)/step, 1);            % model time [day] for every global average point

ind = 1;
tot_t_ind = length(modeltime);

N_lat = length(lat);
N_lon = length(lon);


for t_ind = [1:step:tot_t_ind]
    fprintf(['Condensing to global average time series.\nCurrent time index: ', ...
        num2str(t_ind), ' out of ', num2str(tot_t_ind) ,'. \n\n']);
    
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
    % have to scale by the area (i.e. cos(latitude))
    for i=1:3
       if i == 1
           arr = squeeze(Tnx(:,:,z_want));
       elseif i == 2
           arr = squeeze(denx(:,:,z_want));
       else
           arr = squeeze(mbarx(:,:,z_want));
       end
        % take mean over lon and lat (i.e. rows and cols)
        globe_avg(ind,i) = GlobeAreaAvg(arr, lon, lat);
    end
    
    mtime(ind) = double(modeltime(1,t_ind))+double(modeltime(2,t_ind))/24+double(modeltime(3,t_ind))/1440;
    ind = ind+1;
   
end

%% now create struct and use this in other programs
t_series = struct('t', mtime, 'T', globe_avg(:,1), 'Den', globe_avg(:,2), ...
    'mbar', globe_avg(:,3), 'Zlvl', z_want);
%%
clearvars -except t_series
run startup.m
