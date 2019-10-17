function [x_avg] = GlobeAreaAvg(x, lon, lat)
% this function weights the average of the value x by the solid angle on a
% sphere specific to your lat/lon location
% x = matrix of lon, lat values
% lon = longitude points (deg)
% lat = latitude points (deg)

N_lat = length(lat);        
N_lon = length(lon);

% make sure in terms of radians
dphi = (lon(2) - lon(1)) * pi/180;
dtheta = (lat(2) - lat(1)) * pi/180;

% transpose to make x = x(lon, lat)
x_size = size(x);
if x_size(1) ~= N_lon
   x = x';  
   fprintf('Transposing x array so that x = x(lon, lat).\n\n')
end


sum = 0;
for i = 1:N_lon
    for j = 1:N_lat
        dOmega = dphi*dtheta*cos(lat(j) * pi/180);      % solid angle for specific location [sr]
        sum = sum + x(i,j)*dOmega/(4*pi);   % work with radians
    end 
end

x_avg = sum/(N_lat*N_lon);
end

