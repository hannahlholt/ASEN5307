function [] = PLOT_TotalGas(zp, Z, lon, omega, omegaGrad, Tot_Mdiv, lon_want, lat_want, plotname, saveFig)

% % ###################################################################################
% % ###################################################################################
% 
% % since the geopotential altitude is different for every longitude value,
% % let your y value be the average of each column 

zp = zp/1000;               % we want to plot in km
y_avg = mean(zp, 1);        % row vector with average of each column.
x = lon;
[X, Y] = meshgrid(x, y_avg);        % mesh grid used for every subplot

ytop = y_avg(end-2);                % For N2, do end-5??
ybot = 100;
logic = y_avg >= ybot & y_avg <= ytop;

X = X(logic,:);              % needed to do this b/c had trouble plotting
Y = Y(logic,:);

num_cont = 300;
c = redblue(300);

thisfig = figure('units','normalized','outerposition',[0 1 .7 1]);     % left start, bottom start, side length, tall length

% ----------------
subplot(221)
Q = omega(:,logic)';      % omega for every longitude and alt for a specific latitute [m/s]
ztop = max(abs(Q(:)));  
zbot = -ztop;

colormap(c); 
contourf(X, Y, Q, num_cont, 'linecolor', 'none')
hold on;   
plot(lon_want*ones(length(y_avg), 1), y_avg, 'k--', 'Linewidth', 2.5)   
cbar = colorbar();
caxis([zbot ztop]); 
ylim([ybot ytop]);
cbar.Label.String = '\omega [1/s]';
xlabel('Longtitude [Deg]');
ylabel('Geopotential Altitude [km]');
title('TIEGCM Omega Output (ilev)')
grid on;

% ----------------
subplot(222)
Q = omegaGrad(:,logic)';
ztop = max(abs(Q(:)));  
zbot = -ztop;

colormap(c);
contourf(X, Y, Q, num_cont, 'linecolor', 'none')
hold on;   
plot(lon_want*ones(length(y_avg), 1), y_avg, 'k--', 'Linewidth', 2.5) 
cbar = colorbar();
caxis([zbot ztop]); 
ylim([ybot ytop]);
cbar.Label.String = '\partial \omega / \partial z [1/s]';
xlabel('Longtitude [Deg]');
ylabel('Geopotential Altitude [km]');
title('Gradient of Omega (ilev)');
grid on;

% ----------------
subplot(223)
Q = omega(:,logic)' - omegaGrad(:,logic)';
ztop = max(abs(Q(:)));  
zbot = -ztop;

colormap(c);
contourf(X, Y, Q, num_cont, 'linecolor', 'none')
hold on;   
plot(lon_want*ones(length(y_avg), 1), y_avg, 'k--', 'Linewidth', 2.5) 
cbar = colorbar();
caxis([zbot ztop]); 
ylim([ybot ytop]);
cbar.Label.String = '\omega - \partial \omega/ \partial z';

xlabel('Longtitude [Deg]');
ylabel('Geopotential Altitude [km]');
title('\omega - \partial \omega/ \partial z')
grid on;

% % ----------------
subplot(224)
Q = Tot_Mdiv(:,logic)';      % postive divergence = negative sign
ztop = max(abs(Q(:)));   
zbot = -ztop;

colormap(c);
contourf(X, Y, Q, num_cont, 'linecolor', 'none')
hold on;   
plot(lon_want*ones(length(y_avg), 1), y_avg, 'k--', 'Linewidth', 2.5) 
hold on;
for i =(14:4:55)
    plot(x, zp(:,i));
    
    %label the lines
    xpt = 150;
    ypt = double(zp(end, i));
    lbl = ['z = ', num2str(Z(i))];
    text(xpt, ypt, lbl, 'FontSize', 9, 'BackgroundColor', 'white', 'HorizontalAlignment', 'Center', 'Margin', 0.5)
    hold on;   
end
hold off;
cbar = colorbar();
caxis([zbot ztop]); 
ylim([ybot ytop]);
cbar.Label.String = 'Total Gas Mass Flux Divergence [1/s]';
xlabel('Longtitude [Deg]');
ylabel('Geopotential Altitude [km]');
title('Total Mass Flux Divergence (ilev)');
grid on;

sgtitle(['TIEGCM Model Run at UT = 0, Lat = ', num2str(lat_want), ', ', plotname]);

if saveFig ~= '0'
    saveas(thisfig, ['./Figures/Total_Gas/TIEGCM_LonMap_', saveFig, '.png'])  
end

end

