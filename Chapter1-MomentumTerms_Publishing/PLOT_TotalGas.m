function [] = PLOT_TotalGas(res, x_label, zp, Z, x_slc, xline_want, other_want, omega, omegaGrad, Tot_Mdiv, plotname, saveFig)

% % ###################################################################################
% % ###################################################################################
% 
% % since the geopotential altitude is different for every longitude value,
% % let your y value be the average of each column 

zp = zp/1000;               % we want to plot in km
y_avg = mean(zp, 1);        % row vector with average of each column.
x = x_slc;
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
plot(xline_want(1)*ones(length(y_avg), 1), y_avg, 'r--', 'Linewidth', 2.5);
plot(xline_want(2)*ones(length(y_avg), 1), y_avg, 'b--', 'Linewidth', 2.5);
cbar = colorbar();
caxis([zbot ztop]); 
ylim([ybot ytop]);
cbar.Label.String = '\omega [1/s]';
xlabel([x_label, ' [Deg]']);
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
plot(xline_want(1)*ones(length(y_avg), 1), y_avg, 'r--', 'Linewidth', 2.5);
plot(xline_want(2)*ones(length(y_avg), 1), y_avg, 'b--', 'Linewidth', 2.5); 
cbar = colorbar();
caxis([zbot ztop]); 
ylim([ybot ytop]);
cbar.Label.String = '\partial \omega / \partial z [1/s]';
xlabel([x_label, ' [Deg]']);
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
plot(xline_want(1)*ones(length(y_avg), 1), y_avg, 'r--', 'Linewidth', 2.5);
plot(xline_want(2)*ones(length(y_avg), 1), y_avg, 'b--', 'Linewidth', 2.5);  
cbar = colorbar();
caxis([zbot ztop]); 
ylim([ybot ytop]);
cbar.Label.String = '\omega - \partial \omega/ \partial z';

xlabel([x_label, ' [Deg]']);
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
plot(xline_want(1)*ones(length(y_avg), 1), y_avg, 'r--', 'Linewidth', 2.5);
plot(xline_want(2)*ones(length(y_avg), 1), y_avg, 'b--', 'Linewidth', 2.5);  
hold on;

switch res
    case 2.5
        start = 14;
        step = 4;           % move by one scale height
        stop = 55;
    case 5
        start = 8;
        step = 2;
        stop = 27;

end
for j =(start:step:stop)
    plot(x, zp(:,j));

    %label the lines
    xpt = x_slc(end-3);
    ypt = double(zp(end, j));
    lbl = ['z = ', num2str(Z(j))];
    text(xpt, ypt, lbl, 'FontSize', 9, 'BackgroundColor', 'white', 'HorizontalAlignment', 'Center', 'Margin', 0.5)
    hold on;   
end

hold off;
cbar = colorbar();
caxis([zbot ztop]); 
ylim([ybot ytop]);
cbar.Label.String = 'Total Gas Mass Flux Divergence [1/s]';
xlabel([x_label, ' [Deg]']);
ylabel('Geopotential Altitude [km]');
title('Total Mass Flux Divergence (ilev)');
grid on;

if strcmp(x_label,'Latitude')
    sgtitle(['TIEGCM ', num2str(res) ,'deg Res Model Run at UT = 0, Lon = ', num2str(other_want), ', ', plotname]);
else
   sgtitle(['TIEGCM ', num2str(res) ,'deg Res Model Run at UT = 0, Lat = ', num2str(other_want), ', ', plotname]);
end

if saveFig ~= '0'
    saveas(thisfig, ['./Figures/Total_Gas/TIEGCM_LonMap_', saveFig, '.png'])  
end

end

