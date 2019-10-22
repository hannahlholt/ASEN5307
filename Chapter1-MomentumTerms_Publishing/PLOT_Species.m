function [] = PLOT_Species(res, x_label, zp, Z, x_slc, xline_want, other_want, omega, plotname, saveFig, varargin)

% same as PLOT_TotalGas, except plots are done for specific species given
% number of species objects given


% THIS HAS BEEN UPDATE TO DEAL WITH BOTH latitude and longitude slices!!


zp = zp./1000;              % plot in km

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


for i = 1:length(varargin)
    Obj = varargin{1,i};
       
    thisfig = figure('units','normalized','outerposition',[0 1 .7 1]);     % left start, bottom start, side length, tall length

    % ----------------
    subplot(221)
    Q = omega(:,logic)';      % omega for every longitude and alt for a specific latitute [m/s]
    ztop = max(abs(Q(:)));  
    zbot = -ztop;

    colormap(c); 
    contourf(X, Y, Q, num_cont, 'linecolor', 'none')
    hold on;   
    plot(xline_want*ones(length(y_avg), 1), y_avg, 'k--', 'Linewidth', 2.5)   
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
    Q = Obj.H_percent(:,logic)';
    ztop = max(abs(Q(:)));  
    zbot = -ztop;

    colormap(c);
    contourf(X, Y, Q, num_cont, 'linecolor', 'none')
    hold on;   
    plot(xline_want*ones(length(y_avg), 1), y_avg, 'k--', 'Linewidth', 2.5) 
    cbar = colorbar();
    caxis([zbot ztop]); 
    ylim([ybot ytop]);
    cbar.Label.String = 'Scale Height Perc. Difference';
    xlabel([x_label, ' [Deg]']);
    ylabel('Geopotential Altitude [km]');
    title([Obj.name, ' Scale Height Perc. Difference from Diffusive Eq. (ilev)']);
    grid on;

    % ----------------
    subplot(223)
    Q = Obj.Vert_Adv(:,logic)';
    ztop = max(abs(Q(:)));  
    zbot = -ztop;

    colormap(c);
    contourf(X, Y, Q, num_cont, 'linecolor', 'none')
    hold on;   
    plot(xline_want*ones(length(y_avg), 1), y_avg, 'k--', 'Linewidth', 2.5) 
    cbar = colorbar();
    caxis([zbot ztop]); 
    ylim([ybot ytop]);
    cbar.Label.String = 'Vertical Advection [kg/(m s) * m^2]';

    xlabel([x_label, ' [Deg]']);
    ylabel('Geopotential Altitude [km]');
    title([Obj.name, ' Vertical Advection (+ Up) (ilev)']);
    grid on;

    % % ----------------
    subplot(224)
    Q = Obj.Hor_MDiv(:,logic)';     % postive divergence = negative sign         
    ztop = max(abs(Q(:)));   
    zbot = -ztop;

    colormap(c);
    contourf(X, Y, Q, num_cont, 'linecolor', 'none')
    hold on;   
    plot(xline_want*ones(length(y_avg), 1), y_avg, 'k--', 'Linewidth', 2.5) 
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
    cbar.Label.String = 'Mass Flux Divergence [kg/(m s) * m^2]';
    cbar.Label.String = 'Total Divergence [1/s]';
    xlabel([x_label, ' [Deg]']);
    ylabel('Geopotential Altitude [km]');
    title([Obj.name, ' Mass Flux Divergence (ilev)']);
    grid on;

    if strcmp(x_label,'Latitude')
        sgtitle(['TIEGCM ', num2str(res) ,'deg Res Model Run at UT = 0, Lon = ', num2str(other_want), ', ', plotname]);
    else
       sgtitle(['TIEGCM ', num2str(res) ,'deg Res Model Run at UT = 0, Lat = ', num2str(other_want), ', ', plotname]);
    end
        
    if saveFig ~= '0'
        saveas(thisfig, ['./Figures/', Obj.name, '/TIEGCM_LonMap_', saveFig, '.png'])  
    end
    
    
   
    
    
    
end % end 1:length(varargin) (number of species objects)
end

