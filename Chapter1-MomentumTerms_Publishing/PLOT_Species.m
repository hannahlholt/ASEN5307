function [] = PLOT_Species(res, x_label, zp, z, x_slc, xline_want, other_want, factor, omega, omegaGrad, plotname, saveFig, varargin)

% same as PLOT_TotalGas, except plots are done for specific species given
% number of species objects given


% THIS HAS BEEN UPDATED TO DEAL WITH BOTH latitude and longitude slices!!
str1 = '$$ \alpha \: \psi_{He} \nabla_h \bullet \vec{V}_h  [kg m/s]$$';
str2 = '$$ \alpha \omega \frac{\partial \psi_{He}}{\partial Z} [kg m/s]$$';
str3 = '$$ - \alpha \psi_{He} \omega $$';
str4 = '$$ \alpha \psi_{He} \frac{\partial \omega}{\partial Z} $$';


% Find average geopotential altitude for plotting
zp = zp./1000;               % plot in km
zp_avg = mean(zp, 1);        % row vector with average of each column.
x = x_slc;
x_ind = [find(x_slc == xline_want(1)), find(x_slc == xline_want(2)) ];    % find indices for lat or lon dotted line
% [X, Y] = meshgrid(x, zp_avg);     % mesh grid used for every subplot
[X, Y] = meshgrid(x, z);         % mesh grid for plotting on pressure lvls

zp_top = zp_avg(end-2);                % For N2, do end-5??
zp_bot = 100;
z_bot = -6;
z_top = 6;
logic = zp_avg >= zp_bot & zp_avg <= zp_top;  % dont plot everything!! Boundaries are bad.

X = X(logic,:);              % needed to do this b/c had trouble plotting;
Y = Y(logic,:);

num_cont = 300;
c = redblue(300);
num_lines = 15;

for i = 1:length(varargin)
    Obj = varargin{1,i};
       
    fig1 = figure('units','normalized','outerposition',[0 1 .7 1]);     % left start, bottom start, side length, tall length

% HORIZONTAL MASS FLUX DIVERGENCE
% ---------------------------------------------------------------------
    subplot(221)
    Q1 = Obj.Hor_MDiv(:,logic)';     % postive divergence = negative sign         
    Qtop = max(abs(Q1(:)));   
    Qbot = -Qtop;

    colormap(c);
    contourf(X, Y, Q1, num_cont, 'linecolor', 'none')
    hold on; 
    lines = linspace(min(Q1(:)), max(abs(Q1(:))), num_lines);
    contour(X, Y, Q1, lines, 'LineWidth', 0.1, 'linecolor', 'k')
    plot(xline_want(1)*ones(length(zp_avg), 1), z, 'r--', 'Linewidth', 2.5);
    plot(xline_want(2)*ones(length(zp_avg), 1), z, 'b--', 'Linewidth', 2.5);  
    cbar = colorbar();
    caxis([Qbot Qtop]); 
    ylim([z_bot z_top]);
    cbar.Label.Interpreter = 'latex';
    cbar.Label.String = str1; 
    xlabel([x_label, ' [Deg]']);
    ylabel('Pressure Level');
    
    % add another y axis for geopotential altitude
    yyaxis right
    ax = gca;
    Ticks = ax.YAxis(1).TickValues;
    [ind, ~] = find(z == Ticks);
    ax.YTick = linspace(0, 1, length(Ticks));
    ax.YAxis(2).TickLabels = round(zp_avg(ind));
    ylabel('~Z_p [km]')
    title([Obj.name, ' Mass Flux Divergence']);
    grid on;
    
% VERTICAL ADVECTION
% ---------------------------------------------------------------------
    subplot(222)
    Q2 = Obj.Vert_Adv(:,logic)';     % postive divergence = negative sign         
    Qtop = max(abs(Q2(:)));   
    Qbot = -Qtop;

    colormap(c);
    contourf(X, Y, Q2, num_cont, 'linecolor', 'none')
    hold on;   
    lines = linspace(min(Q2(:)), max(abs(Q2(:))), num_lines);
    contour(X, Y, Q2, lines, 'LineWidth', 0.1, 'linecolor', 'k')
    plot(xline_want(1)*ones(length(zp_avg), 1), z, 'r--', 'Linewidth', 2.5);
    plot(xline_want(2)*ones(length(zp_avg), 1), z, 'b--', 'Linewidth', 2.5);  
    cbar = colorbar();
    caxis([Qbot Qtop]); 
    ylim([z_bot z_top]);
    cbar.Label.Interpreter = 'latex';
    cbar.Label.String = str2;
    xlabel([x_label, ' [Deg]']);
    ylabel('Pressure Level');
    
    % add another y axis for geopotential altitude
    yyaxis right
    ax = gca;
    Ticks = ax.YAxis(1).TickValues;
    [ind, ~] = find(z == Ticks);
    ax.YTick = linspace(0, 1, length(Ticks));
    ax.YAxis(2).TickLabels = round(zp_avg(ind));
    ylabel('~Z_p [km]')
    
    title('Vertical Advection');
    grid on;
    
    
% -FACTOR * HE MMR * OMEGA 
% ---------------------------------------------------------------------
    subplot(223)
    Q3 = -factor .* Obj.mmr .* omega;   
    Q3 = Q3(:, logic)';
    Qtop = max(abs(Q3(:)));   
    Qbot = -Qtop;

    colormap(c);
    contourf(X, Y, Q3, num_cont, 'linecolor', 'none')
    hold on;   
    lines = linspace(min(Q3(:)), max(abs(Q3(:))), num_lines);
    contour(X, Y, Q3, lines, 'LineWidth', 0.1, 'linecolor', 'k')
    plot(xline_want(1)*ones(length(zp_avg), 1), z, 'r--', 'Linewidth', 2.5);
    plot(xline_want(2)*ones(length(zp_avg), 1), z, 'b--', 'Linewidth', 2.5);  
    cbar = colorbar();
    caxis([Qbot Qtop]); 
    ylim([z_bot z_top]);
    cbar.Label.Interpreter = 'latex';
    cbar.Label.String = str3;
    xlabel([x_label, ' [Deg]']);
    ylabel('Pressure Level');
    
    % add another y axis for geopotential altitude
    yyaxis right
    ax = gca;
    Ticks = ax.YAxis(1).TickValues;
    [ind, ~] = find(z == Ticks);
    ax.YTick = linspace(0, 1, length(Ticks));
    ax.YAxis(2).TickLabels = round(zp_avg(ind));
    ylabel('~Z_p [km]');
    
    title('Vertical Omega Term');
    grid on;
    
% FACTOR * HE MMR * GRAD(OMEGA)
% ---------------------------------------------------------------------
    subplot(224)
    Q4 = factor.* Obj.mmr .* omegaGrad; 
    Q4 = Q4(:, logic)';
    Qtop = max(abs(Q4(:)));   
    Qbot = -Qtop;

    colormap(c);
    contourf(X, Y, Q4, num_cont, 'linecolor', 'none')
    hold on;
    lines = linspace(min(Q4(:)), max(abs(Q4(:))), num_lines);
    contour(X, Y, Q4, lines, 'LineWidth', 0.1, 'linecolor', 'k')
    plot(xline_want(1)*ones(length(zp_avg), 1), z, 'r--', 'Linewidth', 2.5);
    plot(xline_want(2)*ones(length(zp_avg), 1), z, 'b--', 'Linewidth', 2.5);  
    cbar = colorbar();
    caxis([Qbot Qtop]); 
    ylim([z_bot z_top]);
    cbar.Label.Interpreter = 'latex';
    cbar.Label.String = str4;
    xlabel([x_label, ' [Deg]']);
    ylabel('Pressure Level');
    
    % add another y axis for geopotential altitude
    yyaxis right
    ax = gca;
    Ticks = ax.YAxis(1).TickValues;
    [ind, ~] = find(z == Ticks);
    ax.YTick = linspace(0, 1, length(Ticks));
    ax.YAxis(2).TickLabels = round(zp_avg(ind));
    ylabel('~Z_p [km]')
    title('Vertical Omega Gradient Term');
    grid on;
    
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
  
    if strcmp(x_label,'Latitude')
        sgtitle(['TIEGCM ', num2str(res) ,'deg Res Model Run at UT = 0, Lon = ', num2str(other_want), ', ', plotname]);
    else
       sgtitle(['TIEGCM ', num2str(res) ,'deg Res Model Run at UT = 0, Lat = ', num2str(other_want), ', ', plotname]);
    end
        
    
   %% PLOT HELIUM DENSITY SCALE HEIGHT at every Z Vs Longitude
    fig2 = figure('units','normalized','outerposition',[0 1 .5 .7]);
    
    % ----------------
    Q = Obj.H_percent(:,logic)';
    Qtop = 100; 
    Qbot = -100;
  
    colormap(c);
    contourf(X, Y, Q, num_cont, 'linecolor', 'none')
    hold on; 
    lines = linspace(min(Q(:)), max(abs(Q(:))), 30);
    contour(X, Y, Q, lines, 'LineWidth', 0.1, 'linecolor', 'k')
    plot(xline_want(1)*ones(length(zp_avg), 1), z, 'r--', 'Linewidth', 2.5);
    plot(xline_want(2)*ones(length(zp_avg), 1), z, 'b--', 'Linewidth', 2.5);  
    cbar = colorbar();
    caxis([Qbot Qtop]); 
    ylim([z_bot z_top]);
    cbar.Label.String = 'Scale Height Percent Difference';
    xlabel([x_label, ' [Deg]']);
    ylabel('Pressure Level');
    
    % add another y axis for geopotential altitude
    yyaxis right
    ax = gca;
    Ticks = ax.YAxis(1).TickValues;
    [ind, ~] = find(z == Ticks);
    ax.YTick = linspace(0, 1, length(Ticks));
    ax.YAxis(2).TickLabels = round(zp_avg(ind));
    ylabel('~Z_p [km]')
    
    title([Obj.name, ' Scale Height Difference from Diffusive Eq., UT = 0, Lat = ', num2str(other_want)]);
    grid on;
    
   
    %% PLOT VERTICAL STRUCTURE OF DIVERGENCE AND ADVECTION
    fig3 = figure('units','normalized','outerposition',[0 1 .7 .6]);     % left start, bottom start, side length, tall length
    for j = 1:length(x_ind)
        subplot(1,2,j)
        
        hold on;
        xtop = 1E-13;  
        xbot = -xtop;
        xdis = (xtop - xbot)/2;
        r1 = rectangle('Position',[xbot zp_bot xdis 400]);
        r1.FaceColor = [.85 .85 1];     % very light blue
        r2 = rectangle('Position',[0 zp_bot xdis 400]);
        r2.FaceColor = [1 .9 .9];       % very light red
        plot(NaN, 'color', [.85 .85 1], 'Linewidth', 10);
        plot(NaN, 'color', [1 .9 .9], 'Linewidth', 10);
        
        
        p1 = plot(Q1(:,x_ind(j)), zp(x_ind(j),logic)', 'Linewidth', 2.5);
        p2 = plot(Q2(:,x_ind(j)), zp(x_ind(j),logic)', 'Linewidth', 2.5);
        legend([p1 p2], 'Hor. Divergence', 'Vert. Advection', 'Location', 'best');
        ylim([zp_bot 400]);
        xlim([xbot xtop])
        ylabel('Geopotential Altitude [km]')
        hold off;
        
        grid on;    
        set(gca, 'Layer', 'top');
        
        if j == 1
            title(['Helium Maximum at Lon=', num2str(xline_want(j)), ', Lat=', num2str(other_want)])
        else
            title(['Helium Minimum at Lon=', num2str(xline_want(j)), ', Lat=', num2str(other_want)])
        end
    end
    

    %% PLOT HELIUM DIFFUSIVE AND ACTUAL DENSITY SCALE HEIGHTS
    fig4 = figure('units','normalized','outerposition',[0 1 .7 .6]);
    for j = 1:length(x_ind)
        subplot(1,2,j)
        hold on;
        plot(Obj.H_diff(x_ind(j), :), zp(x_ind(j), :), '--k' , 'Linewidth', 1.5);
        plot(Obj.H_star(x_ind(j), :), zp(x_ind(j), :), 'k' , 'Linewidth', 1.5);
        xlabel('Scale Height [km]');
        ylabel('Geopoential altitude [km]');
        ylim([zp(x_ind(j), 5) 550]);
        legend('H_{\rho}','H^*_{\rho}', 'location', 'best');
        grid on;       
        title(['Helium Scale Heights at Lon=', num2str(xline_want(j)), ', Lat=', num2str(other_want)]);
    
    end
    
    
     if saveFig ~= '0'
        saveas(fig1, ['./Figures/Helium_Paper/', Obj.name, 'LonMap_MomentumTerms_', saveFig, '.svg'])  
        saveas(fig2, ['./Figures/Helium_Paper/', Obj.name, 'LonMap_ScaleHtPercDiff_', saveFig, '.svg']) 
        saveas(fig3, ['./Figures/Helium_Paper/', Obj.name, 'LonMap_1DAltitude_', saveFig, '.svg'])
        saveas(fig4, ['./Figures/Helium_Paper/', Obj.name, '_Density_ScaleHeights_', saveFig, '.svg'])
    end
    
    
    
end % end 1:length(varargin) (number of species objects)
end



% ----- ADDING LINE LABELS -------
% %     switch res
% %         case 2.5
% %             start = 14;
% %             step = 4;           % move by one scale height
% %             stop = 55;
% %         case 5
% %             start = 8;
% %             step = 2;
% %             stop = 27;
% %             
% %     end
% %     for j =(start:step:stop)
% %         plot(x, zp(:,j));
% % 
% %         %label the lines
% %         xpt = x_slc(end-3);
% %         ypt = double(zp(end, j));
% %         lbl = ['z = ', num2str(Z(j))];
% %         text(xpt, ypt, lbl, 'FontSize', 9, 'BackgroundColor', 'white', 'HorizontalAlignment', 'Center', 'Margin', 0.5)
% %         hold on;   
% %     end
%     hold off;

