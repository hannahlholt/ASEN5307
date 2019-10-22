function [H_P, H_T, H_T_He, H_rho_diff, H_rho_star] = ScaleHeight(N2, O2, O1, He, tn, mbar, zp, den, g0, Av, kb)
    
    %-----Finding Pressure Scale Heights for gas constituents from kT/mg-----
    H_P = kb * tn ./ (mbar/Av .* g0);                % mean pressure scale height [m]
    Hp_n2 = kb * tn ./ (N2.mass * g0);               % N2 press. scale height [m]
    Hp_o2 = kb * tn ./ (O2.mass * g0);               % O2 press. scale height [m]
    Hp_o1 = kb * tn ./ (O1.mass * g0);               % O1 press. scale height [m]
    Hp_he = kb * tn ./ (He.mass * g0);               % He press. scale height [m]

    % Calculate Scale Heights Using 3pt Differentiation Technique    
    points = zeros(size(tn));
    H_he = points;
    H_n2 = points;
    H_o1 = points;
    H_o2 = points;
    H_T = points;
    H_mass = points;
    H_rho_star = points; 

    % ----- 2 Dimensional 3 Point Gradient -------
    for l = 1:N2.slcPts
        for z = 1:N2.altPts
            if z == 1   % First Point gradient technique
                coeff1 = (2*zp(l,1)-zp(l,2)-zp(l,3))/((zp(l,1)-...
                    zp(l,2))*(zp(l,1)-zp(l,3)));

                coeff2 = (2*zp(l,1)-zp(l,1)-zp(l,3))/((zp(l,2)-...
                    zp(l,1))*(zp(l,2)-zp(l,3)));

                coeff3 = (2*zp(l,1)-zp(l,1)-zp(l,2))/((zp(l,3)-...
                    zp(l,1))*(zp(l,3)-zp(l,2)));

                H_he(l,1) = -1/He.rho(l,1)*(He.rho(l,1)*coeff1+He.rho(l,2)*coeff2+...
                    He.rho(l,3)*coeff3);
                H_rho_star(l,1) = -1/den(l,1)*(den(l,1)*coeff1+den(l,2)*coeff2+...
                    den(l,3)*coeff3);
                H_n2(l,1) = -1/N2.rho(l,1)*(N2.rho(l,1)*coeff1+N2.rho(l,2)*coeff2+...
                    N2.rho(l,3)*coeff3);
                H_o1(l,1) = -1/O1.rho(l,1)*(O1.rho(l,1)*coeff1+O1.rho(l,2)*coeff2+...
                    O1.rho(l,3)*coeff3);
                H_o2(l,1) = -1/O2.rho(l,1)*(O2.rho(l,1)*coeff1+O2.rho(l,2)*coeff2+...
                    O2.rho(l,3)*coeff3);
                H_T(l,1) = 1/tn(l,1)*(tn(l,1)*coeff1+tn(l,2)*coeff2+...
                    tn(l,3)*coeff3);
                H_mass(l,1) = -1/mbar(l,1)*(mbar(l,1)*coeff1+mbar(l,2)*coeff2+...
                    mbar(l,3)*coeff3);


            elseif z == N2.altPts %Last point gradient technique
                coeff1 = (2*zp(l,z)-zp(l,z-1)-zp(l,z))/((zp(l,z-2)-...
                    zp(l,z-1))*(zp(l,z-2)-zp(l,z)));

                coeff2 = (2*zp(l,z)-zp(l,z-2)-zp(l,z))/((zp(l,z-1)-...
                    zp(l,z-2))*(zp(l,z-1)-zp(l,z)));

                coeff3 = (2*zp(l,z)-zp(l,z-2)-zp(l,z-1))/((zp(l,z)-...
                    zp(l,z-2))*(zp(l,z)-zp(l,z-1)));

                H_he(l,z) = -1/He.rho(l,z)*(He.rho(l,z-2)*coeff1+He.rho(l,z-1)*coeff2+...
                    He.rho(l,z)*coeff3);
                H_rho_star(l,z) = -1/den(l,z)*(den(l,z-2)*coeff1+den(l,z-1)*coeff2+...
                    den(l,z)*coeff3);
                H_n2(l,z) = -1/N2.rho(l,z)*(N2.rho(l,z-2)*coeff1+N2.rho(l,z-1)*coeff2+...
                    N2.rho(l,z)*coeff3);            
                H_o1(l,z) = -1/O1.rho(l,z)*(O1.rho(l,z-2)*coeff1+O1.rho(l,z-1)*coeff2+...
                    O1.rho(l,z)*coeff3); 
                H_o2(l,z) = -1/O2.rho(l,z)*(O2.rho(l,z-2)*coeff1+O2.rho(l,z-1)*coeff2+...
                    O2.rho(l,z)*coeff3);
                H_T(l,z) = 1/tn(l,z)*(tn(l,z-2)*coeff1+tn(l,z-1)*coeff2+...
                    tn(l,z)*coeff3);
                H_mass(l,z) = -1/mbar(l,z)*(mbar(l,z-2)*coeff1+mbar(l,z-1)*coeff2+...
                    mbar(l,z)*coeff3);


            else % Middle Points gradient technique
                coeff1 = (2*zp(l,z)-zp(l,z)-zp(l,z+1))/((zp(l,z-1)-...
                    zp(l,z))*(zp(l,z-1)-zp(l,z+1)));

                coeff2 = (2*zp(l,z)-zp(l,z-1)-zp(l,z+1))/((zp(l,z)-...
                    zp(l,z-1))*(zp(l,z)-zp(l,z+1)));

                coeff3 = (2*zp(l,z)-zp(l,z-1)-zp(l,z))/((zp(l,z+1)-...
                    zp(l,z-1))*(zp(l,z+1)-zp(l,z)));

                H_he(l,z) = -1/He.rho(l,z)*(He.rho(l,z-1)*coeff1+He.rho(l,z)*coeff2+...
                    He.rho(l,z+1)*coeff3);
                H_rho_star(l,z) = -1/den(l,z)*(den(l,z-1)*coeff1+den(l,z)*coeff2+...
                    den(l,z+1)*coeff3);
                H_n2(l,z) = -1/N2.rho(l,z)*(N2.rho(l,z-1)*coeff1+N2.rho(l,z)*coeff2+...
                    N2.rho(l,z+1)*coeff3);
                H_o1(l,z) = -1/O1.rho(l,z)*(O1.rho(l,z-1)*coeff1+O1.rho(l,z)*coeff2+...
                    O1.rho(l,z+1)*coeff3);
                H_o2(l,z) = -1/O2.rho(l,z)*(O2.rho(l,z-1)*coeff1+O2.rho(l,z)*coeff2+...
                    O2.rho(l,z+1)*coeff3);
                H_T(l,z) = 1/tn(l,z)*(tn(l,z-1)*coeff1+tn(l,z)*coeff2+...
                    tn(l,z+1)*coeff3);                               
                H_mass(l,z) = -1/mbar(l,z)*(mbar(l,z-1)*coeff1+mbar(l,z)*coeff2+...
                    mbar(l,z+1)*coeff3);

            end
        end
    end
    
    %-----Get Scale Height from Inverse-----
    H_N2_star = 1./H_n2;
    H_O2_star = 1./H_o2;
    H_O1_star = 1./H_o1;
    H_He_star = 1./H_he;
    
    H_T = 1./H_T;
    H_T_He = H_T/.62;                     % <---- Alpha for helium is -.38
    H_m = 1./H_mass;
    H_rho_diff = (1./H_T + 1./H_P + 1./H_m).^-1;   % Atmospheric diffusive density scale height
    H_rho_star = 1./H_rho_star;


    %----Put Together Diffusive Profiles and Mean Mass Profile-----
    H_N2_diff = (1./H_T + 1./Hp_n2).^-1;        % N2 Diffusive profile
    H_O2_diff = (1./H_T + 1./Hp_o2).^-1;        % O2 diffusie profile
    H_O1_diff = (1./H_T + 1./Hp_o1).^-1;        % O1 Diffusive profile
    H_He_diff = (1./H_T_He + 1./Hp_he).^-1;     % Helium diffusive profile


    % Set all the object scale heights and return general values to main
    N2.setScaleHeights(Hp_n2, H_N2_diff, H_N2_star);
    O2.setScaleHeights(Hp_o2, H_O2_diff, H_O2_star);
    O1.setScaleHeights(Hp_o1, H_O1_diff, H_O1_star);
    He.setScaleHeights(Hp_he, H_He_diff, H_He_star);
     
end

