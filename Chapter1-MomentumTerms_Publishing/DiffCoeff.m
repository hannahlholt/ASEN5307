function [] = DiffCoeff(T, P, N2, O2, O1, He)
    % set all the molecular diffusion coefficients
    N2.setD( Di_TIEGCM('N2', T, P) );      % [m^2/s]
    O2.setD( Di_TIEGCM('O2', T, P) );
    O1.setD( Di_TIEGCM('O1', T, P) );
    He.setD( Di_TIEGCM('He', T, P) );
    
end

