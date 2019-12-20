function [Tot_MDiv] = HorMassFluxDivergence(omegaExpGrad, Z, p0, g0, varargin)
% Sets the horizontal mass flux divergence terms for each object and
% returns the mass flux divergence for the total gas []

    val = p0 / g0 .* omegaExpGrad;
    
    % number of species objects given
    for i = 1:length(varargin)
        Obj = varargin{1,i};
        Obj.setHorMassFluxDiv(val .* Obj.mmr);
    end

    Tot_MDiv = -exp(Z) .* omegaExpGrad;
end
