function [] = VertAdvection(omega, Z, p0, g0, varargin)
% Sets the vertical advection terms for each object

    val = -p0 / g0 .* exp(-Z) .* omega;
  
    % number of species objects given
    for i = 1:length(varargin)
        Obj = varargin{1,i};
        Obj.setVertAdvec(val .* Obj.mmrGrad);
    end

end

