classdef TIEGCMspecies < handle
    % Class definition for each species of the TIEGCM output. All units are
    % SI, i.e. [kg, m, s]
    
    properties (Constant)
        ATOM_UNIT = 1.6605e-27;         % [kg/amu]
    end
    
    properties (SetAccess = immutable)
        % Can only set in Constructor
        slcPts          % slice pts can either be lat or lon
        altPts
        name
        weight
        mass
        mmr
        mmrGrad
        n
        rho        
    end
    
    properties (SetAccess = private)  
        Hp
        H_diff
        H_star 
        H_percent
        Di
        Vert_Adv
        Hor_MDiv
    end
    
    
    methods
        
        % ---- CONSTRUCTOR ------ 
        function obj = TIEGCMspecies(name, amu, MMR, Den, Z)   
            sizeofObj = size(MMR);
            obj.altPts = sizeofObj(end);
            obj.slcPts = sizeofObj(1:end-1);            % takes into account if MMR is 2D or 3D
            
            obj.name = name;                            % name of atom [string]
            obj.weight = amu/1000;                      % weight of atom [kg/mol]
            obj.mass = amu * obj.ATOM_UNIT;                         % atomic mass of species [kg]
            obj.mmr = CONVERT2ILEV(MMR, obj.slcPts, obj.altPts);    % mass mixing ratio on ilevs

            % Number density [1/m^3] on ilevs
            obj.n = obj.mmr .* Den;
            
            % Mass density [kg/m^3] on ilevs
            obj.rho = obj.mmr .* Den; 
            
            % MMR gradient w.r.t pressure coord Z = ln(p/p0) for every
            if length(obj.slcPts) == 1                
                for l = 1:obj.slcPts
                    obj.mmrGrad(l,:) = ThreePtGrad(Z(l,:), obj.mmr(l,:));
                end
            else % we have a 3D matrix
                for lon = 1:obj.slcPts(1)
                    for lat = 1:obj.slcPts(2)
                        obj.mmrGrad(lon,lat,:) = ThreePtGrad(Z(lon,lat,:), obj.mmr(lon,lat,:));
                    end
                end
                
            end
            
        end
        
        % ---- Change The Scale Heights ------  
        % NOTE - ONLY DESIGNED FOR 2D object 
        
        function obj = setScaleHeights(obj, Hp, H_diff, H_star)
            
            test = [size(Hp); size(H_diff); size(H_star)];
            matrx = [obj.slcPts, obj.altPts];   
            
            if (test(1,:) ~=  matrx) | (test(2,:) ~=  matrx) | (test(3,:) ~=  matrx)
                error("Scale height matrix has wrong size.");
            else                
                % Initialize the scale heights
                obj.Hp = Hp;
                obj.H_diff = H_diff;
                obj.H_star = H_star; 
                
                % Calculate Percents Difference from Diffusive Eq.
                obj.H_percent = ((obj.H_star ./ obj.H_diff) - 1) * 100;
            end  
           
            
        end
        
        % set the molecular Diffusion Coefficient 
        function obj = setD(obj, Di)
            if size(Di) ~= [obj.slcPts, obj.altPts]
                error("Diffusion Coefficient Error: Wrong Dimensions");
            else
               obj.Di = Di; 
            end
        end
        
        % set the vertical advection term 
        function obj = setVertAdvec(obj, Vert_Adv)
            if size(Vert_Adv) ~= [obj.slcPts, obj.altPts]
                error("Vertical Advection Error: Wrong Dimensions");
            else
               obj.Vert_Adv = Vert_Adv; 
            end
        end
        
        % set the horizontal divergence term
        function obj = setHorMassFluxDiv(obj, Hor_MDiv)
            if size(Hor_MDiv) ~= [obj.slcPts, obj.altPts]
                error("Horizontal Divergence Error: Wrong Dimensions");
            else
               obj.Hor_MDiv = Hor_MDiv; 
            end
        end
        
    end        
end

