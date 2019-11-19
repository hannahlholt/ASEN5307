function [feat_ilev] = CONVERT2ILEV(feat_lev, slicePts, altPts)
% converts TIEGCM matrix on the midpoint levels to matrix on interface
% levels

% feat_lev = matrix containing features on levs with size of [slicePts, altPts]
% slicePts can be either number of lat points or lon points. whatever you want
% altPts = number of altitude points
    
feat_ilev = zeros( size(feat_lev) );

% you are dealing with a 2D matrix
if length(slicePts) == 1
% find average between elements to get on the interfaces
    if slicePts ~= 0
        for alt = 1:altPts
            if alt ~= 1           
                for slc = 1:slicePts
                    if slc ~= 1           
                        feat_ilev(slc,alt) = 0.5*(feat_lev(slc,alt) + feat_lev(slc-1, alt-1));
                    else       
                        feat_ilev(slc,alt) = 0.5*(feat_lev(slc,alt) + feat_lev(slc,alt-1));        
                    end
                end

            else
                for slc = 1:slicePts
                    if slc ~= 1   
                        feat_ilev(slc,alt) = 0.5*(feat_lev(slc,alt) + feat_lev(slc-1, alt));
                    else       
                        feat_ilev(slc,alt) = feat_lev(slc,alt);          
                    end
                end
            end
        end

    else    % only converting an array
        for alt = 1:altPts
            if alt ~= 1           
                feat_ilev(alt) = 0.5*(feat_lev(alt) + feat_lev(alt-1));
            else
                feat_ilev(alt) = feat_lev(alt);
            end
        end
    end

% you want to put all lat/lon points onto the interfaces (3D matrix)   
% i.e. feat_lev = (lon, lat, alt)
else  
    for alt = 1:altPts
        if alt ~= 1 
            for slc2 = 1:slicePts(2)
                for slc1 = 1:slicePts(1)
                    if slc1 ~= 1  && slc2 ~= 1         
                        feat_ilev(slc1,slc2,alt) = 0.5*(feat_lev(slc1,slc2,alt) + feat_lev(slc1-1,slc2-1,alt-1));

                    elseif slc1 ~= 1 && slc2 == 1
                        feat_ilev(slc1,slc2,alt) = 0.5*(feat_lev(slc1,slc2,alt) + feat_lev(slc1-1,slc2,alt-1));

                    elseif slc1 == 1 && slc2 ~= 1
                        feat_ilev(slc1,slc2,alt) = 0.5*(feat_lev(slc1,slc2,alt) + feat_lev(slc1,slc2-1,alt-1));

                    else 
                        feat_ilev(slc1,slc2,alt) = 0.5*(feat_lev(slc1,slc2,alt) + feat_lev(slc1,slc2,alt-1));  
                    end
                end
            end

        else              
            for slc2 = 1:slicePts(2)
                for slc1 = 1:slicePts(1)
                    if slc1 ~= 1 && slc2 ~= 1
                        feat_ilev(slc1,slc2,alt) = 0.5*(feat_lev(slc1,slc2,alt) + feat_lev(slc1-1,slc2-1,alt));

                    elseif slc1 ~= 1 && slc2 == 1
                        feat_ilev(slc1,slc2,alt) = 0.5*(feat_lev(slc1,slc2,alt) + feat_lev(slc1-1,slc2,alt));

                    elseif slc1 == 1 && slc2 ~= 1
                        feat_ilev(slc1,slc2,alt) = 0.5*(feat_lev(slc1,slc2,alt) + feat_lev(slc1,slc2-1,alt));

                    else      
                        feat_ilev(slc1,slc2,alt) = feat_lev(slc1,slc2,alt);

                    end
                end
            end
        end
    end    
end
    
end



