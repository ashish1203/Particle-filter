function [SpeciesCount] = GetSpeciesCounts(q, inf)
% This loop assumes the population is already sorted from most fit to least fit.
for i = 1 : 500
    if q(1,i) < inf
        SpeciesCount(1,i) = 500 - i;
    else
        SpeciesCount(1,i) = 0;
    end
end
return;
