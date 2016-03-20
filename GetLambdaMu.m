function [lambda, mu] = GetLambdaMu(SpeciesCount, I, E)

% Compute immigration rate and extinction rate for each species count.
% lambda(i) is the immigration rate for individual i.
% mu(i) is the extinction rate for individual i.

for i = 1 : 500
    lambda(1,i) = I *(1 - SpeciesCount(1,i) / 500);
    mu(1,i) = E *SpeciesCount(1,i) / 500;
end
return;
