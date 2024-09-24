function [Ubar, nf, ne] = compute_extended_kv(xi)
% Compute extended kv of specified local knot vector, i.e. add knots to the
% front and end of xi until they have multiplictiy p+1
% INPUT:  xi   - A local knot vector belonging to a B-Spline of degree p, i.e. 
%                the knot vector xi has length p+2
% OUTPUT: Ubar - The extended kv of xi, where the first and last knot of
%                Ubar have multiplicity p+1
%         nf   - The number of knots added to the front to compute the
%                extended knot vector. This number [along with multiplicities
%                inside the knot vector xi] can be used to obtain the
%                single rows needed for the B-Spline specified by xi. 
% AUTHOR:   Robin Goermer {goermer at ifam.uni-hannover.de}

p = length(xi) - 2;
count = 1; k = 1;
while ( xi(k) == xi(k+1) )
    count = count + 1;
    k = k + 1;
end % while
nf = p - count + 1;

% Assume size(xi) = [1, m]

k = length(xi); count = 1;
while (xi(k) == xi(k-1))
    count = count + 1;
    k = k - 1;
end % while
ne = p - count + 1;

Ubar = [xi(1)*ones(1, nf) xi xi(end)*ones(1, ne)];

end