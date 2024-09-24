function [C, nt] = bezier_extraction_full(xi)
% Compute complete Bezier extraction operators for the B-Spline basis
% obtained by adding knots to the front and end of xi until the knot vector
% is p-open. 
% INPUT:  xi - A local knot vector belonging to a B-Spline of degree p, i.e. 
%              the knot vector xi has length p+2
% OUTPUT: C  - A (p+1) x (p+1) x nb Matrix, where each sub-matrix C(:, :, j)
%              corresponds to the extraction operator of the (p+1) splines
%              having support on element j. The row C(i, :, j) is thus the
%              linear combination of Bernstein Polynomials over element j
%              for B-Spline i in this element. 
%         nt - The number of knots added to the front to compute the
%              extended knot vector. This number [along with multiplicities
%              inside the knot vector xi] can be used to obtain the
%              single rows needed for the B-Spline specified by xi. 
% AUTHOR:   Robin Goermer {goermer@ifam.uni-hannover.de}


p = length(xi) - 2;

[U, nt] = compute_extended_kv(xi);
m = length(U);

C = eye(p+1, p+1);
% C = zeros(1, p+1);
% C(1, nt+1) = 1;

nb = 1; a = p+1; b = a+1;

while (b < m)
    C(:, :, end+1) = eye(p+1, p+1);
    
    i = b;
    while (b < m && U(b+1) == U(b))
        b = b+1;
    end % while ( b )
    mult = b - i + 1;
    
    if (mult < p)
        numer = U(b) - U(a);
        alphas = zeros(1, p-mult);
        for j = p : -1 : mult+1
            alphas(j-mult) = numer/ ( U(a+j) - U(a) );
        end % for ( j )
        r = p - mult;
        
        % Update coefficients
        for j = 1:r
            save = r- j + 1;
            s = mult + j;
            for k = p+1 : -1 : s+1
                alpha = alphas(k - s);
                C(:, k, nb) = alpha * C(:, k, nb) + (1 - alpha) * C(: , k-1, nb);
            end % for( k )
            
            if (b < m)
                C(save:j+save, save, nb+1) = C(p-j+1:p+1, p+1, nb);
            end % if ( b < m )
        end % for ( j )
        
        
        nb = nb + 1;
        if b < m
            a = b;
            b = b + 1;
        end % if ( b < m )
    end % if ( mult )
    

end % while ( b )

end % function;

