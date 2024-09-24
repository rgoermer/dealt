function [C, nf] = bezier_extraction( xi, knot_in)
% Compute Bezier extraction operator rows specified by the knot vector xi,
% where the knots in knot_in are inserted into the extended knot vector of
% xi.
% INPUT:  xi - A local knot vector belonging to a B-Spline of degree p, i.e. 
%              the knot vector xi has length p+2
%    knot_in - A non-decreasing sequence of knots to be inserted into the
%              extended knot vector of xi
% OUTPUT: C  - A (nb+m) x (p+1) Matrix, where each row corresponds to the
%              extraction operator for the B-Spline specified by the local
%              knot vector xi on element j, where m = length(knot_in).
%         nf - The number of knots added to the front to compute the
%              extended knot vector. 
% AUTHOR:   Robin Goermer {goermer at ifam.uni-hannover.de}
% Last Edited: 17th of March, 2023.
% Changelog: 17th of March -- Update for nb in special case fixed


% Ensure the knot vector is non-decreasing:
xi = sort(xi); 

p = length(xi) - 2;
m = length(knot_in);

[Ubar, nf, ne] = compute_extended_kv(xi);

a = p+1;
b = a+1;

C = zeros(1, p+1);
C(1, nf+1) = 1;
mbar = ne + 2 + nf + m;

si = 1;
nb = 1;
while (b < mbar)
    % Initialize extraction operator
    C(end+1, :) = zeros(1, p+1);
    
    add = 0;
    
    % Check if additional knot needs to be inserted:
    if (si <= m && Ubar(b) >= knot_in(si))
        % Yup
        mult = 0;
        add = 1;
        
        Ubar = [Ubar(1:b-1) knot_in(si) Ubar(b:end)];
        si = si+1;
    else
        % Nope
        i = b;
        while (b < mbar && Ubar(b+1) == Ubar(b))
            b = b + 1;
        end
        mult = b - i + 1;
    end % if ( si )
    
    % Where to put the 1 depends on the multiplicity of the current knot
    C(end, nf + 2 - nb + (si-1) - mult - add) = 1; 

    if (mult < p)
        numer = Ubar(b) - Ubar(a);
        alphas = zeros(1, p - mult);
        for j = p:-1:mult+1
            alphas(j - mult) = numer / (Ubar(a + j + add) - Ubar(a));
        end % for ( j )
        
        r = p - mult;
        
        %Update Matrix coefficients:
        for j = 1 : r
            save = r - j + 1;
            s = mult + j;
            for k = p+1 : -1 :s+1
                alpha = alphas(k - s);
                C(nb, k) = alpha * C(nb, k) + (1 - alpha)*C(nb, k-1);
            end % for ( k )
            
            if b < mbar
                % Update overlapping coefficients
                C(nb + 1, save) = C(nb, p + 1);
            end % if ( b )
        end % for ( j )
        
        % Advance to the next element:
        nb = nb + 1;
        
        % Update coefficients for next iteration
        if (b < mbar)
            a = b;
            b = b + 1;
        end % if ( b )
        
    else
        
        % Special case: mult = p => There are only two elements given by the
        % right-most, resp. left-most, bezier splines
        C(nb + 1, 1) = 1; 
        a = b; 
        b = b+1;
        
        % Advance to next element
        nb = nb + 1;
        
    end % if ( mult )
end % while
% C(end, :) = [];

end % function bezier_extraction
