function C = visualize_bezier_extraction(xi, knot_in)
%% Introduction
% Visualize bezier extraction for TSplines by using a local knot vector xi
% and a vector of intermediate knots to be inserted knot_in. This script
% will generate 3 subplots:
% 1. A plot of the B-Spline basis functions for the extended knot vector of
%    xi, where the B-Spline of interest [specified by xi] is the (nt+1)st 
%    spline in the plot
% 2. The different parts of xi over the knot intervals [together with the
%    knots from knot_in] as linear combination of Bernstein Polynomials
%    over the corresponding intervals / elements
% 3. For each knot interval the Bernstein Polynomials 
% Each plot is set to have the same (x,y)-ticks. The y-ticks are simply 
% [0, 1] and the x-ticks are constructed from the distinct knots in xi and
% knot_in, to stress the dependency of the decomposition on each knot
% interval.
%
% This Work is based on "Isogeometric finite element data structures based
% on Bezier extraction of T-Splines", by M. A. Scott, M. J. Borden, C. V.
% Verhoosel, T. W. Sederberg,  and T. J. R. Hughes, in International
% Journal for Numerical Methods in Engineering 2011, DOI: 10.1002/nme.3167
%
% The original version of Algorithm 1 in that paper had minor typos that 
% were corrected, a missing initialization of the extraction operator rows
% as well as a missing special case, that may occur, though this case may
% only be subject to MATLAB. 
%
% The results of Figure 18 are given at the end of this script. 
%
% Note, that we dismissed the construction of knot spans in this
% application. 
%
% Author: Robin Goermer {goermer at ifam.uni-hannover.de}
% Last Edited: 18th of January, 2023.
% Changelog: 18th of January -- Used Bernstein polynomials on reference
%                               element [-1, 1] to compute bezier
%                               extraction
%            10th of February -- Changed reference element to [0, 1] and
%                                put linear combination for derivatives up
%                                to the k-th derivative

%% User Parameters
% The user can either use their own knot vectors, or make use of the few
% examples we have listed here 
% xi = [0 1 2 3 4]; knot_in = [0.5, 3.5];
% xi = [0 0 1 2 2]; knot_in = [0.5, 1.75];
% xi = [0 1 2 2 2]; knot_in = [0.5, 1.75];
% xi = [0 1 1 1 1]; knot_in = [1/3 2/3];
% xi = [0 1 1 1 1]; knot_in = [];
% xi = [0.5 0.75 1 1  ]; knot_in =[];
% xi = [0 0 0.25 0.5  ]; knot_in = [0.375];
% xi = [0 0 0 0 1]; knot_in = [1/3 1/2 2/3];
% xi = [0.3125 0.375 0.4375 0.46875]; knot_in = [];
% xi = [0 1 1 1 2]; knot_in = [0.5, 3.5];
% xi = [0 0 1 1 2]; knot_in = [0.5, 3.5];
% xi = [0 1 2 3 3]; knot_in = [0.5, 3.5];
% xi = [0 1 1 2 2 3 3 4]; knot_in = [0.5, 3.5];
% xi = [0 0 1 1 1 1 2 2]; knot_in = [0.5, 3.5];

% How many intermediate evaluation points we use for any knot interval.
% Note: If this is constant, and the spacing is between intermediate
% evaluation points is the same, we can use the bernstein values on the
% reference to populate the values on the physical cell.
segments_length = 25;


%% The actual script to run through -- Changes are made at own risk

% Define necessary parameters
p = length(xi) - 2;
[U, nt] = compute_extended_kv(xi);

% Single out the distinct knots
single = xi; 
single(single(2:end) - single(1:end-1) == 0) = [];
single = sort([knot_in single]);

% Construct the knot vector for the bernstein polynomials over each
% interval
Ubar = augknt(single, p+1, p);

% Construct a list of evaluation points, where each interval consists of
% previously defined segment_length points.
tau = linspace(single(1), single(2), segments_length)';
for i = 2:length(single) - 1
    tau = [tau 
            linspace(single(i) + 1e-15, single(i+1), segments_length)' ]; 
end

% Calculate the basis functions in both cases
colmat = spcol(U, p+1, brk2knt(tau, p)); 
colmatbar = spcol(Ubar, p+1, tau);

subplot(3,1,1)
plot(tau, colmat(1:p:end, :));
str = '[';
for i = 1:length(xi)
    str = [str ' ' num2str(xi(i))];
end % for
str = [str ' ]'];
title(['Basis Functions of extended kv from \Xi = ' str ', with nt = ' num2str(nt)])
set(gca, 'xtick', single)
set(gca, 'ytick', [0 1])
set(gca, 'xlim',[min(xi), max(xi)]);

text(sum(xi) / length(xi), min(max(colmat(:, nt+1, 1)) + 0.05, 0.9), 'B_{\Xi}(\xi)')
grid on

subplot(3,1,2)
plot(tau, colmat(1:p:end, nt+1))
hold on

C = bezier_extraction(xi, knot_in);
[n, ~] = size(C); 

ref_knts = augknt([-1 1], p+1, p);
colmat_ref = spcol(ref_knts, p+1, ...
                brk2knt(linspace(-1, 1, segments_length), p));

% Get the linear combination of Bernstein polynomials in each interval
B = zeros(length(tau), p);

for j = 1 : n
    % Evaluation indices of Bersntein Polynomials
    indx = (1 : length(tau)/n) + (j-1)*length(tau)/n;
    
    % Indices of Bernstein Polynomials to be used
    indy = (1:(p+1)) + (j-1)*p;
    
    for d = 1 : p
        B(indx, d) = (single(j+1) - single(j))^(-d+1) * (2^(d-1)) * ...
                        sum(C(j, :) .* colmat_ref(d:p:end, :), 2); 
    end
    
    
    plot(tau(indx), B(indx, 1))
end 

% Calculate the errors in each derivative: 
err = zeros(1, p); 
for d = 1:p
    err(d) = norm(colmat(d:p:end, nt+1) - B(:, d)); 
end

hold off
title(['Visual Bezier-Decomposition of initial Spline - error given by err = ' num2str(err(1))])
set(gca, 'xtick', single)
set(gca, 'ytick', [0 1])
set(gca, 'ylim',[0, 1]);
set(gca, 'xlim',[min(xi), max(xi)]);
grid on

subplot(3,1,3)
plot(tau, colmatbar)
str = '[';
for i = 1:length(knot_in)
    str = [str ' ' num2str(knot_in(i))];
end % for
str = [str ' ]'];
title(['Bernstein Polynomials after inserting ' str])
set(gca, 'xtick', single)
set(gca, 'ytick', [0 1])
set(gca, 'xlim',[min(xi), max(xi)]);

grid on

%% Results of Paper, Figure 18
% 
% % The knot vectors of T-Spline 40 are given by xi_x, xi_y. The Bezier
% % decomposition over Element 13 is given. Element 13 is given by 
% %   e_{13} = [1.5, 2] x [1.5, 2]
% % Since the value 1.5 is not present in either xi_x, or xi_y, those
% % knots need to be inserted
% xi_x = [1 2 3 4 4]; in_x = [1.5];
% xi_y = [0 1 2 3 4]; in_y = [1.5];
% 
% % Next compute the extraction rows for the B-Spline given by the local knot
% % vectors xi_x, resp. xi_y
% c_x = bezier_extraction(xi_x, in_x);
% c_y = bezier_extraction(xi_y, in_y);
% 
% % Locally, we have now the extraction operator rows over the knot intervals 
% %   [1, 1.5], [1.5, 2], [2, 3], and [3, 4]
% % stored in c_x and the extraction operator rows over the knot intervals
% %   [0, 1], [1, 1.5], [1.5, 2], [2, 3], and [3, 4]
% % stored in c_y. Using the tensor product for the extraction operators, it
% % follows, that C = c_y(2, :) \tensor c_x(2, :) is the bezier extraction
% % operator over element e_{13}
% 
% % C = c_y(2, :)' * c_x(2, :); % Note, that this is C^T
% C = zeros(size(c_y, 2), size(c_x, 2));
% for y = 1:size(c_y, 2)
%     for x = 1:size(c_x, 2)
%         C(x, y) = c_y(2, y) * c_x(2, x);
%     end
% end
% C
% % 
% % This result fits the results stated in the paper, up to rounding errors

end % function 
