function outline_linear_elasticity

kv_x = [0 0 0 1 1 1];
kv_y = [0 0 0 1 1 1]; 

Nx = 100;
Ny = 100;
xi_x = linspace(0, 1, Nx); 
xi_y = linspace(0, 1, Ny); 

p_x = 2;
p_y = 2; 

B_x = spcol(kv_x, p_x + 1, xi_x);
B_y = spcol(kv_y, p_y + 1, xi_y);

P = [+0.0 +0.0; ...
    +0.5 +0.0; ...
    +1.0 +0.0; ...
    +0.0 +0.5; ...
    +0.5 +0.5; ...
    +1.0 +0.5; ...
    +0.0 +1.0; ...
    +0.5 +1.0; ...
    +1.0 +1.0  ...
    ];

% P = [-1.0 -1.0; ...
%      +0.0 -1.0; ...
%      +1.0 -1.0; ...
%      -1.0 +0.0; ...
%      +0.0 +0.0; ...
%      +1.0 +0.0; ...
%      -1.0 +1.0; ...
%      +0.0 +1.0; ...
%      +1.0 +1.0  ...
%     ];

 
splines_x = size(B_x, 2);
splines_y = size(B_y, 2);

B = zeros(Nx, Ny, splines_x * splines_y); 

ind = 1;
for i = 1:splines_y
    for j = 1:splines_x
        B(:, :, ind) = B_y(:, i) * B_x(:, j)'; 
        ind = ind + 1;
    end
end

% Check if splines are correct ...
% [X, Y] = meshgrid(xi_x, xi_y);
% for i = 1:splines_x*splines_y
%     subplot(splines_x, splines_y, i)
%     contourf(X, Y, B(:, :, i), 'EdgeColor', 'none');
%     title(['i = ' num2str(i)])
% end

Phi_x = zeros(Nx, Ny);
Phi_y = zeros(Nx, Ny);
for i = 1:splines_x * splines_y
    Phi_x = Phi_x + P(i, 1) * B(:, :, i);
    Phi_y = Phi_y + P(i, 2) * B(:, :, i);
end 

u =@(x, y) [sin(pi * x) * sin(pi * y); -cos(pi * x / 2) * cos(pi * y / 2)];

X = Phi_x;
Y = Phi_y;
Z = zeros(size(Phi_x));
U = zeros(size(Phi_x));
V = zeros(size(Phi_y));

for i = 1:Nx
    for j = 1:Ny
        uh = 1/(pi * pi) * u(Phi_x(i, j), Phi_y(i, j));
        U(i, j) = uh(1);
        V(i, j) = uh(2);
        Z(i, j) = norm(uh);
    end
end
cmap = load('smooth-cool-warm.dat') / 255;

figure(2)
colormap(cmap(:, 2:4));
subplot(2, 1, 1)
contourf(X, Y, U, 'EdgeColor', 'none');
colorbar
subplot(2, 1, 2)
colormap(cmap(:, 2:4));
contourf(X, Y, V, 'EdgeColor', 'none');
colorbar





end