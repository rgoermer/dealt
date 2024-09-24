u = @(x, y) 1./(x.^2 + y.^2).^(1/10);

x = linspace(-0.25, 1); 
y = linspace(0, 1);
[X, Y] = meshgrid(x, y);
Z = zeros(size(X));
for i = 1:100
    for j = 1:100
        Z(j, i) = u(x(i), y(j));
    end
end

surf(X, Y, Z);
colorbar