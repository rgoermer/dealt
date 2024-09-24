kvx = [0, 0, 0, 1, 1, 1];
kvy = kvx;
kvz = [0 0 1 1];
px = 2;
py = px;
pz = 1;

N=50;
eval_x = linspace(min(kvx), max(kvx), N);
eval_y = linspace(min(kvy), max(kvy), N);
eval_z = 1; %linspace(min(kvz), max(kvz), N);
n_evals = length(eval_x) * length(eval_y) * length(eval_z);

Bx = spcol(kvx, px + 1, eval_x);
By = spcol(kvy, py + 1, eval_y);
Bz = spcol(kvz, pz + 1, eval_z);
n_splines = size(Bx, 2) * size(By, 2) * size(Bz, 2);

r = 0.9; R = 1;
P = zeros(4, 18);
Pbar  = [ 1 1 0 1 1 0 0 0 0 ; ...
          0 0 0 1 1 0 1 1 0 ; ...
          0 1 1 0 1 1 0 1 1 ];
      
w = [1, 1/sqrt(2), 1, 1/sqrt(2), 0.5, 1/sqrt(2), 1, 1/sqrt(2), 1, ...
     1, 1/sqrt(2), 1, 1/sqrt(2), 0.5, 1/sqrt(2), 1, 1/sqrt(2), 1];
P(1:3, 1:9  ) = r * Pbar(:, :); % [1 4 7 2 5 8 3 6 9] 
P(1:3, 10:18) = R * Pbar(:, :);
P(1:3, :) = w .* P(1:3, :); 
P(4, :)   = w;

B = zeros(n_evals, n_splines );

spline = 1;
for z = 1:size(Bz, 2)
    for y = 1:size(By, 2)
        for x = 1:size(Bx, 2)
            eval = 1;
            for m = 1 : length(eval_z)
                for l = 1 : length(eval_y)
                    for k = 1 : length(eval_x)
                        B(eval, spline) = Bx(k, x) * By(l, y) * Bz(m, z);
                        eval = eval + 1;
                    end
                end
            end
            spline = spline + 1;
        end
    end
end


IPF = zeros(n_evals, 3);
W   = zeros(n_evals, 1);
% W = (w * B')';
for i = 1:n_splines
    W = W + P(4, i) * B(:, i); 
    IPF(:, 1) = IPF(:, 1) + P(1, i) * B(:, i);
    IPF(:, 2) = IPF(:, 2) + P(2, i) * B(:, i);
    IPF(:, 3) = IPF(:, 3) + P(3, i) * B(:, i);
end

if (min(abs(W)) < 1e-16)
    error('Division by zero');
end


IPF(:, 1) = IPF(:, 1) ./ W;
IPF(:, 2) = IPF(:, 2) ./ W;
IPF(:, 3) = IPF(:, 3) ./ W;

% IPF = IPF ./ W;

IPFx = reshape(IPF(:, 1), [length(eval_x), length(eval_y), length(eval_z)]);
IPFy = reshape(IPF(:, 2), [length(eval_x), length(eval_y), length(eval_z)]);
IPFz = reshape(IPF(:, 3), [length(eval_x), length(eval_y), length(eval_z)]);

plot3(IPF(:, 1), IPF(:, 2), IPF(:, 3), '.'); hold on
xlabel('x');
ylabel('y');
zlabel('z');
% for z = 1:length(eval_z)
%     surf(IPFx(:, :, z), IPFy(:, :, z), IPFz(:, :, z), 'EdgeColor', 'none'); hold on;
% end
Pbar = P';
plot3(Pbar(:, 1), Pbar(:, 2), Pbar(:, 3), 'x')
hold off;