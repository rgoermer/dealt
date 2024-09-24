function z = full_lshape_geometry(x, y)
z = zeros(2, 1);
    if 0 <= x && x < 0.5
        z(1) = -1 - 4 * x^2 * (1 + y) + 4 * x * (1 + y);
        z(2) = y;
    else 
        z(1) = y;
        z(2) = -1 - 4 * x^2 * (1 + y) + 4 * x * (1 + y);
    end
end