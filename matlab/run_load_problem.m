%% Poster output
o = 2; l = 9;
data = {'grid', 'physical_grid', 'physical_bezier_grid', 'solution', 'spline', 'knot_vector'};
data = {'solution'};
problem = '../out/poisson_neumann_adaptive/';
out_path = '~/Dokumente/Poster/bezier_grid';
load_problem(data, problem, out_path, o, l)



%% Standard FEM solution
o = 2;
data = {'solution'};
problem = '../out/poisson_neumann_standard_no_nc/';
for l = 1:4
    out_path = ['~/Dokumente/Poster/standard_l' num2str(l) '_'];
    load_standard_problem(data, problem, out_path, o, l, 1)
end

%% Solution
data = {'real_solution'};
out_path = '~/Dokumente/Poster/';
load_standard_problem(data, problem, out_path, o, l, 2)


%% Control Grid

data = {'real_solution', 'control_grid'};
out_path = './img_out/squished_rectangle_';
load_standard_problem(data, '', out_path, 0, 0)

%% Grid and sparsity pattern: benchmark

o = 4; l = 10;
data = {'physical_grid', 'sparsity_pattern'};
problem = '../out/poisson_neumann_adaptive/';
out_path = './img_out/squished_rectangle_';
load_problem(data, problem, [out_path 'iga_fem_o' num2str(o) '_l' num2str(l) '_'], o, l, 1)

l = 5;
data = {'physical_grid', 'sparsity_pattern'};
problem = '../out/poisson_neumann_standard_no_nc/';
load_standard_problem(data, problem, [out_path 'standard_fem_o' num2str(o) '_l' num2str(l) '_'], o, l, 3)

%% Grid and sparsity pattern: lshape

o = 3; l = 10; 
data = {'physical_grid'};
problem = '../out/lshape/lshape_ac/';
out_path = ['./img_out/lshape_ac_'];

load_problem(data, problem, out_path, o, l, 1)

%% Poisson benchmark 3D case
o = 2; l = 5; 
data = {'solution'};
problem = '../out/poisson_benchmark_3d/';
out_path = '';

for l = 3 : 3 : 12
    load_problem(data, problem, out_path, o, l, l+1, 1, 3);
end

%% Linear Elasticity

o = 2; l = 10; dim = 2;
out_path = './img_out/linear_elasticity_inhomogeneous';
data = {'physical_grid'};
problem = ['../out/linear_elasticity_inhomogeneous/' num2str(dim) 'd/'];
load_problem(data, problem, out_path, o, l, 2, dim)
%% 


data = {'displacement', 'control_grid', 'real_solution'};
load_standard_problem(data, '', out_path, 0, 0, 2, dim)

%% Inhomogeneous Linear elasticity
o = 2; l = 0; dim = 3; 
% out_path = ['./img_out/linear_elasticity_inhomogeneous/' num2str(dim) 'd/o' num2str(o) '_'];
out_path = '';
data = {'solution'};
problem = ['../out/linear_elasticity_inhomogeneous/' num2str(dim) 'd/'];
for l = 0:3:9
    load_problem(data, problem, out_path, o, l, l+1, dim);
end

%% LShape domain

o = 4; l = 8; dim = 2;
out_path = ['./img_out/lshape_standard_o' num2str(o) '_l' num2str(l) '_'];
data = {'physical_mesh'};
problem = '../out/lshape/lshape_ac/';
load_problem(data, problem, out_path, o, l, 1);


















