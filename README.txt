This implementation of AS T-splines depends on deal.II v9.3.3 or higher. 

If you use this code for your computations, please do cite our work: https://doi.org/10.1007/s00366-024-02002-1

After installation of deal.ii, configure CMakeLists.txt to find your version of deal.ii. 
The code provided in src/ and inc/ then can be compiled with your standard make commands 
for deal.ii applications.

0. cmake .
1. make

To run the examples:

2. ./tsplines -n <problem> -r <levels> -o <degree> 

where <problem> is [to this extent] either of the following benchmark problems
	benchmark 	-- Solves Laplace's equation on the unit square in 2D with a 
			   distorted inner domain using heuristic refinement approach
			   for a given solution
	poisson_nc 	-- Solves Poisson's equation with given solution on a curved domain,
			   boundary is also curved in this case. The suffixes 
				_uniform
				_adaptive
			   may be added to fix refinement strategy to the respective case.
			   Further, the suffix 
				_standard
			   may be added to solve the same problem using standard Q_p 
			   elements.
	lshape 		-- solves Poisson's equation on the Lshape domain [-1, 1]^2\[-1, 0]^2
			   for a given solution. Suffixes have to be added for a proper run:
				_ac
				_uc
				_ab
				_ub
			   where 'a' stands for adaptive refinement, 'u' stands for uniform 
			   refinement, 'c' refers to the classical L-shape problem, and 
			   'b' to a benchmark problem on the L-shape domain with homogeneous 
			   Dirichlet boundary conditions. Each suffix may be preceeded with 
				_standard
			   to obtain the corresponding solution using standard Q_p elements. 
  linear_elasticity_2d	-- solves a standard 2D linear elasticity equation on the unit square
			   without curved boundaries and no distortions, to demonstrate 
			   vector-valued Tsplines as finite elements. 
  linear_elasticity_3d	-- solves a standard 3D linear elasticity equation on the unit square
			   without curved boundaries and no distortions, to demonstrate 
			   vector-valued Tsplines as finite elements. This is not tested yet.

Other examples may be found in src/examples.cc that demonstrate further usage of TS_Triangulation. 

The variables <levels> and <degrees> are numeric values to describe the maximum level, resp. the
degree elevation imposed on T-splines or degree used on standard Q_p elements. 

With this, a valid command would be 

./tsplines -n lshape_ac -r 10 -o 1

To obtain the solution on the Lshape with a maximum of 10 refinements and a degree elevation of 1,
which gives you TSplines of degree 2 in this case, for the classical LShape problem. 

On the other hand 

./tsplines -n lshape_ac_standard -r 10 -o 1

returns the solution on the Lshape with a maximum of 10 refinements and a degree of 1, which gives
you the hat-functions on the reference element, for the classical Lshape problem. 

Since there is no DoFHandler involved with a TS_Triangulation, we have not yet found a way to obtain
a good way for the outputs. Thus, interesting data-fields will be written in the out/ directory. Each
command above should create their own respective path for the outputs. To load the .dat files in each 
output path, please checkout the matlab files. 

