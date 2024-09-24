#ifndef EXAMPLES_H_
#define EXAMPLES_H_

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/base/tensor.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_out.h>

#include <iostream>
#include <vector>

#include <ts_triangulation.h>
#include <ts_laplace.h>
#include <poisson.h>
#include <poisson_benchmark.h>
#include <poisson_neumann.h>
#include <poisson_neumann_standard.h>
#include <poisson_3d.h>
#include <lshape.h>
#include <lshape_standard.h>
#include <linear_elasticity.h>
#include <linear_elasticity_inhomogeneous.h>
#include <mwe.h>


/*
 * Provides some example for BSpline Basis 
 */
using namespace dealii;

Vector<double> other_kv(double a, double b, unsigned int N);

std::vector< Point<3, double> > cp_disk();
std::vector< std::vector<double> > kv_disk();

void print_example(std::string ex = "disc", int ref = 10, int order = 0);





#endif
