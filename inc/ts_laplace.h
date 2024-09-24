#ifndef TS_LAPLACE
#define TS_LAPLACE


#include <memory>

#include <utilities.h>
#include <ts_triangulation.h>
#include <tspline_function.h>
#include <isoparametric_function.h>

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
 
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_bernstein.h>
 
#include <deal.II/dofs/dof_tools.h>
 
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
 
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
 
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/solver_relaxation.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/lac/affine_constraints.h>
 
#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>

#include <utility>

using namespace dealt;

template<int spacedim>
class TS_Laplace 
{
// Member Variables 
private: 
  int ref; 

  TS_Triangulation<2, spacedim>    tria; 
  DoFHandler<2>             dof_handler;
  FE_Bernstein<2>           fe;

  AffineConstraints<double> constraints;

  SparsityPattern           sparsity_pattern;
  SparseMatrix<double>      system_matrix;
  
  Vector<double>            solution;
  Vector<double>            system_rhs; 

  IsoparametricFunction<2, spacedim> parametric_mapping;  

  
// Public functions:         
public: 
  TS_Laplace(int ref, int order = 0); 
  void run(); 

// Private functions
private:
  void get_IPF_data(std::vector< std::vector<double> >& kv,
                    std::vector< Point<spacedim + 1> >& cps,
                    std::vector<unsigned int>& deg
                    );

  IPF_Data<2, spacedim> get_IPF_data();

  void refine_grid();
  void setup_system(); 
  void assemble_system();
  void solve();
  void output_results();
  void output_system();
};



#endif
