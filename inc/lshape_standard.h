#ifndef LShapeStandard
#define LShapeStandard

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_postprocessor.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/base/smartpointer.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/numerics/solution_transfer.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/numerics/error_estimator.h>


#include <lshape.h>

#include <cmath>
#include <array>
#include <fstream>
#include <iostream>
#include <filesystem>


namespace LShape { 
  using namespace dealt;

  class LShapePostprocessor : public DataPostprocessor<2> {
  private: 
    LShape_SOL<2> sol_fcn;

  public: 
    LShapePostprocessor(Specifier type) 
      : DataPostprocessor<2>()
      , sol_fcn(type) {}

    virtual void  evaluate_scalar_field(
      const DataPostprocessorInputs::Scalar<2>  &inputs, 
            std::vector< Vector<double> >       &outputs
    ) const override;

    virtual std::vector< std::string > get_names(
    ) const override;

    UpdateFlags get_needed_update_flags(
    ) const override;
  };


  class LShape_Benchmark_Standard {
  public:
    Specifier type;
    Strategy  strategy;
  private: 
    Triangulation<2>          triangulation;
    FE_Q<2>                   fe;
    DoFHandler<2>             dof_handler;
    
    AffineConstraints<double> constraints;  
    
    SparsityPattern       sparsity_pattern;
    SparseMatrix<double>  system_matrix;
    
    Vector<double> solution;
    Vector<double> system_rhs;

    LShape_RHS<2>           rhs_fcn;
    LShape_SOL<2>           sol_fcn;
    LShape_NC1<2>           neumann_bc1;
    LShape_NC2<2>           neumann_bc2;
    LShape_NC3<2>           neumann_bc3;
    LShape_NC4<2>           neumann_bc4;

    OutputSetup   problem_out;

    unsigned int cycle = 0;
    double H1          = 1;
    double L2          = 1;

  public:
    LShape_Benchmark_Standard(
      const unsigned int order    = 1,
      const Specifier    type     = Specifier::classic,
      const Strategy     strategy = Strategy::adaptive
    );
    void run(
      const unsigned int ref = 10
    );
    
  private:
    void make_grid();
    void setup_system();
    void assemble_system();
    void compute_h1_error();
    void solve();
    void refine_grid();   
    void process_results();

    
  };


  
} // namespace LShape

#endif
