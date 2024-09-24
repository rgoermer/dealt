#ifndef POISSON_BENCHMARK
#define POISSON_BENCHMARK


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
#include <deal.II/base/convergence_table.h>
 
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
#include <filesystem>

#include <utility>

using namespace dealt;
namespace Poisson_Benchmark {
  class Geometry_PF : public Function<2> {
  public: 
    Geometry_PF() : Function<2>(1) {}
    ~Geometry_PF() = default;

    virtual double value(
      const Point<2>&         p,
      const unsigned int      component = 0
    ) const override;

    virtual Tensor<1, 2> gradient(
      const Point<2>&         p,
      const unsigned int      component = 0
    ) const override;

    virtual SymmetricTensor<2, 2> hessian(
      const Point<2>&         p,
      const unsigned int      component = 0
    ) const override;
  };
  
  class Geometry_PB : public Function<2> {
  public: 
    Geometry_PB() : Function<2>(1) {}
    ~Geometry_PB() = default;

    virtual double value(
      const Point<2>&         p,
      const unsigned int      component = 0
    ) const override;

    virtual Tensor<1, 2> gradient(
      const Point<2>&         p,
      const unsigned int      component = 0
    ) const override;

    virtual SymmetricTensor<2, 2> hessian(
      const Point<2>&         p,
      const unsigned int      component = 0
    ) const override;
  };
          

  class Poisson_RHS : public Function<2> {
  public:
    Poisson_RHS() : Function<2>(1) {}
    ~Poisson_RHS() = default;

    virtual double value(
      const Point<2>&    p,
      const unsigned int component = 0
    ) const override;
  };

  class Poisson_SOL : public Function<2> {
  public:
    Poisson_SOL() : Function<2>(1) {}
    ~Poisson_SOL() = default;

    virtual double value(
      const Point<2>&     p,
      const unsigned int  component = 0
    ) const override;

    virtual Tensor<1, 2> gradient(
      const Point<2>&     p,
      const unsigned int  component = 0
    ) const override;

  };

  class Poisson_Benchmark
  {
  // Member Variables
  private:
    unsigned int ref;
    unsigned int order;
    unsigned int offset;

    IPF_Data<2, 2>            data;

    TS_Triangulation<2, 2>    tria;

    SparsityPattern           sparsity_pattern;
    SparseMatrix<double>      system_matrix;

    Vector<double>            solution;
    Vector<double>            system_rhs;

    // Vector<double>            l2_error;
    // Vector<double>            h1_error;
    // Vector<double>            mesh_size;
    // Vector<double>            dofs_per_level;

    OutputSetup               problem_out;
    double                    H1 = 1; 
    double                    L2 = 1;

    std::map<CellId, double>  cell_errors;

    IsoparametricFunction<2, 2> parametric_mapping;

    Poisson_RHS               rhs_fcn;
    Poisson_SOL               sol_fcn;

    bool uniform              = false;
    unsigned int              cycle = 0;
  // Public functions:
  public:
    Poisson_Benchmark(int ref, int order = 0);
    void run();

  // Private functions
  private:
    void get_IPF_data(std::vector< std::vector<double> >& kv,
                      std::vector< Point<2 + 1> >& cps,
                      std::vector<unsigned int>& deg
                      );

    IPF_Data<2, 2> get_IPF_data();

    std::filesystem::path directory_parent; // /path/to/out/poisson 
    std::filesystem::path directory_degree; // /path/to/out/poisson/op/
    std::filesystem::path directory_svg;    // /path/to/out/poisson/op/00svg
  
    void refine_grid();
    void setup_system();
    void assemble_system();
    void solve();
    void compute_h1_error();
    void compute_cell_errors();
    void output_system();

    void print_grid(
      const std::string& name
    ) const ;
  
    void estimate_and_mark();
  };

}


#endif
