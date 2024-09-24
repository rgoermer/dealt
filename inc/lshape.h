#ifndef LSHAPE_H_
#define LSHAPE_H_

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
#include <deal.II/base/auto_derivative_function.h>
#include <deal.II/base/function_derivative.h>

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
#include <deal.II/base/convergence_table.h>

#include <deal.II/lac/affine_constraints.h>
 
#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>

#include <utility>

namespace LShape {
  using namespace dealt;

  enum Specifier {classic, benchmark};
  enum Strategy {adaptive, uniform};
  enum Boundary {dirichleth, // 0 
                  neumann_1, // 1
                  neumann_2, // 2
                  neumann_3, // 3 
                  neumann_4, // 4
                  neumann_diag, // 5
                  none};     // 6

  template<int spacedim = 2>
  class LShape_RHS_benchmark : public Function<spacedim> {
  public:
    LShape_RHS_benchmark() : Function<spacedim>(1){}

    virtual double value(
      const Point<spacedim>& p,
      const unsigned int     component = 0
    ) const override;
  };

  template<int spacedim = 2>
  class LShape_RHS_classic : public Function<spacedim> {
  public:
    LShape_RHS_classic() : Function<spacedim>(1){}

    virtual double value(
      const Point<spacedim>& p,
      const unsigned int     component = 0
    ) const override;
  };

  template<int spacedim = 2>
  class LShape_SOL_benchmark : public AutoDerivativeFunction<spacedim> {
  public:
    LShape_SOL_benchmark() :
      AutoDerivativeFunction<spacedim>(1e-16){  }

    virtual double value(
      const Point<spacedim>&       p,
      const unsigned int    component = 0
    ) const override;

    virtual Tensor<1, spacedim> gradient(
      const Point<spacedim>&    p,
      const unsigned int component = 0
    ) const override;
  };

  template<int spacedim>
  class LShape_SOL_classic : public AutoDerivativeFunction<spacedim> {
  public:
    LShape_SOL_classic() :
      AutoDerivativeFunction<spacedim>(1e-16){  }

    virtual double value(
      const Point<spacedim>&       p,
      const unsigned int    component = 0
    ) const override;

    virtual Tensor<1, spacedim> gradient(
      const Point<spacedim>&    p,
      const unsigned int component = 0
    ) const override;
  };

  template<int dim = 2, int spacedim = dim>
  struct DataGenerator {
    DataGenerator(
        const unsigned int n_elements = 2,
        const bool no_c0_edges = true
    );

    IPF_Data<dim, spacedim> data;
  };

  template<int spacedim = 2>
  class LShape_SOL : public AutoDerivativeFunction<spacedim> {
    Specifier type;
  public:
    LShape_SOL(Specifier type) :
        AutoDerivativeFunction<spacedim>(1e-16),
        type(type){}

    virtual double value(
      const Point<spacedim>& p,
      const unsigned int     component = 0
    ) const override ;

    virtual Tensor<1, spacedim> gradient(
      const Point<spacedim>& p,
      const unsigned int     component = 0
    ) const override ;
  };

  template<int spacedim = 2>
  class LShape_NC1 : public Function<spacedim> {
    Specifier type;
  public:
    LShape_NC1(Specifier other_type) :
        Function<spacedim>(1),
        type(other_type) {};

    virtual double value(
      const Point<spacedim>& p,
      const unsigned int     component = 0
    ) const override;
  };

  template<int spacedim = 2>
  class LShape_NC2 : public Function<spacedim> {
    Specifier type;
  public:
    LShape_NC2(Specifier other_type) :
        Function<spacedim>(1),
        type(other_type) {};

    virtual double value(
      const Point<spacedim>& p,
      const unsigned int     component = 0
    ) const override;
  };

  template<int spacedim = 2>
  class LShape_NC3 : public Function<spacedim> {
    Specifier type;
  public:
    LShape_NC3(Specifier other_type) :
        Function<spacedim>(1),
        type(other_type) {};

    virtual double value(
      const Point<spacedim>& p,
      const unsigned int     component = 0
    ) const override;
  };

  template<int spacedim = 2>
  class LShape_NC4 : public Function<spacedim> {
    Specifier type;
  public:
    LShape_NC4(Specifier other_type) :
        Function<spacedim>(1),
        type(other_type) {};

    virtual double value(
      const Point<spacedim>& p,
      const unsigned int     component = 0
    ) const override;
  };

  template<int spacedim = 2>
  class LShape_NC_Diag : public Function<spacedim> {
    Specifier type;
  public:
    LShape_NC_Diag(Specifier other_type) :
        Function<spacedim>(1),
        type(other_type) {};

    virtual double value(
      const Point<spacedim>& p,
      const unsigned int     component = 0
    ) const override;
  };



  template<int spacedim = 2>
  class LShape_RHS : public Function<spacedim> {
    Specifier type;
  public:
    LShape_RHS(Specifier type) :
        Function<spacedim>(),
        type(type) {}

    virtual double value(
      const Point<spacedim>&  p,
      const unsigned int      component = 0
    ) const override;
  };


  template<int dim>
  class LShape_Problem {
  public:
    Specifier type;
    Strategy  strategy;

  private:
    unsigned int              ref;
    unsigned int              order;
    unsigned int              offset;

    IPF_Data<dim, dim>         data;
    TS_Triangulation<dim, dim> tria;

    SparsityPattern           sparsity_pattern;
    SparseMatrix<double>      system_matrix;

    Vector<double>            solution;
    Vector<double>            system_rhs;


    OutputSetup               problem_out;

    double H1 = 1.;
    double L2 = 1.;

    std::map<CellId, double>  cell_errors;

    LShape_RHS<dim>           rhs_fcn;
    LShape_SOL<dim>           sol_fcn;
    LShape_NC1<dim>           neumann_bc1;
    LShape_NC2<dim>           neumann_bc2;
    LShape_NC_Diag<dim>       neumann_bc_diag;


    unsigned int cycle;

  // Public functions to be called
  public:
    LShape_Problem(
        const unsigned int ref,
        const unsigned int order = 0,
        const Specifier    type = Specifier::classic,
        const Strategy     strategy = Strategy::adaptive
        );

    LShape_Problem(
        const LShape_Problem& other
        ) = delete;

    LShape_Problem() = delete;

    void run();

  private:
    void setup_system();
    void assemble_system();
    void solve();
    void compute_h1_error();
    void output_system();

    void print_error();
    void print_grid(const std::string& name);

    void estimate_and_refine();
    void set_boundary_ids();
    void print_solution();
    void print_numeric_solution();
    void print_numeric_solution(const double x);
  };
  
}


#endif
