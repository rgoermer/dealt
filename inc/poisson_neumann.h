/*
 * poisson_neumann.h
 *
 *  Created on: 29.03.2023
 *      Author: goermer
 */

#ifndef INC_POISSON_NEUMANN_H_
#define INC_POISSON_NEUMANN_H_


#include <memory>

#include <utilities.h>
#include <ts_triangulation.h>
#include <tspline_function.h>
#include <isoparametric_function.h>

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_bernstein.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/solver_relaxation.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/affine_constraints.h>

#include <fstream>
#include <iostream>

#include <utility>

#include <residual_error_estimators.h>

namespace Poisson_Neumann {
  using namespace dealt;

  class IPF_Mapping : public Function<2> {
  public:
    IPF_Mapping() : Function<2>(2) {};

    virtual double value(
      const Point<2>&    p,
      const unsigned int component
    ) const override;

    virtual Tensor<1, 2> gradient(
      const Point<2>&    p,
      const unsigned int component
    ) const override; 

    virtual SymmetricTensor<2, 2> hessian(
      const Point<2>&    p,
      const unsigned int component
    ) const override;
  };

  class IPF_Inverse : public Function<2> {
  public:
    IPF_Inverse() : Function<2>(2) {};

    virtual double value(
      const Point<2>&    p,
      const unsigned int component
    ) const override;

    virtual Tensor<1, 2> gradient(
      const Point<2>&    p,
      const unsigned int component
    ) const override; 

    virtual SymmetricTensor<2, 2> hessian(
      const Point<2>&    p,
      const unsigned int component
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

  class Poisson_NC : public Function<2> {
  public:
    Poisson_NC() :
        Function<2>(1) {}
    ~Poisson_NC() = default;

    virtual double value(
        const Point<2>&     p,
        const unsigned int component = 0
    ) const override;
  };

  enum RefinementStrategy{
    Adaptive, 
    Uniform
  };

  class Poisson_Benchmark
  {
    enum Boundary {
      None, 
      Dirichlet,
      Neumann
    };

  // Member Variables
  private:
    unsigned int ref;
    unsigned int order;

    IPF_Data<2, 2>            data;

    TS_Triangulation<2, 2>    tria;

    SparsityPattern           sparsity_pattern;
    // TrilinosWrappers::SparseMatrix system_matrix;
    SparseMatrix<double>      system_matrix;

    Vector<double>            solution;
    Vector<double>            system_rhs;

    OutputSetup problem_out;
    double      H1 = 1;
    double      L2 = 1;
    unsigned int              cycle;

    std::map<CellId, double>  cell_errors;

    Poisson_RHS               rhs_fcn;
    Poisson_SOL               sol_fcn;
    Poisson_NC                neumann_bc;

    RefinementStrategy        strategy; 
  // Public functions:
  public:
    Poisson_Benchmark(
      int ref, 
      int order = 0,
      const RefinementStrategy strategy = RefinementStrategy::Adaptive
    );
    void run();

  // Private functions
  private:
    IPF_Data<2, 2> get_IPF_data();

    void setup_system();
    void assemble_system();
    void solve();
    void compute_h1_error();
    void output_system();
    void print_error();
    void estimate_and_mark();
  };


}

#endif /* INC_POISSON_NEUMANN_H_ */
