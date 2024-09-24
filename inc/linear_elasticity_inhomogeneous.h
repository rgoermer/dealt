#ifndef LINEAR_ELASTICITY_INHOMOGENEOUS_H_
#define LINEAR_ELASTICITY_INHOMOGENEOUS_H_

#include <memory>

#include <utilities.h>
#include <ts_triangulation.h>
#include <tspline_function.h>
#include <isoparametric_function.h>
#include <linear_elasticity.h>
#include <get_cross_section.h>


#include <deal.II/base/convergence_table.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/function.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/multithread_info.h>

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

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_generator.h>

#include <fstream>
#include <iostream>

#include <utility>

namespace Linear_Elasticity {
  template<int dim> 
  class InhomogeneousProblem_SOL : public Elasticity_SOL<dim> {
  public: 
    InhomogeneousProblem_SOL (
    ) : Elasticity_SOL<dim>() {}

    virtual double value(
      const Point<dim>      &eval,
      const unsigned int     component = 0
    ) const override ;

    virtual void vector_value(
      const Point<dim>      &eval,
            Vector<double>  &out
    ) const override ;

    void tensor_value(
      const Point<dim>      &eval,
            Tensor<1, dim>  &out
    ) const;

    virtual void value_list(
      const std::vector< Point<dim> >   &evals,
            std::vector< double >       &out,
      const unsigned int                 component = 0
    ) const override ;

    virtual void vector_value_list(
      const std::vector< Point<dim> >     &evals,
            std::vector< Vector<double> > &out
    ) const override ;

    virtual Tensor<1, dim> gradient(
      const Point<dim>    &eval,
      const unsigned int   component
    ) const override ;

    virtual void vector_gradient(
      const Point<dim>                    &eval,
            std::vector< Tensor<1, dim> > &out
    ) const override;

    virtual void gradient_list(
      const std::vector< Point<dim> >         &evals,
            std::vector< Tensor<1, dim> >     &out,
      const unsigned int                       component = 0
    ) const override;

    virtual void vector_gradients(
      const std::vector< Point<dim> >                     &evals,
            std::vector< std::vector< Tensor<1, dim> > >  &out
    ) const override;

    double divergence(
      const Point<dim>  &eval
    ) const override ;

    void divergence_list(
      const std::vector< Point<dim> >   &evals,
            std::vector< double >       &out
    ) const override ;

    // Implement \nabla (\nabla \cdot u)
    double gradient_divergence(
      const Point<dim>    &eval,
      const unsigned int   component
    ) const override ;

    void gradient_divergence_list(
      const std::vector< Point<dim> >     &evals,
            std::vector< double >         &out,
      const unsigned int                   component
    ) const override ;

    void gradient_divergence_tensor(
      const Point<dim>        &eval,
            Tensor< 1, dim>   &out
    ) const override ;

    void gradient_divergence_tensor_list(
      const std::vector< Point<dim> >       &evals,
            std::vector< Tensor<1, dim> >   &out
    ) const override ;

    // Implement (\nabla \cdot \nabla) u
    double divergence_gradient(
      const Point<dim>      &eval,
      const unsigned int     component
    ) const override ;

    void divergence_gradient_list(
      const std::vector< Point<dim> >     &evals,
            std::vector< double >         &out,
      const unsigned int                   component
    ) const override ;

    void divergence_gradient_tensor(
      const Point<dim>        &eval,
            Tensor< 1, dim>   &out
    ) const override ;

    void divergence_gradient_tensor_list(
      const std::vector< Point<dim> >       &evals,
            std::vector< Tensor<1, dim> >   &out
    ) const override ;

    Tensor<1, dim> stress(
      const Point<dim>      &eval,
      const unsigned int    &component
    ) const override ;

    Tensor<2, dim> stress(
      const Point<dim>      &eval
    ) const override ;
  };

  template<int dim> 
  class InhomogeneousProblem {
  private: 
    IPF_Data<dim, dim>    data;

    unsigned int          ref;
    unsigned int          order;
    unsigned int          offset;

    TS_Triangulation<dim, dim> tria;

    SparsityPattern       sparsity_pattern; 
    SparseMatrix<double>  system_matrix;

    Vector<double>        solution;
    Vector<double>        system_rhs;

    double                H1 = 1.;
    double                L2 = 1.;
    OutputSetup           problem_out;

    Strategy              strategy;
    
    InhomogeneousProblem_SOL<dim>       sol_fcn;
    ElasticProblem_RHS<dim>             rhs_fcn; 
    Functions::ConstantFunction<dim>    lambda;
    Functions::ConstantFunction<dim>    mu;

    struct AssemblyScratchData {
      AssemblyScratchData(
              TS_Triangulation<dim, dim>    *tria,
        const std::vector<unsigned int>     &degrees
      );
      AssemblyScratchData(
        const AssemblyScratchData& scratch_data
      );

      TS_Triangulation<dim, dim>*    tria;
      TSValues<dim, dim, dim>       ts_values; 
    };

    struct AssemblyCopyData {
      FullMatrix<double>                    cell_matrix;
      Vector<double>                        cell_rhs;
      std::vector<types::global_dof_index>  local_dof_indices;
    };

  public: 
    InhomogeneousProblem (
      const unsigned int ref,
      const unsigned int degree,
      const Strategy     strategy = Strategy::Adaptive
    );

    void run();

  private: 
    void setup_system();
    void assemble_system();
    void local_assemble_system(
      const typename Triangulation<dim, dim>::active_cell_iterator    &cell,
      AssemblyScratchData                                             &scratch,
      AssemblyCopyData                                                &copy_data
    );
    void copy_local_to_global(
      const AssemblyCopyData    &copy_data
    ); 

    void solve();
    void estimate_mark_refine();
    void output_system();
    void output_prior_apply_boundary_values(
      const std::map<types::global_dof_index, double>& boundary_values
    );
    void compute_h1_error();

    void apply_boundary_values(
      const std::map<types::global_dof_index, double> &boundary_values,
      SparseMatrix<double>                            &matrix,
      Vector<double>                                  &solution,
      Vector<double>                                  &right_hand_side,
      const bool                                       eliminate_columns = true
    );

  }; 

} // namespace Linearelasticity

#endif
