#ifndef LINEAR_ELASTICITY_H_
#define LINEAR_ELASTICITY_H_

#include <memory>

#include <utilities.h>
#include <ts_triangulation.h>
#include <tspline_function.h>
#include <isoparametric_function.h>


#include <deal.II/base/convergence_table.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/function.h>

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

#include <fstream>
#include <iostream>

#include <utility>

namespace Linear_Elasticity {
  using namespace dealii;
  using namespace dealt;

  template<int dim>
  class LambdaFunction : public Function<dim> {
  public: 
    LambdaFunction() : Function<dim>(1) {}

    virtual double value(
      const Point<dim>      &eval, 
      const unsigned int    component = 0
    ) const;
  };

  template<int dim>
  class MuFunction : public Function<dim> {
  public: 
    MuFunction() : Function<dim>(1) {}

    virtual double value(
      const Point<dim>      &eval, 
      const unsigned int    component = 0
    ) const;
  };

  template<int dim>
  class Elasticity_SOL : public Function<dim> {
  public:
    Elasticity_SOL() : Function<dim>(dim) {}

    virtual double divergence(
      const Point<dim>  &eval
    ) const = 0;

    virtual void divergence_list(
      const std::vector< Point<dim> >   &evals,
            std::vector< double >       &out
    ) const = 0;

    // Implement \nabla (\nabla \cdot u)
    virtual double gradient_divergence(
      const Point<dim>    &eval,
      const unsigned int   component
    ) const = 0;

    virtual void gradient_divergence_list(
      const std::vector< Point<dim> >     &evals,
            std::vector< double >         &out,
      const unsigned int                   component
    ) const = 0;

    virtual void gradient_divergence_tensor(
      const Point<dim>        &eval,
            Tensor< 1, dim>   &out
    ) const = 0;

    virtual void gradient_divergence_tensor_list(
      const std::vector< Point<dim> >       &evals,
            std::vector< Tensor<1, dim> >   &out
    ) const = 0;

    // Implement (\nabla \cdot \nabla) u
    virtual double divergence_gradient(
      const Point<dim>      &eval,
      const unsigned int     component
    ) const = 0;

    virtual void divergence_gradient_list(
      const std::vector< Point<dim> >     &evals,
            std::vector< double >         &out,
      const unsigned int                   component
    ) const = 0;

    virtual void divergence_gradient_tensor(
      const Point<dim>        &eval,
            Tensor< 1, dim>   &out
    ) const = 0;

    virtual void divergence_gradient_tensor_list(
      const std::vector< Point<dim> >       &evals,
            std::vector< Tensor<1, dim> >   &out
    ) const = 0;

    virtual Tensor<1, dim> stress(
      const Point<dim>      &eval,
      const unsigned int    &component
    ) const = 0;

    virtual Tensor<2, dim> stress(
      const Point<dim>      &eval
    ) const = 0;

  }; // class Elasticity_SOL

  template<int dim> 
  class ElasticProblem_SOL : public Elasticity_SOL<dim> {
  public: 
    ElasticProblem_SOL(
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
    ) const override;

    void divergence_list(
      const std::vector< Point<dim> >   &evals,
            std::vector< double >       &out
    ) const override;

    // Implement \nabla (\nabla \cdot u)
    double gradient_divergence(
      const Point<dim>    &eval,
      const unsigned int   component
    ) const override;

    void gradient_divergence_list(
      const std::vector< Point<dim> >     &evals,
            std::vector< double >         &out,
      const unsigned int                   component
    ) const override;

    void gradient_divergence_tensor(
      const Point<dim>        &eval,
            Tensor< 1, dim>   &out
    ) const override;

    void gradient_divergence_tensor_list(
      const std::vector< Point<dim> >       &evals,
            std::vector< Tensor<1, dim> >   &out
    ) const override;

    // Implement (\nabla \cdot \nabla) u
    double divergence_gradient(
      const Point<dim>      &eval,
      const unsigned int     component
    ) const override;

    void divergence_gradient_list(
      const std::vector< Point<dim> >     &evals,
            std::vector< double >         &out,
      const unsigned int                   component
    ) const override;

    void divergence_gradient_tensor(
      const Point<dim>        &eval,
            Tensor< 1, dim>   &out
    ) const override;

    void divergence_gradient_tensor_list(
      const std::vector< Point<dim> >       &evals,
            std::vector< Tensor<1, dim> >   &out
    ) const override;

    Tensor<1, dim> stress(
      const Point<dim>      &eval,
      const unsigned int    &component
    ) const override;

    Tensor<2, dim> stress(
      const Point<dim>      &eval
    ) const override;
  };

  template<int dim> 
  class ElasticProblem_RHS : public Function<dim> {
    Elasticity_SOL<dim>*                 sol_fcn;
    Functions::ConstantFunction<dim>     lambda;
    Functions::ConstantFunction<dim>     mu;
  public: 
    ElasticProblem_RHS(
      Elasticity_SOL<dim>* sol_fcn
    ) : Function<dim>(dim) 
      , sol_fcn(sol_fcn)
      , lambda(1.)
      , mu(1.) {}

    virtual double value(
      const Point<dim>      &eval,
      const unsigned int     component = 0
    ) const override ;

    virtual void vector_value(
      const Point<dim>      &eval,
            Vector<double>  &out
    ) const override ;

    virtual void value_list(
      const std::vector< Point<dim> >   &evals,
            std::vector< double >       &out,
      const unsigned int                 component = 0
    ) const override ;

    virtual void vector_value_list(
      const std::vector< Point<dim> >     &evals,
            std::vector< Vector<double> > &out
    ) const override ;

    void tensor_value(
      const Point<dim>      &eval,
            Tensor<1, dim>  &out
    ) const;

    void tensor_value_list(
      const std::vector< Point<dim> >     &evals,
            std::vector< Tensor<1, dim> > &out
    ) const;
  };

  template<int dim>
  class ElasticProblem_NC_Y0 : public Function<dim> {
    IsoparametricFunction<dim, dim> geometry;
    ElasticProblem_SOL<dim>         solution;
  public: 
    ElasticProblem_NC_Y0(const IsoparametricFunction<dim, dim>& geometry) 
      : Function<dim>(dim)
      , geometry(geometry)
      , solution() {}

    virtual double value(
      const Point<dim>      &eval,
      const unsigned int     component = 0
    ) const override ;


    void tensor_value(
      const Point<dim>      &eval,
            Tensor<1, dim>  &out
    ) const;

  };

  template<int dim>
  class ElasticProblem_NC_X0 : public Function<dim> {
    IsoparametricFunction<dim, dim> geometry;
    ElasticProblem_SOL<dim>         solution;
  public: 
    ElasticProblem_NC_X0(const IsoparametricFunction<dim, dim>& geometry) 
      : Function<dim>(dim)
      , geometry(geometry)
      , solution() {}

    virtual double value(
      const Point<dim>      &eval,
      const unsigned int     component = 0
    ) const override ;

    void tensor_value(
      const Point<dim>      &eval,
            Tensor<1, dim>  &out
    ) const;

  };

  template<int dim>
  class ElasticProblem_NC_Z : public Function<dim> {
    IsoparametricFunction<dim, dim> geometry;
    ElasticProblem_SOL<dim>         solution;
  public: 
    ElasticProblem_NC_Z(
      const IsoparametricFunction<dim, dim>& geometry
    ) : Function<dim>(dim)
      , geometry(geometry)
      , solution() {}

    virtual double value(
      const Point<dim>      &eval,
      const unsigned int     component = 0
    ) const override ;

    void tensor_value(
      const Point<dim>      &eval,
            Tensor<1, dim>  &out
    ) const;

  };

  enum Boundaries {
    None,
    Dirichleth,
    Dirichleth0,
    Neumann_Y0,
    Neumann_X0,
    Neumann_Z,
  };
   
  enum Strategy {
    Adaptive, 
    Uniform
  };

  template<int dim> 
  class ElasticProblem {
  private: 
    IPF_Data<dim, dim>    data;

    unsigned int          ref;
    unsigned int          order;

    TS_Triangulation<dim, dim> tria;

    SparsityPattern       sparsity_pattern; 
    SparseMatrix<double>  system_matrix;

    Vector<double>        solution;
    Vector<double>        system_rhs;

    double                H1 = 1.;
    double                L2 = 1.;
    OutputSetup           problem_out;

    Strategy              strategy;
    
    ElasticProblem_SOL<dim>    sol_fcn;
    ElasticProblem_RHS<dim>    rhs_fcn; 
    ElasticProblem_NC_Y0<dim>  neumann_y0;
    ElasticProblem_NC_X0<dim>  neumann_x0;
    ElasticProblem_NC_Z<dim>   neumann_z;
    Functions::ConstantFunction<dim>      lambda;
    Functions::ConstantFunction<dim>      mu;
  public: 
    ElasticProblem(
      const unsigned int ref,
      const unsigned int degree,
      const Strategy     strategy = Strategy::Adaptive
    );

    void run();

  private: 
    void setup_system();
    void assemble_system();
    void solve();
    void estimate_mark_refine();
    void output_system();
    void compute_h1_error();
    void right_hand_side(
      const std::vector< Point<dim> >       &points,
            std::vector< Tensor<1, dim> >   &values
    ) const;

    void right_hand_side(
      const Point<dim>      &point,
            Tensor<1, dim>  &value
    ) const;

  }; 

  template<int dim>
  struct DataGenerator {
    IPF_Data<dim>   data;

    DataGenerator(const bool distortion = false);
  };

} // Linear_Elasticity namespace

#endif
