#ifndef POISSON_3D
#define POISSON_3D

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

#include <deal.II/lac/affine_constraints.h>
 
#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>

#include <utility>

namespace Poisson_Singularity_3D {
  using namespace dealt;

  enum Strategy {adaptive, uniform};
  enum Boundary {none, 
                  dirichleth_x0, 
                  dirichleth_x1, 
                  dirichleth_y0,
                  dirichleth_y1, 
                  neumann_z0, 
                  neumann_z1};

  enum Problem {simple, singularity};

  struct SingularityDataGenerator {
    SingularityDataGenerator ();

    IPF_Data<3, 3> data;
  };

  class IPF_mapping : public Function<3> {
  public: 
    IPF_mapping() : Function<3>(3) {}

    virtual double value(
      const Point<3>&       p,
      const unsigned int    component = 0
    ) const override; 

    virtual void vector_value(
      const Point<3>&       p,
      Vector< double >&     values
    ) const override; 

    Point<3> point_value(
      const Point<3>&   p
    ) const;

    Tensor<1, 3> tensor_value(
      const Tensor<1, 3>& p
    ) const; 

    virtual Tensor<1, 3> gradient(
      const Point<3>&       p,
      const unsigned int    component = 0
    ) const override;

    virtual void vector_gradient(
      const Point<3>&       p,
      std::vector< Tensor<1, 3> >& gradients
    ) const override;

    Tensor<2, 3> tensor_gradient(
      const Tensor<1, 3>&       p
    ) const;

    Tensor<1, 3> normal(
      const Point<3>&           p,
      const unsigned int        f
    ) const;

    virtual SymmetricTensor<2, 3> hessian(
      const Point<3>&       p,
      const unsigned int    component = 0
    ) const override;
  };

  class IPF_mapping_simplified : public Function<3> {
  public: 
    IPF_mapping_simplified() : Function<3>(3) {}

    virtual double value(
      const Point<3>&       p,
      const unsigned int    component = 0
    ) const override; 

    virtual void vector_value(
      const Point<3>&       p,
      Vector< double >&     values
    ) const override; 

    Point<3> point_value(
      const Point<3>&   p
    ) const;

    Tensor<1, 3> tensor_value(
      const Tensor<1, 3>& p
    ) const; 

    virtual Tensor<1, 3> gradient(
      const Point<3>&       p,
      const unsigned int    component = 0
    ) const override;

    virtual void vector_gradient(
      const Point<3>&       p,
      std::vector< Tensor<1, 3> >& gradients
    ) const override;

    Tensor<2, 3> tensor_gradient(
      const Tensor<1, 3>&       p
    ) const;

    Tensor<1, 3> normal(
      const Point<3>&           p,
      const unsigned int        f
    ) const;

    virtual SymmetricTensor<2, 3> hessian(
      const Point<3>&       p,
      const unsigned int    component = 0
    ) const override;
  };

  class IPF_inverse : public AutoDerivativeFunction<3> {
  public: 
    IPF_inverse() : AutoDerivativeFunction<3>(1e-12, 3) {}

    virtual void vector_value(
      const Point<3>&       p,
      Vector< double >&     values
    ) const override; 


    virtual void vector_hessian(
      const Point<3>&       p, 
      std::vector< SymmetricTensor< 2,3 > >& hessians
    ) const override;
  };

  class Singularity_RHS : public Function<3> {
    static Problem prob;
  public:
    Singularity_RHS() : Function<3>(1) {} 

    virtual double value(
      const Point<3>&       p,
      const unsigned int    component = 0
    ) const override; 

  }; // Singularity_RHS

  class Singularity_SOL : public Function<3> {
    static Problem prob;
  public: 
    Singularity_SOL() : Function<3>(1) {} 

    virtual double value(
      const Point<3>&       p,
      const unsigned int    component = 0
    ) const override; 

    virtual Tensor<1, 3> gradient(
      const Point<3>&       p,
      const unsigned int    component = 0
    ) const override;

  }; // Singularity_SOL

  class Singularity_NumSol : public Function<3> {
  private:
    Vector<double>                    solution; 
    std::vector< TSplineFunction<3, 3> > base; 

  public: 
    Singularity_NumSol(
      const Vector<double>&                        solution,
      const std::vector< TSplineFunction<3, 3> >&  base
    ) : Function<3>(1), solution(solution), base(base) {}
    
    virtual double value(
      const Point<3>&       p,
      const unsigned int    component = 0
    ) const override; 

    virtual void value_list(
      const std::vector< Point< 3 > >&  points,
      std::vector<double>&              values,
      const unsigned int                component = 0
    ) const override;

    virtual Tensor<1, 3> gradient(
      const Point<3>&       p,
      const unsigned int    component = 0
    ) const override;
  }; // Singularity_NumSol


  class Singularity_NC : public Function<3> {
  public:
    Singularity_NC() : Function<3>(1) {} 

    virtual double value(
      const Point<3>&       p,
      const unsigned int    component = 0
    ) const override; 

  }; // Singularity_NC

  class Singularity_NC_Z0 : public Function<3> {
    IPF_mapping phi; 
  public:
    Singularity_NC_Z0() : Function<3>(1), phi() {} 

    virtual double value(
      const Point<3>&       p,
      const unsigned int    component = 0
    ) const override; 

    Tensor<1, 3> normal(
      const Point<3>& p
    ) const;
  }; // Singularity_NC

  class Singularity_NC_Z1 : public Function<3> {
    IPF_mapping phi;
  public:
    Singularity_NC_Z1() : Function<3>(1), phi() {} 

    virtual double value(
      const Point<3>&       p,
      const unsigned int    component = 0
    ) const override; 

  }; // Singularity_NC

  class Poisson_Singularity {
  // typenames 
  private: 
    using active_cell_iterator  = TriaActiveIterator<CellAccessor<3, 3>>;
    using cell_iterator         = TriaIterator<CellAccessor<3, 3>>;

  // Member declarations
  private: 
    // Define polynomial degrees and maximum refinement
    unsigned int ref;
    unsigned int order; 
    unsigned int offset = 0;

    // Store the grid used
    IPF_Data<3, 3>            data;
    TS_Triangulation<3, 3>    tria;

    // Define the numeric system to be solved and 
    // the solution vector
    SparsityPattern           sparsity_pattern;
    SparseMatrix<double>      system_matrix;
    Vector<double>            system_rhs;
    Vector<double>            solution; 


    // Define vectors for important outputs
    Vector<double>            l2_error;
    Vector<double>            h1_error;
    Vector<double>            mesh_size;
    Vector<double>            dofs_per_level;

    // save error for each cell
    std::map< active_cell_iterator,
              double >        cell_errors;

    // Data used for the system, rhs, bc, and 
    // solution function for h1 and l2 errors
    Singularity_RHS           rhs_fcn; 
    Singularity_SOL           sol_fcn;
    Singularity_NC            neumann_bc;
    Singularity_NC_Z0         neumann_bc_z0;
    Singularity_NC_Z1         neumann_bc_z1;

    // Output strings repeatedly use the same 
    // strings, so we save it 
    std::string               out_init;

    struct AssemblyScratchData {
      TSValues<3>                 ts_values;
      TSFaceValues<3>             face_values;

      std::vector<double>         rhs_values;
      std::vector<double>         face_boundary_values;

      Singularity_RHS             rhs_fcn;

      AssemblyScratchData();
      AssemblyScratchData(
        const AssemblyScratchData& scratch_data
      );
    };

    struct AssemblyCopyData{
      FullMatrix<double>                   cell_matrix;
      Vector<double>                       cell_rhs;
      std::vector<types::global_dof_index> local_dof_indices;
    };


  // Private functions. These functions are 
  // self explanatory
  private:
    void setup_system();
    void assemble_system();
    void local_assemble_system(
      const TriaIterator<CellAccessor<3, 3>>& cell,
      AssemblyScratchData&                      scratch,
      AssemblyCopyData&                         copy_data
    );
    void copy_local_to_global(
      const AssemblyCopyData& copy_data
    );
    void solve();
    void compute_h1_error();
    void output_system();

    void print_grid(
        const bool bezier = false
    ) const ;

    void estimate_and_refine();
    void set_boundary_ids();

    void print_solution();
    void print_numeric_solution();
    void print_numeric_solution(const double x);

    void print_pointwise_error() const; 

    void print_surface_z0() const; 
    void print_reference_values(const double z) const;

  // Public functions
  public: 
    Poisson_Singularity (
      unsigned int ref,
      unsigned int order = 0
    );

    void run();
  }; 
} // namespace Poisson_Neumann

#endif
