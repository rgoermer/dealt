#ifndef POISSON_BENCHMARK_H_ 
#define POISSON_BENCHMARK_H_ 


#include <memory>

#include <utilities.h>
#include <get_cross_section.h>
#include <ts_triangulation.h>
#include <tspline_function.h>
#include <isoparametric_function.h>
#include <poisson.h>

#include <deal.II/base/function.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/convergence_table.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
 
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
#include <deal.II/lac/affine_constraints.h>
 
#include <fstream>
#include <iostream>
#include <filesystem>

#include <utility>

using namespace dealt; 
namespace Poisson_Benchmark {
  class Poisson_Benchmark_3D {
  enum Boundary_IDs {
    Dirichlet, 
    Neumann
  };
  private: 
    // Member variables
    unsigned int            ref; 
    unsigned int            order; 

    IPF_Data<3, 3>          data; 
    TS_Triangulation<3, 3>  tria;

    SparsityPattern         sparsity_pattern;
    SparseMatrix<double>    system_matrix; 

    Vector<double>          solution; 
    Vector<double>          system_rhs;

    double                  residual = 1;
    double                  max_residual = 1;
    OutputSetup             problem_out; 
    

    Functions::ConstantFunction<3>     rhs_fcn;
    Functions::ConstantFunction<3>     nc_fcn; 
    Functions::ConstantFunction<3>     dirichlet_fcn;
    Functions::ConstantFunction<3>     eps_fcn;

    unsigned int            cycle = 0;
    unsigned int            offset;

    struct AssemblyScratchData {
      TS_Triangulation<3, 3>*    tria;
      TSValues<3>                ts_values; 
      TSFaceValues<3>            face_values;

      AssemblyScratchData(
              TS_Triangulation<3, 3>        *tria,
        const std::vector<unsigned int>     &degrees
      ) : tria(tria)
        , ts_values(this->tria, 
                    degrees, 
                    update_values |
                    update_gradients | 
                    update_JxW_values |
                    update_quadrature_points)
        , face_values(this->tria,
                      degrees, 
                      update_values | 
                      update_JxW_values |
                      update_quadrature_points) {}

      AssemblyScratchData(
        const AssemblyScratchData& scratch_data
      ) : tria(scratch_data.tria)
        , ts_values(scratch_data.ts_values) 
        , face_values(scratch_data.face_values){}

    };

    struct AssemblyCopyData {
      FullMatrix<double>                    cell_matrix;
      Vector<double>                        cell_rhs;
      std::vector<types::global_dof_index>  local_dof_indices;
    };

  public: 
    Poisson_Benchmark_3D(
      int ref = 10,
      int order = 0
    ); // done
    
    void run(); // done

  private: 
    void setup_system(); // done
    void assemble_system(); // done
    void local_assemble_system(
      const typename Triangulation<3, 3>::active_cell_iterator        &cell,
      AssemblyScratchData                                             &scratch,
      AssemblyCopyData                                                &copy_data
    ); // done 
    void copy_local_to_global(
      const AssemblyCopyData    &copy_data
    );  // done 

    void solve(); // done
    void estimate_mark_refine(); // done
    void output_system();

    IPF_Data<3, 3> get_data(
      const double r,
      const double R
    );
  };
} // Poisson_Benchmark





#endif
