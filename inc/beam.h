#ifndef BEAM_H_ 
#define BEAM_H_ 

#include <memory>
#include <chrono>

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

#include <deal.II/fe/fe_bernstein.h>
#include <deal.II/fe/fe_system.h>

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

namespace Classical {
  using namespace dealii;
  using namespace dealt;

  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::duration;
  using std::chrono::milliseconds;

  class ClassicalProblemRHS : public Function<3> {
  public:
    ClassicalProblemRHS() : Function<3>(3) {}
    
    virtual double value(
      const Point<3>     &eval, 
      const unsigned int  component
    ) const override;
  }; // class ClassicalProblemRHS

  class ClassicalProblemNeumann : public Function<3> {
  public:
    ClassicalProblemNeumann() : Function<3>(3) {}
    
    virtual double value(
      const Point<3>     &eval, 
      const unsigned int  component
    ) const override;
  }; // class ClassicalProblemRHS

  enum Boundaries {
    None, 
    Dirichlet,
    Neumann,
    HomogeneousNeumann
  }; 

  enum Type {
    Slab, 
    Beam
  };

  struct DataGenerator {
    IPF_Data<3> data;
    DataGenerator(const Type problem);

    private: 
      void generate_beam_data();
      void generate_slab_data();
  };


  class ClassicalProblem {
  private: 
    IPF_Data<3, 3>        data;

    unsigned int          ref;
    unsigned int          order;
    unsigned int          offset;

    TS_Triangulation<3>   tria;

    SparsityPattern       sparsity_pattern; 
    SparseMatrix<double>  system_matrix;

    Vector<double>        solution;
    Vector<double>        system_rhs;

    double                residual = 1;
    double                max_residual = 1;
    OutputSetup           problem_out;

    dealii::ConvergenceTable  run_times; 

    ClassicalProblemRHS                    rhs_fcn; 
    ClassicalProblemNeumann                nc_fcn;

    Functions::ConstantFunction<3>    lambda;
    Functions::ConstantFunction<3>    mu;

    struct AssemblyScratchData {
      AssemblyScratchData(
              TS_Triangulation<3, 3>        *tria,
        const std::vector<unsigned int>     &degrees
      );
      AssemblyScratchData(
        const AssemblyScratchData& scratch_data
      );

      TS_Triangulation<3, 3>*    tria;
      TSValues<3, 3, 3>          ts_values; 
      TSFaceValues<3, 3, 3>      face_values;
    };

    struct AssemblyCopyData {
      FullMatrix<double>                    cell_matrix;
      Vector<double>                        cell_rhs;
      std::vector<types::global_dof_index>  local_dof_indices;
    };

  public: 
    ClassicalProblem(
      const unsigned int ref,
      const unsigned int degree,
      const Type         problem
    );

    void run();

  private: 
    void setup_system();
    void assemble_system();
    void local_assemble_system(
      const typename Triangulation<3, 3>::active_cell_iterator        &cell,
      AssemblyScratchData                                             &scratch,
      AssemblyCopyData                                                &copy_data
    );
    void copy_local_to_global(
      const AssemblyCopyData    &copy_data
    ); 

    void solve();
    void estimate_mark_refine();
    void compute_residual();
    void output_system();

    void prepare_assembly_and_measure_time();
    const std::string convert_time(
      const std::chrono::duration<double>& time
    ) const;

  }; // class ClassicalProblem

} // namespace Classical 

#endif
