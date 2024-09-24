#ifndef PoissonNeumannStandard
#define PoissonNeumannStandard

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
#include <deal.II/fe/mapping_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/base/smartpointer.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/numerics/solution_transfer.h>


// step-6: adapative refinement
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/numerics/error_estimator.h>

#include <poisson_neumann.h>


#include <cmath>
#include <array>
#include <fstream>
#include <iostream>
#include <filesystem>

namespace Poisson_Neumann { 
  using namespace dealt;

  class Poisson_SOL2 : public  Function<2> {
  public: 
    Poisson_SOL2() : Function<2>(1){}
    ~Poisson_SOL2() = default; 

    virtual double value(
      const Point<2>&     p,
      const unsigned int  component = 0
    ) const override;

    virtual Tensor<1, 2> gradient(
      const Point<2>&     p,
      const unsigned int  component = 0
    ) const override;
  }; 

  class Poisson_RHS2 : public Function<2> {
  public:
    Poisson_RHS2() : Function<2>(1) {}
    ~Poisson_RHS2() = default;

    virtual double value(
      const Point<2>&    p,
      const unsigned int component = 0
    ) const override;
  };

  class SquishedGeometry : public ChartManifold<2, 2>{
  public:
  	virtual Point<2> pull_back(
      const Point<2> &space_point
    ) const override;
  	virtual Point<2> push_forward(
      const Point<2> &chart_point
    ) const override;
  	virtual DerivativeForm<1, 2, 2> push_forward_gradient (
      const Point<2> &chart_point
    ) const override;
  	
  	virtual std::unique_ptr<Manifold<2, 2>> clone(
    ) const override;
  };

  class Poisson_Benchmark_Standard {
  	enum Boundary {
      None, 
  		Dirichleth, 
  		Neumann
  	};
  private: 
  	Triangulation<2> 					triangulation;
  	FE_Q<2>                   fe;
  	DoFHandler<2> 	 					dof_handler;
    MappingManifold<2>        map;
    // MappingQ<2>               map;
  	
  	AffineConstraints<double> constraints;	
  	
  	SparsityPattern 	    sparsity_pattern;
  	SparseMatrix<double>  system_matrix;
  	
  	Vector<double> solution;
  	Vector<double> system_rhs;

    Poisson_SOL       sol_fcn;
    Poisson_RHS       rhs_fcn;
    Poisson_NC        nc_fcn;

    double H1 = 1;
    double L2 = 1;

    OutputSetup problem_out;
    // std::filesystem::path directory_parent; 
    // std::filesystem::path directory_degree; 
    // std::filesystem::path directory_svg;    
    // std::filesystem::path directory_vtg;    
  	
  public:
  	Poisson_Benchmark_Standard(const unsigned int& p);
  	void run(
        const unsigned int ref
    );
  	
  private:
  	void make_grid();
  	void setup_system();
  	void assemble_system();
  	void solve();
  	void refine_grid();		
    void compute_h1_error();
  	void process_results(
      const unsigned int   cycle
    );
  	void print_numerical_solution() const ;
    void print_mesh_file() const ;

    void test_face_quadrature_points() const;
  	
  };


  
} // namespace Poisson_Neumann

#endif 
