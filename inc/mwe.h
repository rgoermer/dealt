#ifndef MWE_H_
#define MWE_H_

#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/data_out.h>

#include <iostream>

template<unsigned int dim>
class MWE {
  	Triangulation<dim> 					triangulation;
  	DoFHandler<dim> 	 					dof_handler;
  	FE_Q<dim>                   fe;

  	Vector<double> solution;
  	
public:
    MWE(const unsigned int p) 
      : dof_handler(triangulation) 
      , fe(p)
    {
      // Make a triangulation
    	GridGenerator::hyper_cube(triangulation);

      // Setup system: 
      dof_handler.distribute_dofs(fe);
      solution.reinit(dof_handler.n_dofs());

      // Refining 
      triangulation.refine_global(3); 

  	  DataOut<dim> data_out;
  	  data_out.attach_dof_handler(dof_handler);
  	  data_out.add_data_vector(solution, "solution");
  	  data_out.build_patches();

      std::string vtu_name = "mwe.vtu";
      std::ofstream vtu_out(vtu_name);
  	  data_out.write_vtu(vtu_out);
    }
}; 

#endif
