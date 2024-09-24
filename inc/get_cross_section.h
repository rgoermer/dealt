#include <deal.II/base/point.h>

#include <deal.II/lac/vector.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out.h>

#include <utilities.h>
#include <ts_triangulation.h>
#include <isoparametric_function.h>

#ifndef CROSS_SECTION_H_
#define CROSS_SECTION_H_

namespace dealt {
  namespace GridPartialOut {
    template<int spacedim>
    void output_grid_to_vtk(
      Triangulation<2, spacedim>*       tria, 
      const std::string&                name,
      const IsoparametricManifold<spacedim>* geometry = NULL
    ); // output_grid_to_vtk

    template<int spacedim>
    void output_grid(
      Triangulation<2, spacedim>*       tria, 
      const std::string&                name,
      const IsoparametricManifold<spacedim>* geometry = NULL
    ); // output_grid()


    void get_cross_section(
      TS_TriangulationBase<3, 3>* tria, 
      const unsigned int          direction,
      const double                value,
      const bool                  physical,
      const std::string&          name
    );  // get_cross_section()
  } // GridPartialOut
} // namespace dealt

#endif
