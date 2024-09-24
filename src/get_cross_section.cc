#include <get_cross_section.h>

namespace dealt {
  namespace GridPartialOut { 
    template<int spacedim>
    void output_grid_to_vtk(
            Triangulation<2, spacedim>* tria, 
      const std::string&                name,
      const IsoparametricManifold<spacedim>* geometry
    ) { 
      Triangulation<2, spacedim> other;
      other.copy_triangulation(*tria);

      if (geometry != NULL){
        Assert(spacedim == 3, ExcImpossibleInDim(spacedim));
        GridTools::transform(
          [&geometry](const Point<spacedim>& p) {
            return geometry->push_forward(p);
          }, 
          other);
      }

      DataOut<2, spacedim>  data_out; 
      data_out.attach_triangulation(other); 
      data_out.build_patches();

      std::ofstream vtk_out(name + ".vtk"); 
      data_out.write_vtk(vtk_out);
    }

    template<int spacedim>
    void output_grid(
      Triangulation<2, spacedim>*       , /* tria */
      const std::string&                , /* name */
      const IsoparametricManifold<spacedim>*  /* geometry */
    ) {
      Assert(spacedim == 2 || spacedim == 3, 
              ExcNotImplemented());
    }

    template</* dim == 2 */>
    void output_grid(
      Triangulation<2, 2>*        tria,
      const std::string&          name,
      const IsoparametricManifold<2>*  /* geometry */
    ) { 
      // Print the grid as svg together with vtk
      output_grid_to_vtk(tria, name);

      GridOutFlags::Svg svg_flags;
      svg_flags.label_level_number      = true;
      svg_flags.label_cell_index        = true;
      svg_flags.label_boundary_id       = true;
      svg_flags.boundary_line_thickness = 1;
      svg_flags.line_thickness          = 1;
      svg_flags.coloring    = GridOutFlags::Svg::Coloring::none;
      svg_flags.background  = GridOutFlags::Svg::Background::transparent;


      std::ofstream svg_out(name + ".svg");
      GridOut       grid_out;
      grid_out.set_flags(svg_flags);
      grid_out.write_svg(*tria, svg_out);
    }

    template</* dim == 3 */>
    void output_grid(
      Triangulation<2, 3>*        tria,
      const std::string&          name,
      const IsoparametricManifold<3>*  geometry
    ) {
      output_grid_to_vtk(tria, name, geometry);
      return; 
    }

    void get_cross_section(
      TS_TriangulationBase<3, 3>* tria, 
      const unsigned int          direction,
      const double                value,
      const bool                  physical,
      const std::string&          name
    ) {
      std::vector< bool >               treated_vertices(tria->n_vertices(), false);
      std::vector< int >                global_to_slice(tria->n_vertices(), -1);
      std::vector< Point<3> >           vertices; 
      std::vector< CellData<2> >        cells;

      // Retrieve cell data at specified slice
      unsigned int face         = 2 * direction;
      unsigned int n_vertices_at_slice = 0;
      for (const auto& cell : tria -> active_cell_iterators()) {
        // Does the slice cut through the current cell?
        if ((cell->vertex(0)).operator()(direction) <= value && 
            (cell->vertex(7)).operator()(direction) >  value) {
          // If so, map the vertices of the correspoding face
          // to the slice at value.
          CellData<2>             cell_data;
          for (unsigned int vertex = 0; vertex<4; vertex++) {
            if (!treated_vertices[cell->face(face)->vertex_index(vertex)]){
              Point<3> v = cell->face(face)->vertex(vertex);
              v(direction) = value;

              vertices.push_back( v );
              treated_vertices[cell->face(face)->vertex_index(vertex)] = true;

              global_to_slice[cell->face(face)->vertex_index(vertex)] = n_vertices_at_slice++;
            }
            Assert(global_to_slice[cell->face(face)->vertex_index(vertex)] != -1, 
                            ExcInternalError());
            cell_data.vertices[vertex] = global_to_slice[cell->face(face)->vertex_index(vertex)];
          }
          cells.push_back(cell_data);
        }
      }

      if (vertices.size() == 0)
        return; // Nothing to do! 


      std::string           add_to_name;
      TableIndices<2>       indices;
      std::string           value_string = std::to_string(value);
      value_string.erase(value_string.find_last_not_of('0') + 1, std::string::npos);
      value_string.erase(value_string.find_last_not_of('.') + 1, std::string::npos);
      if (direction == 0) {
        indices[0] = 1;
        indices[1] = 2;
        add_to_name = "_yz_plane_at_x" + value_string;
      } else if (direction == 1) {
        indices[0] = 0;
        indices[1] = 2;
        add_to_name = "_xz_plane_at_y" + value_string;
      } else {
        indices[0] = 0;
        indices[1] = 1;
        add_to_name = "_xy_plane_at_z" + value_string;
      }
      // With the given data, we can create a Triangulation
      if (physical) {
        Triangulation<2, 3> slice;
        slice.create_triangulation(vertices, cells, SubCellData());
        const IsoparametricManifold<3> pf(tria->get_IPF());      

        output_grid(&slice, name + add_to_name, &pf);
      } else {
        Triangulation<2, 2>   slice; 
        std::vector<Point<2>> reduced_vertices(vertices.size(), Point<2>());

        for (unsigned int n = 0; n < vertices.size(); n++){
          reduced_vertices[n](0) = vertices[n](indices[0]);
          reduced_vertices[n](1) = vertices[n](indices[1]);
        }

        slice.create_triangulation(reduced_vertices, cells, SubCellData());
        output_grid(&slice, name + add_to_name); 
      }
    } // get_cross_section()
  } // namespace GridPartialOut
} // namespace dealt
