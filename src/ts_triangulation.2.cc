#include <ts_triangulation.h>

namespace dealt {
  // ==========================================================================
  //            Specialization TS_Triangulation<2, spacedim>
  // ==========================================================================

// Constsructors:
  template<int spacedim>
  TS_Triangulation<2, spacedim>::TS_Triangulation(
    const std::vector< std::vector<double> >&      knot_vectors , 
    const std::vector< unsigned int >&             degree ,
    const std::vector< Point<spacedim+1> >& wcps 
  ) : TS_TriangulationBase<2, spacedim>(){
    AssertIndexRange(dimension, space_dimension+1);
    this -> create_triangulation(knot_vectors, degree, wcps);
  }

  template<int spacedim>
  TS_Triangulation<2, spacedim>::TS_Triangulation(
    const IPF_Data<dimension, space_dimension>& data
  ) : TS_TriangulationBase<2, spacedim>(){
    AssertIndexRange(dimension, space_dimension+1);
    this -> create_triangulation(data);
  }
  
  template<int spacedim>
  TS_Triangulation<2, spacedim>::TS_Triangulation(
  ) : TS_TriangulationBase<2, spacedim>(){
    AssertIndexRange(dimension, space_dimension+1);
  }

  // We consider create_triangulation() as a constructor
  template<int spacedim>
  void TS_Triangulation<2, spacedim>::create_triangulation(
      const IPF_Data<dimension, space_dimension>& data
  ) {
    this -> create_triangulation(data.kv, data.deg, data.cps);
  } // create_triangulation [1 / 2]

  template<int spacedim>
  void TS_Triangulation<2, spacedim>::create_triangulation(
      const std::vector< std::vector< double > >&   knot_vectors,
      const std::vector< unsigned int >&            degree,
      const std::vector< Point<spacedim+1> >&       wcps
  ) {
    Assert(knot_vectors.size() == dimension, ExcInternalError("Number of knot vectors does not match dimension of the parametric space. (possibly not enough knot vectors)"));

#ifdef DEBUG
    debug_info_start();
    std::string s("Initializing TS_Triangulation with ");
    s += std::to_string(dimension);
    s += " knot vector(s):";
    std::cout << s << std::endl << std::endl;
#endif

    // Initialize variables
          unsigned int n_cells;
    const unsigned int nvc = GeometryInfo<dimension>::vertices_per_cell;
    std::vector< Point<dimension> > vertices;
    std::vector< std::array<int, nvc> > cell_vertices;

    std::vector< int > sizes(dimension);

    std::vector<double> kx = knot_vectors[0];
    std::vector<double> ky = knot_vectors[1];

    const unsigned int px = degree[0];
    const unsigned int py = degree[1];


    this->ensure_open_knot_vector(kx, px);
    this->ensure_open_knot_vector(ky, py);

    const unsigned int nx = kx.size(), ny = ky.size();
  
    // extract singleton knots from given knot vectors to generate
    // the bezier knot vector
    std::vector<double> kx_bezier({kx[0]});
    std::vector<double> ky_bezier({ky[0]});
    std::vector<unsigned int> mult_x;
    std::vector<unsigned int> mult_y;
    unsigned int count_x = 1, count_y = 1;
    for (unsigned int x = 1; x < nx; x++){
      if (kx[x-1] == kx[x])
        count_x++;
      else {
        mult_x.push_back(count_x);
        kx_bezier.push_back(kx[x]);
        count_x = 1;
      }
    } // for ( x )
    mult_x.push_back(count_x);
  
    for (unsigned int y = 1; y < ny; y++){
      if (ky[y-1] == ky[y])
        count_y++;
      else {
        mult_y.push_back(count_y);
        ky_bezier.push_back(ky[y]);
        count_y = 1;
      }
    } // for ( y )
    mult_y.push_back(count_y);

    this->base_knots = {kx_bezier, ky_bezier};
    this->multiplicities    = {mult_x, mult_y};

  #ifdef DEBUG
      std::cout << indent() << "Knot vector(s) are valid and p-open.\n";
      std::cout << indent() << "Initializing vertices for grid...\n";
  #endif


    {
      Point<dimension> lower(kx[0], ky[0]);
      Point<dimension> upper(kx[nx - 1], ky[ny - 1]);

      this->kv_lower = lower;
      this->kv_upper = upper;
    }

    Assert(wcps.size() == (nx-px-1)*(ny-py-1), ExcInternalError());

    sizes[0] = nx;
    sizes[1] = ny;

    unsigned int nx_bezier = kx_bezier.size();
    unsigned int ny_bezier = ky_bezier.size();

    for (unsigned int y = 0; y < ny_bezier; y++){
      for (unsigned int x = 0; x < nx_bezier; x++){
        Point<dimension> vertex(kx_bezier[x], ky_bezier[y]);
        vertices.push_back(vertex);

  // #ifdef DEBUG_TS
  //         std::cout << indent(2) << "Input vertex: ";
  //         printPoint<dim, double>(vertex, y*nx_bezier + x);
  // #endif

       if (x < nx_bezier - 1 && y < ny_bezier - 1){
          std::array<int, nvc> cell;
          cell[0] = y*nx_bezier + x;
          cell[1] = y*nx_bezier + x + 1;
          cell[2] = (y+1)*nx_bezier + x;
          cell[3] = (y+1)*nx_bezier + x + 1;
          cell_vertices.push_back(cell);
        }
      }
    }

    TSpline::setSplineData(degree);

    // Init Tsplines:
    unsigned int n = 0;
    for (unsigned int y = 0; y < ny - py - 1; y++){
      std::vector< double > knots_y(py + 2);
      for (unsigned int i = 0; i < py + 2; i++) knots_y[i] = ky[i+y];

      for (unsigned int x = 0; x < nx - px - 1; x++){
        std::vector<double> knots_x(px + 2);
        for (unsigned int i = 0; i < px + 2; i++) knots_x[i] = kx[i+x];

        std::vector<std::vector<double>> kv = {knots_x, knots_y};
        this->active_splines.push_back(std::make_shared<TSpline>(kv, wcps[n++]));
      }
    }
    

    this->p = degree;

    n_cells = cell_vertices.size();
    std::vector<CellData<dimension>> cells(n_cells, CellData<dimension>());
    for (unsigned int i = 0; i < n_cells; i++){
      for (unsigned int j = 0; j < nvc; j++){
        cells[i].vertices[j] = cell_vertices[i][j];
      }
    }

  #ifdef DEBUG
    std::cout << indent() << "Constructing grid...\n";
  #endif
    try {
      this->Triangulation<dimension>::create_triangulation(vertices, cells, SubCellData());
    } catch (...) { }
  #ifdef DEBUG
    std::cout << indent() << "... done!\n";
  #endif

  #ifdef DEBUG
    std::cout << indent(2) << "Initializing IPF from Spline Basis ... " << std::endl;
  #endif

    IsoparametricFunction<dimension, space_dimension> other_IPF(this->active_splines);

  #ifdef DEBUG
    std::cout << indent(2) << "... done!" << std::endl;
  #endif

    this->IPF = other_IPF;

    this->setup_mof();
    this -> n_splines = this->active_splines.size();
    for (unsigned int i = 0; i < this->n_splines; i++)
      this->active_splines[i] -> set_level_index(i);

    //set_boundary_dofs();

    // Lastly, define the level offset, depending on the direction of the
    // longest edge:
    active_face_iterator le;
    double el = 0;
    for (const auto& f : this -> active_face_iterators()){
      const Point<dimension>& lower = f -> vertex(0);
      const Point<dimension>& upper = f -> vertex( GeometryInfo<dimension>::vertices_per_face - 1);
      const Point<dimension>& lengths = upper + (-1.)*lower;
      if (lengths(0) > el || lengths(1) > el) {
        le = f;
        el = std::max(lengths(0), lengths(1));
      }
    }

    // figure out which direction of face is the longest
    const Point<dimension>& lengths = (-1.) * le -> vertex(0)
            + le->vertex(GeometryInfo<dimension>::vertices_per_face - 1);
    if (lengths(0) < lengths(1))
      this -> level_offset = 1; // cut along y first
    else
      this -> level_offset = 0;
  } // create_triangulation [2 / 2] 

  template<int spacedim>
  const FullMatrix<double>
    TS_Triangulation<2, spacedim>::get_bezier_coefficients(
      const active_cell_iterator& cell
  ) const {
    const auto&                        operators = this->extraction_operators.at(cell);
    const std::vector< unsigned int >& p         = this->p;
    const unsigned int n_fcns = (p[0] + 1) * (p[1] + 1);
    Assert(n_fcns == operators.size(), ExcInternalError());
    FullMatrix<double> out(n_fcns, n_fcns);

    unsigned int i = 0;
    for (const auto& fcn_row : operators){
      unsigned int ind = 0;
      for (unsigned int y = 0; y < p[1] + 1; y++)
        for (unsigned int x = 0; x < p[0] + 1; x++)
          out(i, ind++) = fcn_row[1][y] *
                            fcn_row[0][x];

      i++;
    } // for ( fcn_row )
    return out;
  } // get_bezier_coefficients() [1 / 2] [for dimension = 2]

  template<int spacedim>
  const std::vector< FullMatrix<double> > 
    TS_Triangulation<2, spacedim>::get_bezier_coefficients(
      const active_cell_iterator& cell,
      const unsigned int face
  ) const {
#ifdef DEBUG
    const auto& face_operators = this->face_operators;
    Assert(face_operators.size() != 0, ExcInternalError());
#endif
    const std::vector<unsigned int>& p = this->p;
    if (!cell -> face(face) -> has_children()){
      // Since we simply copied the data upon construction for
      // faces without children, we simply return the coefficients
      // for the cell
      return {get_bezier_coefficients(cell)};
    } else {
#ifdef DEBUG
      Assert(face_operators.at(cell).size() != 0, ExcInternalError());
#endif
      const auto& operators0 =
              this->face_operators.at(cell).at(cell -> face(face) -> child(0) -> index());
      const auto& operators1 =
              this->face_operators.at(cell).at(cell -> face(face) -> child(1) -> index());
      const unsigned int n_fcns = (p[0] + 1) * (p[1] + 1);
  
      Assert(n_fcns == operators0.size(), ExcInternalError());
      Assert(n_fcns == operators1.size(), ExcInternalError());
      std::vector<FullMatrix<double>> out(2, FullMatrix<double>(n_fcns, n_fcns));

      unsigned int i = 0;
      for (unsigned int j = 0; j < n_fcns; j++){
        const std::array< Vector<double>, dimension>& fcn_row0 =
                operators0[j];
        const std::array< Vector<double>, dimension>& fcn_row1 =
                operators1[j];
        unsigned int ind = 0;
        for (unsigned int y = 0; y < p[1] + 1; y++){
          for (unsigned int x = 0; x < p[0] + 1; x++){
            out[0](i, ind) = fcn_row0[1][y] *
                               fcn_row0[0][x];
            out[1](i, ind) = fcn_row1[1][y] *
                               fcn_row1[0][x];
            ind++;
          }
        }
        i++;
      } // for ( j )
      return out;
    }
  } // get_bezier_coefficients() [2 / 2]

  template<int spacedim>
  void TS_Triangulation<2, spacedim>::get_coarse_neighborhood(
    std::vector<active_cell_iterator>& coarse_neighborhood
  ) {
    const unsigned int nvc = GeometryInfo<dimension>::vertices_per_cell;
    const std::vector<unsigned int>& p = this->p;

    // Iterate over all cells [twice, actually] to generate the coarse neighborhood
    for (const auto& marked_cell : this -> active_cell_iterators()){
      // Is the cell marked for refinement?
      if (!marked_cell -> refine_flag_set() || marked_cell -> has_children())
        continue; // no, continue


      // safety check:
      Assert(RefinementCase<dimension>::cut_axis((this->level_offset + marked_cell -> level()) % dimension) 
              == marked_cell->refine_flag_set(), ExcInvalidState());

      coarse_neighborhood.push_back( marked_cell );
      // for each marked cell, check the criterion for coarse neighborhood,
      // this means, we iterate over all active cells -- again

      // First, get the hadamard product for the open environment
      Point<dimension> boundary;

      const Point<dimension>& lengths = (-1.) * 
                                        marked_cell-> vertex(0) + 
                                        marked_cell -> vertex(nvc - 1);
      int off;
      if ((marked_cell -> level() - this->level_offset) % dimension == 0)
        off = 1;
      else
        off = -1;

      boundary(0) = lengths(0) * (p[0] - off*(p[0]%2) + 1.)/2.;
      boundary(1) = lengths(1) * 1
                               * (p[1] + off*(p[1]%2) + 1.)/2.;

      for (auto cell : this->active_cell_iterators()){
        // some checks can be done pre emptively:
        // 1. is the element already finer than the marked cell?
        // 2. is the cell already marked for refinement?
        // The second condition is more to prevent filling the array for the coarse neighborhood
        // with lots and lots of duplicates.
        if ((cell -> level() !=  marked_cell -> level() - 1) || cell -> refine_flag_set())
          continue;

        // check if the cell is in the open environment:
        Point<dimension> dist_n;
        this->vector_dist(marked_cell, cell, dist_n);

        const bool is_open_env = (dist_n(0) <= boundary(0) &&
                                  dist_n(1) <= boundary(1) );


        if (!is_open_env)
          continue; // Cell is disqualified for refinement as it is not sufficiently close to the cell marked for ref

        const RefinementCase<dimension> ref_case = RefinementCase<dimension>::
                                                   cut_axis((this->level_offset + cell->level()) % dimension);
        cell -> set_refine_flag(ref_case);
        coarse_neighborhood.push_back(cell);

      } // for cell
    } // for marked_cell

    std::sort(coarse_neighborhood.begin(), coarse_neighborhood.end(),
                    [](const active_cell_iterator& c1, const active_cell_iterator& c2){
                      if ( (c1 -> level()) != (c2 -> level()))
                        return c1 -> level() < c2 -> level();
                      else if ( (c1 -> level() % 2) != (c2 -> level() % 2) )
                        return (c1 -> level() % 2) < (c2 -> level() % 2);
                      else
                        return ((c1 -> center()).distance(Point<2>()) < (c2 -> center()).distance(Point<2>()));
                    });

  } // get_coarse_neighborhood() 

  template<int spacedim>
  void TS_Triangulation<2, spacedim>::find_bezier_elements(
  ) {
    // get some references from the base class
    const unsigned int                          nfc             = GeometryInfo<dimension>::faces_per_cell;
    const std::vector<unsigned int>&            p               = this->p;
    const std::map<unsigned int, unsigned int>& mof             = this->mof;
          std::vector< cell_iterator >&         bezier_elements = this->bezier_elements;

    bezier_elements = {};
    for (const auto& cell : this->active_cell_iterators()){
      unsigned int face_no = 0;
      for (; face_no < nfc &&
              !(cell -> face(face_no) -> has_children()); face_no++);

      if (face_no == nfc)
        continue;


      // the hanging interface is on cell -> face(face_no), hence
      // the neighbor on the other side of this face exists and is valid.
      cell_iterator cn_cell = cell;
      const unsigned int neighbor_dir = cell -> neighbor_face_no(face_no);
      unsigned int count = 0;
      std::vector< cell_iterator > local_bezier_elements;
      while (count < (p[neighbor_dir/2]+1)/2  &&
              cn_cell.state() == IteratorState::valid &&
              !(cn_cell -> has_children())
      ){
        // Using the mof-array, we can take knot-multiplicities into account,
        // by using the multiplicity of the face which determines the neighbor
        try {
          count += mof.at(cn_cell -> face(neighbor_dir) -> index());
        } catch (...) {
          std::cout << "Couldn't access mof at face with index " << cn_cell -> face(neighbor_dir) -> index();
          std::cout << std::endl;
          throw ExcInternalError();
        }
        local_bezier_elements.push_back(cn_cell);
        cn_cell = cn_cell -> neighbor(neighbor_dir);
      }

      bezier_elements.insert(bezier_elements.end(),
                              local_bezier_elements.begin(),
                              local_bezier_elements.end());
    } // for ( cell )

    // In case we marked cells as bezier elements, that are already refined,
    // remove them from the list
    auto it = bezier_elements.begin();
    for (; it != bezier_elements.end(); ){
      if ((*it) -> has_children())
        bezier_elements.erase(it);
      else
        ++it;
    } // for ( it )
  }// find_bezier_elements

  // Explicit instantiations
  template class TS_Triangulation<2, 2>;
  template class TS_Triangulation<2, 3>;

} // namespace dealt

