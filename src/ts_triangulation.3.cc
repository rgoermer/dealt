#include <ts_triangulation.h>

namespace dealt {
  // ==========================================================================
  //            Specialization TS_Triangulation<3, spacedim>
  // ==========================================================================
  

  template<int spacedim>
  TS_Triangulation<3, spacedim>::TS_Triangulation(
  ) : TS_TriangulationBase<3, spacedim>(){
    AssertIndexRange(dimension, space_dimension+1);
  }

  template<int spacedim>
  TS_Triangulation<3, spacedim>::TS_Triangulation(
      const std::vector< std::vector< double > >&    knot_vectors,
      const std::vector< unsigned int >&             degree,
      const std::vector< Point<spacedim+1> >&        wcps
  ) : TS_TriangulationBase<3, spacedim>() {
    AssertIndexRange(dimension, space_dimension+1);
    this -> create_triangulation(knot_vectors, degree, wcps);
  } // constructor from raw data

  template<int spacedim>
  TS_Triangulation<3, spacedim>::TS_Triangulation(
      const IPF_Data<dimension, space_dimension>& data
  ) : TS_TriangulationBase<3, spacedim>() {
    AssertIndexRange(dimension, space_dimension+1);
    this -> create_triangulation(data);
  } // constructor from IPF_Data


  template<int spacedim>
  void TS_Triangulation<3, spacedim>::create_triangulation(
      const IPF_Data<dimension, space_dimension>& data
  ) {
    this -> create_triangulation(data.kv, data.deg, data.cps);
  } // create_triangulation [1 / 2]

  template<int spacedim>
  void TS_Triangulation<3, spacedim>::create_triangulation(
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

          unsigned int n_cells;
    const unsigned int nvc = GeometryInfo<dimension>::vertices_per_cell;
    std::vector< Point<dimension> > vertices;
    std::vector< std::array<int, nvc> > cell_vertices;

    std::vector< unsigned int > sizes(dimension);

    {
      std::vector<double> kx = knot_vectors[0];
      std::vector<double> ky = knot_vectors[1];
      std::vector<double> kz = knot_vectors[2];

      const unsigned int px = degree[0];
      const unsigned int py = degree[1];
      const unsigned int pz = degree[2];

      this->ensure_open_knot_vector(kx, px);
      this->ensure_open_knot_vector(ky, py);
      this->ensure_open_knot_vector(kz, pz);

      unsigned int nx = kx.size();
      unsigned int ny = ky.size();
      unsigned int nz = kz.size();

      std::vector<double> kx_bezier({kx[0]});
      std::vector<double> ky_bezier({ky[0]});
      std::vector<double> kz_bezier({kz[0]});
      std::vector<unsigned int> mult_x;
      std::vector<unsigned int> mult_y;
      std::vector<unsigned int> mult_z;
      unsigned int count_x = 1, count_y = 1, count_z = 1;
      for (unsigned int x = 1; x < nx; x++){
        if (kx[x-1] == kx[x])
          count_x++;
        else {
          mult_x.push_back(count_x);
          kx_bezier.push_back(kx[x]);
          count_x = 1;
        }
      }
      mult_x.push_back(count_x);

      for (unsigned int y = 1; y < ny; y++){
        if (ky[y-1] == ky[y])
          count_y++;
        else {
          mult_y.push_back(count_y);
          ky_bezier.push_back(ky[y]);
          count_y = 1;
        }
      }
      mult_y.push_back(count_y);
      
      for (unsigned int z = 1; z < nz; z++){
        if (kz[z-1] == kz[z])
          count_z++;
        else {
          mult_z.push_back(count_z);
          kz_bezier.push_back(kz[z]);
          count_z = 1;
        }
      }
      mult_z.push_back(count_z);

      this->base_knots = {kx_bezier, ky_bezier, kz_bezier};
      this->multiplicities    = {mult_x, mult_y, mult_z};


  #ifdef DEBUG
      std::cout << indent() << "Knot vector(s) are valid and p-open.\n";
      std::cout << indent() << "Initializing control points..." << std::endl;
  #endif


      {
        Point<dimension> lower(kx[0], ky[0], kz[0]);
        Point<dimension> upper(kx[nx - 1], ky[ny - 1], kz[nz - 1]);

        this->kv_lower = lower;
        this->kv_upper = upper;
      }

      Assert(wcps.size() == (nx-px-1)*(ny-py-1)*(nz-pz-1), ExcInternalError());

      sizes[0] = nx;
      sizes[1] = ny;
      sizes[2] = nz;

      unsigned int nx_bezier = kx_bezier.size();
      unsigned int ny_bezier = ky_bezier.size();
      unsigned int nz_bezier = kz_bezier.size();
      for (unsigned int z = 0; z < nz_bezier; z++){
        for (unsigned int y = 0; y < ny_bezier; y++){
          for (unsigned int x = 0; x < nx_bezier; x++){
            Point<dimension> vertex(kx_bezier[x], ky_bezier[y], kz_bezier[z]);
            vertices.push_back(vertex);

  //#ifdef DEBUG_TS
  //          std::cout << indent(2) << "Input vertex: ";
  //          printPoint<dim, double>(vertex, z*nx_bezier*ny_bezier + y*nx_bezier + x);
  //#endif

            if (x < nx_bezier - 1 && y < ny_bezier - 1 && z < nz_bezier - 1){
              std::array<int, nvc> cell;
              cell[0] = z*nx_bezier*ny_bezier + y*nx_bezier + x;
              cell[1] = z*nx_bezier*ny_bezier + y*nx_bezier + x + 1;
              cell[2] = z*nx_bezier*ny_bezier + (y+1)*nx_bezier + x;
              cell[3] = z*nx_bezier*ny_bezier + (y+1)*nx_bezier + x + 1;
              cell[4] = (z+1)*nx_bezier*ny_bezier + y*nx_bezier + x;
              cell[5] = (z+1)*nx_bezier*ny_bezier + y*nx_bezier + x + 1;
              cell[6] = (z+1)*nx_bezier*ny_bezier +(y+1)*nx_bezier + x;
              cell[7] = (z+1)*nx_bezier*ny_bezier +(y+1)*nx_bezier + x + 1;
              cell_vertices.push_back(cell);
            }
          }
        }
      }

      TSpline::setSplineData(degree);

      // Init Tsplines:
      unsigned int n = 0;
      for (unsigned int z = 0; z < nz - pz - 1; z++){
        std::vector< double > knots_z(pz + 2);
        for (unsigned int i = 0; i < pz + 2; i++) knots_z[i] = kz[i+z];

        for (unsigned int y = 0; y < ny - py - 1; y++){
          std::vector< double > knots_y(py + 2);
          for (unsigned int i = 0; i < py + 2; i++) knots_y[i] = ky[i+y];

          for (unsigned int x = 0; x < nx - px - 1; x++){
            std::vector<double> knots_x(px+2);
            for (unsigned int i = 0; i < px + 2; i++) knots_x[i] = kx[i+x];

            std::vector<std::vector<double>> kv = {knots_x, knots_y, knots_z};
            this->active_splines.push_back(std::make_shared<TSpline>(kv, wcps[n++]));
          }
        }
      }

    } 
    this->p = degree;

    n_cells = cell_vertices.size();
    std::vector<CellData<dimension>> cells(n_cells, CellData<dimension>());
    for (unsigned int i = 0; i < n_cells; i++){
      for (unsigned int j = 0; j < nvc ; j++){
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

    // set_boundary_dofs();

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
    if (lengths(0) < lengths(1)) {
      if (lengths(1) < lengths(2))
        this -> level_offset = 2; // cut along z axis first
      else 
        this -> level_offset = 1; // cut along y first
    } else {
      if (lengths(0) < lengths(2))
        this -> level_offset = 2;
      else 
        this -> level_offset = 0;
    }
  } //create_triangulation [2 / 2] 

  template<int spacedim>
  void TS_Triangulation<3, spacedim>::find_bezier_elements(
  ) {
    // get some references from the base class
    const unsigned int                          nfc             = GeometryInfo<dimension>::faces_per_cell;
    const std::vector<unsigned int>&            p               = this->p;
          std::map<unsigned int, unsigned int>& mof             = this->mof;
          std::vector< cell_iterator >&         bezier_elements = this->bezier_elements;

    bezier_elements = {};
    const uint8_t& cut_x = RefinementCase<3>::cut_x;
    const uint8_t& cut_y = RefinementCase<3>::cut_y;
    for (const auto& cell : this -> active_cell_iterators()){
      //get face numbers, where hanging interfaces are on the
      //associated cell:
      std::vector< unsigned int > face_nos;
      for (unsigned int face_no = 0; face_no < nfc; face_no++)
        if (cell -> face(face_no) -> has_children()) {
          face_nos.push_back(face_no);
          Assert(!(static_cast<uint8_t>(cell->face(face_no)->refinement_case()) & 
                    RefinementCase<dimension-1>::isotropic_refinement), 
                  ExcInternalError());
        }

      if (face_nos.size() == 0)
        continue;

      // get local bezier elements around hanging interfaces:
      std::vector< cell_iterator > local_bezier_elements;
      for (unsigned int face_no : face_nos){
        // First Step: Determine the orthogonal direction [pointing direction is given by face_no]
  
        // Since face face_no has children, the neighbor neighbor(face_no) will be a coarser element,
        // that was once refined. From that, we can deduce the necessary directions.
        // Note, that cell -> neighbor(face_no) is a valid CellAccessor, otherwise, there could not
        // be a hanging edge on that face of the cell
        Assert(static_cast<uint8_t>(cell->neighbor(face_no)->refinement_case()) 
                        & RefinementCase<3>::cut_x || 
               static_cast<uint8_t>(cell->neighbor(face_no)->refinement_case()) 
                        & RefinementCase<3>::cut_y ||
               static_cast<uint8_t>(cell->neighbor(face_no)->refinement_case()) 
                        & RefinementCase<3>::cut_z, 
                dealii::StandardExceptions::ExcMessage(
                "Only anisotropic cuts in a single direction are allowed. \n"
                "When using dealt's autonomous flagging for refinement, this \n"
                "shouldn't happen. If you have flagged cells manually, \n"
                "ensure that refinement is only set in a single direction, \n"
                "i.e. RefinementCase<dim>::cut_x, RefinementCase<dim>::cut_y \n"
                "or RefinementCase<dim>::cut_z in an alternating, level-\n"
                "dependent manner."
                        ));
        Assert(cell->neighbor(face_no)->level() == 1 + cell->level(),
                        ExcInternalError());


        const std::uint8_t ref_case =
                static_cast<std::uint8_t>(
                  cell->neighbor(face_no)
                      ->refinement_case()
                );
        // The few cases that may occur for rdir are listed as follows. Here, rdir is the remaining direction
        // for the T-junction extensions, the orthogonal direction odir will not be used.
        const unsigned int pdir = (cell->neighbor_face_no(face_no))/2;
        const unsigned int rdir = (ref_case & cut_x ? (
                                    pdir == 2 ? 1 : 2
                                  ) : (
                                    ref_case & cut_y ? (
                                      pdir == 2 ? 0 : 2
                                    ) : (
                                      pdir == 1 ? 0 : 1)
                                    )
                                  );


        // The next step is to "walk these directions" using neighbor relationships in the grid
        // We will have two big while-blocks around two smaller while blocks, each quadrant
        // corresponds to a direction we walk along. Again, we use the information in mof to
        // account for knot multiplicities

        // cn_cell will be moved along the grid using neighboring conditions
        cell_iterator cn_cell = cell;

        // pcount counts the "knots" in the pointing direction pdir using the mof array,
        // rcount counts the "knots" in the remaining direction rdir using the mof array.
        unsigned int pcount = 0, rcount = 0;
  
        // For the T-Junction extension, only the T-Junction knot vector is important.
        // Walking the scheme explained above (moving one along pdir and collecting the cells
        // along rdir) we might end up with elements that are out of the T-Junction extension.
        // Hence, we save the maximum distances we walked at the first cell along rdir
        // and break the loops at later points, whenever we exceed that limits
        double rmin, rmax;

        // begin by moving along rdir and collecting bezier elements in this direction
        unsigned int neighbor_rdir = 2*rdir, rcount_max = (p[rdir] + 2 + p[rdir]%2)/2;
        while (cn_cell.state() == IteratorState::valid
                && rcount < rcount_max){
          // Add the amount of knots at the specific face to the counter
          rcount += mof[cn_cell -> face(neighbor_rdir) -> index()];

          // collect the current cell
          if (!cn_cell->has_children()){
            local_bezier_elements.push_back(cn_cell);
            Assert(RefinementCase<3>::cut_axis((this->level_offset+cn_cell->level())%3) 
                    & RefinementCase<3>::cut_axis((this->level_offset+cell->neighbor(face_no)->level())%3),
                   ExcInternalError());
          }

          // and move the cell along rdir, if possible, i.e. if the
          // neighbor actually exists
          if ( (cn_cell -> neighbor(neighbor_rdir)).state() == IteratorState::valid )
            cn_cell = cn_cell -> neighbor(neighbor_rdir);
          else
            break;
        } // while ( cn_cell )
        // Initialize rmin with the last cell
        rmin = (cn_cell -> face(neighbor_rdir) -> vertex(0)).operator()(rdir);

        // reset for opposite direction.
        // Note, that the cell already was pushed into the local_bezier_elements,
        // thus we start with an offset, i.e. the neighbor of the cell
        neighbor_rdir++;
        if (cell -> neighbor(neighbor_rdir).state() == IteratorState::valid){ // special case
          // get the offset cell
          cn_cell = cell -> neighbor(neighbor_rdir);

          // and take the number of faces from the origin cell into account to
          // enter the loop only if necessary
          rcount = mof[cell -> face(neighbor_rdir) -> index()];
          while (cn_cell.state() == IteratorState::valid
                    && rcount < rcount_max){
            // Update counter
            rcount += mof[cn_cell -> face(neighbor_rdir) -> index()];

            // Add element
            if (!cn_cell->has_children()){
              local_bezier_elements.push_back(cn_cell);
              Assert(RefinementCase<3>::cut_axis((this->level_offset+cn_cell->level())%3) 
                      & RefinementCase<3>::cut_axis((this->level_offset+cell->neighbor(face_no)->level())%3),
                     ExcInternalError());
            }

            // check if neighbor is actually valid, and move along
            if ( (cn_cell -> neighbor(neighbor_rdir)).state() == IteratorState::valid )
              cn_cell = cn_cell -> neighbor(neighbor_rdir);
            else
              break;
          } // while ( cn_cell )
          rmax = (cn_cell -> face(neighbor_rdir) -> vertex(0)).operator()(rdir);
        } else {
          rmax = (cell -> face(neighbor_rdir) ->vertex(0)).operator()(rdir);
        }

        // using rmin and rmax, we can now collect the bezier elements row-by-row
        // along pdir without going too far
        unsigned int neighbor_pdir = cell -> neighbor_face_no(face_no);
        unsigned int pcount_max    = (p[pdir] + 1)/2;
        cn_cell                    = cell -> neighbor(neighbor_pdir);
        pcount                     = mof[cell -> face(neighbor_pdir) -> index()];
        while (pcount < pcount_max
                && cn_cell.state() == IteratorState::valid){
          pcount += mof[cn_cell -> face(neighbor_pdir) -> index()];

          // as before, move along rdir and collect elements
          auto cn_rcell = cn_cell;
          neighbor_rdir--; // it was increased before ...

          // first loop:
          rcount = 0;
          while (cn_rcell.state() == IteratorState::valid){
            rcount += mof[cn_rcell -> face(neighbor_rdir) -> index()];
            if (!cn_cell->has_children()){
              local_bezier_elements.push_back(cn_cell);
              Assert(RefinementCase<3>::cut_axis((this->level_offset+cn_cell->level())%3) 
                      & RefinementCase<3>::cut_axis((this->level_offset+cell->neighbor(face_no)->level())%3),
                     ExcInternalError());
            }
            if (rcount < rcount_max &&
                (cn_rcell -> face(neighbor_rdir) -> vertex(0)).operator()(rdir) > rmin)
              cn_rcell = cn_rcell -> neighbor(neighbor_rdir);
            else
              break;
          }

          //second loop:
          neighbor_rdir++;
          cn_rcell = cn_cell -> neighbor(neighbor_rdir);
          rcount   = mof[cell -> face(neighbor_rdir) -> index()];
          while (cn_rcell.state() == IteratorState::valid){
            rcount += mof[cn_rcell -> face(neighbor_rdir) -> index()];
            if (!cn_cell->has_children()){
              local_bezier_elements.push_back(cn_cell);
              Assert(RefinementCase<3>::cut_axis((this->level_offset+cn_cell->level())%3) 
                      & RefinementCase<3>::cut_axis((this->level_offset+cell->neighbor(face_no)->level())%3),
                     ExcInternalError());
            }
            if (rcount < rcount_max &&
                (cn_rcell -> face(neighbor_rdir) -> vertex(0)).operator()(rdir) < rmax )
              cn_rcell = cn_rcell -> neighbor(neighbor_rdir);
            else
              break;
          }

          if (pcount < pcount_max)
            cn_cell = cn_cell -> neighbor(neighbor_pdir);
          else
            break;
        } // while ( cn_cell )
      } // for ( face_no )

      // sort by index and level
      std::sort(local_bezier_elements.begin(), local_bezier_elements.end(),
                      [](const auto& c1, const auto& c2) {
                        if (c1 -> level() != c2 -> level())
                          return c1 -> level() < c2 -> level();
                        else
                          return c1 -> index() < c2 -> index();
                      });

      // // count multiplicities and remove elements, whose count are
      // // less then face_nos.size(), since these elements are not in the
      // // TJunction extension of all hanging interfaces of the associated
      // // cell
      // typename
      // std::vector< cell_iterator >::iterator it = local_bezier_elements.begin();
      // unsigned int it_count = 1;

      // for (; it != local_bezier_elements.end() - 1; ){
      //   Assert(it->state() == IteratorState::valid, ExcInternalError());
      //   Assert((it+1)->state() == IteratorState::valid, ExcInternalError());
      //   if ((*it) == (*(it+1))){
      //     it_count++; ++it;
      //   } else {
      //     if (it_count == face_nos.size()){
      //       it_count = 1; ++it; // nothing to be done ...
      //     } else {
      //       Assert(it_count < face_nos.size(), ExcInternalError());
      //       if (it_count == 1)
      //         it = local_bezier_elements.erase(it);
      //       else 
      //         it = local_bezier_elements.erase(it - it_count + 1, it);

      //       it_count = 1;
      //     } // if ( count )
      //   } // if ( it )
      // } // for
      // Assert(it != local_bezier_elements.end(), ExcInternalError());
      // if (it_count != face_nos.size())
      //   if (it_count == 1)
      //     local_bezier_elements.erase(it);
      //   else 
      //     local_bezier_elements.erase(it-it_count+1, it+1);
      

      // All the remaining cells now have multiplicity face_nos.size()
      // within the array, thus we remove redundancies ...
      // it = local_bezier_elements.begin();
      typename
      std::vector< cell_iterator >::iterator it = local_bezier_elements.begin();
      for (; it != local_bezier_elements.end() - 1; ){
        if ((*it) == (*(it+1))){
          local_bezier_elements.erase(it);
        } else {
          ++it;
        } // if ( id )
      } // for

      // Then, check if any of the remaining cells happen to be already
      // refined, and remove them in this case
      it = local_bezier_elements.begin();
      for ( ; it != local_bezier_elements.end() ; ){
        if ((*it) -> has_children())
          local_bezier_elements.erase(it);
        else
          ++it;
      } // for

      // ... and put the elements into our vector of bezier elements!
      bezier_elements.insert(bezier_elements.end(),
                              local_bezier_elements.begin(),
                              local_bezier_elements.end());

    } // for ( cell )

    // In case we marked cells as bezier elements, that are already refined
    auto it = bezier_elements.begin();
    for (; it != bezier_elements.end(); ){
      if ((*it) -> has_children())
        bezier_elements.erase(it);
      else
        ++it;
    } // for ( it )
  } // find_bezier_elements

  template<int spacedim>
  void TS_Triangulation<3, spacedim>::get_coarse_neighborhood(
      std::vector<active_cell_iterator>& coarse_neighborhood
  ) {
    const std::vector< unsigned int >&  p            = this -> p;
    const unsigned int&                 nvc          = GeometryInfo<dimension>::vertices_per_cell;
    const unsigned int&                 level_offset = this->level_offset;
    // Iterate over all cells [twice, actually] to generate the coarse neighborhood
    for (const auto& marked_cell : this -> active_cell_iterators()){
      // Is the cell marked for refinement?
      if (!marked_cell -> refine_flag_set() || marked_cell -> has_children())
        continue; // no, continue

      // safety check:
      Assert(RefinementCase<dimension>::cut_axis((level_offset+marked_cell->level())%dimension) 
                  == marked_cell -> refine_flag_set(), ExcInvalidState());

      coarse_neighborhood.push_back( marked_cell );
      // for each marked cell, check the criterion for coarse neighborhood,
      // this means, we iterate over alle active cells -- again

      // First, get the hadamard product for the open environment
      Point<dimension> boundary;

      const Point<dimension>& lengths = (-1.) * (marked_cell-> vertex(0))
                                + marked_cell -> vertex(nvc - 1);

      // This definition is different from the 2D definition
      boundary(0) = lengths(0) * (p.at(0) + 1.5);
      boundary(1) = lengths(1) * (p.at(1) + 1.5);
      boundary(2) = lengths(2) * (p.at(2) + 1.5);
      
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
        const bool is_open_env = dist_n(0) <= boundary(0) &&
                                 dist_n(1) <= boundary(1) &&
                                 dist_n(2) <= boundary(2);

        if (!is_open_env)
          continue; // Cell is disqualified for refinement as it is not sufficiently close to the cell marked for ref
        
        const RefinementCase<dimension> ref_case = RefinementCase<dimension>::cut_axis(
                        (level_offset + cell->level()) % dimension);
        cell -> set_refine_flag(ref_case);
        coarse_neighborhood.push_back(cell);
      } // for cell
    } // for marked_cell

    std::sort(coarse_neighborhood.begin(), coarse_neighborhood.end(),
                    [](const active_cell_iterator& c1, const active_cell_iterator& c2){
                      if ( (c1 -> level()) != (c2 -> level()))
                        return c1 -> level() < c2 -> level();
                      else if ( (c1 -> level() % 3) != (c2 -> level() % 3) )
                        return (c1 -> level() % 3) < (c2 -> level() % 3);
                      else
                        return ((c1 -> center()).distance(Point<3>()) < (c2 -> center()).distance(Point<3>()));
                    });
  } // get_coarse_neighborhood

  template<int spacedim>
  const FullMatrix<double> 
    TS_Triangulation<3, spacedim>::get_bezier_coefficients(
      const active_cell_iterator& cell
  ) const {
    const auto&                        operators = this->extraction_operators.at(cell);
    const std::vector< unsigned int >& p         = this->p;
    const unsigned int n_fcns = (p[0] + 1) * (p[1] + 1) * (p[2] + 1);

    if ( n_fcns != operators.size() ){
      std::string message = "An internal check noticed that the number of T-splines on cell "
                          + std::to_string(cell->level()) + "." + std::to_string(cell -> level()) 
                          + " \ngiven by " + std::to_string(operators.size()) + " differs from the "
                          + "assumed amunt of T-splines on any cell given as " + std::to_string(n_fcns)
                          + ".\n\n"
                          + "There is not very much you can do if you encounter this error\n"
                          + " since it indicates an error in deal.T, not in your own program. \n"
                          + "Try to come up with the smallest possible program that still \n"
                          + "demonstrates the error and contact us with it to obtain help.";
      throw ExcMessage(message);
    }
    // Assert(n_fcns == operators.size(), ExcInternalError());

    FullMatrix<double> out(n_fcns, n_fcns);

    unsigned int i = 0;
    for (const auto& fcn_row : operators){
      unsigned int ind = 0;
      for (unsigned int z = 0; z < p[2] + 1; z++)
        for (unsigned int y = 0; y < p[1] + 1; y++)
          for (unsigned int x = 0; x < p[0] + 1; x++)
            out(i, ind++) = fcn_row[2][z] *
                              fcn_row[1][y] *
                              fcn_row[0][x];

      i++;
    } // for ( fcn_row )
    return out;
  } // get_bezier_coefficients with dim = 3


  template<int spacedim>
  const std::vector< FullMatrix<double> > 
    TS_Triangulation<3, spacedim>::get_bezier_coefficients(
      const active_cell_iterator& cell,
      const unsigned int face
  ) const {
    const std::vector<unsigned int>& p = this->p;
    if ( cell -> face(face) -> has_children()){
      const auto& operators0 =
              this->face_operators.at(cell).at(cell -> face(face) -> child(0) -> index());
      const auto& operators1 =
              this->face_operators.at(cell).at(cell -> face(face) -> child(1) -> index());
      const unsigned int n_fcns = (p[0] + 1) * (p[1] + 1) * (p[2] + 1);

      Assert(n_fcns == operators0.size(), ExcInternalError());
      Assert(n_fcns == operators1.size(), ExcInternalError());

      std::vector< FullMatrix<double> > out(2, FullMatrix<double>(n_fcns, n_fcns));
      
      unsigned int i = 0;
      for (unsigned int j = 0; j < n_fcns; j++) {
        std::array< Vector<double>, dimension> fcn_row0 =
                operators0[j];
        std::array< Vector<double>, dimension> fcn_row1 =
                operators1[j];
        unsigned int ind = 0;
        for (unsigned int z = 0; z < p[2] + 1; z++) {
          for (unsigned int y = 0; y < p[1] + 1; y++) {
            for (unsigned int x = 0; x < p[0] + 1; x++) {
              out[0](i, ind) = fcn_row0[2][z] *
                                fcn_row0[1][y] *
                                fcn_row0[0][x];
              out[1](i, ind) = fcn_row1[2][z] *
                                fcn_row1[1][y] *
                                fcn_row1[0][x];
              ind++;
            }
          }
        }

        i++;
      } // for ( fcn_row )
      return out;
    } else {
      // We simply copied the necessary data upon creation from
      // the cell, hence, we return that data
      return {get_bezier_coefficients(cell)};
    }
  } // get_bezier_coefficients with dim = 3
  
  // Explicit instantiations
  template class TS_Triangulation<3, 3>;

} // namespace dealt
