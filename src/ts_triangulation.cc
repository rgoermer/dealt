/*
 * TS_Triangulation.cc
 *
 *  Created on: Mar 27, 2020
 *      Author: goermer
 */

#include <ts_triangulation.h>

namespace dealt {


  // =======================================================================
  // =======================================================================
  //                  TS_TriangulationBase implementation
  // =======================================================================
  // =======================================================================
  //
  // =============================================================
  //        Getter / Setter section
  // =============================================================
  template<int dim, int spacedim>
  const IsoparametricFunction<dim, spacedim>& TS_TriangulationBase<dim, spacedim>::get_IPF(
  ) const {
    return IPF;
  } // get_IPF()

  template<int dim, int spacedim>
  const unsigned int& 
      TS_TriangulationBase<dim, spacedim>::n_active_splines(
  ) const {
    return n_splines;
  } // n_active_splines()

  template<int dim, int spacedim>
  const std::vector< unsigned int >& 
      TS_TriangulationBase<dim, spacedim>::get_degree(
  ) const {
    return p;
  } // get_degree() [1 / 2]
  
  template<int dim, int spacedim>
  const unsigned int& 
      TS_TriangulationBase<dim, spacedim>::get_degree(
    const unsigned int d
  ) const {
    return p.at(d);
  } // get_degree() [2 / 2]

  template<int dim, int spacedim>
  auto TS_TriangulationBase<dim, spacedim>::get_IEN_array(
  ) const -> const std::map<
          const active_cell_iterator, 
          std::vector< types::global_dof_index > 
                     >& {
#ifdef DEBUG
    Assert(is_bezier_mesh, ExcInvalidState());
#else 
    if (!is_bezier_mesh)
      throw ExcInvalidState();
#endif
    return IEN_array;
  } // get_IEN_array() [1 / 4]

  template<int dim, int spacedim>
  const std::vector< types::global_dof_index >& 
      TS_TriangulationBase<dim, spacedim>::get_IEN_array(
    const active_cell_iterator& cell
  ) const {
#ifdef DEBUG
    Assert(is_bezier_mesh, ExcInvalidState());
    return IEN_array.at(cell);
#else
    if (!is_bezier_mesh)
      throw ExcInvalidState();

    return IEN_array.at(cell);
#endif
  } // get_IEN_array() [2 / 4]

  template<int dim, int spacedim>
  auto TS_TriangulationBase<dim, spacedim>::get_IEN_array(
    const unsigned int n_components
  ) const -> std::map<
          const active_cell_iterator, 
          std::vector< types::global_dof_index > 
                     > {
#ifdef DEBUG
    Assert(is_bezier_mesh, ExcInvalidState());
#else 
    if (!is_bezier_mesh)
      throw ExcInvalidState();
#endif
    std::map<
      const active_cell_iterator,
            std::vector< types::global_dof_index >
            > out_map; 
    for (const auto& [cell, arr] : IEN_array){
      std::vector< types::global_dof_index > cell_array;
      for (const unsigned int i : arr)
        for (unsigned int n = 0; n < n_components; n++)
          cell_array.push_back(n_components * i + n);
      out_map[cell] = cell_array;
    }
    return out_map;
  } // get_IEN_array() [3 / 4]

  template<int dim, int spacedim>
  std::vector< types::global_dof_index > 
      TS_TriangulationBase<dim, spacedim>::get_IEN_array(
    const active_cell_iterator& cell,
    const unsigned int n_components
  ) const {
#ifdef DEBUG
    Assert(is_bezier_mesh, ExcInvalidState());
#else
    if (!is_bezier_mesh)
      throw ExcInvalidState();
#endif
    const std::vector< types::global_dof_index >& cell_array = IEN_array.at(cell);
          std::vector< types::global_dof_index >  out;
    for (const unsigned int i : cell_array)
      for (unsigned int n = 0; n < n_components; n++)
        out.push_back(n_components * i + n);

    return out;
  } // get_IEN_array() [4 / 4]

  template<int dim, int spacedim>
  auto TS_TriangulationBase<dim, spacedim>::get_splines(
  ) const -> const std::vector< ts_ptr >& {
    return active_splines;
  } // get_splines() [1 / 2]

  template<int dim, int spacedim>
  auto TS_TriangulationBase<dim, spacedim>::get_splines(
    const  active_cell_iterator& cell
  ) const -> const std::vector< ts_ptr > {
    const std::vector< unsigned int >& splines = get_IEN_array(cell);
    std::vector< ts_ptr > out;
    unsigned int n_splines = splines.size();
    for (unsigned int j = 0; j < n_splines; j++)
      out.push_back(active_splines[splines[j]]);

    return out;
  } // get_splines() [2 / 2]

  template<int dim, int spacedim>
  std::pair< const Point<dim>&, const Point<dim>&> 
      TS_TriangulationBase<dim, spacedim>::get_bounding_box(
  ) const {
    return std::make_pair<const Point<dim>&, const Point<dim>& >(kv_lower, kv_upper);
  } // get_bounding_box()


  template<int dim, int spacedim>
  const std::map<
      types::boundary_id,
      std::vector< types::global_dof_index >
    >& TS_TriangulationBase<dim, spacedim>::get_boundary_dofs(
  ) const {
    return boundary_dofs;
  } // get_boundary_dofs() [1 / 2]

  template<int dim, int spacedim>
  std::map<
      types::boundary_id,
      std::vector< types::global_dof_index >
    > TS_TriangulationBase<dim, spacedim>::get_boundary_dofs(
    const unsigned int n_components
  ) const {
    std::map<
      types::boundary_id,
      std::vector< types::global_dof_index >
      > bdry_dofs;
    for (const auto& [b_id, dofs] : boundary_dofs){
      std::vector< types::global_dof_index > populated_dofs;
      for (const auto& dof : dofs)
        for (unsigned int n = 0; n < n_components; n++)
          populated_dofs.push_back(n_components * dof + n);

      bdry_dofs[b_id] = populated_dofs;
    } 
    return bdry_dofs;
  } // get_boundary_dofs() [2 / 2]

  template<int dim, int spacedim>
  auto TS_TriangulationBase<dim, spacedim>::get_bezier_elements(
  ) const -> const std::vector< cell_iterator >& {
    return bezier_elements;
  } // get_bezier_elements()

  template<int dim, int spacedim>
  auto TS_TriangulationBase<dim, spacedim>::get_control_points(
      const active_cell_iterator& cell
  ) const -> const std::vector< Point<space_dimension+1> > {
    const std::vector< unsigned int >& ts_indices = get_IEN_array(cell);
    std::vector< Point<spacedim+1> > out(ts_indices.size());

    for (unsigned int i = 0; i < ts_indices.size(); i++)
      out[i] = active_splines[ts_indices[i]] -> get_cp();

    return out;
  } // get_control_points()
  
  template<int dim, int spacedim>
  unsigned int 
    TS_TriangulationBase<dim, spacedim>::get_max_nnz_supports(
  ) const {
    unsigned int deg = p[0]; // Assert all degrees are equal

  #ifdef DEBUG
      Assert(deg == p[1], ExcMessage("The maximum number of non-zero support intersections is only given if the polynomial degree is equal in every direction."));
  #endif

    return (deg*(deg+1)*5 + 1 < n_splines) ? deg*(deg+1)*5 + 1 : n_splines;
  } // get_max_nnz_supports()

  template<int dim, int spacedim>
  std::vector< unsigned int > 
    TS_TriangulationBase<dim, spacedim>::get_max_entries_per_row(
    const unsigned int n_components
  ) const {
    const unsigned int n_dofs = n_components * this->n_active_splines();
    const auto& IEN_array = this->get_IEN_array(n_components);
    
    std::vector<unsigned int> out(n_dofs, 0);
    std::vector<std::vector<bool>> 
          pattern(n_dofs, std::vector<bool>(n_dofs, false));

    for (const auto& [_, arr] : IEN_array) 
      for (const unsigned int i : arr)
        for (const unsigned int j : arr)
          if (!pattern[i][j]){
            pattern[i][j] = true;
            out[i] += 1;
          }


    return out;
  }

  template<int dim, int spacedim>
  const std::map<unsigned int, unsigned int>& 
      TS_TriangulationBase<dim, spacedim>::get_mof(
  ) const {
    return mof;
  } // get_mof()

  // =============================================================
  //        Assembly section
  // =============================================================

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::set_refine_flags(
      const std::vector< cell_iterator >& mark
  ) {
    for (const auto& cell : mark) {
#ifdef DEBUG
      bool marked_cell_belongs_to_triangulation = false;
      for (const auto& c : this -> active_cell_iterators()) {
        marked_cell_belongs_to_triangulation = (c == cell); 
        if (marked_cell_belongs_to_triangulation)
          break;
      }
      Assert(marked_cell_belongs_to_triangulation, ExcMessage("Cells of different Triangulation given!"));
#endif
      cell -> set_refine_flag(
              RefinementCase<dim>::cut_axis((level_offset + cell -> level()) % dim )
          );
    }
  } // set_refine_flags()

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::prepare_assembly(
  ) {
    this->set_boundary_dofs();
    this->refine_bezier_elements();
    this->compute_extraction_operators();
  } // TS_TriangulationBase<dim, spacedim>::prepare_assembly();

  template<int dim, int spacedim>
  template<int soldim>
  void TS_TriangulationBase<dim, spacedim>::project_boundary_values(
    const std::map< 
              types::boundary_id, 
              const Function<spacedim, double>* 
          >& boundary_functions,
    const std::vector<unsigned int>& n_gauss_points,
          std::map< 
              types::global_dof_index,
              double 
          >& boundary_values
  ) const {
    Assert(soldim == 1 || soldim == spacedim, ExcDimensionMismatch2(soldim, 1, spacedim));

    // Prepare Degrees of Freedom, i.e. 
    // count total amount of DoFs on boundaries specified by
    // boundary_functions and store DoFs seperately
    unsigned int n_boundary_dofs = 0; 
    std::map<types::global_dof_index, unsigned int/*local_dof_index*/>
            global_to_local;
    std::map<unsigned int/*local_dof_index*/, types::global_dof_index>
            local_to_global;
    { 
      const auto& boundary_dofs = this -> get_boundary_dofs(soldim);
      for (const auto& [b_id, _] : boundary_functions){
        // Make sure that there are DoFs assigned to the specified 
        // boundary id
        Assert(boundary_dofs.find(b_id) != boundary_dofs.end(), 
                ExcMessage("You are trying to access boundary DoFs through\n boundary_functions, which was not specified at\n the definition of the boundary DoFs!"));
        for (const auto& dof : boundary_dofs.at(b_id)){
          // On vertices or edges (3D), there may be DoFs belonging to 
          // two types of boundaries. If the dof was already added to 
          // the set of boundary dofs, we skip it.
          if (global_to_local.find(dof) != global_to_local.end()) continue;
          else global_to_local[dof] = n_boundary_dofs++; 
        }
      }
      // local_to_global is needed later to move the solution vector to
      // the outbout boundary_values
      for (const auto& [global, local] : global_to_local)
        local_to_global[local] = global;
    }

    // std::cout << "local_to_global: ";
    // for (const auto& [_, dof] : local_to_global)
    //   printf("%2i ", dof);
    // std::cout << std::endl;

    // We only compute face values for the projector
    TSFaceValues<dim, spacedim, soldim> face_values(
      this, 
      n_gauss_points, 
      update_values | 
      update_JxW_values |
      update_quadrature_points
    );

    // Setup system matrix
    SparsityPattern sparsity_pattern(
      n_boundary_dofs, 
      n_boundary_dofs,
      n_boundary_dofs
    );
    const auto& IEN_array = this -> get_IEN_array(soldim); 
    for (const auto& [_, arr] : IEN_array){
      for (const unsigned int n : arr){
        for (const unsigned int m : arr){
          if (global_to_local.find(n) != global_to_local.end() 
                && global_to_local.find(m) != global_to_local.end())
            sparsity_pattern.add(global_to_local[n], global_to_local[m]); 
        }
      }
    }
    sparsity_pattern.compress();
    SparseMatrix<double> system_matrix(sparsity_pattern);
    Vector<double>       system_rhs(n_boundary_dofs);

    // As in the assembly store values in a cell-matrix
    FullMatrix<double> cell_matrix(face_values.n_dofs_per_cell(),
                                    face_values.n_dofs_per_cell());
    Vector<double>     cell_rhs(face_values.n_dofs_per_cell());

    // We need to find u_h \in V_h, s.t.
    //  \int_\Gamma \varphi_i u_h = \sum_{\Gamma_k} \int_{\Gamma_k} f_k \varphi_i
    // for all \varphi_i. This yields a system of equations to be solved, as we have 
    //  u_h = \sum_{j} c_j \varphi_j
    for (const auto& cell : this -> active_cell_iterators()){
      if (cell -> at_boundary()){
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; f++){
          const auto& face = cell -> face(f);
          const auto& s = boundary_functions.find(face -> boundary_id());
          if (face -> at_boundary() && s != boundary_functions.end()){
            face_values.reinit(cell, f);
            cell_matrix = 0;
            cell_rhs    = 0;
            const auto& bdry_fcn = s -> second;
            for (const unsigned int q : face_values.quadrature_point_indices()){
              // get the function value
              const double dx         = face_values.JxW(q);
              for (const unsigned int i : face_values.dof_indices()){
                // Get the index of the non-zero component. For soldim = 1, this is just 0.
                const unsigned int comp_i =
                  face_values.system_to_component_index(i).first; 
                const double bdry_value = 
                  bdry_fcn -> value(face_values.quadrature_point(q), comp_i);

                AssertIsFinite(bdry_value);
                AssertIsFinite(dx);
                AssertIsFinite(face_values.shape_value(i, q));
                for (const unsigned int j : face_values.dof_indices()) {
                  const unsigned int comp_j =
                    face_values.system_to_component_index(j).first;
                  if (comp_i == comp_j)
                    cell_matrix(i, j) += face_values.shape_value(i, q) *
                                         face_values.shape_value(j, q) * 
                                         dx;
                } // for ( j )
                cell_rhs(i) += face_values.shape_value(i, q) * 
                               bdry_value * 
                               dx;
              } // for ( i )
            } // for ( q )

            // At the end of q-loop, move entries to the system matrix. 
            // 1. Get the local IEN_array
            const auto& local_IEN_array = this -> get_IEN_array(cell, soldim); 
            
            // 2. Translate it to the local boundary dofs
            for (unsigned int n = 0; n < face_values.n_dofs_per_cell(); n++){
              if (global_to_local.find(local_IEN_array[n]) == global_to_local.end()) 
                continue; 

              const unsigned int& i = global_to_local[local_IEN_array[n]]; 
              for (unsigned int m = 0; m < face_values.n_dofs_per_cell(); m++){
                if (global_to_local.find(local_IEN_array[m]) != global_to_local.end()){
                  const unsigned int& j = global_to_local[local_IEN_array[m]]; 
                  system_matrix.add( i, j, cell_matrix(n, m) );  
                }
              }

              system_rhs(i) += cell_rhs(n);
            } // for ( n )
          } // if ( face -> at_boundary ?)
        } // for ( f )
      } // if ( cell -> at_boundary? )
    } // for ( cell )
    

    // Now, the system is ready to be solved: 
    // 1. Define a solver
    SolverControl             solver_control( 750 * n_boundary_dofs, 1e-15 );
    SolverCG<Vector<double>>  solver(solver_control);

    // 2. Use a diagonal preconditioner
    PreconditionJacobi<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix); 

    // 3. Define a helper vector for the solution
    Vector<double>    solution(n_boundary_dofs);

    // 4. Solve the system
    solver.solve(system_matrix, solution, system_rhs, preconditioner); 

    // 5. Move values to output
    for(unsigned int local_dof = 0;
          local_dof < n_boundary_dofs; 
          local_dof++)
      boundary_values[local_to_global[local_dof]] = solution(local_dof);


  } // project_boundary_values()


  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::
    poisson_residual_error_estimate(
      const std::vector< unsigned int >&          n_gauss_points,
      const Function<space_dimension>*            rhs_fcn,
      const Function<space_dimension>*            sigma,
      const std::map< 
                    types::boundary_id,
              const Function<space_dimension>* 
            >&                                    neumann_bc,
      const Vector<double>&                       solution,
      std::map< cell_iterator,
                double >&                         residuals
  ) const {
    TSValues<dimension, space_dimension> ts_values(
      this, n_gauss_points,
      update_values |
      update_gradients |
      update_quadrature_points |
      update_JxW_values | 
      update_hessians
    );
    TSFaceValues<dimension, space_dimension> face_values(
      this, n_gauss_points,
      update_values |
      update_gradients |
      update_quadrature_points |
      update_JxW_values  |
      update_normal_vectors);
    const unsigned int nvc = GeometryInfo<dimension>::vertices_per_cell;
    const unsigned int nvf = GeometryInfo<dimension>::vertices_per_face;

    std::map< unsigned int, double > face_residuals;
    for (const auto& cell : this -> active_cell_iterators()) {
      std::vector< unsigned int > local_dof_indices = get_IEN_array(cell);

      // Compute the jumps of faces at C0 continuity
      for (unsigned int f = 0;
              f < GeometryInfo<dimension>::faces_per_cell;
              f++) {
        // Get orientation of current face:
        const auto& face = cell -> face(f);
        const Point<dimension>& c = (-1.) * 
                                    face -> vertex(0) +
                                    face -> vertex(nvf - 1);
        unsigned int orientation = 0;
        for ( ; orientation < dimension && c(orientation) != 0; orientation++);

        if (face -> has_children() ) {
          // This is merely a special case, as refined faces are considered
          // in-active and hence their index is not represented in mof
          if (mof.at(face -> child(0) -> index()) == p[orientation]) {
            // C0 Edge detected:
            face_values.reinit(cell, f);
            const unsigned int ppf = 2 * face_values.n_quadrature_points_per_face();
            for (const unsigned int q_index : face_values.quadrature_point_indices()){
              double g = 0;
              for (const unsigned int i : face_values.dof_indices())
                g += solution(local_dof_indices[i]) *
                        face_values.shape_grad(i, q_index) *
                        face_values.normal_vector(q_index);

              const unsigned int ch = q_index > ppf ? 1 : 0;
              face_residuals[face -> child(ch) -> index()] +=
                      0.25 * g * g * face_values.JxW(q_index);
            } // for ( q_index )
          } // otherwise it is atleast a C1 edge, and cannot be at the boundary
        } else if (mof.at(face -> index()) == p[orientation]){
          // C0 edge detected
          face_values.reinit(cell, f);

          // Since this face is not refined, we can simply compute
          // the jump terms along it abd store the values at the
          // corresponding place in face_residuals
          for (const unsigned int q_index : face_values.quadrature_point_indices()){
            double g = 0;
            for (const unsigned int i : face_values.dof_indices())
              g += solution(local_dof_indices[i]) *
                      face_values.shape_grad(i, q_index) *
                      face_values.normal_vector(q_index);

            face_residuals[face -> index()] += 0.25 * g * g * face_values.JxW(q_index);
          } // for ( q_index )

        } else if (face -> at_boundary()) {
          const auto& bc = neumann_bc.find(face -> boundary_id());
          if (bc == neumann_bc.end())
            continue;

          face_values.reinit(cell, f);
          // If the face is at the boundary, compute the jump from
          // the face's boundary id
          for (const unsigned int q_index : face_values.quadrature_point_indices()){
            double g = neumann_bc.at(face -> boundary_id())
                          -> value(face_values.quadrature_point(q_index));
            for (const unsigned int i : face_values.dof_indices())
              g -= solution(local_dof_indices[i]) *
                      face_values.shape_grad(i, q_index) *
                      face_values.normal_vector(q_index);

            face_residuals[face -> index()] += g * g * face_values.JxW(q_index);
          } // for ( q_index )
        } // if ( ... )
      } // for ( face )
    } // for ( cell )

    // After every Jump on every face is computed,
    // we compute the residuals on the cells
    for (const auto& cell : this -> active_cell_iterators()){
      double& local_residual = residuals[cell];
      const std::vector< unsigned int >& local_dof_indices = get_IEN_array(cell);
      ts_values.reinit(cell);
      for ( const unsigned int q_index : ts_values.quadrature_point_indices()){
        double slu = 0;
        double gsgu = 0;

        const Point<space_dimension>& Q = ts_values.quadrature_point(q_index);
        const Tensor<1, space_dimension>& sigma_grad = sigma -> gradient(Q);
        const double sigma_val = sigma -> value(Q);
        for (const unsigned int i : ts_values.dof_indices() ){
          double l = 0;
          for (unsigned int d = 0; d < space_dimension; d++)
            l += ts_values.shape_hessian(i, q_index)[d][d];

          slu   += solution(local_dof_indices[i]) * sigma_val * l;
          gsgu  += solution(local_dof_indices[i]) *
                      sigma_grad *
                      ts_values.shape_grad(i, q_index);
        } // for ( i )

        double g = (slu + gsgu + rhs_fcn -> value(Q));
        local_residual += g * g * ts_values.JxW(q_index);
      } // for ( q_index )
      const double cell_width = (cell->vertex(0)).distance(
                                cell->vertex(nvc - 1));

      local_residual *= cell_width * cell_width;

      // And then also add the values from the faces:
      for (unsigned int f = 0;
              f < GeometryInfo<dimension>::faces_per_cell;
              f++) {
        if (cell -> face(f) -> has_children()){
          const double child0_width = 
             (cell -> face(f) -> child(0) -> vertex(0)).distance(
              cell -> face(f) -> child(0) -> vertex(nvf - 1));
          const double child1_width = 
             (cell -> face(f) -> child(1) -> vertex(0)).distance(
              cell -> face(f) -> child(1) -> vertex(nvf - 1));
          local_residual += child0_width *
                  face_residuals[cell -> face(f) -> child(0) -> index()];
          local_residual += child1_width *
                  face_residuals[cell -> face(f) -> child(1) -> index()];
        } else {
          const double face_width = 
              (cell -> face(f) -> vertex(0)).distance(
               cell -> face(f) -> vertex(nvf - 1));
          local_residual += face_width *
                  face_residuals[cell -> face(f) -> index()];
        }
      } // for ( f )

      local_residual = std::sqrt(local_residual);
    } // for ( cell )
  } // poisson_residual_error_estiamte() [1 / 4]

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::
    poisson_residual_error_estimate(
      const std::vector< unsigned int >&          n_gauss_points,
      const Function<space_dimension>*            rhs_fcn,
      const Function<space_dimension>*            sigma,
      const Function<space_dimension>*            a,
      const std::map< 
                    types::boundary_id,
              const Function<space_dimension>* 
            >&                                    neumann_bc,
      const Vector<double>&                       solution,
            Vector<double>&                       residuals
  ) const {
    TSValues<dimension, space_dimension> ts_values(
      this, n_gauss_points,
      update_values |
      update_gradients |
      update_quadrature_points |
      update_JxW_values | 
      update_hessians
    );
    TSFaceValues<dimension, space_dimension> face_values(
      this, n_gauss_points,
      update_values |
      update_gradients |
      update_quadrature_points |
      update_JxW_values  |
      update_normal_vectors);
    const unsigned int nvc = GeometryInfo<dimension>::vertices_per_cell;
    const unsigned int nvf = GeometryInfo<dimension>::vertices_per_face;

    std::map< unsigned int, double > face_residuals;
    for (const auto& cell : this -> active_cell_iterators()) {
      std::vector< unsigned int > local_dof_indices = get_IEN_array(cell);

      // Compute the jumps of faces at C0 continuity
      for (unsigned int f = 0;
              f < GeometryInfo<dimension>::faces_per_cell;
              f++) {
        // Get orientation of current face:
        const auto& face = cell -> face(f);
        const Point<dimension>& c = (-1.) * 
                                    face -> vertex(0) +
                                    face -> vertex(nvf - 1);
        unsigned int orientation = 0;
        for ( ; orientation < dimension && c(orientation) != 0; orientation++);

        if (face -> has_children() ) {
          // This is merely a special case, as refined faces are considered
          // in-active and hence their index is not represented in mof
          if (mof.at(face -> child(0) -> index()) == p[orientation]) {
            // C0 Edge detected:
            face_values.reinit(cell, f);
            const unsigned int ppf = 2 * face_values.n_quadrature_points_per_face();
            for (const unsigned int q_index : face_values.quadrature_point_indices()){
              double g = 0;
              for (const unsigned int i : face_values.dof_indices())
                g += solution(local_dof_indices[i]) *
                        face_values.shape_grad(i, q_index) *
                        face_values.normal_vector(q_index);

              const unsigned int ch = q_index > ppf ? 1 : 0;
              face_residuals[face -> child(ch) -> index()] +=
                      0.25 * g * g * face_values.JxW(q_index);
            } // for ( q_index )
          } // otherwise it is atleast a C1 edge, and cannot be at the boundary
        } else if (mof.at(face -> index()) == p[orientation]){
          // C0 edge detected
          face_values.reinit(cell, f);

          // Since this face is not refined, we can simply compute
          // the jump terms along it abd store the values at the
          // corresponding place in face_residuals
          for (const unsigned int q_index : face_values.quadrature_point_indices()){
            double g = 0;
            for (const unsigned int i : face_values.dof_indices())
              g += solution(local_dof_indices[i]) *
                      face_values.shape_grad(i, q_index) *
                      face_values.normal_vector(q_index);

            face_residuals[face -> index()] += 0.25 * g * g * face_values.JxW(q_index);
          } // for ( q_index )

        } else if (face -> at_boundary()) {
          const auto& bc = neumann_bc.find(face -> boundary_id());
          if (bc == neumann_bc.end())
            continue;

          face_values.reinit(cell, f);
          // If the face is at the boundary, compute the jump from
          // the face's boundary id
          for (const unsigned int q_index : face_values.quadrature_point_indices()){
            double g = neumann_bc.at(face -> boundary_id())
                          -> value(face_values.quadrature_point(q_index));
            for (const unsigned int i : face_values.dof_indices())
              g -= solution(local_dof_indices[i]) *
                      face_values.shape_grad(i, q_index) *
                      face_values.normal_vector(q_index);

            face_residuals[face -> index()] += g * g * face_values.JxW(q_index);
          } // for ( q_index )
        } // if ( ... )
      } // for ( face )
    } // for ( cell )

    // After every Jump on every face is computed,
    // we compute the residuals on the cells
    unsigned int kk = 0;
    for (const auto& cell : this -> active_cell_iterators()){
      double& local_residual = residuals(kk++);
      const std::vector< unsigned int >& local_dof_indices = get_IEN_array(cell);
      ts_values.reinit(cell);
      for ( const unsigned int q_index : ts_values.quadrature_point_indices()){
        double slu  = 0;
        double gsgu = 0;
        double au   = 0;

        const Point<space_dimension>& Q = ts_values.quadrature_point(q_index);
        const Tensor<1, space_dimension>& sigma_grad = sigma -> gradient(Q);
        const double sigma_val = sigma -> value(Q);
        const double a_val     = a -> value(Q);
        for (const unsigned int i : ts_values.dof_indices() ){
          double l = 0;
          for (unsigned int d = 0; d < space_dimension; d++)
            l += ts_values.shape_hessian(i, q_index)[d][d];

          slu   += solution(local_dof_indices[i]) * sigma_val * l;
          gsgu  += solution(local_dof_indices[i]) *
                      sigma_grad *
                      ts_values.shape_grad(i, q_index);
          au    += solution(local_dof_indices[i]) * a_val;
        } // for ( i )

        double g = (slu + gsgu - au + rhs_fcn -> value(Q));
        local_residual += g * g * ts_values.JxW(q_index);
      } // for ( q_index )
      const double cell_width = (cell->vertex(0)).distance(
                                cell->vertex(nvc - 1));

      local_residual *= cell_width * cell_width;

      // And then also add the values from the faces:
      for (unsigned int f = 0;
              f < GeometryInfo<dimension>::faces_per_cell;
              f++) {
        if (cell -> face(f) -> has_children()){
          const double child0_width = 
             (cell -> face(f) -> child(0) -> vertex(0)).distance(
              cell -> face(f) -> child(0) -> vertex(nvf - 1));
          const double child1_width = 
             (cell -> face(f) -> child(1) -> vertex(0)).distance(
              cell -> face(f) -> child(1) -> vertex(nvf - 1));
          local_residual += child0_width *
                  face_residuals[cell -> face(f) -> child(0) -> index()];
          local_residual += child1_width *
                  face_residuals[cell -> face(f) -> child(1) -> index()];
        } else {
          const double face_width = 
              (cell -> face(f) -> vertex(0)).distance(
               cell -> face(f) -> vertex(nvf - 1));
          local_residual += face_width *
                  face_residuals[cell -> face(f) -> index()];
        }
      } // for ( f )

      local_residual = std::sqrt(local_residual);
    } // for ( cell )
  } // poisson_residual_error_estiamte() [1 / 4]

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::
    poisson_residual_error_estimate(
      const std::vector< unsigned int >&          n_gauss_points,
      const Function<space_dimension>*            rhs_fcn,
      const Functions::ConstantFunction<space_dimension>*    sigma,
      const Function<space_dimension>*            a,
      const std::map< 
                    types::boundary_id,
              const Function<space_dimension>* 
            >&                                    neumann_bc,
      const Vector<double>&                       solution,
            Vector<double>&                       residuals
  ) const {
    TSValues<dimension, space_dimension> ts_values(
      this, n_gauss_points,
      update_values |
      update_gradients |
      update_quadrature_points |
      update_JxW_values | 
      update_hessians
    );
    TSFaceValues<dimension, space_dimension> face_values(
      this, n_gauss_points,
      update_values |
      update_gradients |
      update_quadrature_points |
      update_JxW_values  |
      update_normal_vectors);
    const unsigned int nvc = GeometryInfo<dimension>::vertices_per_cell;
    const unsigned int nvf = GeometryInfo<dimension>::vertices_per_face;

    std::map< unsigned int, double > face_residuals;
    for (const auto& cell : this -> active_cell_iterators()) {
      std::vector< unsigned int > local_dof_indices = get_IEN_array(cell);

      // Compute the jumps of faces at C0 continuity
      for (unsigned int f = 0;
              f < GeometryInfo<dimension>::faces_per_cell;
              f++) {
        // Get orientation of current face:
        const auto& face = cell -> face(f);
        const Point<dimension>& c = (-1.) * 
                                    face -> vertex(0) +
                                    face -> vertex(nvf - 1);
        unsigned int orientation = 0;
        for ( ; orientation < dimension && c(orientation) != 0; orientation++);

        if (face -> has_children() ) {
          // This is merely a special case, as refined faces are considered
          // in-active and hence their index is not represented in mof
          if (mof.at(face -> child(0) -> index()) == p[orientation]) {
            // C0 Edge detected:
            face_values.reinit(cell, f);
            const unsigned int ppf = 2 * face_values.n_quadrature_points_per_face();
            for (const unsigned int q_index : face_values.quadrature_point_indices()){
              double g = 0;
              for (const unsigned int i : face_values.dof_indices())
                g += solution(local_dof_indices[i]) *
                        face_values.shape_grad(i, q_index) *
                        face_values.normal_vector(q_index);

              const unsigned int ch = q_index > ppf ? 1 : 0;
              face_residuals[face -> child(ch) -> index()] +=
                      0.25 * g * g * face_values.JxW(q_index);
            } // for ( q_index )
          } // otherwise it is atleast a C1 edge, and cannot be at the boundary
        } else if (mof.at(face -> index()) == p[orientation]){
          // C0 edge detected
          face_values.reinit(cell, f);

          // Since this face is not refined, we can simply compute
          // the jump terms along it abd store the values at the
          // corresponding place in face_residuals
          for (const unsigned int q_index : face_values.quadrature_point_indices()){
            double g = 0;
            for (const unsigned int i : face_values.dof_indices())
              g += solution(local_dof_indices[i]) *
                      face_values.shape_grad(i, q_index) *
                      face_values.normal_vector(q_index);

            face_residuals[face -> index()] += 0.25 * g * g * face_values.JxW(q_index);
          } // for ( q_index )

        } else if (face -> at_boundary()) {
          const auto& bc = neumann_bc.find(face -> boundary_id());
          if (bc == neumann_bc.end())
            continue;

          face_values.reinit(cell, f);
          // If the face is at the boundary, compute the jump from
          // the face's boundary id
          for (const unsigned int q_index : face_values.quadrature_point_indices()){
            double g = neumann_bc.at(face -> boundary_id())
                          -> value(face_values.quadrature_point(q_index));
            for (const unsigned int i : face_values.dof_indices())
              g -= solution(local_dof_indices[i]) *
                      face_values.shape_grad(i, q_index) *
                      face_values.normal_vector(q_index);

            face_residuals[face -> index()] += g * g * face_values.JxW(q_index);
          } // for ( q_index )
        } // if ( ... )
      } // for ( face )
    } // for ( cell )

    // After every Jump on every face is computed,
    // we compute the residuals on the cells
    const double sigma_val = sigma -> value(Point<space_dimension>());
    unsigned int kk = 0;
    for (const auto& cell : this -> active_cell_iterators()){
      double& local_residual = residuals(kk++);
      const std::vector< unsigned int >& local_dof_indices = get_IEN_array(cell);
      ts_values.reinit(cell);
      for ( const unsigned int q_index : ts_values.quadrature_point_indices()){
        double slu  = 0;
        double au   = 0;

        const Point<space_dimension>& Q = ts_values.quadrature_point(q_index);
        const double a_val     = a -> value(Q);
        for (const unsigned int i : ts_values.dof_indices() ){
          double l = 0;
          for (unsigned int d = 0; d < space_dimension; d++)
            l += ts_values.shape_hessian(i, q_index)[d][d];

          slu   += solution(local_dof_indices[i]) * sigma_val * l;
          au    += solution(local_dof_indices[i]) * a_val;
        } // for ( i )

        double g = (slu - au + rhs_fcn -> value(Q));
        local_residual += g * g * ts_values.JxW(q_index);
      } // for ( q_index )
      const double cell_width = (cell->vertex(0)).distance(
                                cell->vertex(nvc - 1));

      local_residual *= cell_width * cell_width;

      // And then also add the values from the faces:
      for (unsigned int f = 0;
              f < GeometryInfo<dimension>::faces_per_cell;
              f++) {
        if (cell -> face(f) -> has_children()){
          const double child0_width = 
             (cell -> face(f) -> child(0) -> vertex(0)).distance(
              cell -> face(f) -> child(0) -> vertex(nvf - 1));
          const double child1_width = 
             (cell -> face(f) -> child(1) -> vertex(0)).distance(
              cell -> face(f) -> child(1) -> vertex(nvf - 1));
          local_residual += child0_width *
                  face_residuals[cell -> face(f) -> child(0) -> index()];
          local_residual += child1_width *
                  face_residuals[cell -> face(f) -> child(1) -> index()];
        } else {
          const double face_width = 
              (cell -> face(f) -> vertex(0)).distance(
               cell -> face(f) -> vertex(nvf - 1));
          local_residual += face_width *
                  face_residuals[cell -> face(f) -> index()];
        }
      } // for ( f )

      local_residual = std::sqrt(local_residual);
    } // for ( cell )
  } // poisson_residual_error_estiamte() [1 / 4]

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::
    poisson_residual_error_estimate(
      const std::vector< unsigned int >&          n_gauss_points,
      const Function<space_dimension>*            rhs_fcn,
      const Function<space_dimension>*            sigma,
      const Functions::ConstantFunction<space_dimension>*    a,
      const std::map< 
                    types::boundary_id,
              const Function<space_dimension>* 
            >&                                    neumann_bc,
      const Vector<double>&                       solution,
            Vector<double>&                       residuals
  ) const {
    TSValues<dimension, space_dimension> ts_values(
      this, n_gauss_points,
      update_values |
      update_gradients |
      update_quadrature_points |
      update_JxW_values | 
      update_hessians
    );
    TSFaceValues<dimension, space_dimension> face_values(
      this, n_gauss_points,
      update_values |
      update_gradients |
      update_quadrature_points |
      update_JxW_values  |
      update_normal_vectors);
    const unsigned int nvc = GeometryInfo<dimension>::vertices_per_cell;
    const unsigned int nvf = GeometryInfo<dimension>::vertices_per_face;

    std::map< unsigned int, double > face_residuals;
    for (const auto& cell : this -> active_cell_iterators()) {
      std::vector< unsigned int > local_dof_indices = get_IEN_array(cell);

      // Compute the jumps of faces at C0 continuity
      for (unsigned int f = 0;
              f < GeometryInfo<dimension>::faces_per_cell;
              f++) {
        // Get orientation of current face:
        const auto& face = cell -> face(f);
        const Point<dimension>& c = (-1.) * 
                                    face -> vertex(0) +
                                    face -> vertex(nvf - 1);
        unsigned int orientation = 0;
        for ( ; orientation < dimension && c(orientation) != 0; orientation++);

        if (face -> has_children() ) {
          // This is merely a special case, as refined faces are considered
          // in-active and hence their index is not represented in mof
          if (mof.at(face -> child(0) -> index()) == p[orientation]) {
            // C0 Edge detected:
            face_values.reinit(cell, f);
            const unsigned int ppf = 2 * face_values.n_quadrature_points_per_face();
            for (const unsigned int q_index : face_values.quadrature_point_indices()){
              double g = 0;
              for (const unsigned int i : face_values.dof_indices())
                g += solution(local_dof_indices[i]) *
                        face_values.shape_grad(i, q_index) *
                        face_values.normal_vector(q_index);

              const unsigned int ch = q_index > ppf ? 1 : 0;
              face_residuals[face -> child(ch) -> index()] +=
                      0.25 * g * g * face_values.JxW(q_index);
            } // for ( q_index )
          } // otherwise it is atleast a C1 edge, and cannot be at the boundary
        } else if (mof.at(face -> index()) == p[orientation]){
          // C0 edge detected
          face_values.reinit(cell, f);

          // Since this face is not refined, we can simply compute
          // the jump terms along it abd store the values at the
          // corresponding place in face_residuals
          for (const unsigned int q_index : face_values.quadrature_point_indices()){
            double g = 0;
            for (const unsigned int i : face_values.dof_indices())
              g += solution(local_dof_indices[i]) *
                      face_values.shape_grad(i, q_index) *
                      face_values.normal_vector(q_index);

            face_residuals[face -> index()] += 0.25 * g * g * face_values.JxW(q_index);
          } // for ( q_index )

        } else if (face -> at_boundary()) {
          const auto& bc = neumann_bc.find(face -> boundary_id());
          if (bc == neumann_bc.end())
            continue;

          face_values.reinit(cell, f);
          // If the face is at the boundary, compute the jump from
          // the face's boundary id
          for (const unsigned int q_index : face_values.quadrature_point_indices()){
            double g = neumann_bc.at(face -> boundary_id())
                          -> value(face_values.quadrature_point(q_index));
            for (const unsigned int i : face_values.dof_indices())
              g -= solution(local_dof_indices[i]) *
                      face_values.shape_grad(i, q_index) *
                      face_values.normal_vector(q_index);

            face_residuals[face -> index()] += g * g * face_values.JxW(q_index);
          } // for ( q_index )
        } // if ( ... )
      } // for ( face )
    } // for ( cell )

    // After every Jump on every face is computed,
    // we compute the residuals on the cells
    const double a_val     = a -> value(Point<space_dimension>());
    unsigned int kk = 0;
    for (const auto& cell : this -> active_cell_iterators()){
      double& local_residual = residuals(kk++);
      const std::vector< unsigned int >& local_dof_indices = get_IEN_array(cell);
      ts_values.reinit(cell);
      for ( const unsigned int q_index : ts_values.quadrature_point_indices()){
        double slu  = 0;
        double au   = 0;
        double gsgu = 0;

        const Point<space_dimension>& Q = ts_values.quadrature_point(q_index);
        const double sigma_val = sigma -> value(Q);
        for (const unsigned int i : ts_values.dof_indices() ){
          double l = 0;
          for (unsigned int d = 0; d < space_dimension; d++)
            l += ts_values.shape_hessian(i, q_index)[d][d];

          const Tensor<1, space_dimension>& sigma_grad = sigma -> gradient(Q);
          slu   += solution(local_dof_indices[i]) * sigma_val * l;
          gsgu  += solution(local_dof_indices[i]) *
                      sigma_grad *
                      ts_values.shape_grad(i, q_index);
          au    += solution(local_dof_indices[i]) * a_val;
        } // for ( i )

        double g = (slu + gsgu - au + rhs_fcn -> value(Q));
        local_residual += g * g * ts_values.JxW(q_index);
      } // for ( q_index )
      const double cell_width = (cell->vertex(0)).distance(
                                cell->vertex(nvc - 1));

      local_residual *= cell_width * cell_width;

      // And then also add the values from the faces:
      for (unsigned int f = 0;
              f < GeometryInfo<dimension>::faces_per_cell;
              f++) {
        if (cell -> face(f) -> has_children()){
          const double child0_width = 
             (cell -> face(f) -> child(0) -> vertex(0)).distance(
              cell -> face(f) -> child(0) -> vertex(nvf - 1));
          const double child1_width = 
             (cell -> face(f) -> child(1) -> vertex(0)).distance(
              cell -> face(f) -> child(1) -> vertex(nvf - 1));
          local_residual += child0_width *
                  face_residuals[cell -> face(f) -> child(0) -> index()];
          local_residual += child1_width *
                  face_residuals[cell -> face(f) -> child(1) -> index()];
        } else {
          const double face_width = 
              (cell -> face(f) -> vertex(0)).distance(
               cell -> face(f) -> vertex(nvf - 1));
          local_residual += face_width *
                  face_residuals[cell -> face(f) -> index()];
        }
      } // for ( f )

      local_residual = std::sqrt(local_residual);
    } // for ( cell )
  } // poisson_residual_error_estiamte() [1 / 4]

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::
    poisson_residual_error_estimate(
      const std::vector< unsigned int >&          n_gauss_points,
      const Function<space_dimension>*            rhs_fcn,
      const Functions::ConstantFunction<space_dimension>*    sigma,
      const Functions::ConstantFunction<space_dimension>*    a,
      const std::map< 
                    types::boundary_id,
              const Function<space_dimension>* 
            >&                                    neumann_bc,
      const Vector<double>&                       solution,
            Vector<double>&                       residuals
  ) const {
    TSValues<dimension, space_dimension> ts_values(
      this, n_gauss_points,
      update_values |
      update_gradients |
      update_quadrature_points |
      update_JxW_values | 
      update_hessians
    );
    TSFaceValues<dimension, space_dimension> face_values(
      this, n_gauss_points,
      update_values |
      update_gradients |
      update_quadrature_points |
      update_JxW_values  |
      update_normal_vectors);
    const unsigned int nvc = GeometryInfo<dimension>::vertices_per_cell;
    const unsigned int nvf = GeometryInfo<dimension>::vertices_per_face;

    std::map< unsigned int, double > face_residuals;
    for (const auto& cell : this -> active_cell_iterators()) {
      std::vector< unsigned int > local_dof_indices = get_IEN_array(cell);

      // Compute the jumps of faces at C0 continuity
      for (unsigned int f = 0;
              f < GeometryInfo<dimension>::faces_per_cell;
              f++) {
        // Get orientation of current face:
        const auto& face = cell -> face(f);
        const Point<dimension>& c = (-1.) * 
                                    face -> vertex(0) +
                                    face -> vertex(nvf - 1);
        unsigned int orientation = 0;
        for ( ; orientation < dimension && c(orientation) != 0; orientation++);

        if (face -> has_children() ) {
          // This is merely a special case, as refined faces are considered
          // in-active and hence their index is not represented in mof
          if (mof.at(face -> child(0) -> index()) == p[orientation]) {
            // C0 Edge detected:
            face_values.reinit(cell, f);
            const unsigned int ppf = 2 * face_values.n_quadrature_points_per_face();
            for (const unsigned int q_index : face_values.quadrature_point_indices()){
              double g = 0;
              for (const unsigned int i : face_values.dof_indices())
                g += solution(local_dof_indices[i]) *
                        face_values.shape_grad(i, q_index) *
                        face_values.normal_vector(q_index);

              const unsigned int ch = q_index > ppf ? 1 : 0;
              face_residuals[face -> child(ch) -> index()] +=
                      0.25 * g * g * face_values.JxW(q_index);
            } // for ( q_index )
          } // otherwise it is atleast a C1 edge, and cannot be at the boundary
        } else if (mof.at(face -> index()) == p[orientation]){
          // C0 edge detected
          face_values.reinit(cell, f);

          // Since this face is not refined, we can simply compute
          // the jump terms along it abd store the values at the
          // corresponding place in face_residuals
          for (const unsigned int q_index : face_values.quadrature_point_indices()){
            double g = 0;
            for (const unsigned int i : face_values.dof_indices())
              g += solution(local_dof_indices[i]) *
                      face_values.shape_grad(i, q_index) *
                      face_values.normal_vector(q_index);

            face_residuals[face -> index()] += 0.25 * g * g * face_values.JxW(q_index);
          } // for ( q_index )

        } else if (face -> at_boundary()) {
          const auto& bc = neumann_bc.find(face -> boundary_id());
          if (bc == neumann_bc.end())
            continue;

          face_values.reinit(cell, f);
          // If the face is at the boundary, compute the jump from
          // the face's boundary id
          for (const unsigned int q_index : face_values.quadrature_point_indices()){
            double g = neumann_bc.at(face -> boundary_id())
                          -> value(face_values.quadrature_point(q_index));
            for (const unsigned int i : face_values.dof_indices())
              g -= solution(local_dof_indices[i]) *
                      face_values.shape_grad(i, q_index) *
                      face_values.normal_vector(q_index);

            face_residuals[face -> index()] += g * g * face_values.JxW(q_index);
          } // for ( q_index )
        } // if ( ... )
      } // for ( face )
    } // for ( cell )

    // After every Jump on every face is computed,
    // we compute the residuals on the cells
    const double a_val     = a -> value(Point<space_dimension>());
    const double sigma_val = sigma -> value(Point<space_dimension>());
    unsigned int kk = 0;
    for (const auto& cell : this -> active_cell_iterators()){
      double& local_residual = residuals(kk++);
      const std::vector< unsigned int >& local_dof_indices = get_IEN_array(cell);
      ts_values.reinit(cell);
      for ( const unsigned int q_index : ts_values.quadrature_point_indices()){
        double slu  = 0;
        double au   = 0;

        const Point<space_dimension>& Q = ts_values.quadrature_point(q_index);
        for (const unsigned int i : ts_values.dof_indices() ){
          double l = 0;
          for (unsigned int d = 0; d < space_dimension; d++)
            l += ts_values.shape_hessian(i, q_index)[d][d];

          slu   += solution(local_dof_indices[i]) * sigma_val * l;
          au    += solution(local_dof_indices[i]) * a_val;
        } // for ( i )

        double g = (slu - au + rhs_fcn -> value(Q));
        local_residual += g * g * ts_values.JxW(q_index);
      } // for ( q_index )
      const double cell_width = (cell->vertex(0)).distance(
                                cell->vertex(nvc - 1));

      local_residual *= cell_width * cell_width;

      // And then also add the values from the faces:
      for (unsigned int f = 0;
              f < GeometryInfo<dimension>::faces_per_cell;
              f++) {
        if (cell -> face(f) -> has_children()){
          const double child0_width = 
             (cell -> face(f) -> child(0) -> vertex(0)).distance(
              cell -> face(f) -> child(0) -> vertex(nvf - 1));
          const double child1_width = 
             (cell -> face(f) -> child(1) -> vertex(0)).distance(
              cell -> face(f) -> child(1) -> vertex(nvf - 1));
          local_residual += child0_width *
                  face_residuals[cell -> face(f) -> child(0) -> index()];
          local_residual += child1_width *
                  face_residuals[cell -> face(f) -> child(1) -> index()];
        } else {
          const double face_width = 
              (cell -> face(f) -> vertex(0)).distance(
               cell -> face(f) -> vertex(nvf - 1));
          local_residual += face_width *
                  face_residuals[cell -> face(f) -> index()];
        }
      } // for ( f )

      local_residual = std::sqrt(local_residual);
    } // for ( cell )
  } // poisson_residual_error_estiamte() [1 / 4]

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::
    poisson_residual_error_estimate(
      const std::vector< unsigned int >&          n_gauss_points,
      const Function<space_dimension>*            rhs_fcn,
      const Function<space_dimension>*            sigma,
      const std::map< 
                    types::boundary_id,
                    Function<space_dimension>* 
            >&                                    neumann_bc,
      const Vector<double>&                       solution,
      std::map< cell_iterator,
                double >&                         residuals
  ) const {
    std::map<
            types::boundary_id,
      const Function<space_dimension>*
      > const_bc;
    for (const auto& [id, fcn] : neumann_bc)
      const_bc[id] = fcn;
    poisson_residual_error_estimate(n_gauss_points, rhs_fcn, sigma, const_bc, solution, residuals);
  } // poisson_residual_error_estimate() [1.5 / 4]

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::
    poisson_residual_error_estimate(
      const std::vector< unsigned int >&         n_gauss_points,
      const Function<space_dimension>*           rhs_fcn,
      const Function<space_dimension>*           sigma,
      const Vector<double>&                      solution,
      std::map< cell_iterator,
                double >&                        residuals
  ) const {
    TSValues<dimension, space_dimension> ts_values(
      this, n_gauss_points,
      update_values |
      update_gradients |
      update_quadrature_points |
      update_JxW_values | 
      update_hessians
    );
    const unsigned int nvc = GeometryInfo<dimension>::vertices_per_cell;

    for (const auto& cell : this -> active_cell_iterators()){
      double& local_residual = residuals[cell];
      const std::vector< unsigned int >& local_dof_indices = get_IEN_array(cell);
      ts_values.reinit(cell);
      for ( const unsigned int q_index : ts_values.quadrature_point_indices()){
        double slu = 0;
        double gsgu = 0;

        const Point<space_dimension>& Q = ts_values.quadrature_point(q_index);
        const Tensor<1, space_dimension>& sigma_grad = sigma -> gradient(Q);
        const double sigma_val = sigma -> value(Q);
        for (const unsigned int i : ts_values.dof_indices() ){
          double l = 0;
          for (unsigned int d = 0; d < space_dimension; d++)
            l += ts_values.shape_hessian(i, q_index)[d][d];

          slu   += solution(local_dof_indices[i]) * sigma_val * l;
          gsgu  += solution(local_dof_indices[i]) *
                      sigma_grad *
                      ts_values.shape_grad(i, q_index);
        } // for ( i )

        double g = (slu + gsgu + rhs_fcn -> value(Q));
        local_residual += g * g * ts_values.JxW(q_index);
      } // for ( q_index )
      const double cell_width = (cell->vertex(0)).distance(
                              cell->vertex(nvc - 1));

      local_residual *= cell_width * cell_width;
      local_residual = std::sqrt(local_residual);
    } // for ( cell )
  } // poisson_residual_error_estiamte() [2 / 4]

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::
    poisson_residual_error_estimate(
      const std::vector< unsigned int >&         n_gauss_points,
      const Function<space_dimension>*           rhs_fcn,
      const std::map< 
                    types::boundary_id,
              const Function<space_dimension>* 
            >&                                   neumann_bc,
      const Vector<double>&                      solution,
      std::map< cell_iterator,
                double >&                        residuals
  ) const { // sigma = 1
    TSValues<dimension, space_dimension> ts_values(
      this, n_gauss_points,
      update_values |
      update_gradients |
      update_quadrature_points |
      update_JxW_values | 
      update_hessians
    );
    TSFaceValues<dimension, space_dimension> face_values(
      this, n_gauss_points,
      update_values |
      update_gradients |
      update_quadrature_points |
      update_JxW_values | 
      update_normal_vectors );
    const unsigned int nvc = GeometryInfo<dimension>::vertices_per_cell;
    const unsigned int nvf = GeometryInfo<dimension>::vertices_per_face;

    std::map< unsigned int, double > face_residuals;
    for (const auto& cell : this -> active_cell_iterators()) {
      std::vector< unsigned int > local_dof_indices = get_IEN_array(cell);

      // Compute the jumps of faces at C0 continuity
      for (unsigned int f = 0;
              f < GeometryInfo<dimension>::faces_per_cell;
              f++) {
        // Get orientation of current face:
        const auto& face = cell -> face(f);
        const Point<dimension>& c = (-1.) * 
                                    face -> vertex(0) +
                                    face -> vertex(nvf - 1);
        unsigned int orientation = 0;
        for ( ; orientation < dimension && c(orientation) != 0; orientation++);

#ifdef DEBUG
        if (!face -> has_children()){
          if (mof.find(face->index()) == mof.end()) {
            std::cout << "Out of range error occured while accessing elements of mof" << std::endl;
            std::cout << "Size of mof: " << mof.size() << std::endl; 
            std::cout << "Number of active faces: " << this->n_active_faces() << std::endl;
            std::cout << "Couldn't access face with index " << face -> index() << std::endl;
            std::cout << "Indices stored are: " << std::endl;
            for (const auto& [ind, mult] : mof)
              std::cout << "[" << ind << ", " << mult << "]" << std::endl;
            std::cout << "face has bounds (" << face->vertex(0) << ") x ("
                      << face->vertex(GeometryInfo<dimension>::vertices_per_face-1) << ")"
                      << std::endl;
            std::cout << std::boolalpha << "has_children():" << face -> has_children() << std::endl;
            std::cout << "Cells with faces: " << std::endl;
            for (const auto& cell : this->active_cell_iterators()){
              std::cout << "Cell: " << cell -> level() 
                        << "." 
                        << cell -> index() 
                        << " faces: ";
              for (const auto& face : cell -> face_iterators())
                std::cout << face->index() << " ";
              std::cout << std::endl;
            } // for cell

            if (dimension == 2) {
              GridOutFlags::Svg svg_flags;
              svg_flags.coloring = GridOutFlags::Svg::Coloring::none;
              svg_flags.label_level_number  = true;
              svg_flags.label_cell_index    = true;

              GridOut       grid_out;
              grid_out.set_flags(svg_flags);
  
              std::string name = "log/mof_error_grid.svg";
              try {
                std::ofstream out(name);
                grid_out.write_svg(*this, out);
              } catch (const dealii::StandardExceptions::ExcIO& e){
                name = "mof_error_grid.svg"; 
                std::ofstream out(name);
                grid_out.write_svg(*this, out);
              }
              std::cout << "Grid printed to " << name << std::endl;
            }
            throw ExcInternalError();
          }
        }
#endif
        if (face -> has_children() ) {
          // This is merely a special case, as refined faces are considered
          // in-active and hence their index is not represented in mof
          if (mof.at(face -> child(0) -> index()) == p.at(orientation)) {
            // C0 Edge detected:
            face_values.reinit(cell, f);
            const unsigned int ppf = 2 * face_values.n_quadrature_points_per_face();
            for (const unsigned int q_index : face_values.quadrature_point_indices()){
              double g = 0;
              for (const unsigned int i : face_values.dof_indices())
                g += solution(local_dof_indices[i]) *
                        face_values.shape_grad(i, q_index) *
                        face_values.normal_vector(q_index);

              const unsigned int ch = q_index > ppf ? 1 : 0;
              face_residuals[face -> child(ch) -> index()] +=
                      0.25 * g * g * face_values.JxW(q_index);
            } // for ( q_index )
          } // otherwise it is atleast a C1 edge, and cannot be at the boundary
        } else if (mof.at(face -> index()) == p.at(orientation)){
          // C0 edge detected
          face_values.reinit(cell, f);

          // Since this face is not refined, we can simply compute
          // the jump terms along it abd store the values at the
          // corresponding place in face_residuals
          for (const unsigned int q_index : face_values.quadrature_point_indices()){
            double g = 0;
            for (const unsigned int i : face_values.dof_indices())
              g += solution(local_dof_indices[i]) *
                      face_values.shape_grad(i, q_index) *
                      face_values.normal_vector(q_index);

            face_residuals[face -> index()] += 0.25 * g * g * face_values.JxW(q_index);
          } // for ( q_index )

        } else if (face -> at_boundary()) {
          const auto& bc = neumann_bc.find(face -> boundary_id());
          if (bc == neumann_bc.end())
            continue;

          face_values.reinit(cell, f);
          // If the face is at the boundary, compute the jump from
          // the face's boundary id
          for (const unsigned int q_index : face_values.quadrature_point_indices()){
            double g = neumann_bc.at(face -> boundary_id())
                          -> value(face_values.quadrature_point(q_index));
            for (const unsigned int i : face_values.dof_indices()){
              g -= solution(local_dof_indices[i]) *
                      face_values.shape_grad(i, q_index) *
                      face_values.normal_vector(q_index);
            }

            face_residuals[face -> index()] += g * g * face_values.JxW(q_index);
          } // for ( q_index )
        } // if ( ... )
      } // for ( face )
    } // for ( cell )

    // After every Jump on every face is computed,
    // we compute the residuals on the cells
    for (const auto& cell : this -> active_cell_iterators()){
      double& local_residual = residuals[cell];
      const std::vector< unsigned int >& local_dof_indices = get_IEN_array(cell);
      ts_values.reinit(cell);
      for ( const unsigned int q_index : ts_values.quadrature_point_indices()){
        const Point<space_dimension>& Q = ts_values.quadrature_point(q_index);
        double g = rhs_fcn -> value(Q);
        for (const unsigned int i : ts_values.dof_indices() ){
          double l = 0;
          for (unsigned int d = 0; d < space_dimension; d++)
            l += ts_values.shape_hessian(i, q_index)[d][d];

          g += solution(local_dof_indices[i]) * l;
        } // for ( i )

        local_residual += g * g * ts_values.JxW(q_index);
      } // for ( q_index )
      const double cell_width = (cell->vertex(0)).distance(
                              cell->vertex(nvc - 1));

      local_residual *= cell_width * cell_width;

      // And then also add the values from the faces:
      for (unsigned int f = 0;
              f < GeometryInfo<dimension>::faces_per_cell;
              f++) {
        if (cell -> face(f) -> has_children()){
          const double child0_width = 
             (cell -> face(f) -> child(0) -> vertex(0)).distance(
              cell -> face(f) -> child(0) -> vertex(nvf - 1));
          const double child1_width = 
             (cell -> face(f) -> child(1) -> vertex(0)).distance(
              cell -> face(f) -> child(1) -> vertex(nvf - 1));
          local_residual += child0_width *
                  face_residuals[cell -> face(f) -> child(0) -> index()];
          local_residual += child1_width *
                  face_residuals[cell -> face(f) -> child(1) -> index()];
        } else {
          const double face_width = 
              (cell -> face(f) -> vertex(0)).distance(
               cell -> face(f) -> vertex(nvf - 1));
          local_residual += face_width *
                  face_residuals[cell -> face(f) -> index()];
        }
      } // for ( f )

      local_residual = std::sqrt(local_residual);
    } // for ( cell )
  } // poisson_residual_error_estimate() [ 3 / 4]

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::
    poisson_residual_error_estimate(
      const std::vector< unsigned int >&         n_gauss_points,
      const Function<space_dimension>*           rhs_fcn,
      const std::map< 
                    types::boundary_id,
                    Function<space_dimension>* 
            >&                                   neumann_bc,
      const Vector<double>&                      solution,
      std::map< cell_iterator,
                double >&                        residuals
  ) const { // sigma = 1
    std::map<
            types::boundary_id,
      const Function<space_dimension>*
      > const_bc;
    for (const auto& [id, fcn] : neumann_bc)
      const_bc[id] = fcn;
    poisson_residual_error_estimate(n_gauss_points, rhs_fcn, const_bc, solution, residuals);
  }

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::
    poisson_residual_error_estimate(
      const std::vector< unsigned int >&         n_gauss_points,
      const Function<space_dimension>*           rhs_fcn,
      const Vector<double>&                      solution,
      std::map< cell_iterator,
                double >&                        residuals
  ) const { // sigma = 1, g = 0
    TSValues<dimension, space_dimension> ts_values(
      this, n_gauss_points,
      update_values |
      update_gradients |
      update_quadrature_points |
      update_JxW_values | 
      update_hessians
    );
    const unsigned int nvc = GeometryInfo<dimension>::vertices_per_cell;

    // After every Jump on every face is computed,
    // we compute the residuals on the cells
    for (const auto& cell : this -> active_cell_iterators()){
      double& local_residual = residuals[cell];
      const std::vector< unsigned int >& local_dof_indices = get_IEN_array(cell);
      ts_values.reinit(cell);
      for ( const unsigned int q_index : ts_values.quadrature_point_indices()){
        const Point<space_dimension>& Q = ts_values.quadrature_point(q_index);
        double g = rhs_fcn -> value(Q);
        for (const unsigned int i : ts_values.dof_indices() ){
          double l = 0;
          for (unsigned int d = 0; d < space_dimension; d++)
            l += ts_values.shape_hessian(i, q_index)[d][d];

          g += solution(local_dof_indices[i]) * l;
        } // for ( i )

        local_residual += g * g * ts_values.JxW(q_index);
      } // for ( q_index )
      const double cell_width = (cell->vertex(0)).distance(
                                  cell->vertex(nvc - 1));

      local_residual = std::sqrt(local_residual * cell_width);
    } // for ( cell )
  } // poisson_residual_error_estimate() [ 4 / 4]


  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::
    poisson_residual_error_estimate(
      const std::vector< unsigned int >&          n_gauss_points,
      const Function<space_dimension>*            rhs_fcn,
      const Function<space_dimension>*            sigma,
      const std::map< 
                    types::boundary_id,
              const Function<space_dimension>* 
            >&                                    neumann_bc,
      const Vector<double>&                       solution,
            Vector<double>&                       residuals
  ) const {
    Assert(residuals.size() == this->n_active_cells(),
            ExcDimensionMismatch(residuals.size(), this->n_active_cells()));
    // Pass arguments to previous function
    std::map< cell_iterator, double > res;
    poisson_residual_error_estimate(
      n_gauss_points, 
      rhs_fcn,
      sigma,
      neumann_bc, 
      solution,
      res
    );

    unsigned int i = 0; 
    for (const auto& [_, e] : res)
      residuals(i++) = e; 
  }

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::
    poisson_residual_error_estimate(
      const std::vector< unsigned int >&          n_gauss_points,
      const Function<space_dimension>*            rhs_fcn,
      const Function<space_dimension>*            sigma,
      const std::map< 
                    types::boundary_id,
                    Function<space_dimension>* 
            >&                                    neumann_bc,
      const Vector<double>&                       solution,
            Vector<double>&                       residuals
  ) const {
    Assert(residuals.size() == this->n_active_cells(),
            ExcDimensionMismatch(residuals.size(), this->n_active_cells()));
    // Pass arguments to previous function
    std::map< cell_iterator, double > res;
    poisson_residual_error_estimate(
      n_gauss_points, 
      rhs_fcn,
      sigma,
      neumann_bc, 
      solution,
      res
    );

    unsigned int i = 0; 
    for (const auto& [_, e] : res)
      residuals(i++) = e; 
  }

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::
    poisson_residual_error_estimate(
      const std::vector< unsigned int >&          n_gauss_points,
      const Function<space_dimension>*            rhs_fcn,
      const Function<space_dimension>*            sigma,
      const Vector<double>&                       solution,
            Vector<double>&                       residuals
  ) const {
    Assert(residuals.size() == this->n_active_cells(),
            ExcDimensionMismatch(residuals.size(), this->n_active_cells()));
    // Pass arguments to previous function
    std::map< cell_iterator, double > res;
    poisson_residual_error_estimate(
      n_gauss_points, 
      rhs_fcn,
      sigma,
      solution,
      res
    );

    unsigned int i = 0; 
    for (const auto& [_, e] : res)
      residuals(i++) = e; 
  }

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::
    poisson_residual_error_estimate(
      const std::vector< unsigned int >&          n_gauss_points,
      const Function<space_dimension>*            rhs_fcn,
      const std::map< 
                    types::boundary_id,
              const Function<space_dimension>* 
            >&                                   neumann_bc,
      const Vector<double>&                       solution,
            Vector<double>&                           residuals
  ) const {
    Assert(residuals.size() == this->n_active_cells(),
            ExcDimensionMismatch(residuals.size(), this->n_active_cells()));
    // Pass arguments to previous function
    std::map< cell_iterator, double > res;
    poisson_residual_error_estimate(
      n_gauss_points, 
      rhs_fcn,
      neumann_bc,
      solution,
      res
    );

    unsigned int i = 0; 
    for (const auto& [_, e] : res)
      residuals(i++) = e; 
  }

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::
    poisson_residual_error_estimate(
      const std::vector< unsigned int >&          n_gauss_points,
      const Function<space_dimension>*            rhs_fcn,
      const std::map< 
                    types::boundary_id,
                    Function<space_dimension>* 
            >&                                   neumann_bc,
      const Vector<double>&                       solution,
            Vector<double>&                       residuals
  ) const {
    Assert(residuals.size() == this->n_active_cells(),
            ExcDimensionMismatch(residuals.size(), this->n_active_cells()));
    // Pass arguments to previous function
    std::map< cell_iterator, double > res;
    poisson_residual_error_estimate(
      n_gauss_points, 
      rhs_fcn,
      neumann_bc,
      solution,
      res
    );

    unsigned int i = 0; 
    for (const auto& [_, e] : res)
      residuals(i++) = e; 
  }

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::
    poisson_residual_error_estimate(
      const std::vector< unsigned int >&          n_gauss_points,
      const Function<space_dimension>*            rhs_fcn,
      const Vector<double>&                       solution,
            Vector<double>&                       residuals
  ) const {
    Assert(residuals.size() == this->n_active_cells(),
            ExcDimensionMismatch(residuals.size(), this->n_active_cells()));
    // Pass arguments to previous function
    std::map< cell_iterator, double > res;
    poisson_residual_error_estimate(
      n_gauss_points, 
      rhs_fcn,
      solution,
      res
    );

    unsigned int i = 0; 
    for (const auto& [_, e] : res)
      residuals(i++) = e; 
  }

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::
    linear_elasticity_residual_error_estimate(
      const std::vector< unsigned int >&          n_gauss_points,
      const Function<space_dimension>*            rhs_fcn,
      const Function<space_dimension>*            lambda_ptr,
      const Function<space_dimension>*            mu_ptr,
      const std::map< 
                    types::boundary_id,
              const Function<space_dimension>* 
            >&                                    neumann_bc,
      const Vector<double>&                       solution,
      std::map< cell_iterator,
                double >&                         residuals
  ) const {
    Assert(rhs_fcn -> n_components == space_dimension, 
            ExcDimensionMismatch(rhs_fcn -> n_components, space_dimension));
    Assert(lambda_ptr -> n_components == 1, 
            ExcDimensionMismatch(lambda_ptr -> n_components, space_dimension));
    Assert(mu_ptr -> n_components == 1, 
            ExcDimensionMismatch(mu_ptr -> n_components, space_dimension));
#ifdef DEBUG
    for (auto bc_fcn_it = neumann_bc.begin(); bc_fcn_it != neumann_bc.end(); ++bc_fcn_it){
      Assert(bc_fcn_it -> second -> n_components == rhs_fcn -> n_components, 
                ExcDimensionMismatch(bc_fcn_it -> second -> n_components, rhs_fcn -> n_components));
    }
#endif
    const unsigned int &nfc = GeometryInfo<dimension>::faces_per_cell;
    const unsigned int &nvf = GeometryInfo<dimension>::vertices_per_face;
    const unsigned int &n_components = space_dimension;
                
    TSValues<dimension, space_dimension, space_dimension> ts_values(
      this, n_gauss_points,
      update_values |
      update_gradients |
      update_quadrature_points |
      update_JxW_values | 
      update_hessians
    );
    TSFaceValues<dimension, space_dimension, space_dimension> face_values(
      this, n_gauss_points,
      update_values |
      update_gradients |
      update_quadrature_points |
      update_JxW_values  |
      update_normal_vectors);
    for (const auto& cell : this -> active_cell_iterators()) {
      double& cell_residual = residuals[cell];
      cell_residual = 0.;
      const std::vector< unsigned int > local_IEN_array
        = get_IEN_array(cell, n_components);
      ts_values.reinit(cell); 
      
      Tensor<1, space_dimension>    div_c_u;
      for (const unsigned int q : ts_values.quadrature_point_indices()){
        Tensor<1, space_dimension> 
              tmp1, // \nabla \lambda (\nabla \cdot u)
              tmp2, // (\nabla \cdot mu + \mu(\nabla \cdot \nabla)) u
              tmp3; // \nabla \cdot \mu (nabla u)^T
        const Point<space_dimension>  Q = ts_values.quadrature_point(q);
        const double                  lambda_q = lambda_ptr -> value(Q);
        const double                  mu_q = mu_ptr -> value(Q);

        // Initialize tmp1 with the first term
        //    (\nabla \lambda) (\nabla \cdot u)
        {
          // Compute \nabla \cdot u
          tmp1 = lambda_ptr -> gradient(Q);
          double nabla_dot_u = 0;
          for (const unsigned int i : ts_values.dof_indices()) {
            const double        &c_i = solution(local_IEN_array[i]);
            const unsigned int  &comp_i = ts_values.system_to_component_index(i).first;
            nabla_dot_u += c_i * ts_values.shape_grad(i, q)[comp_i];
          }

          // compute (\nabla \lambda) (\nabla \cdot u)
          tmp1 *= nabla_dot_u;
        } 

        {
          const Tensor<1, space_dimension>& mu_grad = mu_ptr -> gradient(Q);
                double                      nabla_dot_mu = 0.;
          for (unsigned int d = 0; d < space_dimension; d++)
            nabla_dot_mu += mu_grad[d];

          // Initialize tmp2 with the first term 
          //      (\nabla \cdot \mu) u
          // Initialize tmp3 with the first term
          //      (\nabla u)^T \nabla \mu
          Tensor<1, space_dimension>  shape_u;
          for (const unsigned int i : ts_values.dof_indices()){
            const double        &c_i = solution(local_IEN_array[i]);
            const unsigned int  &comp_i = ts_values.system_to_component_index(i).first;
            shape_u[comp_i] += c_i * ts_values.shape_value(i, q);

            const DerivativeForm<1, space_dimension, space_dimension> hess = 
                    DerivativeForm<1, space_dimension, space_dimension>(
                        ts_values.shape_hessian(i, q)
                    ).transpose();
            for (unsigned int d = 0; d < space_dimension; d++)
              tmp3[d] += c_i * hess[d] * mu_grad;
          }

          tmp2 = nabla_dot_mu * shape_u;
        } 

        for (const unsigned int i : ts_values.dof_indices()) {
          const double c_i = solution(local_IEN_array[i]);
          const Tensor<2, space_dimension>& hess = 
            ts_values.shape_hessian(i, q);
          const unsigned int comp_i = 
            ts_values.system_to_component_index(i).first;

          // \lambda \nabla (\nabla \cdot u)
          tmp1 += c_i * lambda_q * hess[comp_i];

          // \mu (\nabla \cdot \nabla u)
          for (unsigned int d = 0; d < space_dimension; d++) 
            tmp2[comp_i] += mu_q * hess[d][d];

          // \mu \nabla \cdot (\nabla u)^T
          tmp3 += c_i * mu_q * hess[comp_i];
        } // for ( i )

        Tensor<1, space_dimension> integrand =
                tmp1 + tmp2 + tmp3;
        for (unsigned int d = 0; d < space_dimension; d++)
          integrand[d] += rhs_fcn -> value(Q, d);

        cell_residual += integrand * integrand * ts_values.JxW(q);
      } // for ( q )
      cell_residual *= cell -> diameter() * cell -> diameter();

      for (unsigned int f = 0; f < nfc; f++){
        double face_residuals = 0.;
        const auto& face = cell -> face(f);
        if (face -> at_boundary()){
          const auto& bdry_fcn_ptr = neumann_bc.find(face -> boundary_id()); 
          
          // Is the boundary condition present for neumann? 
          if (bdry_fcn_ptr == neumann_bc.end())
            continue; // no, continue

          const Function<space_dimension>* bdry_fcn = bdry_fcn_ptr -> second;
          Assert(bdry_fcn -> n_components == space_dimension , 
                  ExcDimensionMismatch(bdry_fcn -> n_components, space_dimension));
          face_values.reinit(cell, f); 
          for (const unsigned int q : face_values.quadrature_point_indices()){
            const Point<space_dimension>     &Q      = face_values.quadrature_point(q);
                  Tensor<2, space_dimension>  function_jacobian;
                  SymmetricTensor<2, space_dimension> stress;
                  Tensor<2, space_dimension>  strain;
                  Tensor<1, space_dimension>  normal_dot_strain;
            const double                lambda_q = lambda_ptr -> value(Q);
            const double                mu_q = mu_ptr -> value(Q);

            // Initialize normal_dot_strain with rhs
            for (unsigned int d = 0; d < space_dimension; d++)
              normal_dot_strain[d] = bdry_fcn -> value(Q, d);
            
            // Initialize the strain tensor
            for (const unsigned int dof : face_values.dof_indices()) {
              const double c_i = solution(local_IEN_array[dof]);
              const unsigned int comp_i =
                face_values.system_to_component_index(dof).first;
              function_jacobian[comp_i] += c_i * face_values.shape_grad(dof, q);
            } // for ( i )

            // Get the stress tensor from the jacobian
            stress = symmetrize(function_jacobian);
            for (unsigned int i = 0; i < space_dimension; i++) {
              for (unsigned int j = 0; j < space_dimension; j++) {
                for (unsigned int k = 0; k < space_dimension; k++) {
                  for (unsigned int l = 0; l < space_dimension; l++) {
                    if (i == j && k == l)
                      strain[i][j] += lambda_q * stress[k][l];
                    if (i == k && j == l)
                      strain[i][j] += mu_q * stress[k][l];
                    if (i == l && j == k)
                      strain[i][j] += mu_q * stress[k][l];
                  } // for ( l )
                } // for ( k )
              } // for ( j )
            } // for ( i )

            for (unsigned int d = 0; d < space_dimension; d++)
              normal_dot_strain[d] -= strain[d] * face_values.normal_vector(q);

            face_residuals += normal_dot_strain * normal_dot_strain * face_values.JxW(q);
          } // for ( q )
          cell_residual += face_residuals * face -> diameter();
        } else {
          // In this case, we have to check the multiplicty of thecurrent face
          // to check whether this is a C0 edge or not. 
          const Point<dimension>& c = (-1.) * face -> vertex(0) 
                                            + face -> vertex(nvf - 1);
                unsigned int orientation = 0;
          for (; orientation < dimension && c(orientation) != 0; orientation++);
          // face is not at the boundary and hence an interior face. 

          std::vector< unsigned int > neighboring_faces;
          std::vector< cell_iterator > neighbors;

          if ( face -> has_children() ) {
            // Contnue loop if face has multiplicity less then p, hence it is 
            // not a C0 edge
            if (mof.at(face -> child(0) -> index()) < p[orientation])
              continue; 
            neighboring_faces = {f, cell -> neighbor_face_no(f), cell->neighbor_face_no(f)};
            neighbors = {cell, cell -> neighbor(f) -> child(0), cell -> neighbor(f) -> child(1)};
          } else {
            // Contnue loop if face has multiplicity less then p, hence it is 
            // not a C0 edge
            if (mof.at(face -> index()) < p[orientation])
              continue; 
            neighboring_faces = {f, cell -> neighbor_face_no(f)};
            neighbors = {cell, cell -> neighbor(f)};
          } // if ( face -> has_children() ) 

          for (unsigned int e = 0; e < neighbors.size(); e++){
            face_values.reinit(neighbors[e], neighboring_faces[e]);

            // essentially perform operations from above 
            for (const unsigned int q : face_values.quadrature_point_indices()){
              const Point<space_dimension>     &Q      = face_values.quadrature_point(q);
                    Tensor<2, space_dimension>  function_jacobian;
                    SymmetricTensor<2, space_dimension> stress;
                    Tensor<2, space_dimension>  strain;
                    Tensor<1, space_dimension>  normal_dot_strain;
              const double                lambda_q = lambda_ptr -> value(Q);
              const double                mu_q = mu_ptr -> value(Q);
              
              // Initialize the strain tensor
              for (const unsigned int dof : face_values.dof_indices()) {
                const double c_i = solution(local_IEN_array[dof]);
                const unsigned int comp_i =
                  face_values.system_to_component_index(dof).first;
                function_jacobian[comp_i] += c_i * face_values.shape_grad(dof, q);
              } // for ( i )

              // Get the stress tensor from the jacobian
              stress = symmetrize(function_jacobian);
              for (unsigned int i = 0; i < space_dimension; i++) {
                for (unsigned int j = 0; j < space_dimension; j++) {
                  for (unsigned int k = 0; k < space_dimension; k++) {
                    for (unsigned int l = 0; l < space_dimension; l++) {
                      if (i == j && k == l)
                        strain[i][j] += lambda_q * stress[k][l];
                      if (i == k && j == l)
                        strain[i][j] += mu_q * stress[k][l];
                      if (i == l && j == k)
                        strain[i][j] += mu_q * stress[k][l];
                    } // for ( l )
                  } // for ( k )
                } // for ( j )
              } // for ( i )

              for (unsigned int d = 0; d < space_dimension; d++)
                normal_dot_strain[d] += strain[d] * face_values.normal_vector(q);

              face_residuals += 0.25 * 
                                normal_dot_strain * 
                                normal_dot_strain * 
                                face_values.JxW(q);
            } // for ( q )
          } // for ( e )
          cell_residual += face_residuals * face->diameter();
        } // if ( boundary )
      } // for ( f )

      cell_residual = std::sqrt(cell_residual);
    } // for ( cell )
  } // TS_TriangulationBase<dim, spacedim>::linear_elasticity_residual_error_estimate() [1 / X]

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::
    linear_elasticity_residual_error_estimate(
      const std::vector< unsigned int >&          n_gauss_points,
      const Function<space_dimension>*            rhs_fcn,
      const Function<space_dimension>*            lambda_ptr,
      const Function<space_dimension>*            mu_ptr,
      const std::map< 
                    types::boundary_id,
              const Function<space_dimension>* 
            >&                                    neumann_bc,
      const Vector<double>&                       solution,
            Vector<double>&                       residuals
  ) const {
    Assert(residuals.size() == this-> n_active_cells(), 
              ExcDimensionMismatch(residuals.size(), this->n_active_cells()));

    std::map< cell_iterator, double >   map_residuals;
    linear_elasticity_residual_error_estimate(
      n_gauss_points,
      rhs_fcn,
      lambda_ptr,
      mu_ptr,
      neumann_bc,
      solution,
      map_residuals
    );

    // Move results to different output
    unsigned int i = 0; 
    for (const auto& [_, res] : map_residuals)
      residuals(i++) = res;
  } // TS_Triangulation<dim, spacedim>::linear_elasticity_residual_error_estimate() [2 / X]

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::
    linear_elasticity_residual_error_estimate(
      const std::vector< unsigned int >&          n_gauss_points,
      const Function<space_dimension>*            rhs_fcn,
      const Function<space_dimension>*            lambda_ptr,
      const Function<space_dimension>*            mu_ptr,
      const std::map< 
                    types::boundary_id,
                    Function<space_dimension>* 
            >&                                    neumann_bc,
      const Vector<double>&                       solution,
            Vector<double>&                       residuals
  ) const {
    Assert(residuals.size() == this-> n_active_cells(), 
              ExcDimensionMismatch(residuals.size(), this->n_active_cells()));

    std::map< cell_iterator, double >   map_residuals;
    std::map<
            types::boundary_id,
      const Function<space_dimension>*
      > const_bc;
    for (const auto& [id, fcn] : neumann_bc)
      const_bc[id] = fcn;
    linear_elasticity_residual_error_estimate(
      n_gauss_points,
      rhs_fcn,
      lambda_ptr,
      mu_ptr,
      const_bc,
      solution,
      map_residuals
    );

    // Move results to different output
    unsigned int i = 0; 
    for (const auto& [_, res] : map_residuals)
      residuals(i++) = res;
  } // TS_Triangulation<dim, spacedim>::linear_elasticity_residual_error_estimate() [3 / X]

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::
    linear_elasticity_residual_error_estimate(
      const std::vector< unsigned int >&          n_gauss_points,
      const Function<space_dimension>*            rhs_fcn,
      const Function<space_dimension>*            lambda_ptr,
      const Function<space_dimension>*            mu_ptr,
      const Vector<double>&                       solution,
            Vector<double>&                       residuals
  ) const {
    const std::map< types::boundary_id, Function<space_dimension>* > empty_bc;

    linear_elasticity_residual_error_estimate(
      n_gauss_points,
      rhs_fcn,
      lambda_ptr,
      mu_ptr,
      empty_bc, 
      solution,
      residuals
    );

  } // TS_TriangulationBase<dim, spacedim>::linear_elasticity_residual_error_estimate() [1 / X]

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::
    linear_elasticity_residual_error_estimate(
      const std::vector< unsigned int >&          n_gauss_points,
      const Function<space_dimension>*            rhs_fcn,
      const Functions::ConstantFunction<space_dimension>*    lambda_ptr,
      const Functions::ConstantFunction<space_dimension>*    mu_ptr,
      const std::map< 
                    types::boundary_id,
              const Function<space_dimension>* 
            >&                                    neumann_bc,
      const Vector<double>&                       solution,
      std::map< cell_iterator,
                double >&                         residuals
  ) const {
    Assert(rhs_fcn -> n_components == space_dimension, 
            ExcDimensionMismatch(rhs_fcn -> n_components, space_dimension));
    Assert(lambda_ptr -> n_components == 1, 
            ExcDimensionMismatch(lambda_ptr -> n_components, space_dimension));
    Assert(mu_ptr -> n_components == 1, 
            ExcDimensionMismatch(mu_ptr -> n_components, space_dimension));
#ifdef DEBUG
    for (auto bc_fcn_it = neumann_bc.begin(); bc_fcn_it != neumann_bc.end(); ++bc_fcn_it){
      Assert(bc_fcn_it -> second -> n_components == rhs_fcn -> n_components, 
                ExcDimensionMismatch(bc_fcn_it -> second -> n_components, rhs_fcn -> n_components));
    }
#endif
    const unsigned int &nfc = GeometryInfo<dimension>::faces_per_cell;
    const unsigned int &nvf = GeometryInfo<dimension>::vertices_per_face;
    const unsigned int &n_components = space_dimension;
                
    TSValues<dimension, space_dimension, space_dimension> ts_values(
      this, n_gauss_points,
      update_values |
      update_gradients |
      update_quadrature_points |
      update_JxW_values | 
      update_hessians
    );
    TSFaceValues<dimension, space_dimension, space_dimension> face_values(
      this, n_gauss_points,
      update_values |
      update_gradients |
      update_quadrature_points |
      update_JxW_values  |
      update_normal_vectors);
    for (const auto& cell : this -> active_cell_iterators()) {
      double& cell_residual = residuals[cell];
      cell_residual = 0.;
      const std::vector< unsigned int > local_IEN_array
        = get_IEN_array(cell, n_components);
      ts_values.reinit(cell); 
      
      Tensor<1, space_dimension>    div_c_u;
      for (const unsigned int q : ts_values.quadrature_point_indices()){
        Tensor<1, space_dimension> 
              tmp1, // \nabla \lambda (\nabla \cdot u)
              tmp2, // (\nabla \cdot mu + \mu(\nabla \cdot \nabla)) u
              tmp3; // \nabla \cdot \mu (nabla u)^T
        const Point<space_dimension>  Q = ts_values.quadrature_point(q);
        const double                  lambda_q = lambda_ptr -> value(Q);
        const double                  mu_q = mu_ptr -> value(Q);

        for (const unsigned int i : ts_values.dof_indices()) {
          const double c_i = solution(local_IEN_array[i]);
          const Tensor<2, space_dimension>& hess = 
            ts_values.shape_hessian(i, q);
          const unsigned int comp_i = 
            ts_values.system_to_component_index(i).first;

          // \lambda \nabla (\nabla \cdot u)
          tmp1 += c_i * lambda_q * hess[comp_i];

          // \mu (\nabla \cdot \nabla u)
          for (unsigned int d = 0; d < space_dimension; d++) 
            tmp2[comp_i] += mu_q * hess[d][d];

          // \mu \nabla \cdot (\nabla u)^T
          tmp3 += c_i * mu_q * hess[comp_i];
        } // for ( i )

        Tensor<1, space_dimension> integrand =
                tmp1 + tmp2 + tmp3;
        for (unsigned int d = 0; d < space_dimension; d++)
          integrand[d] += rhs_fcn -> value(Q, d);

        cell_residual += integrand * integrand * ts_values.JxW(q);
      } // for ( q )
      cell_residual *= cell -> diameter() * cell -> diameter();

      for (unsigned int f = 0; f < nfc; f++){
        double face_residuals = 0.;
        const auto& face = cell -> face(f);
        if (face -> at_boundary()){
          const auto& bdry_fcn_ptr = neumann_bc.find(face -> boundary_id()); 
          
          // Is the boundary condition present for neumann? 
          if (bdry_fcn_ptr == neumann_bc.end())
            continue; // no, continue

          const Function<space_dimension>* bdry_fcn = bdry_fcn_ptr -> second;
          Assert(bdry_fcn -> n_components == space_dimension , 
                  ExcDimensionMismatch(bdry_fcn -> n_components, space_dimension));
          face_values.reinit(cell, f); 
          for (const unsigned int q : face_values.quadrature_point_indices()){
            const Point<space_dimension>     &Q      = face_values.quadrature_point(q);
                  Tensor<2, space_dimension>  function_jacobian;
                  SymmetricTensor<2, space_dimension> stress;
                  Tensor<2, space_dimension>  strain;
                  Tensor<1, space_dimension>  normal_dot_strain;
            const double                lambda_q = lambda_ptr -> value(Q);
            const double                mu_q = mu_ptr -> value(Q);

            // Initialize normal_dot_strain with rhs
            for (unsigned int d = 0; d < space_dimension; d++)
              normal_dot_strain[d] = bdry_fcn -> value(Q, d);
            
            // Initialize the strain tensor
            for (const unsigned int dof : face_values.dof_indices()) {
              const double c_i = solution(local_IEN_array[dof]);
              const unsigned int comp_i =
                face_values.system_to_component_index(dof).first;
              function_jacobian[comp_i] += c_i * face_values.shape_grad(dof, q);
            } // for ( i )

            // Get the stress tensor from the jacobian
            stress = symmetrize(function_jacobian);
            for (unsigned int i = 0; i < space_dimension; i++) {
              for (unsigned int j = 0; j < space_dimension; j++) {
                for (unsigned int k = 0; k < space_dimension; k++) {
                  for (unsigned int l = 0; l < space_dimension; l++) {
                    if (i == j && k == l)
                      strain[i][j] += lambda_q * stress[k][l];
                    if (i == k && j == l)
                      strain[i][j] += mu_q * stress[k][l];
                    if (i == l && j == k)
                      strain[i][j] += mu_q * stress[k][l];
                  } // for ( l )
                } // for ( k )
              } // for ( j )
            } // for ( i )

            for (unsigned int d = 0; d < space_dimension; d++)
              normal_dot_strain[d] -= strain[d] * face_values.normal_vector(q);

            face_residuals += normal_dot_strain * normal_dot_strain * face_values.JxW(q);
          } // for ( q )
          cell_residual += face_residuals * face -> diameter();
        } else {
          // In this case, we have to check the multiplicty of thecurrent face
          // to check whether this is a C0 edge or not. 
          const Point<dimension>& c = (-1.) * face -> vertex(0) 
                                            + face -> vertex(nvf - 1);
                unsigned int orientation = 0;
          for (; orientation < dimension && c(orientation) != 0; orientation++);
          // face is not at the boundary and hence an interior face. 

          std::vector< unsigned int > neighboring_faces;
          std::vector< cell_iterator > neighbors;

          if ( face -> has_children() ) {
            // Contnue loop if face has multiplicity less then p, hence it is 
            // not a C0 edge
            if (mof.at(face -> child(0) -> index()) < p[orientation])
              continue; 
            neighboring_faces = {f, cell -> neighbor_face_no(f), cell->neighbor_face_no(f)};
            neighbors = {cell, cell -> neighbor(f) -> child(0), cell -> neighbor(f) -> child(1)};
          } else {
            // Contnue loop if face has multiplicity less then p, hence it is 
            // not a C0 edge
            if (mof.at(face -> index()) < p[orientation])
              continue; 
            neighboring_faces = {f, cell -> neighbor_face_no(f)};
            neighbors = {cell, cell -> neighbor(f)};
          } // if ( face -> has_children() ) 

          for (unsigned int e = 0; e < neighbors.size(); e++){
            face_values.reinit(neighbors[e], neighboring_faces[e]);

            // essentially perform operations from above 
            for (const unsigned int q : face_values.quadrature_point_indices()){
              const Point<space_dimension>     &Q      = face_values.quadrature_point(q);
                    Tensor<2, space_dimension>  function_jacobian;
                    SymmetricTensor<2, space_dimension> stress;
                    Tensor<2, space_dimension>  strain;
                    Tensor<1, space_dimension>  normal_dot_strain;
              const double                lambda_q = lambda_ptr -> value(Q);
              const double                mu_q = mu_ptr -> value(Q);
              
              // Initialize the strain tensor
              for (const unsigned int dof : face_values.dof_indices()) {
                const double c_i = solution(local_IEN_array[dof]);
                const unsigned int comp_i =
                  face_values.system_to_component_index(dof).first;
                function_jacobian[comp_i] += c_i * face_values.shape_grad(dof, q);
              } // for ( i )

              // Get the stress tensor from the jacobian
              stress = symmetrize(function_jacobian);
              for (unsigned int i = 0; i < space_dimension; i++) {
                for (unsigned int j = 0; j < space_dimension; j++) {
                  for (unsigned int k = 0; k < space_dimension; k++) {
                    for (unsigned int l = 0; l < space_dimension; l++) {
                      if (i == j && k == l)
                        strain[i][j] += lambda_q * stress[k][l];
                      if (i == k && j == l)
                        strain[i][j] += mu_q * stress[k][l];
                      if (i == l && j == k)
                        strain[i][j] += mu_q * stress[k][l];
                    } // for ( l )
                  } // for ( k )
                } // for ( j )
              } // for ( i )

              for (unsigned int d = 0; d < space_dimension; d++)
                normal_dot_strain[d] += strain[d] * face_values.normal_vector(q);

              face_residuals += 0.25 * 
                                normal_dot_strain * 
                                normal_dot_strain * 
                                face_values.JxW(q);
            } // for ( q )
          } // for ( e )
          cell_residual += face_residuals * face->diameter();
        } // if ( boundary )
      } // for ( f )

      cell_residual = std::sqrt(cell_residual);
    } // for ( cell )
  } // TS_TriangulationBase<dim, spacedim>::linear_elasticity_residual_error_estimate() [1 / X]

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::
    linear_elasticity_residual_error_estimate(
      const std::vector< unsigned int >&          n_gauss_points,
      const Function<space_dimension>*            rhs_fcn,
      const Functions::ConstantFunction<space_dimension>*    lambda_ptr,
      const Functions::ConstantFunction<space_dimension>*    mu_ptr,
      const std::map< 
                    types::boundary_id,
              const Function<space_dimension>* 
            >&                                    neumann_bc,
      const Vector<double>&                       solution,
            Vector<double>&                       residuals
  ) const {
    Assert(residuals.size() == this-> n_active_cells(), 
              ExcDimensionMismatch(residuals.size(), this->n_active_cells()));

    std::map< cell_iterator, double >   map_residuals;
    linear_elasticity_residual_error_estimate(
      n_gauss_points,
      rhs_fcn,
      lambda_ptr,
      mu_ptr,
      neumann_bc,
      solution,
      map_residuals
    );

    // Move results to different output
    unsigned int i = 0; 
    for (const auto& [_, res] : map_residuals)
      residuals(i++) = res;
  } // TS_Triangulation<dim, spacedim>::linear_elasticity_residual_error_estimate() [2 / X]

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::
    linear_elasticity_residual_error_estimate(
      const std::vector< unsigned int >&          n_gauss_points,
      const Function<space_dimension>*            rhs_fcn,
      const Functions::ConstantFunction<space_dimension>*    lambda_ptr,
      const Functions::ConstantFunction<space_dimension>*    mu_ptr,
      const std::map< 
                    types::boundary_id,
                    Function<space_dimension>* 
            >&                                    neumann_bc,
      const Vector<double>&                       solution,
            Vector<double>&                       residuals
  ) const {
    Assert(residuals.size() == this-> n_active_cells(), 
              ExcDimensionMismatch(residuals.size(), this->n_active_cells()));

    std::map< cell_iterator, double >   map_residuals;
    std::map<
            types::boundary_id,
      const Function<space_dimension>*
      > const_bc;
    for (const auto& [id, fcn] : neumann_bc)
      const_bc[id] = fcn;
    linear_elasticity_residual_error_estimate(
      n_gauss_points,
      rhs_fcn,
      lambda_ptr,
      mu_ptr,
      const_bc,
      solution,
      map_residuals
    );

    // Move results to different output
    unsigned int i = 0; 
    for (const auto& [_, res] : map_residuals)
      residuals(i++) = res;
  } // TS_Triangulation<dim, spacedim>::linear_elasticity_residual_error_estimate() [3 / X]

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::
    linear_elasticity_residual_error_estimate(
      const std::vector< unsigned int >&          n_gauss_points,
      const Function<space_dimension>*            rhs_fcn,
      const Functions::ConstantFunction<space_dimension>*    lambda_ptr,
      const Functions::ConstantFunction<space_dimension>*    mu_ptr,
      const Vector<double>&                       solution,
            Vector<double>&                       residuals
  ) const {
    const std::map< types::boundary_id, Function<space_dimension>* > empty_bc;

    linear_elasticity_residual_error_estimate(
      n_gauss_points,
      rhs_fcn,
      lambda_ptr,
      mu_ptr,
      empty_bc, 
      solution,
      residuals
    );

  } // TS_TriangulationBase<dim, spacedim>::linear_elasticity_residual_error_estimate() [1 / X]

  // =======================================================
  //        Printer section
  // =======================================================
  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::printIPF(
      const std::vector< Point<dim> >& Points,
      const std::string& out_name,
      const unsigned int precision,
      const bool print_splines, 
      const bool print_IPF 
  ) const {
    if (print_splines){
      std::cout << indent(1) << "Printing values of splines ... " << std::endl;
      FullMatrix<double> spline_values(Points.size(), this -> n_active_splines());


      for (unsigned int i = 0; i < Points.size(); i++)
        for (unsigned int j = 0; j < this -> n_active_splines(); j++)
          spline_values(i, j) = active_splines[j] -> value(Points[i]);

      std::string name = out_name + "_splines.dat";
      
      std::filebuf f;
      f.open(name.c_str(), std::ios::out);
      std::ostream out(&f);


      try {
        spline_values.print_formatted(out, precision, true, 1, "0");
        std::cout << indent(2) << "Output printed to " << name << std::endl;
      } catch (const dealii::StandardExceptions::ExcIO& e) {
        const std::size_t found = out_name.find_last_of("/");
        const std::string path = out_name.substr(0, found);
        std::cout << e.what();
        std::cout << "Additional Info:" << std::endl;
        std::cout << "\n\n";
        std::cout << "================================================================================== " << std::endl;
        std::cout << "Output failed, cannot write to " 
                  << path 
                  << ". Ensure you have permissions to write in the specified folder!"
                  << std::endl;
        std::cout << "================================================================================== " << std::endl;
        throw;
      }
      f.close();
      std::cout << indent(1) << "... done!" << std::endl;
    }

    if (print_IPF){
      std::cout << indent(1) << "Printing value of IPF ... " << std::endl;

      IPF.print(Points, out_name);

      std::cout << "... done! " << std::endl;
    }
  } // printIPF [1 / 2]
  
  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::printIPF(
      const unsigned int n_components,
      const std::vector< Point<dim> >& Points,
      const std::string& out_name,
      const unsigned int precision,
      const bool print_splines, 
      const bool print_IPF 
  ) const {
    if (print_splines){
      std::cout << indent(1) << "Printing values of splines ... " << std::endl;
      const unsigned int n_pnts = Points.size();
      const unsigned int n_splines = n_components * this -> n_active_splines();

      std::vector< FullMatrix<double> > spline_values(
        n_components,
        FullMatrix<double>(n_pnts, n_splines)
      );

      for (unsigned int i = 0; i < n_pnts; i++) {
        for (unsigned int j = 0; j < n_splines; j++) {
          const unsigned int base_j = j / n_components;
          const unsigned int comp_j = j % n_components;
          spline_values[comp_j](i, j) = active_splines[base_j] -> value(Points[i]);
        } // for ( j )
      } // for ( i )

      try {
        for (unsigned int n = 0; n < n_components; n++) {
          std::string name = out_name + "_d" + std::to_string(n) + "_splines.dat";
          
          std::filebuf f;
          f.open(name.c_str(), std::ios::out);
          std::ostream out(&f);
          spline_values[n].print_formatted(out, precision, true, 1, "0");
          f.close();
          std::cout << indent(2) << "Output printed to " << name << std::endl;
        } 
      } catch (const dealii::StandardExceptions::ExcIO& e) {
        const std::size_t found = out_name.find_last_of("/");
        const std::string path = out_name.substr(0, found);
        std::cout << e.what();
        std::cout << "Additional Info:" << std::endl;
        std::cout << "\n\n";
        std::cout << "================================================================================== " << std::endl;
        std::cout << "Output failed, cannot write to " 
                  << path 
                  << ". Ensure you have permissions to write in the specified folder!"
                  << std::endl;
        std::cout << "================================================================================== " << std::endl;
        throw;
      }
      std::cout << indent(1) << "... done!" << std::endl;
    }

    if (print_IPF){
      std::cout << indent(1) << "Printing value of IPF ... " << std::endl;

      IPF.print(Points, out_name);

      std::cout << "... done! " << std::endl;
    }
  } // TS_TriangulationBase<dim, spacedim>::printIPF()

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::printIPF(
      const std::string& out_name,
      const unsigned int precision,
      const bool print_splines,
      const bool print_IPF 
  ) const {
    Assert(is_bezier_mesh, ExcInvalidState());
#ifndef DEBUG
      if (!is_bezier_mesh)
        throw ExcInvalidState();
#endif

    //std::cout << "====================================================================" << std::endl;
    //std::cout << "================= Printing IPF and/or Splines ======================" << std::endl;
    //std::cout << "====================================================================" << std::endl;


    std::vector< unsigned int > n_gauss_points = p; 
    for (unsigned int d = 0; d < dimension; d++)
      n_gauss_points[d] += 1;
    TSValues<dimension, space_dimension> ts_values(
          this, n_gauss_points,
          update_values |
          update_gradients |
          update_quadrature_points |
          update_JxW_values);
    const unsigned int n_cells                    = this -> n_active_cells();
    const unsigned int n_quadrature_points       = ts_values.n_quadrature_points_per_cell();
    const std::vector< unsigned int>& n_q_points = ts_values.n_q_points_per_cell();
    if (print_splines && print_IPF){
      FullMatrix<double> B1(n_quadrature_points * n_cells, n_splines);
      FullMatrix<double> B2(n_quadrature_points * n_cells, space_dimension);
      FullMatrix<double> b1(n_quadrature_points, n_splines);
      FullMatrix<double> b2(n_quadrature_points, space_dimension);
      unsigned int index = 0;
      for (const auto& cell : this -> active_cell_iterators()){
        const auto& loc_to_glob = this -> get_IEN_array(cell);
        unsigned int row = 0;
        ts_values.reinit(cell);

        // For each quadrature point store the values of quadrature point
        for (const unsigned int q : ts_values.quadrature_point_indices()){
          for (const unsigned int j : ts_values.dof_indices() )
            b1(row, loc_to_glob[j]) = ts_values.shape_value(j, q);

          for (unsigned int d = 0; d < spacedim; d++)
            b2(row, d) = ts_values.quadrature_point(q)(d);

        row++;
        }

        B1.fill(b1, index, 0);
        B2.fill(b2, index, 0);

        index += n_quadrature_points;
      } // for ( cell )

      std::string name1 = out_name
                          + "_splines.dat";
      std::string name2 = out_name
                          + "_IPF"
                          + std::to_string(spacedim)
                          + "d.dat";
      std::filebuf f1, f2;
      f1.open(name1.c_str(), std::ios::out);
      f2.open(name2.c_str(), std::ios::out);

      std::ostream out1(&f1);
      std::ostream out2(&f2);

      try {
        B1.print_formatted(out1, precision, true, 1, "0");
        B2.print_formatted(out2, precision, true, 1, "0");
        std::cout << indent(2) << "Splines printed to " << name1 << std::endl;
        std::cout << indent(2) << "IPF printed to " << name2 << std::endl;
      } catch (const dealii::StandardExceptions::ExcIO& e) {
        const std::size_t found = out_name.find_last_of("/");
        const std::string path = out_name.substr(0, found);
        std::cout << e.what();
        std::cout << "Additional Info:" << std::endl;
        std::cout << "\n\n";
        std::cout << "================================================================================== " << std::endl;
        std::cout << "Output failed, cannot write to " 
                  << path 
                  << ". Ensure you have permissions to write in the specified folder!"
                  << std::endl;
        std::cout << "================================================================================== " << std::endl;
        throw;
      }
      f1.close();
      f2.close();

      std::cout << indent(1) << "... done!" << std::endl;
    } else if (print_splines) {
      FullMatrix<double> B(n_quadrature_points * n_cells, n_splines);
      FullMatrix<double> b(n_quadrature_points , n_splines);
      unsigned int index = 0;
      for (const auto& cell : this -> active_cell_iterators()){
        unsigned int row = 0;
        ts_values.reinit(cell);

        // For each quadrature point store the values of quadrature point
        for (const unsigned int q : ts_values.quadrature_point_indices()){
          for (const unsigned int j : ts_values.dof_indices() )
            b(row, j) = ts_values.shape_value(j, q);

          row++;
        }

        B.fill(b, index, 0);
        index += n_quadrature_points;
      } // for ( cell )
      
      std::string name = out_name
                         + "_splines.dat";
      std::filebuf f;
      f.open(name.c_str(), std::ios::out);
      std::ostream out(&f);


      try {
        B.print_formatted(out, precision, true, 1, "0");
        std::cout << indent(2) << "Output printed to " << name << std::endl;
      } catch (const dealii::StandardExceptions::ExcIO& e) {
        const std::size_t found = out_name.find_last_of("/");
        const std::string path = out_name.substr(0, found);
        std::cout << e.what();
        std::cout << "Additional Info:" << std::endl;
        std::cout << "\n\n";
        std::cout << "================================================================================== " << std::endl;
        std::cout << "Output failed, cannot write to " 
                  << path 
                  << ". Ensure you have permissions to write in the specified folder!"
                  << std::endl;
        std::cout << "================================================================================== " << std::endl;
        throw;
      }
      f.close();
      std::cout << indent(1) << "... done!" << std::endl;
    } else if (print_IPF) {
      FullMatrix<double> B(n_quadrature_points * n_cells, space_dimension);
      FullMatrix<double> b(n_quadrature_points, space_dimension);
      unsigned int index = 0;
      for (const auto& cell : this -> active_cell_iterators()){
        unsigned int row = 0;
        ts_values.reinit(cell);

        // For each quadrature point store the values of quadrature point
        for (const unsigned int q : ts_values.quadrature_point_indices()){
          for (unsigned int d = 0; d < space_dimension; d++)
            b(row, d) = ts_values.quadrature_point(q)(d);

          row++;
        }

        B.fill(b, index, 0);

        index += n_quadrature_points;
      } // for ( cell )

      std::string name = out_name
                         + "_IPF_"
                         + std::to_string(spacedim)
                         + "d.dat";
      std::filebuf f;
      f.open(name.c_str(), std::ios::out);
      std::ostream out(&f);

      try {
        B.print_formatted(out, precision, true, 1, "0");
        std::cout << indent(2) << "Output printed to " << name << std::endl;
      } catch (const dealii::StandardExceptions::ExcIO& e) {
        const std::size_t found = out_name.find_last_of("/");
        const std::string path = out_name.substr(0, found);
        std::cout << e.what();
        std::cout << "Additional Info:" << std::endl;
        std::cout << "\n\n";
        std::cout << "================================================================================== " << std::endl;
        std::cout << "Output failed, cannot write to " 
                  << path 
                  << ". Ensure you have permissions to write in the specified folder!"
                  << std::endl;
        std::cout << "================================================================================== " << std::endl;
        throw;
      }
      f.close();
      std::cout << indent(1) << "... done!" << std::endl;
    } // if
    Vector<double> data(4);
    data(0) = n_cells;
    data(3) = n_splines;
    for (unsigned int d = 0; d < dimension; d++)
      data(d+1) = n_q_points[d];

    std::string data_name = out_name + "_data.dat";
    std::filebuf dat;
    dat.open(data_name.c_str(), std::ios::out);

    std::ostream data_out(&dat);

    data.print(data_out, 5);
    dat.close();

  } // printIPF [2 / 2]

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::print_IPF_wireframe(
      const std::string& out_name,
      const unsigned int precision,
      const unsigned int n_intermediate_points
  ) const {
    // get each line:
    std::vector< line_iterator > lines;
    const unsigned int nlc = GeometryInfo<dimension>::lines_per_cell;
    for (const auto& cell : this -> active_cell_iterators())
      for (unsigned int n = 0; n < nlc; n++)
        if ( cell -> line(n) -> used() )
          lines.push_back(cell -> line(n));

    const unsigned int n_lines = lines.size();
          unsigned int n = 0;
    FullMatrix<double> B(n_lines*(n_intermediate_points), 2 * space_dimension);
    FullMatrix<double> b(n_intermediate_points-1, 2 * space_dimension);
    for (unsigned int j = 0; j < n_lines; j++){
      const auto& l = lines[j];
      const Point<dimension>& p0 = l -> vertex(0);
      const Point<dimension>& p1 = l -> vertex(1);
      for (unsigned int i = 0; i < n_intermediate_points-1; i++){
        const Point<dimension>& lower = (p0 + i*( p1 + (-1.)*p0 )/(n_intermediate_points-1));
        const Point<dimension>& upper = (p0 + (i+1)*( p1 + (-1.)*p0 )/(n_intermediate_points-1));

        Vector<double> IPF_lower, IPF_upper;
        IPF.vector_value(lower, IPF_lower);
        IPF.vector_value(upper, IPF_upper);

        for (unsigned int d = 0; d < space_dimension; d++){
          b(i, d) = IPF_lower(d);
          b(i, d + spacedim) = IPF_upper(d);
        }
      }
      B.fill(b, n, 0);
      n += n_intermediate_points-1;
    }

    std::filebuf f;
    std::string name = out_name + "_physical";
    if (is_bezier_mesh)
      name += "_bezier_grid.dat";
    else 
      name += "_grid.dat";
    f.open(name.c_str(), std::ios::out);
    std::ostream out(&f);

    try {
      B.print_formatted(out, precision, true, 1, "0");
      std::cout << indent(2) << "Output printed to " << name << std::endl;
    } catch (const dealii::StandardExceptions::ExcIO& e) {
      const std::size_t found = out_name.find_last_of("/");
      const std::string path = out_name.substr(0, found);
      std::cout << e.what();
      std::cout << "Additional Info:" << std::endl;
      std::cout << "\n\n";
      std::cout << "================================================================================== " << std::endl;
      std::cout << "Output failed, cannot write to " 
                << path 
                << ". Ensure you have permissions to write in the specified folder!"
                << std::endl;
      std::cout << "================================================================================== " << std::endl;
      throw;
    }
    f.close();
  } // print_IPF_wireframe

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::print_IPF_surface_wireframe(
      const std::string& out_name,
      const unsigned int precision ,
      const unsigned int n_intermediate_points 
  ) const {
    // get each line:
    std::vector< line_iterator > lines;
    const unsigned int nlc = GeometryInfo<dimension>::lines_per_cell;
    for (const auto& cell : this -> active_cell_iterators())
      for (unsigned int n = 0; n < nlc; n++)
        if ( cell -> line(n) -> used() && 
              cell -> line(n) -> at_boundary() )
          lines.push_back(cell -> line(n));

    const unsigned int n_lines = lines.size();
          unsigned int n = 0;
    FullMatrix<double> B(n_lines*(n_intermediate_points), 2*space_dimension);
    FullMatrix<double> b(n_intermediate_points, 2*space_dimension);
    for (unsigned int j = 0; j < n_lines; j++){
      b = 0;
      const auto& l = lines[j];
      const Point<dimension>& p0 = l -> vertex(0);
      const Point<dimension>& p1 = l -> vertex(1);
      for (unsigned int i = 0; i < n_intermediate_points-1; i++){
        const Point<dimension>& lower = (p0 + i*( p1 + (-1.)*p0 )/(n_intermediate_points-1));
        const Point<dimension>& upper = (p0 + (i+1)*( p1 + (-1.)*p0 )/(n_intermediate_points-1));

        Vector<double> IPF_lower, IPF_upper;
        IPF.vector_value(lower, IPF_lower);
        IPF.vector_value(upper, IPF_upper);

        for (unsigned int d = 0; d < space_dimension; d++){
          b(i, d) = IPF_lower(d);
          b(i, d + spacedim) = IPF_upper(d);
        }
      }
      B.fill(b, n, 0);
      n += n_intermediate_points;
    } // for ( j )

    std::filebuf f;
    std::string name = out_name + "_parametric_surface";
    if (is_bezier_mesh)
      name += "_bezier_grid.dat";
    else 
      name += "_grid.dat";

    f.open(name.c_str(), std::ios::out);
    std::ostream out(&f);

    try {
      B.print_formatted(out, precision, true, 1, "0");
      std::cout << indent(2) << "Output printed to " << name << std::endl;
    } catch (const dealii::StandardExceptions::ExcIO& e) {
      const std::size_t found = out_name.find_last_of("/");
      const std::string path = out_name.substr(0, found);
      std::cout << e.what();
      std::cout << "Additional Info:" << std::endl;
      std::cout << "\n\n";
      std::cout << "================================================================================== " << std::endl;
      std::cout << "Output failed, cannot write to " 
                << path 
                << ". Ensure you have permissions to write in the specified folder!"
                << std::endl;
      std::cout << "================================================================================== " << std::endl;
      throw;
    }
    f.close();
  } // print_IPF_durface_wireframe

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::print_IPF_vertices(
      const std::string& out_name,
      const unsigned int precision
  ) const {
    const unsigned int nac = this -> n_active_cells();
    const unsigned int nvc = GeometryInfo<dimension>::vertices_per_cell;
    FullMatrix<double> B(nac, space_dimension * nvc);
    FullMatrix<double> b(1, space_dimension * nvc);

    unsigned int row = 0;
    for (const auto& cell : this -> active_cell_iterators()){
      unsigned int col = 0;
      for (unsigned int n = 0; n < nvc; n++){
        const Point<dimension>& p = cell -> vertex(n);
        Vector<double> val_p(space_dimension);

        IPF.vector_value(p, val_p);
        for (int d = 0; d < spacedim; d++)
          b(0, col++) = val_p(d);
      }
      B.fill(b, row++);
    }

    std::filebuf f;
    std::string name = out_name + "_parametric_grid_" +std::to_string(spacedim) + "d_vertices.dat";
    f.open(name.c_str(), std::ios::out);
    std::ostream out(&f);

    try {
      B.print_formatted(out, precision, true, 1, "0");
      std::cout << indent(2) << "Output printed to " << name << std::endl;
    } catch (const dealii::StandardExceptions::ExcIO& e) {
      const std::size_t found = out_name.find_last_of("/");
      const std::string path = out_name.substr(0, found);
      std::cout << e.what();
      std::cout << "Additional Info:" << std::endl;
      std::cout << "\n\n";
      std::cout << "================================================================================== " << std::endl;
      std::cout << "Output failed, cannot write to " 
                << path 
                << ". Ensure you have permissions to write in the specified folder!"
                << std::endl;
      std::cout << "================================================================================== " << std::endl;
      throw;
    }
    f.close();
          
  } // print_IPF_vertices

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::print_grid(
      const std::string& out_name,
      const unsigned int precision
  ) const {
    const unsigned int nlc = GeometryInfo<dimension>::lines_per_cell;
    std::vector< line_iterator > lines;
    for (const auto& cell : this -> active_cell_iterators())
      for (unsigned int n = 0; n < nlc; n++)
        if (cell -> line(n) -> used())
          lines.push_back(cell -> line(n));

    const unsigned int n_lines = lines.size();
    FullMatrix<double> B(n_lines, 4);
    for (unsigned int i = 0; i < n_lines; i++){
      const auto& l = lines[i];
      const Point<dimension>& P0 = l -> vertex(0);
      const Point<dimension>& P1 = l -> vertex(1);

      B(i, 0) = P0(0);
      B(i, 1) = P0(1);
      B(i, 2) = P1(0);
      B(i, 3) = P1(1);
    } // for ( i )

    std::filebuf f;
    //std::string name = "out/" + add + "grid_level_" + std::to_string(this->n_levels() - 1) + ".dat";
    std::string name = out_name + "_grid.dat";
    f.open(name.c_str(), std::ios::out);
    std::ostream out(&f);

    try {
      B.print_formatted(out, precision, true, 1, "0");
      std::cout << indent(2) << "Output printed to " << name << std::endl;
    } catch (const dealii::StandardExceptions::ExcIO& e) {
      const std::size_t found = out_name.find_last_of("/");
      const std::string path = out_name.substr(0, found);
      std::cout << e.what();
      std::cout << "Additional Info:" << std::endl;
      std::cout << "\n\n";
      std::cout << "================================================================================== " << std::endl;
      std::cout << "Output failed, cannot write to " 
                << path 
                << ". Ensure you have permissions to write in the specified folder!"
                << std::endl;
      std::cout << "================================================================================== " << std::endl;
      throw;
    }
    f.close();
  } // print_grid()

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::print_bezier_grid(
      const std::string& out_name,
      const unsigned int precision
  ) const {
    Assert(is_bezier_mesh, ExcInvalidState());
#ifndef DEBUG
    if (!is_bezier_mesh)
      throw ExcInvalidState();
#endif
    print_grid(out_name + "_bezier", precision);
  } // print_bezier_grid()

  template<int dim, int spacedim>
  const std::string TS_TriangulationBase<dim, spacedim>::kv_to_string(
      const std::vector<double>& kx
  ) const {
    std::string s = "";
    for(const double & x : kx)
      s += std::to_string(x)+"  ";
    s+="\n";
    return s;
  } // kv_to_string()

  template<int dim, int spacedim>
  template<int soldim>
  void TS_TriangulationBase<dim, spacedim>::generate_mesh_file(
        const std::string&  out_name, 
        const bool          parametric,
        const unsigned int  precision,
        const std::map< unsigned int, Point<soldim>>&
          vertex_values
  ) const { 
    if (soldim != 0)
      Assert(soldim == 1 || soldim == spacedim, ExcDimensionMismatch2(soldim, 1, spacedim));

    // Declare output names
    const std::string& bezier_string  = is_bezier_mesh ? "_bezier" : "";
    const std::string& par_phy_string = parametric ? "_parametric" : "_physical";
    const std::string& cells = out_name + par_phy_string + bezier_string + "_cell_list.dat";
    const std::string& verts = out_name + par_phy_string + bezier_string + "_vertex_list.dat";
    const std::string& mesh = out_name + par_phy_string + bezier_string + "_mesh.txt";

    // Declare variables
    const unsigned int nvc = GeometryInfo<dimension>::vertices_per_cell;
    const unsigned int solution_dimension = soldim;

    Assert(parametric 
            || soldim == 0 
            || soldim == 1 
            || soldim == spacedim, 
            ExcMessage("Invalid combination of arguments!"));


    // Declare container to store necessary data
    FullMatrix<double> cell_list(this->n_active_cells(), std::pow(2, dimension)); 

    // Declare container for vertices
    std::map<unsigned int, Point<dimension>> vert_map; 

    // Fill the containers with necessary data
    unsigned int ind = 0;
    for (const auto& cell : this -> active_cell_iterators()) {
      for (unsigned int v = 0; v < nvc; v++){
        cell_list(ind, v) = cell -> vertex_index(v);
        vert_map[cell->vertex_index(v)] = cell -> vertex(v);
      }
      ind++;
    }

    // Declare amount of vertices and the list for output
    const unsigned int n_vertices = vert_map.size();
    const unsigned int final_dimension = parametric ? dimension 
                                          : (vertex_values.empty() ? space_dimension 
                                              : space_dimension + solution_dimension);
    if (!vertex_values.empty())
      Assert(vert_map.size() == vertex_values.size(), 
            ExcDimensionMismatch(vert_map.size(), vertex_values.size()));

    FullMatrix<double> vert_list(
         n_vertices, 
         final_dimension
    ); 

      

    if (parametric) {
      // fill the list from the map
      for (const auto& [ind, v] : vert_map)
        for (unsigned int d = 0; d < dimension; d++)
          vert_list(ind, d) = v(d);
    } else {
      // fill the list from the mapped vertices
      for (const auto& [ind, v] : vert_map){
        const Point< spacedim >& mapped_v = IPF.point_value(v);
        for (unsigned int d = 0; d < space_dimension; d++) 
          vert_list(ind, d) = mapped_v(d);

        for (unsigned int d = 0; d < solution_dimension; d++)
          vert_list(ind, d + space_dimension) = vertex_values.at(ind)(d);
      }
    } // if (parametric)


    std::ofstream mesh_out(mesh, std::ios::out | std::ios::trunc);
    // header: 
    mesh_out << std::pow(2, dimension) 
             << " " 
             << this->n_active_cells() 
             << " " 
             << n_vertices 
             << " "
             << solution_dimension
             << " " 
             << final_dimension
             << std::endl;

    // Cell data: 
    for (unsigned int i = 0; i < this -> n_active_cells(); i++){
      for (unsigned int d = 0; d < nvc; d++)
        mesh_out << cell_list(i, d) << " ";

      mesh_out << std::endl;
    }

    // vertex data:
    for (unsigned i = 0; i < n_vertices; i++){
      for (unsigned int d = 0; d < final_dimension; d++)
        mesh_out << vert_list(i, d) << " ";

      mesh_out << std::endl;
    }

    mesh_out.close();


    // output lists to files
    std::filebuf cell_file, vert_file;
    cell_file.open(cells.c_str(), std::ios::out);
    vert_file.open(verts.c_str(), std::ios::out);

    std::ostream cells_out(&cell_file), verts_out(&vert_file);
    cell_list.print_formatted(cells_out, precision, true, 1, "0");
    vert_list.print_formatted(verts_out, precision, true, 1, "0");
  } // generate_mesh_file

  // =======================================================
  //        Manipulating the grid
  // =======================================================
  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::refine_fixed_number(
    const Vector<double>&       residuals, 
    const double                percentile
  ) {
    Assert(residuals.size() == this->n_active_cells(),
            ExcDimensionMismatch(residuals.size(), this->n_active_cells()));
    Assert(is_bezier_mesh, ExcInvalidState());

    unsigned int i = 0; 
    std::vector< cell_iterator > mark;
    std::map< double, 
              active_cell_iterator, 
              std::greater<double> 
            > cell_residuals;
    for (const auto& cell : this->active_cell_iterators())
      cell_residuals[residuals(i++)] = cell;
    
    const unsigned int top_marked = std::ceil(percentile * this -> n_active_cells());
          auto         map_iterator = cell_residuals.begin();
    while (mark.size() < top_marked && 
            map_iterator != cell_residuals.end() ) {
      if (map_iterator -> second -> level() == 0) {
        mark.push_back(map_iterator -> second);
      } else if (std::find(bezier_elements.begin(), 
                           bezier_elements.end(), 
                           map_iterator -> second -> parent() )
                    != bezier_elements.end() 
      ) {
        mark.push_back(map_iterator -> second -> parent());
      } else {
        mark.push_back(map_iterator -> second);
      }
      map_iterator++;
    }


    // prepare refinement: 
    this->coarsen_bezier_elements();
    this->set_refine_flags(mark);

    // execute refinement on marked cells
    this->execute_coarsening_and_refinement();
  } // TS_TriangulationBase<dim, spacedim>::refine_fixed_number() [1 / 2]

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::refine_fixed_number(
      const std::map< 
              active_cell_iterator, 
              double>&        residuals,
      const unsigned int      percentile
  ) {
    Vector<double> res(this -> n_active_cells());
    unsigned int   i = 0;
    for (const auto& [_, r] : residuals)
      res(i++) = r; 

    refine_fixed_number(res, percentile);
  } // TS_TriangulationBase<dim, spacedim>::refine_fixed_number() [2 / 2]

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::execute_coarsening_and_refinement(
  ) {
  #ifdef DEBUG
    std::cout << indent(2) << "calculating coarse neighborhood of marked cells ... " << std::endl;
  #endif
    Assert(!is_bezier_mesh, ExcInvalidState());
#ifndef DEBUG
    if (is_bezier_mesh)
      throw ExcInvalidState();
#endif

    // get coarse neighborhood
    std::vector<active_cell_iterator> coarse_neighborhood = {};
    get_coarse_neighborhood(coarse_neighborhood);

    // In case no cells were marked for refinement, and this function is called by accident, do nothing
    if (coarse_neighborhood.size() == 0){
  #ifdef DEBUG
      std::cout << indent(2) << "... no cells marked for refinement, returning." << std::endl;
  #endif
      return;
    }

    // Iteratively obtain coarse neighborhood
    unsigned int n_cells_cn = 0;
    while (n_cells_cn != coarse_neighborhood.size()){
      n_cells_cn = coarse_neighborhood.size();
      coarse_neighborhood = {};
      get_coarse_neighborhood(coarse_neighborhood);
    }

    Assert(coarse_neighborhood.size() > 0, ExcInternalError());

  #ifdef DEBUG
    std::cout << indent(2) << "... done! Collecting affected splines ... " << std::endl;
  #endif

    // before refinement, collect the TSplines to be updated
    std::vector< ts_ptr > to_be_updated;
    get_affected_splines( coarse_neighborhood, to_be_updated);

    Assert(to_be_updated.size() > 0, ExcInternalError());

  #ifdef DEBUG
    std::cout << indent(2) << "... done! Performing refinement of grid..." << std::endl;
  #endif
  //  {
  //    GridOutFlags::Svg svg_flags;
  //    svg_flags.coloring = GridOutFlags::Svg::Coloring::none;
  //    svg_flags.label_level_number = true;
  //    svg_flags.label_cell_index = true;
  //
  //    std::string name = "out/grid_br" + std::to_string(this->n_levels() - 1) + ".svg";
  //    std::ofstream out(name);
  //    GridOut       grid_out;
  //    grid_out.set_flags(svg_flags);
  //
  //    grid_out.write_svg(*this, out);
  //  }

    // refine coarse neighborhood until it is empty, i.e. there are no
    // coarser elements in the environment of a marked cell.
    this->Triangulation<dimension, dimension>::execute_coarsening_and_refinement();
  //  {
  //    GridOutFlags::Svg svg_flags;
  //    svg_flags.coloring = GridOutFlags::Svg::Coloring::none;
  //    svg_flags.label_level_number = true;
  //    svg_flags.label_cell_index = true;
  //
  //    std::string name = "out/grid_ar" + std::to_string(this->n_levels() - 1) + ".svg";
  //    std::ofstream out(name);
  //    GridOut       grid_out;
  //    grid_out.set_flags(svg_flags);
  //
  //    grid_out.write_svg(*this, out);
  //  }


  #ifdef DEBUG
    std::cout << indent(2) << "... done! Generating new TSplines from new faces ... " << std::endl;
  #endif


    // And update mof and splines simultaneously
    for (const auto& cell : coarse_neighborhood) {
      Assert(cell -> has_children(), ExcInternalError());
      // Get the inner faces of the refined cell
      const std::vector<face_iterator>& inner_faces = get_inner_faces(cell);


      // For each inner face we perform the update for our tsplines
      // Note: at this point, it is arbitrary in which order we refine, this still needs to be adjusted above!
      for (const auto& face : inner_faces){
        // Get orthogonal direction and the sest of splines affected by this refinement:
        unsigned int dr;
        std::vector<std::vector< ts_ptr >> selected_face =
                get_affected_splines(face, to_be_updated, dr);

        // save the value for knot insertion:
        const double in = (face -> vertex(0)).operator()(dr);
        for (auto& selected : selected_face){
          // Some safety check, see get_affected_splines for details
          Assert(selected.size() == p[dr] + 1, ExcInternalError());

  #ifdef DEBUG
          for (const auto& ts : selected)
            Assert(ts -> is_active(), ExcInternalError());
  #endif

          // collect the global knot vector over all TSplines, set corresponding spline
          // to be in-active, and collect the corresponding CPs
          std::vector<double>        global_kv = selected[0] -> get_local_kv(dr);

          std::vector<ControlPoint>  global_cp(p[dr] + 1);
          global_cp[0] = selected[0] -> get_cp();
          selected[0] -> set_active(false);
          for (unsigned int n = 1; n < p[dr]+1; n++){
            const auto& selected_kv = selected[n] -> get_local_kv(dr);
            merge_kv({selected_kv[p[dr] + 1]}, global_kv, dr);
            global_cp[n] = selected[n] -> get_cp();
            selected[n] -> set_active(false); // de-activate spline ...
          }

          // Check if evrything went as planned:
#ifdef DEBUG
          unsigned int n_kv = global_kv.size();
          Assert(n_kv == 2*(p[dr] + 1), ExcInternalError());
#endif

          // Generate the new control points by linear combination of old points
          std::vector<ControlPoint> new_cp(p[dr] + 2);
          new_cp[0]         = global_cp[0];       // the first one will not change
          new_cp[p[dr] + 1] = global_cp[p[dr]];   // the last one will not change
          for (unsigned int k = 1; k < p[dr] + 1 /* = n_kv - (p[dr] + 2)  - 1*/; k++){
            double alpha   = (in - global_kv[k])/(global_kv[k+p[dr]] - global_kv[k]);
            new_cp[k]    = alpha*global_cp[k] + (1. - alpha)*global_cp[k-1];
          }

          // add the new knot into the kv
          merge_kv({in}, global_kv, dr);

          // initialize new splines: The two following for-loops are identical in its body,
          // we can simply generate the new splines with its corresponding father
          unsigned int kappa = p[dr]/2;
          std::vector< ts_ptr > new_ts(p[dr] + 2);
          for (unsigned int k = 0; k < p[dr] + 2 ; k++){
            bool offset = (k > kappa);
            std::vector<std::vector<double>> kv = selected[k - offset] -> get_kv();
            for (unsigned int i = 0; i < p[dr] + 2; i++) kv[dr][i] = global_kv[i+k];

            new_ts[k] = std::make_shared<TSpline>(
                    kv,
                    new_cp[k],
                    selected[k - offset]);
          }

          // Special case oocurs for different degrees:
          // If p[dr] is even, we generate two new TSplines who share the same father
          // If p[dr] is odd, we generate one new TSpline with no father. In this case,
          // the kv can be initialized by one of its neighbors, e.g. kappa + 1, since the
          // mesh will be analysis suitable at all times.
          if (p[dr] % 2 == 1)
            new_ts[kappa+1] -> clear_father();

          // Add elements to basis:
          // 1. Locally for the next inserted face
          to_be_updated.insert( to_be_updated.end(), new_ts.begin(), new_ts.end() );

          // 2. To the global base of active splines
          this -> active_splines.insert(active_splines.end(), new_ts.begin(), new_ts.end());

          // 3. Put old splines [in-active splines] from the active set to the in-active set
          typename std::vector< ts_ptr >::iterator it = this -> active_splines.begin();
          for (; it != active_splines.end(); ){
            if (!(*it) -> is_active()){
              active_splines.erase(it);
              old_splines.insert(old_splines.end(), *it);
            }
            else
              ++it;
          } // for ( it )

          it = to_be_updated.begin();
          for (; it != to_be_updated.end() ;){
            if (!(*it) -> is_active())
                    to_be_updated.erase(it);
            else
                    ++it;
          } // for( it )
        } // for ts_slice
      } // for face
    } // for cell : coarse_neighborhood

  #ifdef DEBUG
    std::cout << indent(2) << "... done! Finishing refinement process ..." << std::endl;
  #endif

    this -> n_splines = active_splines.size();
    for (unsigned int i = 0; i < n_splines; i++)
      active_splines[i] -> set_level_index(i);

    // Update IPF for evaluation at beyier elements
    this -> IPF = IsoparametricFunction<dimension, space_dimension>(active_splines);

    // Setup new mof map for refined grid
    this -> setup_mof();

    // Prepare the grid to switch between T-mesh and Beyier mesh
    this -> find_bezier_elements();

#ifdef DEBUG
  std::cout << indent(2) << "... done!" << std::endl;
#endif

  } // execute_coarsening_and_refinement()

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::refine_global(
      const unsigned int times
  ) {
    for (unsigned int i = 0; i < times; i++){
      for (auto& cell : this->active_cell_iterators()){
        RefinementCase<dimension> ref = RefinementCase<dimension>::cut_axis(
                        (level_offset + cell -> level()) % dimension);
        cell -> set_refine_flag(ref);
      }
      this -> execute_coarsening_and_refinement();
    }
  } // refine_global ()

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::degree_elevate_global(
      const unsigned int times
  ) {
    if (times == 0)
      return; // nothing to do here

    for (unsigned int d = 0; d < dimension; d++)
      degree_elevate(d, times);
  } // degree_elevate_global()

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::degree_elevate(
      const unsigned int d,
      const unsigned int times
  ) {
  #ifdef DEBUG
    std::cout << indent(1) << "Elevating degree in direction " << d << " by " << times << " ... " << std::endl;
  #endif


    if (times == 0)
      return;

  #ifdef DEBUG
    Assert(this->n_levels() == 1, ExcMessage("The Triangulation has already been refined. Degree elevation is not possible in this case."));
    Assert(d < dimension, ExcLowerRange(d, dimension));
  #else
    if (this -> n_levels() > 1)
      throw ExcInvalidState();

    if (d > dimension)
      throw ExcLowerRange(d, dimension);
  #endif

    const unsigned int n_bezier_knots  = base_knots[d].size();
    // use update of multiplicities to create global kv
    std::vector< double > global_kv_old;
    std::vector< double >::iterator it = global_kv_old.end();
    for (unsigned int i = 0; i < n_bezier_knots; i++){
      global_kv_old.insert(it, multiplicities[d][i], base_knots[d][i]);
      it = global_kv_old.end();
    }
    const unsigned int n_old_knots = global_kv_old.size();


    // Calculate the amount of new splines:
    unsigned int s = 0;
    for (unsigned int i = 1; i < n_bezier_knots + 1; i++) s++; // Get the amount of distinct inner knots

    const unsigned int n_splines_per_chunk = n_old_knots - (p[d] + 1);
    const unsigned int n_spline_chunks     = active_splines.size() / n_splines_per_chunk;
    const unsigned int n_new_splines       = n_splines_per_chunk + (s - 1)*times;

    {
      // We don't need d1, d2 anymore, after we have sorted the TSplines
      // accordingly, so we can just delete them at the end-of-scope
      unsigned int d1, d2;
      switch(d){
      case 0:
        d1 = 1; d2 = (dimension == 2 ? d1 : 2);
        break;
      case 1:
        d1 = 0; d2 = (dimension == 2 ? d1 : 2);
        break;
      case 2:
        d1 = 0; d2 = 1;
        break;
      default:
        d1 = -1; d2 = -1;
        ExcInternalError();
        break;
      }

      // sort the active tsplines accordingly and seperate them in chunks
      std::sort(active_splines.begin(), active_splines.end(),
                      [d1, d2, d](const ts_ptr& t1, const ts_ptr& t2){
                        if (t1 -> get_barycenter(d2) != t2 -> get_barycenter(d2))
                          return t1 -> get_barycenter(d2) < t2 -> get_barycenter(d2);
                        else if (t1 -> get_barycenter(d1) != t2 -> get_barycenter(d1))
                          return t1 -> get_barycenter(d1) < t2 -> get_barycenter(d1);
                        else
                          return t1 -> get_barycenter(d) < t2 -> get_barycenter(d);
                      });
    }

    int n = 0;
    std::vector< std::vector< ts_ptr > > sorted(n_spline_chunks);
    for (unsigned int i = 0; i < n_spline_chunks; i++){
      std::vector< ts_ptr > sorted_i(n_splines_per_chunk);
      for (unsigned int j = 0; j < n_splines_per_chunk; j++)
        sorted_i[j] = active_splines[n++];

      sorted[i] = sorted_i;
    }

    // Construct the global kv after degree elevation:
    // Update multiplicities of each knot
    for (auto& m : multiplicities[d])
      m += times;

    // Update polynomial degree
    p[d] += times;
    TSplineFunction<dimension, space_dimension>::setSplineData(p);

    // use update of multiplicities to create global kv
    std::vector< double > global_kv;
    it = global_kv.end();
    for (unsigned int i = 0; i < n_bezier_knots; i++){
      global_kv.insert(it, multiplicities[d][i], base_knots[d][i]);
      it = global_kv.end();
    }

    Assert(n_new_splines == global_kv.size() - (p[d] + 1), ExcInternalError());

    // Calculate the factorials for binomial coefficients and store them
    // in an array for faster access.
    const unsigned int n1 = p[d] + 1;
    const unsigned int n2 = p[d] + 1;
    std::vector<long double> factorials(n1);
    std::vector<long double> inv_factorials(n2);

    factorials[0] = 1;
    inv_factorials[0] = 1;

    for (unsigned int i = 1; i < n1; i++)
      factorials[i] = factorials[i-1] * i;

    for (unsigned int i = 1; i < n2; i++)
      inv_factorials[i] = inv_factorials[i-1] / ((double)i);

    auto binom = [&factorials, &inv_factorials]
            (const int n, const int k) -> double{
      if (k == 0 || k == n)
        return 1;
      return factorials[n] * inv_factorials[k] * inv_factorials[n - k];
    };

    // the update coefficients for the controlpoints in each chunk will
    // be the same, so we compute them offline:
    std::vector< std::vector< double > > bezier_alphas(p[d] + 1);
    for (int i = 0; i < (int)p[d] + 1; i++){
      int j_min = std::max(0, (int)i - (int)times);
      int j_max = std::min((int)p[d] - (int)times, (int)i);

      // Calculate the coefficients:
      std::vector< double > alpha_i(j_max - j_min + 1);
      long double div = 1. / binom(p[d], i);
      for (int j = j_min; j < j_max + 1; j++)
        alpha_i[j - j_min] = binom(p[d] - times, j) * binom(times , i - j) * div;

      bezier_alphas[i] = alpha_i;
    }

    std::vector< ts_ptr > new_ts;

    // The coefficients for the updates are ready, run the update on the splines chunk wise:
    for (const auto& chunk : sorted){
      // Save the old control points:
      std::vector< ControlPoint > old_cps(n_splines_per_chunk);
      for (unsigned int i = 0; i < n_splines_per_chunk; i++)
        old_cps[i] = chunk[i] -> get_cp();

      // initialize the vector of new control points
      std::vector< ControlPoint > new_cps(n_new_splines);

      // save some variables to be used later
      int m = n_old_knots - 1;
      int ph = p[d];

      // From Algorithm 5.9 we have the bezier coefficients alredy defined above

      // Further variables:
      int mh               = ph;
      int kind             = ph + 1;
      int r                = -1;
      int a       = p[d] - times;
      int b       = a + 1;
      int cind    = 1;
      double ua   = global_kv_old[0];
      new_cps[0]  = old_cps[0];

      // This array contains the control points for one bezier segment
      std::vector< ControlPoint > bpts(p[d] - times + 1);
      for (unsigned int i = 0; i < p[d] - times + 1; i++) bpts[i] = old_cps[i];

      // This array stores controlpoints for the next iteration, i.e.
      // the next bezier segment
      std::vector< ControlPoint > nextbpts(p[d] - times - 1);

      // This array stores the control points for the elevated bezier segment
      std::vector< ControlPoint > ebpts(p[d] + 1);

      // Big loop through old knot vector
      while (b < m){
        unsigned int i = b;
        while (b < m && global_kv_old[b] == global_kv_old[b+1]) b++;

        int mul   = b - i + 1;
        mh        = mh + mul + times;
        double ub = global_kv_old[b];
        int oldr  = r;
        r         = p[d] - times - mul;

        /* Insert knot u(b) r times  */
        double lbz = oldr > 0 ? (oldr + 2)/2 : 1;
        double rbz = r    > 0 ? ph - (r+1)/2 : ph;

        if (r > 0){
          // Insert knot to get Bezier segment:
          double numer = ub - ua;
          std::vector< double > alfs(p[d] - times - mul + 1);
          for (int k = p[d]-times; k > mul; k--)
            alfs[k-mul-1] = numer / (global_kv_old[a+k] - ua);

          for (int j = 1; j < r+1; j++){
            int save = r-j, s = mul+j;
            for (int k = (int)(p[d] - times); k >= s; k--)
              bpts[k] = alfs[k-s]*bpts[k] + (1. - alfs[k - s])*bpts[k-1];

            nextbpts[save] = bpts[p[d] - times];
          }
        } // if (r > 0)

        // Degree elevate Bezier curve
        for (int i = lbz; i < (int)ph + 1; i++){
          ebpts[i] = 0*ebpts[i];
          int j_max = std::min((int)(p[d] - times), i);
          int j_min = std::max(0, (int)i - (int)times);
          for (int j = j_min; j < j_max + 1; j++)
            ebpts[i] = ebpts[i] + bezier_alphas[i][j-j_min]*bpts[j];
        } // for (int i)

        if (oldr > 1){ // Must remove knot ua oldr times
          unsigned int first = kind - 2, last = kind;
          double den = ub - ua;
          double bet = (ub - global_kv.at(kind-1)) / den;
          for (int tr = 1; (int)tr < oldr; tr++){ // knot removal loop
            int i = first, j = last, kj = j - kind + 1;
            // Loop and compute the new control points for
            // one removal step
            while ( j-i > tr){
              if ( i < cind ){
                double alf = (ub - global_kv.at(i)) / (ua - global_kv.at(i));
                new_cps.at(i) = alf * new_cps.at(i) + (1. - alf)*new_cps.at(i-1);
              } // if (i < cind)

              if ( j >= lbz ){
                if ( j - tr <= (int)(kind) - (int)(ph) + (int)(oldr)){
                  double gam =  (ub - global_kv.at(j - tr)) / den;
                  ebpts.at(kj) = gam * ebpts.at(kj) + (1.-gam)*ebpts.at(kj+1);
                } else {
                  ebpts.at(kj) = bet * ebpts.at(kj) + (1.-bet)*ebpts.at(kj+1);
                }// if (j - r <= ...)
              } // if (j >= lbz)
              i++; j--; kj--;
            } // while (j - i)
            first--; last++;
          } // for (int tr)
        } // if (oldr)

        if (a != (int)(p[d] - times))
          kind += p[d] - oldr;

        // Load control points
        for (unsigned int j = lbz; j <= rbz; j++)
          new_cps.at(cind++) = ebpts.at(j);

        // Setup for next iteration through while loop
        if (b < m){
          for (int j = 0; j < r; j++)
            bpts[j] = nextbpts[j];
          for (int j = r; j < (int)(p[d] - times + 1); j++)
            bpts.at(j) = old_cps.at(b - (p[d] - times) + j);
          a = b; b++; ua = ub;

        } // if (b < m)
      } // while ( b )
  
      // by the end of the while loop, the new Control points should be stored in new_cps.
      // The last step is then to generate new TSplines from those:
      const auto& local_kv = chunk[0] -> get_kv();
      for (unsigned int i = 0; i < n_new_splines; i++){
        std::vector< double > local_kv_d(p[d] + 2);
        for (unsigned int j = 0; j < p[d] + 2; j++)
          local_kv_d[j] = global_kv[i + j];

        auto local_kv_i = local_kv;
        local_kv_i[d] = local_kv_d;

        new_ts.push_back(std::make_shared<TSpline>(local_kv_i, new_cps[i]));
      } // for(int i)
    }// for chunk

    // add new Splines to Collection:
    active_splines = new_ts;
  
    // Sort the set of active splines with respect to its barycentric coordinates:
    const unsigned int d2 = (dimension == 2 ? 1 : 2), d1 = 1, d0 = 0;
    std::sort(active_splines.begin(), active_splines.end(),
        [&d2, &d1, &d0](const auto& t1, const auto& t2) {
          if (t1 -> get_barycenter(d2) != t2 -> get_barycenter(d2))
            return t1 -> get_barycenter(d2) < t2 -> get_barycenter(d2);
          else if (t1 -> get_barycenter(d1) != t2 -> get_barycenter(d1))
            return t1 -> get_barycenter(d1) < t2 -> get_barycenter(d1);
          else
            return t1 -> get_barycenter(d0) < t1 -> get_barycenter(d0);
        });

    // and set new level indices:
    this -> n_splines = active_splines.size();
    for (unsigned int i = 0; i < n_splines; i++)
      active_splines[i] -> set_level_index(i);

    // set_boundary_dofs();

    // This needs to be overwritten, since the degrees for TSplines is a static variable
    IsoparametricFunction<dimension, space_dimension> new_ipf(active_splines);
    IPF = new_ipf;
  #ifdef DEBUG
    std::cout << indent(1) << "... done!" << std::endl;
  #endif
  } // degree_elevate()

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::set_boundary_dofs(
  ) {
    boundary_dofs.clear();
    const unsigned int nvf = GeometryInfo<dimension>::vertices_per_face;
    const unsigned int d2 = (dimension == 2 ? 1 : 2), d1 = 1, d0 = 0;
    for (const auto& cell : this -> active_cell_iterators()){
      for (const auto& face : cell -> face_iterators()){
        if (face -> at_boundary()){
          // std::cout << "On cell " << cell -> level() << "." << cell -> index()
          //           << " at face F" << face -> index() << " = ";
          const Point<dimension>& lower = face -> vertex(0);
          const Point<dimension>& upper = face -> vertex(nvf - 1);
          const Point<dimension>& diff = (-1.)*lower + upper;

          // for (unsigned int i = 0; i < dimension; i++){
          //   printf("(%+1.4e, %+1.4e)", lower(i), upper(i));
          //   if (i < dimension - 1)
          //     std::cout << " x ";
          // }
          // std::cout << std::endl;

          // find the orientation of the current face
          unsigned int d = 0;
          for (; diff(d) != 0; d++);

          // and use it to find the spline whose barycenter is closest to that
          // point
          for (const auto& t : active_splines){
            // std::cout << indent(1) << "Checking spline " << t -> get_level_index() << std::endl;
            // std::cout << indent(2) << "Anchor: ";
            bool anchor_on_face = true;
            const std::pair< Point<dimension>, Point<dimension> >& t_anchor = t -> get_anchor();
            // for (unsigned int i = 0; i < dimension; i++){
            //   if (t_anchor.first(i) == t_anchor.second(i))
            //     printf("{%+1.4e}", t_anchor.first(i));
            //   else 
            //     printf("(%+1.4e, %+1.4e)", t_anchor.first(i), t_anchor.second(i));

            //   if (i < dimension - 1)
            //     std::cout << " x ";
            // }
            // std::cout << std::endl << indent(2);
            anchor_on_face = (t_anchor.first(d0) >= lower(d0)
                                && t_anchor.second(d0) <= upper(d0) ) &&
                             (t_anchor.first(d1) >= lower(d1)
                                && t_anchor.second(d1) <= upper(d1) ) &&
                             (t_anchor.first(d2) >= lower(d2)
                                && t_anchor.second(d2) <= upper(d2) );

            // if (anchor_on_face)
            //   std::cout << "anchor on face.";
            // else 
            //   std::cout << "anchor not on face.";

            // std::cout << std::endl;


            // Is the anchor actually on the face?
            // If not, simply continue.
            if (!anchor_on_face)
              continue;

            // Otherwise, check if the spline is actually a boundary spline
            // Note: This is not a redundant check. E.g. a spline with kv
            // {0, 0, 0, 1, 1} x {0, 0, 1, 1} has its anchor technically on the face
            // F = {0} \times (0, 1), however, it is not a spline on the
            // boundary of this triangulation.
            std::vector< int > mult = t -> get_multiplicities(d);
            unsigned int n_unique_knots = mult.size();
            if (mult.at(0) == (int)p[d] + 1 || mult.at(n_unique_knots-1) == (int)p[d] + 1)
              this->boundary_dofs[face -> boundary_id()].push_back(t -> get_level_index());
          } // for ( t )
        } // if ( at_boundaray )
      } // for ( face )
    } // for ( cell )

    for (auto& [_, dofs] : boundary_dofs){
      std::sort(dofs.begin(), dofs.end());
      auto it = dofs.begin();
      for (; it != (dofs.end()-1); )
        if (*it == *(it+1))
          it = dofs.erase(it);
        else
          it++;
    }
  } // set_boundary_dofs()

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::refine_bezier_elements(
  ) {
    Assert(!is_bezier_mesh, ExcInternalError());
#ifdef DEBUG
    for(const auto& cell : this->active_cell_iterators())
      Assert(!(cell->refine_flag_set()), ExcMessage("Before calling refine_bezier_elements(), make sure no other cell is marked for refinement."));
#endif
    
    typename
    std::vector< cell_iterator >::iterator it_bezier = bezier_elements.begin();
    for (; it_bezier != bezier_elements.end(); ){
      if ((*it_bezier) -> has_children()){
        // This case is not handled in the search for bezier-elements
        bezier_elements.erase(it_bezier);
      } else {
        (*it_bezier) -> set_refine_flag(
                        RefinementCase<dimension>::cut_axis(
                                (level_offset + (*it_bezier) -> level()) % dimension)
                                       );
        ++it_bezier;
      }
    }

    // Refine cells without producing new splines:
    this -> Triangulation<dimension>::execute_coarsening_and_refinement();

    // // afterwards, update mof to include new inner faces
    // for (it_bezier = bezier_elements.begin();
    //         it_bezier != bezier_elements.end();
    //         ++it_bezier) {
    //   const auto& cell = *it_bezier;
    //   const auto& ref_case = static_cast<std::uint8_t>(cell -> refinement_case());
    //   const auto& child = cell -> child(0); 
    //   unsigned int face_no = ref_case == 1 /* cut_x */ ? 1 : 
    //                             ref_case == 2 /* cut_y */ ? 3 : 5;
    //   const auto& face = child -> face(face_no); 
    //   mof[face->index()] = 1;

    //   for (const auto& f : cell -> face_iterators()){
    //     if (f -> has_children()){
    //       mof[f -> child(0)->index()] = mof[f->index()];
    //       mof[f -> child(1)->index()] = mof[f->index()];
    //     }
    //   }
    // }

    is_bezier_mesh = true;
    this -> setup_mof();
    this -> compute_IEN_array();
  }

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::coarsen_bezier_elements(
  ){
    Assert(is_bezier_mesh, ExcInvalidState());
#ifndef DEBUG
    if (!is_bezier_mesh)
      throw ExcInvalidState();
#endif

    // Switch mesh indicator to normal grid
    is_bezier_mesh = false;

    // Force recomputation of IEN_array after coarsening
    IEN_array.clear();

    // in case there were no bezier elements to begin with
    if (bezier_elements.size() == 0)
      return;

  #ifdef DEBUG
    for(const auto& cell : this->active_cell_iterators()){
      Assert(!(cell->coarsen_flag_set()), ExcMessage("Before calling coarsen_bezier_elements(), make sure no other cell is marked for coarsening."));
      Assert(!(cell->refine_flag_set()), ExcMessage("Before calling coarsen_bezier_elements(), make sure no other cell is marked for refinement."));
    }
  #endif

    for (auto& bezier_cell : bezier_elements){
      bezier_cell -> child(0) -> set_coarsen_flag();
      bezier_cell -> child(1) -> set_coarsen_flag();
    }

    this -> Triangulation<dimension>::execute_coarsening_and_refinement();
  }

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::compute_extraction_operators(
  ){
    // This function assumes that the Bezier elements are already refined
  #ifdef DEBUG
    Assert(is_bezier_mesh, ExcInvalidState());
  #else
   if (!is_bezier_mesh)
     throw ExcInvalidState();
  #endif

#ifdef DEBUG
    const std::string name_operators = "log/extraction_operators.txt";
    const std::string name_log       = "log/bezier_coefficients.txt";
    std::ofstream log(name_log, std::ios::out | std::ios::trunc);
    std::ofstream out(name_operators, std::ios::out | std::ios::trunc);
#endif


   extraction_operators.clear();
   const unsigned int nvc = GeometryInfo<dimension>::vertices_per_cell;
   const unsigned int nvf = GeometryInfo<dimension>::vertices_per_face;
   const unsigned int nfc = GeometryInfo<dimension>::faces_per_cell;
   for (const auto& bezier_cell : this -> active_cell_iterators()){
#ifdef DEBUG
     out << std::endl;
     out << " ============================================================== "
         << std::endl;
     out << " >>>>>>>>>>>>>>> Level: "
         << this -> n_levels() - this -> level_offset
         << " <<<<<<<<<<<<<<<" << std::endl;
     out << "On cell " << bezier_cell->level() << "." << bezier_cell -> index()
         << " with bounds ("
         << bezier_cell->vertex(0) << ") x (" 
         << bezier_cell->vertex(nvc-1) << ")" << std::endl;
     out << " ============================================================== "
         << std::endl << std::endl;
#endif
     const unsigned int sz = IEN_array.at(bezier_cell).size();
     
     // Get a reference to the local operator
     std::vector< std::array< Vector<double>, dimension> >& local_extraction_operator
             = extraction_operators[bezier_cell];
     local_extraction_operator.resize(sz);

     // Get Splines from IEN_array
     unsigned int i = 0;
     for (const auto& ts_ind : IEN_array[bezier_cell]){
       const auto& ts = active_splines[ts_ind];
       // check if this cell splits the support into new elements:
       const Point<dimension>& lower = bezier_cell -> vertex(0);
       const Point<dimension>& upper = bezier_cell -> vertex(nvc - 1);

#ifdef DEBUG
       out << "Extraction operator for ts " << ts_ind << " given as: " << std::endl;
#endif
       // For each dimension, setup the extraction rows from the current TSpline
       const auto& kv = ts -> get_kv();
       for (unsigned int d = 0; d < dimension; d++){
         std::vector<double> interior_knots;
         for (unsigned int j = 0; j < p[d] + 1; j++){
           // check if either lower or upper is inside a knot span
           if (kv[d][j] < lower(d) && lower(d) < kv[d][j+1]){
             interior_knots.push_back(lower(d));
           } else if (kv[d][j] < upper(d) && upper(d) < kv[d][j+1]){
             interior_knots.push_back(upper(d));
           }
         }

         // Get the row of this T-spline:
         local_extraction_operator[i][d]
             = ts -> compute_be_operator(interior_knots, d, bezier_cell);
#ifdef DEBUG
         out << "kv: " << "[ ";
         for (const auto& knot : kv[d])
           out << knot << " ";
         out << "], with interior knot [";
         for (const auto& k : interior_knots)
           out << k << " ";
         out << "], yields [";
         for (const auto& c : local_extraction_operator[i][d])
           out << c << " ";
         out << "]" << std::endl;
#endif
       }
       i++;
     } // for ( ts_ind )

     // Set the boundary indicators:
     face_operators.clear();
     for (unsigned int f = 0; f < nfc; f++){
       i = 0;
       auto& local_face_operators = face_operators[bezier_cell];
       if (!bezier_cell -> face(f) -> has_children()) {
         // copy local cell operators:
         local_face_operators[bezier_cell -> face(f) -> index()]
             = local_extraction_operator;
       } else {
         // To match Quadrature points for jumps on faces, compute
         // the extraction operators on this cell at specified face
         // with an interior knot in direction of refinement


         // Firstly, get direction of refinement. The refined face comes
         // from the neighbor of this current face.
         const auto& ref_case = static_cast<std::uint8_t>( 
                         bezier_cell -> neighbor(f) -> refinement_case() );
         Assert(ref_case == 1 || ref_case == 2 || ref_case == 4, ExcInternalError());
         const unsigned int dr = (ref_case == 1 ? 0 : (ref_case == 2 ? 1 : 2));
         const unsigned int d1 = (dr == 0 ? 1 : 0);
         const unsigned int d2 = (dimension == 2 ? d1 : (dr == 0 ? 2 : (dr == 1 ? 2 : 1)));
         const Point<dimension>& face_center = 0.5 * (
                                    bezier_cell -> face(f) -> vertex(0)
                                  + bezier_cell -> face(f) -> vertex(nvf - 1)
                                  );

         // Instead of saving the extraction operators for the parent face,
         // we instead split it into two parts, namely its children, and
         // store the respective data for the two children
         const auto& face0 = bezier_cell -> face(f) -> child(0);
         const auto& face1 = bezier_cell -> face(f) -> child(1);

         // set the size of the local extraction operator on this face:
         local_face_operators[face0 -> index()].resize(sz);
         local_face_operators[face1 -> index()].resize(sz);

         // We will proceed as above, however, for the direction of non-refinement,
         // i.e. d1, we can copy the extraction operators of the cell for the two sub-faces
         for (const auto& ts_ind : IEN_array[bezier_cell]){
           const auto& ts = active_splines[ts_ind];
           // check if this cell splits the support into new elements:
           const Point<dimension>& lower = bezier_cell -> vertex(0);
           const Point<dimension>& upper = bezier_cell -> vertex(nvc - 1);

           const auto& kv = ts -> get_kv();
           std::vector<double> interior_knots;
           for (unsigned int j = 0; j < p[dr] + 1; j++){
             // check if either lower or upper is inside a knot span
             if (kv[dr][j] < lower(dr) && lower(dr) < kv[dr][j+1]){
               interior_knots.push_back(lower(dr));
             } else if (kv[dr][j] < upper(dr) && upper(dr) < kv[dr][j+1]){
               interior_knots.push_back(upper(dr));
             }
           }

           interior_knots.push_back(face_center(dr));
           std::sort(interior_knots.begin(), interior_knots.end());

           // Get the row of this T-spline: compute it anew for the two children
           // of the current face
           std::vector< active_cell_iterator > children(2);
           children[0] = bezier_cell -> neighbor(f) -> child(0);
           children[1] = bezier_cell -> neighbor(f) -> child(1);
           std::vector< Vector<double> > ops(2);
           ts -> compute_be_operator(
                   interior_knots,
                   dr,
                   children,
                   ops);

           local_face_operators[face0 -> index()][i][dr] = ops[0];
           local_face_operators[face1 -> index()][i][dr] = ops[1];


           // And copy the data from beefore for the non-refined direction
           local_face_operators[face0 -> index()][i][d1] =
                   local_extraction_operator[i][d1];
           local_face_operators[face1 -> index()][i][d1] =
                   local_extraction_operator[i][d1];

           if (dimension != 2) {
             local_face_operators[face0 -> index()][i][d2] =
                     local_extraction_operator[i][d2];
             local_face_operators[face1 -> index()][i][d2] =
                     local_extraction_operator[i][d2];
           }
           i++;
         } // for ( ts_ind )
       }
      } // for ( f )
    } // for ( bezier_cell )
#ifdef DEBUG
    // Save the matrix to a seperate file after value are computed:
    for (const auto& bezier_cell : this -> active_cell_iterators()){
      log << std::endl;
      log << " ============================================================== "
          << std::endl;
      log << " >>>>>>>>>>>>>>> Level: "
          << this -> n_levels() - level_offset
          << " <<<<<<<<<<<<<<<" << std::endl;
      log << "On cell " << bezier_cell->level() << "." << bezier_cell -> index()
          << " with bounds ("
          << bezier_cell->vertex(0) << ") x (" 
          << bezier_cell->vertex(nvc-1) << ")" << std::endl;
      log << " ============================================================== "
          << std::endl << std::endl;
     
      const FullMatrix<double>& C = get_bezier_coefficients(bezier_cell); 
      const auto& I = get_IEN_array(bezier_cell);
      for (unsigned int i = 0; i < C.n(); i++){
        log << I[i] << ": ";
        for (unsigned int j = 0; j < C.m(); j++){
          log << C(i, j) << " "; 
        }
        log << std::endl;
      }

    } // for (cell)
    out.close();
    log.close();
#endif
  } // compute_extraction_operators()

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::test_bezier(
  ){
#ifdef DEBUG
    std::cout << indent(3) << "Performing test of Bezier mesh ..." << std::endl;

    // Refine bezier elements to obtain bezier mesh
    refine_bezier_elements();

    // Calculate the amount of TSplines with support on a beyier cell
          unsigned int n_splines = 1;
    const unsigned int nvc = GeometryInfo<dimension>::vertices_per_cell;
    for (unsigned int d = 0; d < dimension; d++)
      n_splines *= (p[d] + 1);

    // Perform the test
    bool is_success = true;
    for (const auto& cell : this->active_cell_iterators()){
      if (!is_success)
        break;

      // Count the splines with support on this cell. Here, a tspline
      // has support on a cell, if each vertex of the cell is supported
      unsigned int n_splines_on_cell = 0;
      for (const auto& ts : active_splines){
        bool has_support = true;
        for (unsigned int n = 0; n < nvc && has_support; n++)
          has_support = has_support && (ts -> has_support(cell -> vertex(n)));

        if (has_support)
          n_splines_on_cell++;
      }

      is_success = (n_splines_on_cell == n_splines);
    } // for (cell)

    if (is_success)
      std::cout << indent(3) << "... test succesfull!" << std::endl;
    else
      std::cout << indent(3) << "... test failed!" << std::endl;

    Assert(is_success, ExcInternalError());
    coarsen_bezier_elements();
#endif
  } // test_bezier()
  // =============================================================================
  //                Helper section 
  // =============================================================================

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::compute_IEN_array(
  ){
    Assert(is_bezier_mesh, ExcInvalidState())
#ifndef DEBUG
    if (!is_bezier_mesh)
      throw ExcInvalidState();
#endif
    const std::vector< unsigned int >& p         = this->p;
    const unsigned int splines_per_cell = (p[0] + 1) * (p[1] + 1) * (p[2] + 1);

    // Reset current IEN_array
    const unsigned int nvc = GeometryInfo<dimension>::vertices_per_cell;
    IEN_array.clear();
    for (auto bezier_cell : this -> active_cell_iterators()){
      std::vector<unsigned int> local_splines;
      for (const auto& ts : active_splines){
        // has_support in its current state is a little flawed, but works
        // for the other purposes right now. The problem is, that it returns
        // true on cells, where only one face, or even one point lies in the
        // the support of the TSpline, which is usually fine. However,
        // with bezier grids, each cell is completely contained inside the
        // supports of multiple splines. Thus, we check if each point of a cell
        // lies in the support of the TSplines
        bool support = true;
        for (unsigned int n = 0; n < nvc; n++)
          support = support && (ts -> has_support(bezier_cell -> vertex(n)));

        if (support)
          local_splines.push_back(ts -> get_level_index());
      } // for ( ts )

      if ( splines_per_cell != local_splines.size() ){
        std::string message = "An internal check noticed that the number of T-splines on cell "
                            + std::to_string(bezier_cell->level()) + "." + std::to_string(bezier_cell -> level()) 
                            + " \ngiven by " + std::to_string(local_splines.size()) + " differs from the "
                            + "assumed amunt of T-splines on any cell given as " + std::to_string(splines_per_cell)
                            + ".\n\n"
                            + "There is not very much you can do if you encounter this error\n"
                            + " since it indicates an error in deal.T, not in your own program. \n"
                            + "Try to come up with the smallest possible program that still \n"
                            + "demonstrates the error and contact us with it to obtain help.";
        throw ExcMessage(message);
      }
      Assert(splines_per_cell == local_splines.size(), ExcInternalError());
      IEN_array[bezier_cell] = local_splines;
    } // for ( bezier_cell )
  } // compute_IEN_array()

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::get_affected_splines(
    const std::vector<active_cell_iterator>& coarse_neighborhood,
    std::vector< ts_ptr >& to_be_updated
  ) const {
    const unsigned int nvc = GeometryInfo<dimension>::vertices_per_cell;
    for (const auto& cell : coarse_neighborhood){
      const Point<dimension>& lower = cell -> vertex(0);
      const Point<dimension>& upper = cell -> vertex(nvc - 1);

      // Consider the following sketch
      //
      //                                      +
      //                                 d1  /|
      //                                  ^ / |
      //                                  |/  |
      //                                  +   | d2
      //                        dr <----- |   |/
      //  z                               |   +
      //  ^  y                            |  /
      //  | /                             | /
      //  |/                              |/
      //  +--->x                          +
      //
      // The considered face derives from a cut in the yz-plane. We then need to find
      // the TSplines along the orthogonal direction of the cut (i.e. dr = x) and
      // check if the anchor is contained within the new face. To not consider
      // anchor entities in this part of the code, we use the fact that
      //    floor( (p[d]+1)/2 ) = ceil( (p[d]+1)/2 ), if p[d] mod 2 = 1,
      // and we simply use the diirections of the cut, hence
      //    center_dir = {{d1, d2}} = {{1, 2}}
      //
      // Now, since some cells may be refined in two directions, we save center_dir
      // as a vector<vector<int>>, doing the above strategy to each direction
      // orthogonal to the plane of the cut
      const std::uint8_t ref_case = static_cast<std::uint8_t>(cell -> refine_flag_set());
      const unsigned int dr = (ref_case == 1 ? 0 /* refined in x */: 
                                (ref_case == 2 ? 1 /* refined in y */: 
                                  (ref_case == 4 ? 2 /* refined in z */: 
                                     0 /* will never be reached */)));
      const unsigned int d1 = (ref_case == 1 ? 1 : 0);
      const unsigned int d2 = (dimension == 2 ? d1 : 
                                (ref_case == 4 ? 1 : 2));

      // Find the active splines, that have support on the current cell,
      // and are aligned with the cell in the corresponding direction
      for (auto& ts : active_splines){
        if (!(ts -> has_support(cell)))
          continue;

        const auto& kv = ts -> get_local_kv(dr);
        const double& k1r = kv[0];
        const double& k2r = kv[p[dr] + 1];
        const double& in  = .5*(lower(dr) + upper(dr));
        if ( k1r > in || k2r < in )
          continue;

        // is the anchor of the current TSpline bigger then the inserted face?
        const double& k11 = ts -> get_anchor(d1).first;
        const double& k12 = ts -> get_anchor(d1).second;
        const double& k21 = ts -> get_anchor(d2).first;
        const double& k22 = ts -> get_anchor(d2).second;
        if (     (k11 >= lower(d1) && k12 <= upper(d1)) 
              && (k21 >= lower(d2) && k22 <= upper(d2)))
          to_be_updated.push_back( ts );
      } // for ts
    } // for cell

    // remove redundancies in to_be_updated, by sorting the splines by its index and removing
    // multiple occurances of an index
    std::sort(to_be_updated.begin(), to_be_updated.end(), [](const auto& t1, const auto& t2){
                      return t1 -> index() < t2 -> index();
                    });

    typename std::vector< ts_ptr >::iterator ts_it = to_be_updated.begin();
    for (; ts_it != to_be_updated.end() - 1; )
      if ( (*ts_it) -> index() == (*(ts_it + 1)) -> index() )
        to_be_updated.erase(ts_it);
      else
        ++ts_it;
  } // get_affected_splines [1 / 2]

  template<int dim, int spacedim>
  auto TS_TriangulationBase<dim, spacedim>::get_affected_splines(
    const active_face_iterator&   face,
    const std::vector<ts_ptr>&    to_be_updated,
          unsigned int&           dr
  ) const -> std::vector<std::vector< ts_ptr >> {
    const unsigned int nvf = GeometryInfo<dimension>::vertices_per_face;
    const Point<dimension>& lower = face -> vertex(0);
    const Point<dimension>& upper = face -> vertex(nvf - 1);

    // get the direction orthogonal to the face:
    const Point<dimension>& dif = face -> vertex(nvf-1) + ((-1.)*face -> vertex(0));
    for (dr = 0;  dr < (int)dimension && dif(dr) != 0; dr++);

    const unsigned int d1 = (dr == 0 ? 1 : 0), d2 = (dimension == 2 ? d1 : (dr == 2 ? 1 : 2));

    const double& knot_in = lower(dr);
    Assert(knot_in == upper(dr), ExcInternalError());

    std::vector< ts_ptr > selected;
    for (auto& ts : to_be_updated){
      if ( !(ts -> is_active()) )
        continue;
      // As before we need to check, that the tspline of interest has indeed support on the
      // current face / cell. Further, consider the following example
      //
      //    T3        T4        T5
      //    + - - - - + - - - - +
      //    |         |         |
      //    |   F1    |    F2   |
      //    * = = = = * = = = = *
      //    |         |         |
      //    |         |         |
      //    + - - - - + - - - - +
      //    T0        T1        T2
      //
      // where px = py = odd, i.e. anchors are on vertices. We update the TSplines
      // face by face. For face F1 (the left line marked by = ), we then iterate
      // over all TSplines and check if their local knot vectors fulfill the criterions
      // already mentioned in the function above, i.e.
      //
      //            anchor(Ti) \cap F1 \neq \emptyset
      //
      // This will give us the Splines T0, T1, T3, and T4 which will be updated along
      // with its corresponding control points. For details see the function
      // execute_coarsening_and_refinement(). The selected splines will then
      // be set in_active and new splines T6, T7, T8, and T9, as well as T10 and T11
      // will be generated as below
      //
      //    T8        T9        T5
      //    + - - - - + - - - - +
      //    |         |         |
      //    |T10      |T11 F2   |
      //    + - - - - + = = = = *
      //    |         |         |
      //    |         |         |
      //    + - - - - + - - - - +
      //    T6        T7        T2
      //
      // where The splines T6 ... T9 get the old splines as fathers. We now want to
      // proceed to update all the splines affected when we insert face F2.
      //
      // The TSplines we would select would be T2, T5, T7, T9, and T11. But T7, T9, and
      // T11 already contain the new knot that we want to insert by cutting the cell,
      // and thus need no further updates from face F2.
      //
      // To remove them from selection, we simply check if the knot [that we want to insert]
      // is alredy given with the is_present(double, int) method.
      if ( !ts -> has_support(face) || ts -> is_present(knot_in, dr) )
        continue;

      // is the anchor of the current TSpline bigger then the inserted face?
      const double& k11 = ts -> get_anchor(d1).first;
      const double& k12 = ts -> get_anchor(d1).second;
      const double& k21 = ts -> get_anchor(d2).first;
      const double& k22 = ts -> get_anchor(d2).second;
      if (     (k11 >= lower(d1) && k12 <= upper(d1)) 
            && (k21 >= lower(d2) && k22 <= upper(d2)))
        selected.push_back(ts);
    } // for ts

    // Sort the selected tsplines in a proper manner, depending on d1, d2, and dr
    // for each dr we desire a "consecutive numbering" in order to sort them
    // by dr
    std::sort(selected.begin(), selected.end(),
                    [d1, d2, dr](
                      const auto& t1, 
                      const auto& t2
                    ) {
                      if (t1 -> get_barycenter(d2) != t2 -> get_barycenter(d2))
                        return t1 -> get_barycenter(d2) < t2 -> get_barycenter(d2);
                      else if (t1 -> get_barycenter(d1) != t2 -> get_barycenter(d1))
                        return t1 -> get_barycenter(d1) < t2 -> get_barycenter(d1);
                      else
                        return t1 -> get_barycenter(dr) < t2 -> get_barycenter(dr);
                    });

    // Sort the TSplines
    const unsigned int n_splines = selected.size(), n_out = n_splines/(p[dr] + 1);
          unsigned int n = 0;
    std::vector< std::vector< std::shared_ptr<TSpline> > > out(n_out);
    for (unsigned int i = 0; i < n_out; i++){
        out[i] = std::vector< ts_ptr >(selected.begin() + n, selected.begin() + n + p[dr] + 1);
        n += p[dr] + 1;
    }

    return out;
  } // get_affected_splines()

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::ensure_open_knot_vector(
      std::vector< double >& kv,
      const unsigned int p
  ) const {
    unsigned int k, n = kv.size();
    // Assert that the knot vector is non-decreasing:
#ifdef DEBUG
    for (k=0; k<n-1; k++)
      Assert( kv[k] <= kv[k+1], ExcInvalidState());
#endif
    
    Assert( kv[0] < kv[n-1], ExcInvalidState());

    // First knot:
    for(k=0; kv[k]==kv[k+1]; k++);
    // The first knot has multiplicity k+1.  Aim is multiplicity p+1.
    auto it = kv.begin();
    if(k>p) kv.erase( it, it + k-p-1 ); // erase superfluous entries
    if(k<p) kv.insert( it, p-k-1, kv[0] ); // insert missing entries
  
    // Last knot:
    for(k=n-1; kv[k]==kv[k-1]; k--);
    k = n-1-k;
    // The last knot has multiplicity k+1.  Aim is multiplicity p+1.
    it = kv.end()-1;
    if(k>p) kv.erase( it - (k-p-1), it ); // erase superfluous entries
    if(k<p) kv.insert( it, p-k-1, kv[n-1] ); // insert missing entries
  } // ensure_open_knot_vector()

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::merge_kv(
    const std::vector<double>& kv1,
          std::vector<double>& kv2,
    const unsigned int d
  ) const {
    for (const double& in: kv1){
      std::vector<double>::iterator it     = kv2.begin();
      const auto&                   it_end = kv2.end();

      for ( ; it != it_end && *it < in; ++it);

      unsigned int count = 0;
      for ( ; it != it_end && *it == in; ++it) count++;

      if (count < p[d] + (unsigned int)(in == kv_lower(d) || in == kv_upper(d)))
        kv2.insert(it, in);
    }
  } // merge_kv()

  template<int dim, int spacedim>
  auto TS_TriangulationBase<dim, spacedim>::get_inner_faces(
    const cell_iterator& cell
  ) -> const std::vector< face_iterator > {
          std::vector< face_iterator > inner_faces;
    const auto& childs = cell -> child_iterators();

    const unsigned int nvc = GeometryInfo<dimension>::vertices_per_cell;
    const unsigned int d0 = 0, d1 = 1, d2 = (dimension == 2 ? d1 : 2);
    const Point<dimension>& lower = cell -> vertex(0);
    const Point<dimension>& upper = cell -> vertex(nvc - 1);

    Assert(cell -> n_children() == 2, ExcInternalError());

    // There is only one face of interest, thus only a child is needed
    for (const face_iterator& face : cell -> child(0) -> face_iterators()){
      // Test whether the point is fully inside the cell
      // we cannot use point_inside, as it returns true for points
      // on the boundary ...
      const Point<dimension>& face_center = face->center();
      const bool inner = !(face_center(d0) == lower(d0) || face_center(d0) == upper(d0)) &&
                         !(face_center(d1) == lower(d1) || face_center(d1) == upper(d1)) &&
                         !(face_center(d2) == lower(d2) || face_center(d2) == upper(d2));
      // the newly inserted face generates a knot with multiplicity 1
      if (inner) 
        inner_faces.push_back(face);
    } // for ( face )

    Assert(inner_faces.size() == 1, ExcInternalError());

    return inner_faces;
  } // get_inner_faces()

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::setup_mof(
  ) {
    // clear contents of mof 
    mof.clear();
    const unsigned int nvf = GeometryInfo<dimension>::vertices_per_face;
    for (auto face = this -> begin_face() ;
              face != this -> end_face();
              ++face){

      const Point<dimension>& diff = face -> vertex(0) + (-1)*face->vertex(nvf - 1);
      const unsigned int d = (diff(0) == 0 ? 0 : (diff(1) == 0 ? 1 : 2));

      const double knot = (face -> vertex(0)).operator()(d);
      unsigned int i = 0;
      for (; i < base_knots[d].size() && base_knots[d][i] != knot; i++);

      if (i == base_knots[d].size())
        mof[face -> index()] = 1;
      else 
        mof[face -> index()] = multiplicities[d][i];
    }
  } // setup_mof()

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::vector_dist(
      const cell_iterator& cell,
      const Point<dim>& z,
            Point<dim>& distance
  ) const {
    const Point<dimension>& mid = cell -> center();
    for (unsigned int d = 0 ; d < dimension; d++) 
      distance(d) = std::abs(mid(d) - z(d));
  } // vector_dist() [1 / 2]

  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::vector_dist(
      const cell_iterator& cell1,
      const cell_iterator& cell2,
            Point<dim>& distance
  ) const {
    const Point<dimension>& mid1 = cell1 -> center();
    const Point<dimension>& mid2 = cell2 -> center();
    for (unsigned int d = 0 ; d < dimension; d++) 
      distance(d) = std::abs(mid1(d) - mid2(d));
  } // vector_dist [2 / 2]
  
#include "ts_triangulation.inst.in"

} // namespace dealt
