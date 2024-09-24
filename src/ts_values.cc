#include <ts_triangulation.h>

namespace dealt {
  template<int dim, int spacedim, int soldim>
  TSValuesBase<dim, spacedim, soldim>::TSValuesBase(
      const TS_TriangulationBase<dim, spacedim>* tria,
      const unsigned int                         n_gauss_points,
      const UpdateFlags                          flags,
      const int&                                 face_no
  ) : face_no(face_no)
    , flags(flags)
    , tria(const_cast<TS_TriangulationBase<dimension, space_dimension>*>(tria)) {
    Assert(dimension < 4, ExcNotImplemented());
    Assert(dimension > 1, ExcNotImplemented());
    Assert(dimension <= space_dimension, ExcIndexRange(dimension, 1, space_dimension));
    Assert(solution_dimension == 1 || solution_dimension == space_dimension, 
                  ExcDimensionMismatch2(solution_dimension, 1, space_dimension));

    degrees = tria->get_degree();
    this -> n_gauss_points = std::vector< unsigned int >(dim, n_gauss_points);

    if ( dimension == 2){
      cell_quadrature = QAnisotropic<dimension>(
                          QGauss<1>(n_gauss_points),
                          QGauss<1>(n_gauss_points)
                        );
      face_quadratures[0] = QAnisotropic<dimension>(
                             QTrapezoid<1>(),
                             QGauss<1>(n_gauss_points)
                            );
      face_quadratures[1] = QAnisotropic<dimension>(
                             QGauss<1>(n_gauss_points),
                             QTrapezoid<1>()
                            );
      quadrature_points_per_face[0] = n_gauss_points;
      quadrature_points_per_face[1] = n_gauss_points;
      face_quadrature_points = 4 * quadrature_points_per_face[0]
                                 * quadrature_points_per_face[1];
    } else {
      cell_quadrature = QAnisotropic<dimension>(
                          QGauss<1>(n_gauss_points),
                          QGauss<1>(n_gauss_points),
                          QGauss<1>(n_gauss_points)
                        );
      face_quadratures[0] = QAnisotropic<dimension>(
                             QTrapezoid<1>(),
                             QGauss<1>(n_gauss_points),
                             QGauss<1>(n_gauss_points)
                            );
      face_quadratures[1] = QAnisotropic<dimension>(
                             QGauss<1>(n_gauss_points),
                             QTrapezoid<1>(),
                             QGauss<1>(n_gauss_points)
                            );
      face_quadratures[2] = QAnisotropic<dimension>(
                             QGauss<1>(n_gauss_points),
                             QGauss<1>(n_gauss_points),
                             QTrapezoid<1>()
                            );
      quadrature_points_per_face[0] = n_gauss_points * n_gauss_points;
      quadrature_points_per_face[1] = n_gauss_points * n_gauss_points;
      quadrature_points_per_face[2] = n_gauss_points * n_gauss_points;
      face_quadrature_points = 8 * quadrature_points_per_face[0]
                                 * quadrature_points_per_face[1]
                                 * quadrature_points_per_face[2];
    }

    init_bernstein_tables();
    set_bernstein_tables();
  } // TSValuesBase constructor from isotropic quadrature

  template<int dim, int spacedim, int soldim>
  TSValuesBase<dim, spacedim, soldim>::TSValuesBase(
      const TS_TriangulationBase<dim, spacedim>* tria,
      const std::vector<unsigned int>&           n_gauss_points,
      const UpdateFlags                          flags,
      const int&                                 face_no
  ) : face_no(face_no)
    , flags(flags)
    , n_gauss_points(n_gauss_points)
    , tria(const_cast<TS_TriangulationBase<dimension, space_dimension>*>(tria)) {
    Assert(dimension <= space_dimension, ExcIndexRange(dimension, 1, space_dimension));
    Assert(solution_dimension == 1 || solution_dimension == space_dimension, 
                  ExcDimensionMismatch2(solution_dimension, 1, space_dimension));
    degrees = tria->get_degree();
    if ( dimension == 2){
      cell_quadrature = QAnisotropic<dimension>(
                          QGauss<1>(n_gauss_points[0]),
                          QGauss<1>(n_gauss_points[1])
                        );
      face_quadratures[0] = QAnisotropic<dim>(
                             QTrapezoid<1>(),
                             QGauss<1>(n_gauss_points[1])
                            );
      face_quadratures[1] = QAnisotropic<dim>(
                             QGauss<1>(n_gauss_points[0]),
                             QTrapezoid<1>()
                            );
      quadrature_points_per_face[0] = n_gauss_points[1];
      quadrature_points_per_face[1] = n_gauss_points[0];
      face_quadrature_points = 2 * (quadrature_points_per_face[0] 
                                    + quadrature_points_per_face[1]);
    } else {
      cell_quadrature = QAnisotropic<dim>(
                          QGauss<1>(n_gauss_points[0]),
                          QGauss<1>(n_gauss_points[1]),
                          QGauss<1>(n_gauss_points[2])
                        );
      face_quadratures[0] = QAnisotropic<dim>(
                             QTrapezoid<1>(),
                             QGauss<1>(n_gauss_points[1]),
                             QGauss<1>(n_gauss_points[2])
                            );
      face_quadratures[1] = QAnisotropic<dim>(
                             QGauss<1>(n_gauss_points[0]),
                             QTrapezoid<1>(),
                             QGauss<1>(n_gauss_points[2])
                            );
      face_quadratures[2] = QAnisotropic<dim>(
                             QGauss<1>(n_gauss_points[0]),
                             QGauss<1>(n_gauss_points[1]),
                             QTrapezoid<1>()
                            );
      quadrature_points_per_face[0] = n_gauss_points[1] * n_gauss_points[2];
      quadrature_points_per_face[1] = n_gauss_points[0] * n_gauss_points[2];
      quadrature_points_per_face[2] = n_gauss_points[0] * n_gauss_points[1];
      face_quadrature_points = 2 * (quadrature_points_per_face[0] 
                                    + quadrature_points_per_face[1]
                                    + quadrature_points_per_face[2]);
    }

    init_bernstein_tables();
    set_bernstein_tables();

    // define shape_quadrature_points with the use of face_no:
    if (face_no != -1){
      // We have only quadrature points on some face defined by face_no/2
      shape_quadrature_points = quadrature_points_per_face[face_no/2]; 
    } else {
      // We have quadrature points in the whole cell. In this case
      // this is equivalent to bernstein_quadrature_points
      shape_quadrature_points = bernstein_quadrature_points;
    }
    total_quadrature_points = shape_quadrature_points; // may change later with subfaces
  } // TSValuesBase constructor from anisotropic quadrature

  template<int dim, int spacedim, int soldim>
  void TSValuesBase<dim, spacedim, soldim>::init_bernstein_tables(
  ) {
    TableIndices<2> bernstein_indices;
    bernstein_indices[0] = 1; bernstein_indices[1] = 1;
    for (unsigned int d = 0; d < dimension; d++) {
      bernstein_indices[1] *= (degrees[d] + 1);
      bernstein_indices[0] *= n_gauss_points[d];
    }
    

    // Set dofs_per_cell
    dofs_per_cell = bernstein_indices[1];

    // is a face number given? 
    // if so, then we want to calculate the tables 
    // pre-emptively for the faces. 
    if (face_no != -1) 
      bernstein_indices[0] = face_quadrature_points; 

    // initializa the number of quadrature points used for bernstein evaluation
    unsigned int &n_qpnts = bernstein_quadrature_points; 
    n_qpnts = bernstein_indices[0];

    // Set the values and gradient tables. These values will 
    // always be needed for the integration
    bernstein_values.reinit(bernstein_indices);
    bernstein_grads.reinit(bernstein_indices);
    // value.reinit(shape_indices);
    // grad.reinit(shape_indices);

    // Initialize tables according to UpdateFlags
    if (flags & update_hessians) {
      bernstein_hessians.reinit(bernstein_indices);
    //  hessian.reinit(shape_indices);
    }
    
    if (flags & update_quadrature_points) {
      // Initialize the table of quadrature points
      quadrature_points.reinit(n_qpnts);
      quadrature_weights.reinit(n_qpnts);
    //   mapped_quadrature_points.resize(n_qpnts);
    } 
    
    // if (flags & update_JxW_values) {
    //   dx.resize(n_qpnts);
    // }
    // 
    // if (flags & update_jacobians) {
    //   J.resize(n_qpnts);
    // }
    // 
    // if (flags & update_inverse_jacobians) {
    //   I.resize(n_qpnts);
    // }

    // if (flags & update_jacobian_grads) {
    //   H.resize(n_qpnts);
    //   HI.resize(n_qpnts);
    // }
  } // init_bernstein_tables

  template<int dim, int spacedim, int soldim>
  void TSValuesBase<dim, spacedim, soldim>::init_shape_tables(
    const unsigned int n_subfaces
  ){
    const UpdateFlags& flags = this -> flags;
    TableIndices<2> shape_indices;

    // Reset shape_quadrature_points
    if (face_no != -1)
      shape_quadrature_points = quadrature_points_per_face[face_no/2];

    shape_indices[1] = this->dofs_per_cell;
    shape_indices[0] = shape_quadrature_points;

    unsigned int& n_qpnts = total_quadrature_points;
    n_qpnts = n_subfaces * shape_indices[0];

    // reset values before reinitilizing
    shape_values.reset_values();
    shape_grads.reset_values();

    shape_values.reinit(shape_indices);
    shape_grads.reinit(shape_indices);

    // Initialize tables according to UpdateFlags
    if (flags & update_hessians) {
      shape_hessians.reset_values();
      shape_hessians.reinit(shape_indices);
    }
    
    // Initialize the table of quadrature points
    if (flags & update_quadrature_points) 
      mapped_quadrature_points.resize(n_qpnts);
    
    if (flags & update_JxW_values) 
      dx.resize(n_qpnts);
    
    if (flags & update_jacobians) 
      J.resize(n_qpnts);
    
    if (flags & update_inverse_jacobians)
      I.resize(n_qpnts);

    if (flags & update_jacobian_grads) {
      H.resize(n_qpnts);
      HI.resize(n_qpnts);
    }

    if (flags & update_normal_vectors) 
      normals.resize(n_qpnts);
  } // init_shape_tables()

  template<int dim, int spacedim, int soldim>
  void TSValuesBase<dim, spacedim, soldim>::set_bernstein_tables(
  ) {
    const std::vector<unsigned int>& p = this -> tria -> get_degree();
    // Everything is setup in size and dimension, we can fill these
    // tables now with values [on the reference cell]
    //
    // First: generate the base for AnisotropicPolynomials
    std::vector< std::vector< Polynomials::Polynomial<double> > >
            base(dim);
    for (unsigned int d = 0; d < dimension; d++)
      base[d] = generate_complete_bernstein_basis<double>(p[d]);

    std::vector<Point<dim>>  evals; 
    std::vector< double >    weights; 

    if (face_no == -1){
      evals = cell_quadrature.get_points(); 
      weights = cell_quadrature.get_weights(); 
    } else {
      for (unsigned int d = 0; d < dimension; d++){
        std::vector< Point<dimension> > quadrature_points = face_quadratures[dimension-d-1].get_points();
        std::vector< double >           quadrature_weights = face_quadratures[dimension-d-1].get_weights();
        std::vector< unsigned int >     quadrature_indices(quadrature_points.size()); 
      
        for (unsigned int i =0; i < quadrature_points.size(); i++)
          quadrature_indices[i] = i;
      
        // Sort quadrature points accordingly
        const unsigned int dh = dimension - d - 1; 
        const unsigned int d1 = (dh == 0 ? 1 : 0);
        const unsigned int d2 = (dimension == 2 ? d1 : 
                                  (d1 == 1 || dh == 1 ? 2 : 1));
        // std::cout << "d1 = " << d1 << ", d2 = " << d2 << std::endl;
        std::sort(quadrature_indices.begin(), quadrature_indices.end(), 
                        [&dh, &d1, &d2, &quadrature_points](const unsigned int i, const unsigned int j){
                          if (quadrature_points[i](dh) != quadrature_points[j](dh))
                            return quadrature_points[i](dh) < quadrature_points[j](dh);
                          else if (quadrature_points[i](d2) != quadrature_points[j](d2))
                            return quadrature_points[i](d2) < quadrature_points[j](d2);
                          else 
                            return quadrature_points[i](d1) < quadrature_points[j](d1);
                        });
      
        const std::vector< Point<dimension> >& tmp_quadratures = face_quadratures[dimension-d-1].get_points();
        const std::vector< double >& tmp_weights = face_quadratures[dimension-d-1].get_weights();
        for (unsigned int i = 0; i < quadrature_points.size(); i++){
          quadrature_points[i] = tmp_quadratures[quadrature_indices[i]];
          // Since the quadrature rule is imposed by a trapezoid rule in one
          // direction, we have to multiply by 2 to nullify the 0.5 weight of the
          // trapezoid rule
          quadrature_weights[i] = 2. * tmp_weights[quadrature_indices[i]];
        }
        
        // std::cout << "d = " << dimension-d-1 << std::endl; 
        // for (const auto& P :quadrature_points){
        //   for (unsigned int d = 0; d < dimension; d++)
        //     printf("%1.8f ", P(d));
        //   std::cout << std::endl;
        // }
      
      
        evals.insert(evals.begin(), 
                      quadrature_points.begin(),
                      quadrature_points.end());
        weights.insert(weights.begin(), 
                        quadrature_weights.begin(),
                        quadrature_weights.end());
      }
    }
      // std::cout << "End" << std::endl;
      // for (const auto &P : evals){
      //   for (unsigned int d = 0 ; d < dimension; d++)
      //     printf("%+1.8f ", P(d));
      //   std::cout << std::endl;
      // }
    // Store evaluation points if desired 
    if (flags & update_quadrature_points){
      for (unsigned int i = 0; i < bernstein_quadrature_points; i++){
        quadrature_points[i]  = evals[i];
        quadrature_weights[i] = weights[i];
      }
    }

    // With the base generate anisotropic bernstein polynomials
    AnisotropicPolynomials<dim> abp(base);

    // Check if everything so far is set correctly
    Assert(abp.n() == dofs_per_cell, ExcInternalError());
    Assert(evals.size() == bernstein_quadrature_points, ExcInternalError());
    Assert(weights.size() == bernstein_quadrature_points, ExcInternalError());

    // Use abp to fill table values at Gauss points:
    std::vector< double >               values;
    std::vector< Tensor<1, dimension> > grads;
    std::vector< Tensor<2, dimension> > grad_grads;
    std::vector< Tensor<3, dimension> > third_derivatives;
    std::vector< Tensor<4, dimension> > fourth_derivatives;

    values.resize(dofs_per_cell);
    grads.resize(dofs_per_cell);
    if (flags & update_hessians) {
      grad_grads.resize(dofs_per_cell);
      for (unsigned int j = 0; j < bernstein_quadrature_points; j++) { // Iterate over points
        // Evaluate polynomials:
        abp.evaluate( evals[j],
                      values,
                      grads,
                      grad_grads,
                      third_derivatives,
                      fourth_derivatives );

        for (unsigned int i = 0; i < dofs_per_cell; i++)
          bernstein_values(j, i) = values[i];

        for (unsigned int i = 0; i < dofs_per_cell; i++)
          bernstein_grads(j, i) = grads[i];

        for (unsigned int i = 0; i < dofs_per_cell; i++)
          bernstein_hessians(j, i) = grad_grads[i];

        // third derivatives not implemented
        //if (third_derivatives.size() > 0)
        //  for (unsigned int i = 0; i < n; i++)
        //    bernstein_3rd_on_reference(i, j) = third_derivatives[i];
      }
    } else {
      for (unsigned int j = 0; j < bernstein_quadrature_points; j++) { // Iterate over points
        // Evaluate polynomials:
        abp.evaluate( evals[j],
                      values,
                      grads,
                      grad_grads,
                      third_derivatives,
                      fourth_derivatives );

        for (unsigned int i = 0; i < dofs_per_cell; i++)
          bernstein_values(j, i) = values[i];

        for (unsigned int i = 0; i < dofs_per_cell; i++)
          bernstein_grads(j, i) = grads[i];
      }
    }
  } // set_bernstein_tables()

  template<int dim, int spacedim, int soldim>
  const double& TSValuesBase<dim, spacedim, soldim>::shape_value(
    const unsigned int function_no,
    const unsigned int point_no
  ) const {
    AssertIndexRange(function_no, soldim * dofs_per_cell);
    AssertIndexRange(point_no,    total_quadrature_points);
    Assert(flags & update_values, ExcNotInitialized());
    const unsigned int non_zero_component = system_to_component_index(function_no).second;
    return shape_values(point_no, non_zero_component);
  }

  template<int dim, int spacedim, int soldim>
  const Tensor<1, spacedim>& TSValuesBase<dim, spacedim, soldim>::shape_grad(
    const unsigned int function_no,
    const unsigned int point_no
  ) const {
    AssertIndexRange(function_no, soldim * dofs_per_cell);
    AssertIndexRange(point_no,    total_quadrature_points);
    Assert(flags & update_gradients, ExcNotInitialized());
    const unsigned int non_zero_component = system_to_component_index(function_no).second;
    return shape_grads(point_no, non_zero_component);
  }

  template<int dim, int spacedim, int soldim>
  const Tensor<2, spacedim>& TSValuesBase<dim, spacedim, soldim>::shape_hessian(
    const unsigned int function_no,
    const unsigned int point_no
  ) const {
#ifdef DEBUG
    Assert(flags & update_hessians, ExcNotInitialized());
    AssertIndexRange(function_no, soldim * dofs_per_cell);
    AssertIndexRange(point_no,    total_quadrature_points);
    const unsigned int non_zero_component = system_to_component_index(function_no).second;
    return shape_hessians(point_no, non_zero_component);
#else 
    if (!flags & update_hessians) {
      throw ExcNotInitialized();
      return shape_hessians(0, 0);
    } else {
      const unsigned int non_zero_component = system_to_component_index(function_no).second;
      return shape_hessians(point_no, non_zero_component);
    }
#endif
  }

  template<int dim, int spacedim, int soldim>
  const Tensor<3, spacedim>& 
      TSValuesBase<dim, spacedim, soldim>::shape_3rd_derivative(
    const unsigned int function_no,
    const unsigned int point_no
  ) const {
    Assert(false, ExcNotImplemented());
#ifdef DEBUG
    // As soon as it is implemented:
    Assert(flags & update_3rd_derivatives, ExcNotInitialized());
    AssertIndexRange(function_no, soldim * dofs_per_cell);
    AssertIndexRange(point_no,    total_quadrature_points);
    const unsigned int non_zero_component = system_to_component_index(function_no).second;
    return third_derivative(point_no, non_zero_component);
#else 
    throw ExcNotInitialized();
    if (!flags & update_3rd_derivatives) {
      throw ExcNotInitialized();
      return third_derivative(0, 0); 
    } else {
      const unsigned int non_zero_component = system_to_component_index(function_no).second;
      return third_derivative(point_no, non_zero_component);
    }
#endif
  }

  template<int dim, int spacedim, int soldim>
  const DerivativeForm<1, dim, spacedim>& 
      TSValuesBase<dim, spacedim, soldim>::jacobian(
    const unsigned int point_no
  ) const {
#ifdef DEBUG
    Assert(flags & update_jacobians, ExcNotInitialized());
    AssertIndexRange(point_no, total_quadrature_points);
    return J[point_no];
#else 
    return J.at(point_no);
#endif
  }

  template<int dim, int spacedim, int soldim>
  const DerivativeForm<1, spacedim, dim>& 
      TSValuesBase<dim, spacedim, soldim>::inverse_jacobian(
    const unsigned int point_no
  ) const {
#ifdef DEBUG
    AssertIndexRange(point_no, total_quadrature_points);
    Assert(flags & update_inverse_jacobians, ExcNotInitialized());
    return I[point_no];
#else 
    return I.at(point_no);
#endif
  }

  template<int dim, int spacedim, int soldim>
  const Tensor<1, spacedim, Tensor<2, dim>>& 
      TSValuesBase<dim, spacedim, soldim>::jacobian_grad(
    const unsigned int point_no
  ) const {
#ifdef DEBUG
    AssertIndexRange(point_no, total_quadrature_points);
    Assert(flags & update_jacobian_grads, ExcNotInitialized());
    return H[point_no];
#else
    return H.at(point_no);
#endif
  }

  template<int dim, int spacedim, int soldim>
  const Tensor<1, dim, Tensor<2, spacedim>>& 
      TSValuesBase<dim, spacedim, soldim>::inverse_jacobian_grad(
    const unsigned int point_no
  ) const {
#ifdef DEBUG
    AssertIndexRange(point_no, total_quadrature_points);
    Assert(flags & update_jacobian_grads, ExcNotInitialized());
    return HI[point_no];
#else
    return HI.at(point_no);
#endif
  }

  template<int dim, int spacedim, int soldim>
  const double&
      TSValuesBase<dim, spacedim, soldim>::JxW(
    const unsigned int point_no
  ) const {
#ifdef DEBUG
    AssertIndexRange(point_no, total_quadrature_points);
    Assert(flags & update_JxW_values, ExcNotInitialized());
    return dx[point_no];
#else 
    return dx.at(point_no);
#endif
  }

  template<int dim, int spacedim, int soldim>
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
      TSValuesBase<dim, spacedim, soldim>::dof_indices(
  ) const {
    return {0, solution_dimension * dofs_per_cell};
  }

  template<int dim, int spacedim, int soldim>
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
      TSValuesBase<dim, spacedim, soldim>::quadrature_point_indices(
  ) const {
    return {0, total_quadrature_points};
  }

  template<int dim, int spacedim, int soldim>
  const Point<spacedim>& 
      TSValuesBase<dim, spacedim, soldim>::quadrature_point(
    const unsigned int q
  ) const {
#ifdef DEBUG
    AssertIndexRange(q, total_quadrature_points);
    Assert(flags & update_quadrature_points, ExcNotInitialized());
    return mapped_quadrature_points[q];
#else 
    return mapped_quadrature_points.at(q);
#endif
  }

  template<int dim, int spacedim, int soldim>
  const Point<dim>& 
      TSValuesBase<dim, spacedim, soldim>::quadrature_point_reference(
    const unsigned int q
  ) const {
#ifdef DEBUG
    AssertIndexRange(q, quadrature_points.size(0));
    return quadrature_points(q + offset);
#else 
    return quadrature_points(q + offset);
#endif
  }

  template<int dim, int spacedim, int soldim>
  const unsigned int&
      TSValuesBase<dim, spacedim, soldim>::n_quadrature_points_per_cell(
  ) const {
    return total_quadrature_points;
  }
  
  template<int dim, int spacedim, int soldim>
  const unsigned int&
      TSValuesBase<dim, spacedim, soldim>::n_quadrature_points_per_face(
  ) const {
    Assert(face_no != -1, ExcMessage("Invalid function called for cells"))
    return quadrature_points_per_face[face_no/2];
  }

  template<int dim, int spacedim, int soldim>
  const std::vector<unsigned int>&
      TSValuesBase<dim, spacedim, soldim>::n_q_points_per_cell(
  ) const {
    return n_gauss_points;
  }

  template<int dim, int spacedim, int soldim>
  unsigned int
      TSValuesBase<dim, spacedim, soldim>::n_dofs_per_cell(
  ) const {
    return solution_dimension * dofs_per_cell;
  }

  template<int dim, int spacedim, int soldim>
  std::pair< unsigned int, unsigned int> 
      TSValuesBase<dim, spacedim, soldim>::system_to_component_index(
    const unsigned int index
  ) const {
    return std::make_pair<unsigned int, unsigned int>(
              index % solution_dimension, 
              index / solution_dimension);
  } // system_to_component_index()

  template<int dim, int spacedim, int soldim>
  unsigned int 
      TSValuesBase<dim, spacedim, soldim>::component_to_system_index(
    const unsigned int component,
    const unsigned int index
  ) const {
    Assert(component < solution_dimension, ExcIndexRange(component, 0, solution_dimension));
    
    return dim * index + component ;
  } // component_to_system()

  template<int dim, int spacedim, int soldim>
  const std::vector< Point< spacedim > >& 
      TSValuesBase<dim, spacedim, soldim>::get_quadrature_points(
  ) const {
    return mapped_quadrature_points;
  } // TSValuesBase<dim, spacedim, soldim>::get_quadrature_points()

  template<int dim, int spacedim, int soldim>
  const unsigned int& 
      TSValuesBase<dim, spacedim, soldim>::n_quadrature_points(
  ) const {
    return total_quadrature_points;
  } // TSValuesBase<dim, spacedim, soldim>::n_quadrature_points()

  template<int dim, int spacedim, int soldim>
  TSFaceValuesBase<dim, spacedim, soldim>::TSFaceValuesBase(
    const TS_TriangulationBase<dim, spacedim>*  tria,
    const std::vector<unsigned int>&            n_gauss_points,
    const UpdateFlags                           flags
  ) : TSValuesBase<dim, spacedim, soldim>(tria, n_gauss_points, flags, 0) { 
    // The additional 0 in the constructor of TSValuesBase, simply sets
    // the bernstein tables up to work on faces, i.e. the offset parameter
    // has to be recomputed on each face to get the corresponding quadrature
    // points and values from the bernstein tables.
    this -> init_shape_tables();
  } // TSFaceValuesBase constructor 2

  template<int dim, int spacedim, int soldim>
  TSFaceValuesBase<dim, spacedim, soldim>::TSFaceValuesBase(
    const TS_TriangulationBase<dim, spacedim>*  tria,
    const unsigned int&                         n_gauss_points,
    const UpdateFlags                           flags
  ) : TSValuesBase<dim, spacedim, soldim>(tria, n_gauss_points, flags, 0) { 
    // The additional 0 in the constructor of TSValuesBase, simply sets
    // the bernstein tables up to work on faces, i.e. the offset parameter
    // has to be recomputed on each face to get the corresponding quadrature
    // points and values from the bernstein tables.

    // If we use an anisotropic quadrature, we can initialize
    // the shape tables at this point already. They will have the same size
    // on every face
    anisotropic_quadrature = false;
    this -> init_shape_tables();
  } // TSFaceValuesBase constructor 2

  template<int dim, int spacedim, int soldim>
  const Tensor<1, spacedim>& 
      TSFaceValuesBase<dim, spacedim, soldim>::normal_vector(
    const unsigned int q
  ) const {
    AssertIndexRange(q, this->total_quadrature_points);
    Assert(this->flags & update_normal_vectors, ExcNotInitialized());
    return this->normals[q];
  }

  template<int dim, int spacedim, int soldim>
  void TSValues<dim, spacedim, soldim>::set_grad_and_hessian_tables(
    const active_cell_iterator& cell
  ) {
    const unsigned int& dofs_per_cell = this -> dofs_per_cell;
    const unsigned int& offset        = this -> offset;
    const UpdateFlags& flags          = this -> flags;

    const Table<2, double>&                 bernstein_values    = this -> bernstein_values;
    const Table<2, Tensor<1, dimension> >&  bernstein_grads     = this -> bernstein_grads;
    const Table<2, Tensor<2, dimension> >&  bernstein_hessians  = this -> bernstein_hessians;
    Assert(space_dimension == dimension, ExcNotImplemented());
    if (space_dimension != dimension)
      return;
    // Get coefficients for linear combination
    const FullMatrix<double>&             coeffs = 
            this -> tria -> get_bezier_coefficients(cell);
    // Get weighted control points for geometric mapping
    const std::vector<Point<spacedim+1>>  control_points = 
            this -> tria -> get_control_points(cell);

    // store intermediate values 
    Table<1, double>                 spline_values(dofs_per_cell);
    Table<1, Tensor<1, dimension>>   spline_grads(dofs_per_cell); 
    Table<1, Tensor<2, dimension>>   spline_grad_grads(dofs_per_cell);

    // Split weighted control points
    std::vector< Tensor<1, spacedim >> controls_reduced(dofs_per_cell);
    std::vector< double >              controls_weights(dofs_per_cell);
    for (unsigned int i = 0; i < dofs_per_cell; i++){
      controls_weights[i] = control_points[i](space_dimension);
      for (unsigned int d = 0; d < space_dimension; d++)
        controls_reduced[i][d] = control_points[i](d);
    } // for ( i )

    // Get a reference to full tables for ts values, grads, and grad grads
    Table<2, double>& ts_values =
      this -> shape_values;
    Table<2, Tensor<1, space_dimension>>& ts_grads =
      this -> shape_grads;
    Table<2, Tensor<2, space_dimension>>& ts_grad_grads =
      this -> shape_hessians;
    // For each quadrature point, compute the corresponding entries in the
    // tables.
    
    const Point<dimension>& lengths = cell -> vertex(GeometryInfo<dimension>::vertices_per_cell-1) 
                                      + (-1.) * cell -> vertex(0);
    const double            measure = cell -> measure();
    // Initialize variables to be computed along the way
    Point<spacedim>                                   mapped_quadrature_point ;
    DerivativeForm<1, dimension, space_dimension>     mapping_grad ;
    DerivativeForm<1, space_dimension, dimension>     inverse_mapping_grad ;
    Tensor<1, space_dimension, Tensor<2, dimension>>  mapping_grad_grad ; 
    Tensor<1, dimension, Tensor<2, space_dimension>>  inverse_mapping_grad_grad ; 
    double                                            JxW;
    for (unsigned int pnt = 0; pnt < this -> shape_quadrature_points; pnt++) {
      // Get a reference to corresponding table entries
      const double& quadrature_weight = this -> quadrature_weights[pnt];
      
      // reset temporary variables
      spline_values.reset_values();
      spline_grads.reset_values();
      spline_grad_grads.reset_values();


      { // Scope section to clear helper variables
        // helper variables
        Tensor<1, space_dimension, Tensor<1, dimension> > N1, N2;
        double D1, D2;

        Tensor<1, space_dimension>                        controls_value;
        Tensor<1, dimension>                              weights_grad;
        Tensor<1, space_dimension, Tensor<1, dimension>>  controls_grad;
        Tensor<1, space_dimension, Tensor<2, dimension>>  controls_hessian;
        Tensor<2, dimension>                              weights_hessian;
        double                                            weight = 0;


        for (unsigned int fcn = 0; fcn < dofs_per_cell; fcn++){
          // Get references to the values to be updated
          double&         spline_value            = spline_values(fcn);
          Tensor<1, dimension>& spline_grad       = spline_grads(fcn);
          Tensor<2, dimension>& spline_grad_grad  = spline_grad_grads(fcn);

          // run linear combination
          for (unsigned int fcn_b = 0; fcn_b < dofs_per_cell; fcn_b++){
            spline_value      += coeffs(fcn, fcn_b) * 
                                 bernstein_values(pnt + offset, fcn_b);
            spline_grad       += coeffs(fcn, fcn_b) * 
                                 bernstein_grads(pnt + offset, fcn_b);
            spline_grad_grad  += coeffs(fcn, fcn_b) *
                                 bernstein_hessians(pnt + offset, fcn_b);
          }
          // Update values for reference cell
          for (unsigned int d = 0; d < dimension; d++){
            spline_grad[d]        /= ( lengths(d) );
            for (unsigned int dd = 0; dd < dimension; dd++)
              spline_grad_grad[d][dd] /= ( lengths(d) * lengths(dd) );
          }

          weight           += controls_weights[fcn] * spline_value;
          controls_value   += controls_reduced[fcn] * spline_value;
          weights_grad     += controls_weights[fcn] * spline_grad;
          weights_hessian  += controls_weights[fcn] * spline_grad_grad;
          for ( unsigned int d = 0; d < space_dimension; d++){
            controls_grad[d]    += controls_reduced[fcn][d] * spline_grad;
            controls_hessian[d] += controls_reduced[fcn][d] * spline_grad_grad;
          }
        } // for ( fcn )

        // store the mapped quadrature point:
        for (unsigned int d = 0; d < space_dimension; d++)
          mapped_quadrature_point(d) = controls_value[d] / weight;

        D1 = weight * weight;
        D2 = D1 * D1;
        for (unsigned int d = 0; d < space_dimension; d++){
          N1[d] = controls_value[d] * weights_grad;
          N2[d] = controls_grad[d]  * weight;

          Tensor<2, dimension>& Hk = mapping_grad_grad[d];

          for (unsigned int m = 0; m < dimension; m++){
            for (unsigned int n = m; n < dimension; n++){
              Hk[n][m] = ( ( controls_hessian[d][n][m] * weight
                                 + controls_grad[d][m] * weights_grad[n]
                                 - controls_grad[d][n] * weights_grad[m]
                                 - controls_value[d] * weights_hessian[n][m] ) * D1 
                               - 2. * (controls_grad[d][m] * weight 
                                         - controls_value[d] * weights_grad[m]) * 
                                         weights_grad[n] * weight ) / D2;
              Hk[m][n] = Hk[n][m];
            } // for ( n )
          } // for ( m )
        } // for ( d )

        // Get value of Jacobian at current gauss point
        mapping_grad = DerivativeForm<1, dim, spacedim>( (1. / D1) *
                N2.Tensor<1, spacedim, Tensor<1, dimension> >::operator-=(N1) );
        // At the end of scope, all the helper variables are not needed anymore and can be destroyed
      } // scope 

      // Store the inverse of the jacobian at the current point. The covariant of a 
      // DerivativeForm is given by C = DF * (DF^T * DF)^{-1}. Hence, taking the 
      // transpose we get C^T = (DF^T * DF)^{-1} DF^T which is exactly the left-inverse
      // at the current point pnt. The function covariant_form() simplifies in case 
      // of rectangular DerivativeForms, i.e. if dim == spacedim, and DF^{-T} is returned. 
      const DerivativeForm<1, dimension, space_dimension>& covariant = 
              mapping_grad.covariant_form();
      inverse_mapping_grad = covariant.transpose();

      // Save values to tables:
      JxW  = measure *
              quadrature_weight *
              std::abs(mapping_grad.determinant());

      // Using the Jacobians and the inverse of the jacobians we can setup the hessian 
      // of the transformation in terms of Tensors. Note, that we have 
      // HI = - I . H[I] . I, 
      // with the following data-types
      //    HI   : Tensor<1, dim, Tensor<2, spacedim> >, i.e. a size dim array with rank-2 tensors
      //    I    : DerivativeForm<1, spacedim, dim> ~ Tensor<1, dim, Tensor<1, spacedim>>, i.e. 
      //           a size dim array, where each entry is a Tensor of length spacedim
      //    H    : Tensor<1, spacedim, Tensor<2, dim> >, i.e. a size spacedim array with rank-2
      //           tensors as entries. 
      //    H[I] : Tensor<1, spacedim, Tensor<1, dim, Tensor<1, spacedim> > >, this object refers
      //           mathematically to the application of H to I, it is a rank-3 Tensor where each 
      //           entry refers to a dim by spacedim matrix, and the data type we use is according 
      //           to this
      // It is important to note, that the dots in the multiplication refer to tensor products
      // and not the typical matrix products!
      for (unsigned int k = 0; k < dimension; k++){
        Tensor<2, space_dimension>& HIk = inverse_mapping_grad_grad[k];
        // first build the sum over hessians of forward
        // mapping
        Tensor<2, dimension> temp;
        for (unsigned int d = 0; d < space_dimension; d++)
          temp += -1. * inverse_mapping_grad[k][d] * mapping_grad_grad[d];
        //    mapping_grad_grad             ~ Tensor<1, spacedim, Tensor<2, dim>>
        // -> mapping_grad_grad[d]          ~ Tensor<2, dim>
        //    inverse_mapping_grad          ~ DerivativeForm<1, spacedim, dim>
        // -> inverse_mapping_grad[k]       ~ Tensor<1, spacedim>
        // -> inverse_mapping_grad[k][d] * mapping_grad_grad[d] ~ Tensor<2, dim>
        // Note that this is only allowed if spacedim == dim

        // perform matrix-matrix product
        const DerivativeForm<1, space_dimension, dimension>& Ipnt = inverse_mapping_grad;
        if (dimension == 2) {
          HIk[0][0] = Ipnt[0][0] * (temp[0][0] * Ipnt[0][0] + temp[0][1] * Ipnt[1][0]) 
                    + Ipnt[1][0] * (temp[1][0] * Ipnt[0][0] + temp[1][1] * Ipnt[1][0]);
          HIk[0][1] = Ipnt[0][0] * (temp[0][0] * Ipnt[0][1] + temp[0][1] * Ipnt[1][1]) 
                    + Ipnt[1][0] * (temp[1][0] * Ipnt[0][1] + temp[1][1] * Ipnt[1][1]);
          HIk[1][1] = Ipnt[0][1] * (temp[0][0] * Ipnt[0][1] + temp[0][1] * Ipnt[1][1]) 
                    + Ipnt[1][1] * (temp[1][0] * Ipnt[0][1] + temp[1][1] * Ipnt[1][1]);
          HIk[1][0] = HIk[0][1];
        } else if (dimension == 3) {
          HIk[0][0] = Ipnt[0][0] * (
                        temp[0][0] * Ipnt[0][0] + temp[0][1] * Ipnt[1][0] + temp[0][2] * Ipnt[2][0]
                      ) 
                      + Ipnt[1][0] * (
                        temp[1][0] * Ipnt[0][0] + temp[1][1] * Ipnt[1][0] + temp[1][2] * Ipnt[2][0]
                      ) 
                      + Ipnt[2][0] * (
                        temp[2][0] * Ipnt[0][0] + temp[2][1] * Ipnt[1][0] + temp[2][2] * Ipnt[2][0]
                      );

          HIk[1][0] = Ipnt[0][1] * (
                        temp[0][0] * Ipnt[0][0] + temp[0][1] * Ipnt[1][0] + temp[0][2] * Ipnt[2][0]
                      ) 
                      + Ipnt[1][1] * (
                        temp[1][0] * Ipnt[0][0] + temp[1][1] * Ipnt[1][0] + temp[1][2] * Ipnt[2][0]
                      ) 
                      + Ipnt[2][1] * (
                        temp[2][0] * Ipnt[0][0] + temp[2][1] * Ipnt[1][0] + temp[2][2] * Ipnt[2][0]
                      );

          HIk[2][0] = Ipnt[0][2] * (
                        temp[0][0] * Ipnt[0][0] + temp[0][1] * Ipnt[1][0] + temp[0][2] * Ipnt[2][0]
                      ) 
                      + Ipnt[1][2] * (
                        temp[1][0] * Ipnt[0][0] + temp[1][1] * Ipnt[1][0] + temp[1][2] * Ipnt[2][0]
                      ) 
                      + Ipnt[2][2] * (
                        temp[2][0] * Ipnt[0][0] + temp[2][1] * Ipnt[1][0] + temp[2][2] * Ipnt[2][0]
                      );

          HIk[1][1] = Ipnt[0][1] * (
                        temp[0][0] * Ipnt[0][1] + temp[0][1] * Ipnt[1][1] + temp[0][2] * Ipnt[2][1]
                      ) 
                      + Ipnt[1][1] * (
                        temp[1][0] * Ipnt[0][1] + temp[1][1] * Ipnt[1][1] + temp[1][2] * Ipnt[2][1]
                      ) 
                      + Ipnt[2][1] * (
                        temp[2][0] * Ipnt[0][1] + temp[2][1] * Ipnt[1][1] + temp[2][2] * Ipnt[2][1]
                      );

          HIk[2][1] = Ipnt[0][2] * (
                        temp[0][0] * Ipnt[0][1] + temp[0][1] * Ipnt[1][1] + temp[0][2] * Ipnt[2][1]
                      ) 
                      + Ipnt[1][2] * (
                        temp[1][0] * Ipnt[0][1] + temp[1][1] * Ipnt[1][1] + temp[1][2] * Ipnt[2][1]
                      ) 
                      + Ipnt[2][2] * (
                        temp[2][0] * Ipnt[0][1] + temp[2][1] * Ipnt[1][1] + temp[2][2] * Ipnt[2][1]
                      );

          HIk[2][2] = Ipnt[0][2] * (
                        temp[0][0] * Ipnt[0][2] + temp[0][1] * Ipnt[1][2] + temp[0][2] * Ipnt[2][2]
                      ) 
                      + Ipnt[1][2] * (
                        temp[1][0] * Ipnt[0][2] + temp[1][1] * Ipnt[1][2] + temp[1][2] * Ipnt[2][2]
                      ) 
                      + Ipnt[2][2] * (
                        temp[2][0] * Ipnt[0][2] + temp[2][1] * Ipnt[1][2] + temp[2][2] * Ipnt[2][2]
                      );

          HIk[0][1] = HIk[1][0];
          HIk[0][2] = HIk[2][0];
          HIk[1][2] = HIk[2][1];
        }


        // The above implementation is rather lengthy, simple for-loops
        // may help. However, I am currently unable to find the mistake in 
        // the code below. In theory, there should be 
        //  HIk = I^{-T} temp  I,
        // which we try to resemble in the code below. 
        // for (unsigned int i = 0; i < space_dimension; i++){
        //   for (unsigned int j = 0; j <= i; j++) {
        //     double&                     HIk_ij = HIk[i][j];
        //     const Tensor<1, dimension>& grad_i = covariant[i];
        //     const Tensor<1, dimension>& grad_j = covariant[j];
        //     //    covariant               ~ DerivativeForm<1, dim, spacedim>
        //     // -> covariant[i]            ~ Tensor<1, dim>
        //     //    inverse_mapping_grad    ~ DerivativeForm<1, spacedim, dim>
        //     // -> inverse_mapping_grad[j] ~ Tensor<1, spacedim>
        //     for (unsigned int m = 0; m < dimension; m++)
        //       for (unsigned int n = 0; n < dimension; n++)
        //         HIk_ij += grad_i[n] * temp[n][m] * grad_j[m];

        //     HIk[j][i] = HIk_ij;
        //   } // for ( j )
        // } // for ( i )
      } // for ( k )

      // Using the inverse hessians, we can then compute the hessian of the 
      // shape functions in physical coordinates from the hessians in parametric
      // coordinates and the hessians of the inverse functions of Phi, using the
      // formula
      // ^H_[T o Phi] = I^tr * H_[T] * I + sum (\partial T / \partial x_k) * H_[Phi^{-1}_k] 
      //              = I^tr * H_[T] * I + sum grad_[T]_k * H_[Phi_k]
      for (unsigned int fcn = 0; fcn < dofs_per_cell; fcn++){
        ts_values(pnt, fcn)                          = spline_values(fcn);
        const Tensor<1, dimension>& spline_grad      = spline_grads(fcn);
        const Tensor<2, dimension>& spline_grad_grad = spline_grad_grads(fcn);
        
        // J_Phi * grad T_i
        if (flags & update_gradients)
          ts_grads(pnt, fcn) = apply_transformation(covariant, spline_grad);

        // For the hessian get two intermediate results
        // First, the latter sum of hessians
        Tensor<2, spacedim> sum;
        for (unsigned int d = 0; d < dimension; d++)
          sum += spline_grad[d] * inverse_mapping_grad_grad[d];

        // Second, the first matrix-matrix product
        Tensor<1, spacedim, Tensor<1, dim> > ip; 
        for (unsigned int k = 0; k < space_dimension; k++) {
          for (unsigned int d = 0; d < dimension; d++)
            ip[k][d] = spline_grad_grad[d] * covariant[k];
        }

        // Use the first matrix-matrix product, to compute the first term of ^H_[T o Phi]
        Tensor<2, spacedim>& h = ts_grad_grads(pnt, fcn);
        for (unsigned int i = 0; i < space_dimension; i++)
          for (unsigned int k = 0; k < space_dimension; k++)
            h[i][k] = (double)(covariant[i] * ip[k]);

        // Add the results: 
        h += sum;

//        std::cout << "Hessian: " << std::endl;
//        for (unsigned int i = 0; i < spacedim; i++)
//          std::cout << std::scientific << h[i] << std::endl;
//
//        // Test for symmetry:
//        for (unsigned int i = 0; i < spacedim; i++)
//          for (unsigned int j = 0; j < spacedim; j++)
//            Assert(h[i][j] == h[j][i], ExcInternalError());

      } // for ( fcn )
      // Get a reference to corresponding table entries
      if (flags & update_quadrature_points)
        this -> mapped_quadrature_points[pnt] = mapped_quadrature_point;
      if (flags & update_jacobians)
        this -> J[pnt] = mapping_grad;
      if (flags & update_inverse_jacobians)
        this -> I[pnt] = inverse_mapping_grad; 
      if (flags & update_jacobian_grads){
        this -> H[pnt]  = mapping_grad_grad;
        this -> HI[pnt] = inverse_mapping_grad_grad;
      }
      if (flags & update_JxW_values)
        this -> dx[pnt] = JxW;
    } // for ( pnt )
  } // set_grad_and_hessian_tables

  template<int dim, int spacedim, int soldim>
  void TSValues<dim, spacedim, soldim>::set_grad_tables(
    const active_cell_iterator& cell
  ) {
    const unsigned int& dofs_per_cell = this -> dofs_per_cell;
    const unsigned int& offset        = this -> offset;
    const UpdateFlags& flags          = this -> flags;

    const Table<2, double>&                 bernstein_values    = this -> bernstein_values;
    const Table<2, Tensor<1, dimension> >&  bernstein_grads     = this -> bernstein_grads;

    // Get coefficients for linear combination
    const FullMatrix<double>&             coeffs = 
            this -> tria -> get_bezier_coefficients(cell);
    // Get weighted control points for geometric mapping
    const std::vector<Point<space_dimension+1>>  control_points = 
            this -> tria -> get_control_points(cell);

    // store intermediate values 
    Table<1, double>                 spline_values(dofs_per_cell);
    Table<1, Tensor<1, dimension>>   spline_grads(dofs_per_cell); 

    // Split weighted control points
    std::vector< Tensor<1, spacedim >> controls_reduced(dofs_per_cell);
    std::vector< double >              controls_weights(dofs_per_cell);
    for (unsigned int i = 0; i < dofs_per_cell; i++){
      controls_weights[i] = control_points[i](space_dimension);
      for (unsigned int d = 0; d < space_dimension; d++)
        controls_reduced[i][d] = control_points[i](d);
    } // for ( i )

    // Get a reference to full tables for ts values, grads, and grad grads
    Table<2, double>& ts_values =
      this -> shape_values;
    Table<2, Tensor<1, space_dimension>>& ts_grads =
      this -> shape_grads;
    // For each quadrature point, compute the corresponding entries in the
    // tables.
    
    const Point<dimension>& lengths = cell -> vertex(GeometryInfo<dimension>::vertices_per_cell-1) 
                                      + (-1.) * cell -> vertex(0);
    const double            measure = cell -> measure();
    
    // Define variables to be computed along the way,
    // but allocate space only once
    Point<space_dimension>                        mapped_quadrature_point;
    DerivativeForm<1, dimension, space_dimension> mapping_grad;
    DerivativeForm<1, space_dimension, dimension> inverse_mapping_grad;
    double                                        JxW;
    for (unsigned int pnt = 0; pnt < this -> shape_quadrature_points; pnt++) {
      const double& quadrature_weight = this -> quadrature_weights[pnt];
      // reset temporary variables
      spline_values.reset_values();
      spline_grads.reset_values();

      { // Scope section to clear helper variables
        // helper variables
        Tensor<1, space_dimension, Tensor<1, dimension> > N1, N2;
        double D1;

        Tensor<1, space_dimension>                        controls_value;
        Tensor<1, dimension>                              weights_grad;
        Tensor<1, space_dimension, Tensor<1, dimension>>  controls_grad;
        double                                            weight = 0;


        for (unsigned int fcn = 0; fcn < dofs_per_cell; fcn++){
          // Get references to the values to be updated
          double&               spline_value = spline_values(fcn);
          Tensor<1, dimension>& spline_grad  = spline_grads(fcn);

          // run linear combination
          for (unsigned int fcn_b = 0; fcn_b < dofs_per_cell; fcn_b++){
            spline_value += coeffs(fcn, fcn_b) * 
                            bernstein_values(pnt + offset, fcn_b);
            spline_grad  += coeffs(fcn, fcn_b) * 
                            bernstein_grads(pnt + offset, fcn_b);
          }
          // Update values for reference cell
          for (unsigned int d = 0; d < dimension; d++)
            spline_grad[d] /= ( lengths(d) );

          weight           += controls_weights[fcn] * spline_value;
          controls_value   += controls_reduced[fcn] * spline_value;
          weights_grad     += controls_weights[fcn] * spline_grad;
          for (unsigned int d = 0; d < space_dimension; d++)
            controls_grad[d] += controls_reduced[fcn][d] * spline_grad;
        } // for ( fcn )

        // store the mapped quadrature point:
        for (unsigned int d = 0; d < space_dimension; d++)
          mapped_quadrature_point(d) = controls_value[d] / weight;

        D1 = weight * weight;
        for (unsigned int d = 0; d < space_dimension; d++){
          N1[d] = controls_value[d] * weights_grad;
          N2[d] = controls_grad[d]  * weight;
        } // for ( d )

        // Get value of Jacobian at current gauss point
        mapping_grad = DerivativeForm<1, dim, spacedim>( (1. / D1) *
                N2.Tensor<1, spacedim, Tensor<1, dimension> >::operator-=(N1) );
        // At the end of scope, all the helper variables are not needed anymore and can be destroyed
      } // scope 

      // Store the inverse of the jacobian at the current point. The covariant of a 
      // DerivativeForm is given by C = DF * (DF^T * DF)^{-1}. Hence, taking the 
      // transpose we get C^T = (DF^T * DF)^{-1} DF^T which is exactly the left-inverse
      // at the current point pnt. The function covariant_form() simplifies in case 
      // of rectangular DerivativeForms, i.e. if dim == spacedim, and DF^{-T} is returned. 
      const DerivativeForm<1, dimension, space_dimension>& covariant = 
              mapping_grad.covariant_form();
      inverse_mapping_grad = covariant.transpose();

      // Save values to tables:
      JxW  = measure *
              quadrature_weight *
              std::abs(mapping_grad.determinant());


      for (unsigned int fcn = 0; fcn < dofs_per_cell; fcn++){
        ts_values(pnt, fcn) = spline_values(fcn);
        
        // J_Phi * grad T_i
        if (flags & update_gradients)
          ts_grads(pnt, fcn)  = apply_transformation(covariant, spline_grads(fcn));
      } // for ( fcn )
      if (flags & update_quadrature_points)
        this -> mapped_quadrature_points[pnt] = mapped_quadrature_point;
      if (flags & update_jacobians)
        this -> J[pnt] = mapping_grad;
      if (flags & update_inverse_jacobians)
        this -> I[pnt] = inverse_mapping_grad; 
      if (flags & update_JxW_values)
        this -> dx[pnt] = JxW;
    } // for ( pnt )
  } // set_grad_tables


  template<int dim, int spacedim, int soldim>
  TSValues<dim, spacedim, soldim>::TSValues(
    const TS_TriangulationBase<dim, spacedim>*  tria,
    const unsigned int                          n_gauss_points,
    const UpdateFlags                           flags
  ) : TSValuesBase<dim, spacedim, soldim>(tria, n_gauss_points, flags) {
  // There's nothing TSValue specific yet
    this -> init_shape_tables();
  }

  template<int dim, int spacedim, int soldim>
  TSValues<dim, spacedim, soldim>::TSValues(
    const TS_TriangulationBase<dim, spacedim>*  tria,
    const std::vector<unsigned int>&            n_gauss_points,
    const UpdateFlags                           flags
  ) : TSValuesBase<dim, spacedim, soldim>(tria, n_gauss_points, flags) {
  // There's nothing TSValue specific yet
    this -> init_shape_tables();
  }

  template<int dim, int spacedim, int soldim>
  void TSValues<dim, spacedim, soldim>::reinit(
    const active_cell_iterator& cell
  ) { 
    Assert(this -> flags & update_values, ExcNotInitialized());
    // Assert(this -> flags & update_gradients, ExcNotInitialized());
    if (this->flags & update_hessians)
      set_grad_and_hessian_tables(cell);
    else 
      set_grad_tables(cell);
  }

  template<int dim, int spacedim, int soldim>
  TSFaceValues<dim, spacedim, soldim>::TSFaceValues(
    const TS_TriangulationBase<dim, spacedim>*  tria,
    const unsigned int                          n_gauss_points,
    const UpdateFlags                           flags
  ) : TSFaceValuesBase<dim, spacedim, soldim>(tria, n_gauss_points, flags) {
  }

  template<int dim, int spacedim, int soldim>
  TSFaceValues<dim, spacedim, soldim>::TSFaceValues(
    const TS_TriangulationBase<dim, spacedim>*  tria,
    const std::vector<unsigned int>&            n_gauss_points,
    const UpdateFlags                           flags
  ) : TSFaceValuesBase<dim, spacedim, soldim>(tria, n_gauss_points, flags) {
  // There's nothing TSFaceValues specific yet
  }

  template<int dim, int spacedim, int soldim>
  void TSFaceValues<dim, spacedim, soldim>::reinit(
    const active_cell_iterator& cell,
    const          int&         face_no
  ) { 
             int& old_face = this->face_no; 
    unsigned int& offset   = this->offset; 
             int  f        = 0;
    unsigned int n_subfaces;
    const UpdateFlags& flags = this->flags;
    
    // define the new offset to read the correct bernstein values
    offset = 0;
    while (f < face_no) 
      offset += this -> quadrature_points_per_face[f++/2];
    
    n_subfaces = 1 + (cell -> face(face_no) -> n_children())/2; // either 1 or 2
    // check if we actually need to reinitialize tables. 
    // We do not need to resize shape tables if we have
    // the same type of face, or if we used
    // isotropic quadrature
    if (old_face/2 != face_no/2 && 
          this -> anisotropic_quadrature && 
          n_subfaces == 2)
      this -> init_shape_tables(n_subfaces);

    old_face = face_no; 

    // const Table<1, Point<dimension>>& quadrature_points = this->quadrature_points; 
    // const unsigned int& n = quadrature_points.size(0); 
    // std::cout << "Quadrature points: " << std::endl;
    // for (unsigned int i = 0; i < n; i++){
    //   printf("%3i: ", i);
    //   for (unsigned int d = 0; d < dimension; d++)
    //     printf("% .8f", quadrature_points(i)(d));
    //   std::cout << std::endl; 
    // }
    // std::cout << "End" << std::endl;

    Assert(flags & update_values, ExcNotInitialized());
    // Assert(flags & update_gradients, ExcNotInitialized());
    if (flags & update_hessians)
      for (unsigned int sub_face = 0; sub_face < n_subfaces; sub_face++)
        set_grad_and_hessian_tables(cell, sub_face);
    else 
      for (unsigned int sub_face = 0; sub_face < n_subfaces; sub_face++)
        set_grad_tables(cell, sub_face);
  }


  template<int dim, int spacedim, int soldim>
  void TSFaceValues<dim, spacedim, soldim>::set_grad_and_hessian_tables(
    const active_cell_iterator& cell, 
    const unsigned int sub_face
  ) {
    Assert(space_dimension == dimension, ExcNotImplemented());
    const unsigned int& dofs_per_cell = this -> dofs_per_cell;
    const          int& face_no       = this -> face_no;
    const unsigned int& offset        = this -> offset;
    const UpdateFlags& flags          = this -> flags;

    const Table<2, double>&                 bernstein_values    = this -> bernstein_values;
    const Table<2, Tensor<1, dimension> >&  bernstein_grads     = this -> bernstein_grads;
    const Table<2, Tensor<2, dimension> >&  bernstein_hessians  = this -> bernstein_hessians;

    // Get coefficients for linear combination
    const std::vector< FullMatrix<double> >& coeffs_vector = 
            this -> tria -> get_bezier_coefficients(cell, face_no);
    const FullMatrix<double>& coeffs = 
            coeffs_vector[sub_face];
    // Get weighted control points for geometric mapping
    const std::vector<Point<spacedim+1>>  control_points = 
            this -> tria -> get_control_points(cell);

    // store the sign for the current face
    // Note that this might be unused
    const int s = GeometryInfo<dim>::unit_normal_orientation[face_no];

    // store intermediate values 
    Table<1, double>                 spline_values(dofs_per_cell);
    Table<1, Tensor<1, dimension>>   spline_grads(dofs_per_cell); 
    Table<1, Tensor<2, dimension>>   spline_grad_grads(dofs_per_cell);

    // Split weighted control points
    std::vector< Tensor<1, spacedim >> controls_reduced(dofs_per_cell);
    std::vector< double >              controls_weights(dofs_per_cell);
    for (unsigned int i = 0; i < dofs_per_cell; i++){
      controls_weights[i] = control_points[i](space_dimension);
      for (unsigned int d = 0; d < space_dimension; d++)
        controls_reduced[i][d] = control_points[i](d);
    } // for ( i )

    // Get a reference to full tables for ts values, grads, and grad grads
    Table<2, double>& ts_values =
      this -> shape_values;
    Table<2, Tensor<1, space_dimension>>& ts_grads =
      this -> shape_grads;
    Table<2, Tensor<2, space_dimension>>& ts_grad_grads =
      this -> shape_hessians;
    // For each quadrature point, compute the corresponding entries in the
    // tables.
    

    const Point<dimension>& lengths = cell -> vertex(GeometryInfo<dimension>::vertices_per_cell-1) 
                                      + (-1.) * cell -> vertex(0);
    const double            measure = (cell -> measure()) / lengths(face_no/2);
    // Initialize variables to be computed along the waz
    Point<spacedim>                                   mapped_quadrature_point ;
    DerivativeForm<1, dimension, space_dimension>     mapping_grad ;
    DerivativeForm<1, space_dimension, dimension>     inverse_mapping_grad ;
    Tensor<1, space_dimension, Tensor<2, dimension>>  mapping_grad_grad ; 
    Tensor<1, dimension, Tensor<2, space_dimension>>  inverse_mapping_grad_grad ; 
    Tensor<1, space_dimension>                        normal;
    double                                            JxW;
    for (unsigned int p = 0; p < this -> shape_quadrature_points; p++) {
      const unsigned int pnt = p + sub_face * this->shape_quadrature_points;
      // Get a reference to corresponding table entries
      const double& quadrature_weight = this -> quadrature_weights[pnt + offset];

      // Reset temporary values
      spline_values.reset_values();
      spline_grads.reset_values();
      spline_grad_grads.reset_values();


      { // Scope section to clear helper variables
        // helper variables
        Tensor<1, space_dimension, Tensor<1, dimension> > N1, N2;
        double D1, D2;

        Tensor<1, space_dimension>                        controls_value;
        Tensor<1, dimension>                              weights_grad;
        Tensor<1, space_dimension, Tensor<1, dimension>>  controls_grad;
        Tensor<1, space_dimension, Tensor<2, dimension>>  controls_hessian;
        Tensor<2, dimension>                              weights_hessian;
        double                                            weight = 0;


        for (unsigned int fcn = 0; fcn < dofs_per_cell; fcn++){
          // Get references to the values to be updated
          double&         spline_value            = spline_values(fcn);
          Tensor<1, dimension>& spline_grad       = spline_grads(fcn);
          Tensor<2, dimension>& spline_grad_grad  = spline_grad_grads(fcn);

          // run linear combination
          for (unsigned int fcn_b = 0; fcn_b < dofs_per_cell; fcn_b++){
            spline_value      += coeffs(fcn, fcn_b) * 
                                 bernstein_values(pnt + offset, fcn_b);
            spline_grad       += coeffs(fcn, fcn_b) * 
                                 bernstein_grads(pnt + offset, fcn_b);
            spline_grad_grad  += coeffs(fcn, fcn_b) *
                                 bernstein_hessians(pnt + offset, fcn_b);
          }
          // Update values for reference cell
          for (unsigned int d = 0; d < dimension; d++){
            spline_grad[d]        /= ( lengths(d) );
            for (unsigned int dd = 0; dd < dimension; dd++)
              spline_grad_grad[d][dd] /= ( lengths(d) * lengths(dd) );
          }

          weight           += controls_weights[fcn] * spline_value;
          controls_value   += controls_reduced[fcn] * spline_value;
          weights_grad     += controls_weights[fcn] * spline_grad;
          weights_hessian  += controls_weights[fcn] * spline_grad_grad;
          for ( unsigned int d = 0; d < space_dimension; d++){
            controls_grad[d]    += controls_reduced[fcn][d] * spline_grad;
            controls_hessian[d] += controls_reduced[fcn][d] * spline_grad_grad;
          }
        } // for ( fcn )

        // store the mapped quadrature point:
        for (unsigned int d = 0; d < space_dimension; d++)
          mapped_quadrature_point(d) = controls_value[d] / weight;

        D1 = weight * weight;
        D2 = D1 * D1;
        for (unsigned int d = 0; d < space_dimension; d++){
          N1[d] = controls_value[d] * weights_grad;
          N2[d] = controls_grad[d]  * weight;

          Tensor<2, dimension>& Hk = mapping_grad_grad[d];

          for (unsigned int m = 0; m < dimension; m++){
            for (unsigned int n = m; n < dimension; n++){
              Hk[n][m] = ( ( controls_hessian[d][n][m] * weight
                                 + controls_grad[d][m] * weights_grad[n]
                                 - controls_grad[d][n] * weights_grad[m]
                                 - controls_value[d] * weights_hessian[n][m] ) * D1 
                               - 2. * (controls_grad[d][m] * weight 
                                         - controls_value[d] * weights_grad[m]) * 
                                         weights_grad[n] * weight ) / D2;
              Hk[m][n] = Hk[n][m];
            } // for ( n )
          } // for ( m )
        } // for ( d )

        // Get value of Jacobian at current gauss point
        mapping_grad = DerivativeForm<1, dim, spacedim>( (1. / D1) *
                N2.Tensor<1, spacedim, Tensor<1, dimension> >::operator-=(N1) );
        // At the end of scope, all the helper variables are not needed anymore and can be destroyed
      } // scope 

      // Store the inverse of the jacobian at the current point. The covariant of a 
      // DerivativeForm is given by C = DF * (DF^T * DF)^{-1}. Hence, taking the 
      // transpose we get C^T = (DF^T * DF)^{-1} DF^T which is exactly the left-inverse
      // at the current point pnt. The function covariant_form() simplifies in case 
      // of rectangular DerivativeForms, i.e. if dim == spacedim, and DF^{-T} is returned. 
      const DerivativeForm<1, dimension, space_dimension>& covariant = 
              mapping_grad.covariant_form();
      inverse_mapping_grad = covariant.transpose();

      normal = s * inverse_mapping_grad[face_no/2] /
                   inverse_mapping_grad[face_no/2].norm();

      long double bf;
      if (dimension == 2){
        const unsigned int d = face_no/2 == 0 ? 1 : 0;
        const Tensor<1, space_dimension> tangent = 
                mapping_grad.transpose().operator[](d);
        if (space_dimension == 2) {
          bf = tangent.norm() ;
        } else {
          bf  = cross_product_3d(tangent, normal).norm();
        }
      } else if (dimension == 3) {
        Assert(space_dimension == 3, ExcInternalError());
        const unsigned int d1 = (face_no/2 == 0) ? 1 : 0;
        const unsigned int d2 = (d1 == 1 || face_no/2 == 1) ? 2 : 1;

        const Tensor<1, space_dimension> tangent1 =
                mapping_grad.transpose().operator[](d1);
        const Tensor<1, space_dimension> tangent2 =
                mapping_grad.transpose().operator[](d2);

        bf = cross_product_3d(tangent1, tangent2).norm();
      } else {
        Assert(false, ExcNotImplemented());
      }

      JxW  =  quadrature_weight *     // quadrature weight
              measure *               // transformation from reference to parametric cell
              bf;                     // boundary form dS on physical cell 

      // Using the Jacobians and the inverse of the jacobians we can setup the hessian 
      // of the transformation in terms of Tensors. Note, that we have 
      // HI = - I . H[I] . I, 
      // with the following data-types
      //    HI   : Tensor<1, dim, Tensor<2, spacedim> >, i.e. a size dim array with rank-2 tensors
      //    I    : DerivativeForm<1, spacedim, dim> ~ Tensor<1, dim, Tensor<1, spacedim>>, i.e. 
      //           a size dim array, where each entry is a Tensor of length spacedim
      //    H    : Tensor<1, spacedim, Tensor<2, dim> >, i.e. a size spacedim array with rank-2
      //           tensors as entries. 
      //    H[I] : Tensor<1, spacedim, Tensor<1, dim, Tensor<1, spacedim> > >, this object refers
      //           mathematically to the application of H to I, it is a rank-3 Tensor where each 
      //           entry refers to a dim by spacedim matrix, and the data type we use is according 
      //           to this
      // It is important to note, that the dots in the multiplication refer to tensor products
      // and not the typical matrix products!
      for (unsigned int k = 0; k < dimension; k++){
        Tensor<2, space_dimension>& HIk = inverse_mapping_grad_grad[k];
        // first build the sum over hessians of forward
        // mapping
        Tensor<2, dimension> temp;
        for (unsigned int d = 0; d < space_dimension; d++)
          temp += -1. * inverse_mapping_grad[k][d] * mapping_grad_grad[d];
        //    mapping_grad_grad     ~ Tensor<1, spacedim, Tensor<2, dim>>
        // -> mapping_grad_grad[d]  ~ Tensor<2, dim>
        //    mapping_grad          ~ DerivativeForm<1, dim, spacedim>
        // -> mapping_grad[k]       ~ Tensor<1, dim>
        // -> mapping_grad[k][d] * mapping_grad_grad[d] ~ Tensor<2, dim>

        // perform matrix-matrix product
        const DerivativeForm<1, space_dimension, dimension>& Ipnt = inverse_mapping_grad;
        if (dimension == 2) {
          HIk[0][0] = Ipnt[0][0] * (temp[0][0] * Ipnt[0][0] + temp[0][1] * Ipnt[1][0]) 
                    + Ipnt[1][0] * (temp[1][0] * Ipnt[0][0] + temp[1][1] * Ipnt[1][0]);
          HIk[0][1] = Ipnt[0][0] * (temp[0][0] * Ipnt[0][1] + temp[0][1] * Ipnt[1][1]) 
                    + Ipnt[1][0] * (temp[1][0] * Ipnt[0][1] + temp[1][1] * Ipnt[1][1]);
          HIk[1][1] = Ipnt[0][1] * (temp[0][0] * Ipnt[0][1] + temp[0][1] * Ipnt[1][1]) 
                    + Ipnt[1][1] * (temp[1][0] * Ipnt[0][1] + temp[1][1] * Ipnt[1][1]);
          HIk[1][0] = HIk[0][1];
        } else if (dimension == 3) {
          HIk[0][0] = Ipnt[0][0] * (
                        temp[0][0] * Ipnt[0][0] + temp[0][1] * Ipnt[1][0] + temp[0][2] * Ipnt[2][0]
                      ) 
                      + Ipnt[1][0] * (
                        temp[1][0] * Ipnt[0][0] + temp[1][1] * Ipnt[1][0] + temp[1][2] * Ipnt[2][0]
                      ) 
                      + Ipnt[2][0] * (
                        temp[2][0] * Ipnt[0][0] + temp[2][1] * Ipnt[1][0] + temp[2][2] * Ipnt[2][0]
                      );

          HIk[1][0] = Ipnt[0][1] * (
                        temp[0][0] * Ipnt[0][0] + temp[0][1] * Ipnt[1][0] + temp[0][2] * Ipnt[2][0]
                      ) 
                      + Ipnt[1][1] * (
                        temp[1][0] * Ipnt[0][0] + temp[1][1] * Ipnt[1][0] + temp[1][2] * Ipnt[2][0]
                      ) 
                      + Ipnt[2][1] * (
                        temp[2][0] * Ipnt[0][0] + temp[2][1] * Ipnt[1][0] + temp[2][2] * Ipnt[2][0]
                      );

          HIk[2][0] = Ipnt[0][2] * (
                        temp[0][0] * Ipnt[0][0] + temp[0][1] * Ipnt[1][0] + temp[0][2] * Ipnt[2][0]
                      ) 
                      + Ipnt[1][2] * (
                        temp[1][0] * Ipnt[0][0] + temp[1][1] * Ipnt[1][0] + temp[1][2] * Ipnt[2][0]
                      ) 
                      + Ipnt[2][2] * (
                        temp[2][0] * Ipnt[0][0] + temp[2][1] * Ipnt[1][0] + temp[2][2] * Ipnt[2][0]
                      );

          HIk[1][1] = Ipnt[0][1] * (
                        temp[0][0] * Ipnt[0][1] + temp[0][1] * Ipnt[1][1] + temp[0][2] * Ipnt[2][1]
                      ) 
                      + Ipnt[1][1] * (
                        temp[1][0] * Ipnt[0][1] + temp[1][1] * Ipnt[1][1] + temp[1][2] * Ipnt[2][1]
                      ) 
                      + Ipnt[2][1] * (
                        temp[2][0] * Ipnt[0][1] + temp[2][1] * Ipnt[1][1] + temp[2][2] * Ipnt[2][1]
                      );

          HIk[2][1] = Ipnt[0][2] * (
                        temp[0][0] * Ipnt[0][1] + temp[0][1] * Ipnt[1][1] + temp[0][2] * Ipnt[2][1]
                      ) 
                      + Ipnt[1][2] * (
                        temp[1][0] * Ipnt[0][1] + temp[1][1] * Ipnt[1][1] + temp[1][2] * Ipnt[2][1]
                      ) 
                      + Ipnt[2][2] * (
                        temp[2][0] * Ipnt[0][1] + temp[2][1] * Ipnt[1][1] + temp[2][2] * Ipnt[2][1]
                      );

          HIk[2][2] = Ipnt[0][2] * (
                        temp[0][0] * Ipnt[0][2] + temp[0][1] * Ipnt[1][2] + temp[0][2] * Ipnt[2][2]
                      ) 
                      + Ipnt[1][2] * (
                        temp[1][0] * Ipnt[0][2] + temp[1][1] * Ipnt[1][2] + temp[1][2] * Ipnt[2][2]
                      ) 
                      + Ipnt[2][2] * (
                        temp[2][0] * Ipnt[0][2] + temp[2][1] * Ipnt[1][2] + temp[2][2] * Ipnt[2][2]
                      );

          HIk[0][1] = HIk[1][0];
          HIk[0][2] = HIk[2][0];
          HIk[1][2] = HIk[2][1];
        }

        // The above implementation is rather lengthy, simple for-loops
        // may help. However, I am currently unable to find the mistake in 
        // the code below. In theory, there should be 
        //  HIk = I^{-T} temp  I,
        // which we try to resemble in the code below. 
        // for (unsigned int i = 0; i < space_dimension; i++){
        //   for (unsigned int j = 0; j <= i; j++) {
        //     double&                     HIk_ij = HIk[i][j];
        //     const Tensor<1,       dimension>& grad_i = covariant[i];
        //     const Tensor<1, space_dimension>& grad_j = inverse_mapping_grad[j];
        //     //    covariant               ~ DerivativeForm<1, dim, spacedim>
        //     // -> covariant[i]            ~ Tensor<1, dim>
        //     //    inverse_mapping_grad    ~ DerivativeForm<1, spacedim, dim>
        //     // -> inverse_mapping_grad[j] ~ Tensor<1, spacedim>
        //     for (unsigned int m = 0; m < space_dimension; m++)
        //       for (unsigned int n = 0; n < dimension; n++)
        //         HIk_ij += grad_i[n] * temp[n][m] * grad_j[m];

        //     HIk[j][i] = HIk_ij;
        //   } // for ( j )
        // } // for ( i )
      } // for ( k )

      // Using the inverse hessians, we can then compute the hessian of the 
      // shape functions in physical coordinates from the hessians in parametric
      // coordinates and the hessians of the inverse functions of Phi, using the
      // formula
      // ^H_[T o Phi] = I^tr * H_[T] * I + sum (\partial T / \partial x_k) * H_[Phi^{-1}_k] 
      //              = I^tr * H_[T] * I + sum grad_[T]_k * H_[Phi_k]
      for (unsigned int fcn = 0; fcn < dofs_per_cell; fcn++){
        ts_values(pnt, fcn)                          = spline_values(fcn);
        const Tensor<1, dimension>& spline_grad      = spline_grads(fcn);
        const Tensor<2, dimension>& spline_grad_grad = spline_grad_grads(fcn);
        
        // J_Phi * grad T_i
        if (flags & update_gradients)
          ts_grads(pnt, fcn) = apply_transformation(covariant, spline_grad);

        // For the hessian get two intermediate results
        // First, the latter sum of hessians
        Tensor<2, spacedim> sum;
        for (unsigned int d = 0; d < dimension; d++)
          sum += spline_grad[d] * inverse_mapping_grad_grad[d];

        // Second, the first matrix-matrix product
        Tensor<1, spacedim, Tensor<1, dim> > ip; 
        for (unsigned int k = 0; k < space_dimension; k++) {
          for (unsigned int d = 0; d < dimension; d++)
            ip[k][d] = spline_grad_grad[d] * covariant[k];
        }

        // Use the first matrix-matrix product, to compute the first term of ^H_[T o Phi]
        Tensor<2, spacedim>& h = ts_grad_grads(pnt, fcn);
        for (unsigned int i = 0; i < space_dimension; i++)
          for (unsigned int k = 0; k < space_dimension; k++)
            h[i][k] = (double)(covariant[i] * ip[k]);

        // Add the results: 
        h += sum;

//        std::cout << "Hessian: " << std::endl;
//        for (unsigned int i = 0; i < spacedim; i++)
//          std::cout << std::scientific << h[i] << std::endl;
//
//        // Test for symmetry:
//        for (unsigned int i = 0; i < spacedim; i++)
//          for (unsigned int j = 0; j < spacedim; j++)
//            Assert(h[i][j] == h[j][i], ExcInternalError());

      } // for ( fcn )
      if (flags & update_quadrature_points)
        this -> mapped_quadrature_points[pnt] = mapped_quadrature_point;

      if (flags & update_jacobians)
        this -> J[pnt] = mapping_grad;

      if (flags & update_inverse_jacobians)
        this -> I[pnt] = inverse_mapping_grad; 

      if (flags & update_jacobian_grads){
        this -> H[pnt]  = mapping_grad_grad;
        this -> HI[pnt] = inverse_mapping_grad_grad;
      }

      if (flags & update_JxW_values)
        this -> dx[pnt] = JxW;

      if (flags & update_normal_vectors)
        this -> normals[pnt] = normal; 
    } // for ( pnt )
  } // set_grad_and_hessian_tables

  template<int dim, int spacedim, int soldim>
  void TSFaceValues<dim, spacedim, soldim>::set_grad_tables(
    const active_cell_iterator& cell,
    const unsigned int sub_face
  ) {
    const unsigned int& dofs_per_cell = this -> dofs_per_cell;
    const          int& face_no       = this -> face_no;
    const unsigned int& offset        = this -> offset;
    const UpdateFlags& flags          = this -> flags;

    const Table<2, double>&                 bernstein_values    = this -> bernstein_values;
    const Table<2, Tensor<1, dimension> >&  bernstein_grads     = this -> bernstein_grads;
    // Get coefficients for linear combination
    const std::vector<FullMatrix<double>>& coeffs_vector = 
            this -> tria -> get_bezier_coefficients(
                        cell,
                        face_no
                    );
    const FullMatrix<double>& coeffs = coeffs_vector[sub_face];
    // Get weighted control points for geometric mapping
    const std::vector<Point<spacedim+1>>  control_points = 
            this -> tria -> get_control_points(cell);

    // store the sign for the current face
    // Note that this might be unused
    const int s = GeometryInfo<dim>::unit_normal_orientation[face_no];

    // store intermediate values 
    Table<1, double>                 spline_values(dofs_per_cell);
    Table<1, Tensor<1, dimension>>   spline_grads(dofs_per_cell); 

    // Split weighted control points
    std::vector< Tensor<1, spacedim >> controls_reduced(dofs_per_cell);
    std::vector< double >              controls_weights(dofs_per_cell);
    for (unsigned int i = 0; i < dofs_per_cell; i++){
      controls_weights[i] = control_points[i](space_dimension);
      for (unsigned int d = 0; d < space_dimension; d++)
        controls_reduced[i][d] = control_points[i](d);
    } // for ( i )

    // Get a reference to full tables for ts values, grads, and grad grads
    Table<2, double>& ts_values =
      this -> shape_values;
    Table<2, Tensor<1, space_dimension>>& ts_grads =
      this -> shape_grads;

    // For each quadrature point, compute the corresponding entries in the
    // tables.
    const Point<dimension>& lengths = cell -> vertex(GeometryInfo<dimension>::vertices_per_cell-1) 
                                      + (-1.) * cell -> vertex(0);
    const double            measure = (cell -> measure()) / lengths(face_no/2);
    // Initialize variables to be computed along the way
    Point<spacedim>                                   mapped_quadrature_point ;
    DerivativeForm<1, dimension, space_dimension>     mapping_grad ;
    DerivativeForm<1, space_dimension, dimension>     inverse_mapping_grad ;
    Tensor<1, space_dimension>                        normal;
    double                                            JxW;
    for (unsigned int p = 0; p < this -> shape_quadrature_points; p++) {
      const unsigned int pnt = p + sub_face * this->shape_quadrature_points;
      // Get a reference to corresponding table entries
      const double& quadrature_weight = this -> quadrature_weights[pnt];

      // Reset temporary values
      spline_values.reset_values();
      spline_grads.reset_values();

      { // Scope section to clear helper variables
        // helper variables
        Tensor<1, space_dimension, Tensor<1, dimension> > N1, N2;
        double D1;

        Tensor<1, space_dimension>                        controls_value;
        Tensor<1, dimension>                              weights_grad;
        Tensor<1, space_dimension, Tensor<1, dimension>>  controls_grad;
        double                                            weight = 0;


        for (unsigned int fcn = 0; fcn < dofs_per_cell; fcn++){
          // Get references to the values to be updated
          double&               spline_value = spline_values(fcn);
          Tensor<1, dimension>& spline_grad  = spline_grads(fcn);

          // run linear combination
          for (unsigned int fcn_b = 0; fcn_b < dofs_per_cell; fcn_b++){
            spline_value += coeffs(fcn, fcn_b) * 
                            bernstein_values(pnt + offset, fcn_b);
            spline_grad  += coeffs(fcn, fcn_b) * 
                            bernstein_grads(pnt + offset, fcn_b);
          }
          // Update values for reference cell
          for (unsigned int d = 0; d < dimension; d++)
            spline_grad[d] /= ( lengths(d) );

          weight           += controls_weights[fcn] * spline_value;
          controls_value   += controls_reduced[fcn] * spline_value;
          weights_grad     += controls_weights[fcn] * spline_grad;
          for ( unsigned int d = 0; d < space_dimension; d++)
            controls_grad[d] += controls_reduced[fcn][d] * spline_grad;
        } // for ( fcn )

        // store the mapped quadrature point:
        for (unsigned int d = 0; d < space_dimension; d++)
          mapped_quadrature_point(d) = controls_value[d] / weight;

        D1 = weight * weight;
        for (unsigned int d = 0; d < space_dimension; d++){
          N1[d] = controls_value[d] * weights_grad;
          N2[d] = controls_grad[d]  * weight;
        } // for ( d )

        // Get value of Jacobian at current gauss point
        mapping_grad = DerivativeForm<1, dim, spacedim>( (1. / D1) *
                N2.Tensor<1, spacedim, Tensor<1, dimension> >::operator-=(N1) );
        // At the end of scope, all the helper variables are not needed anymore and can be destroyed
      } // scope 

      // Store the inverse of the jacobian at the current point. The covariant of a 
      // DerivativeForm is given by C = DF * (DF^T * DF)^{-1}. Hence, taking the 
      // transpose we get C^T = (DF^T * DF)^{-1} DF^T which is exactly the left-inverse
      // at the current point pnt. The function covariant_form() simplifies in case 
      // of rectangular DerivativeForms, i.e. if dim == spacedim, and DF^{-T} is returned. 
      const DerivativeForm<1, dimension, space_dimension>& covariant = 
              mapping_grad.covariant_form();
      inverse_mapping_grad = covariant.transpose();

      normal = s * inverse_mapping_grad[face_no/2] /
                   inverse_mapping_grad[face_no/2].norm();

      long double bf;
      if (dimension == 2){
        const unsigned int d = face_no/2 == 0 ? 1 : 0;
        const Tensor<1, space_dimension> tangent = 
                mapping_grad.transpose().operator[](d);
        if (space_dimension == 2) {
          bf = tangent.norm() ;
        } else {
          bf  = cross_product_3d(tangent, normal).norm();
        }
      } else if (dimension == 3) {
        Assert(space_dimension == 3, ExcInternalError());
        const unsigned int d1 = (face_no/2 == 0) ? 1 : 0;
        const unsigned int d2 = (d1 == 1 || face_no/2 == 1) ? 2 : 1;

        const Tensor<1, space_dimension> tangent1 =
                mapping_grad.transpose().operator[](d1);
        const Tensor<1, space_dimension> tangent2 =
                mapping_grad.transpose().operator[](d2);

        bf = cross_product_3d(tangent1, tangent2).norm();
      } else {
        Assert(false, ExcNotImplemented());
      }

      JxW  =  quadrature_weight *     // quadrature weight
              measure *               // transformation from reference to parametric cell
              bf;                     // boundary form dS on physical cell 


      for (unsigned int fcn = 0; fcn < dofs_per_cell; fcn++){
        ts_values(pnt, fcn) = spline_values(fcn);
        
        // J_Phi * grad T_i
        if (flags & update_gradients)
          ts_grads(pnt, fcn)  = apply_transformation(covariant, spline_grads(fcn));
      } // for ( fcn )
      if (flags & update_quadrature_points)
        this -> mapped_quadrature_points[pnt] = mapped_quadrature_point;
      if (flags & update_jacobians)
        this -> J[pnt] = mapping_grad;
      if (flags & update_inverse_jacobians)
        this -> I[pnt] = inverse_mapping_grad; 
      if (flags & update_JxW_values)
        this -> dx[pnt] = JxW;
      if (flags & update_normal_vectors)
        this -> normals[pnt] = normal; 
    } // for ( pnt )
  } // set_grad_tables

  template<int dim, int spacedim, int soldim>
  void TSValues<dim, spacedim, soldim>::test_geometry_mapping(
    const Function<dim>      *phi_ptr    , 
    const Function<spacedim> *inv_phi_ptr, 
    const long double        &TOL,
    const bool                test_second_derivative
  ) { 
    std::cout << "Testing TSValues... " << std::endl;
    bool success = true; 
    bool test_inverse = inv_phi_ptr != NULL;
    const UpdateFlags old_flags = this -> flags;
    const Functions::ZeroFunction<spacedim> dummy; 
    const Function<dim>         &phi      = *phi_ptr;
    const Function<spacedim>    &inv_phi  = (test_inverse ? *inv_phi_ptr : dummy);
    std::vector< unsigned int >  degrees  = this -> tria->get_degree();
    for (unsigned int& p : degrees) p++;
    const unsigned int nvc = GeometryInfo<dimension>::vertices_per_cell;

    UpdateFlags &flags = this -> flags; 
    flags = update_values | update_gradients | update_quadrature_points |
                        update_JxW_values | update_jacobians;
    if (test_inverse)
      flags = flags | update_inverse_jacobians;

    this -> n_gauss_points = degrees; 
    
    // Setup tables for the test
    if (test_second_derivative){
      flags = flags | update_hessians | update_jacobian_grads;
      this -> init_bernstein_tables();
      this -> set_bernstein_tables();
      this -> init_shape_tables(); 
    } else {
      this -> init_bernstein_tables();
      this -> set_bernstein_tables();
      this -> init_shape_tables(); 
    }

    // ToDo: Once clear how exceptions are generated, encase this in a try catch block
    FILE *forward_result, *inverse_result;
    forward_result = fopen("log/test_tsvalues_push_forward.txt", "w");
    inverse_result = fopen("log/test_tsvalues_pull_back.txt", "w");

    fputs(" ============================================================= \n", forward_result);
    fputs(" ============================================================= \n", forward_result);
    fputs("           Results on push forward mapping\n", forward_result);
    fputs(" ============================================================= \n", forward_result);
    fputs(" ============================================================= \n\n\n", forward_result);
    
    if (test_inverse) {
      fputs(" ============================================================= \n", inverse_result);
      fputs(" ============================================================= \n", inverse_result);
      fputs("           Results on pull back mapping\n", inverse_result);
      fputs(" ============================================================= \n", inverse_result);
      fputs(" ============================================================= \n\n\n", inverse_result);
    }

    // 1. Write a summary of geometry values
    // Note that the inverse values need not be checked at this stage
    fputs( " ============================================================= \n", forward_result);
    fputs( "           Results on geometry values\n", forward_result);
    fputs( " ============================================================= \n", forward_result);
    for ( const auto& cell : this->tria->active_cell_iterators() ){
      const std::string cell_str = std::to_string(cell -> level()) + "." + std::to_string(cell -> index());
      const std::string ResultOnCell = " >>>>>>>>>>>>>>> Results on cell "  + cell_str + " <<<<<<<<<<<<<<<\n\n";
      fputs( ResultOnCell.c_str(), forward_result);
                  
      this -> reinit(cell); 

      const Point<dimension>& b = cell -> vertex(0);
      const Point<dimension>& A = cell -> vertex(nvc - 1) + -1. * cell -> vertex(0);

      Point<dimension>        parametric_Q;
      Point<space_dimension>  sh_physical_Q_push_forward;
      Point<space_dimension>  e_push_forward; 
      for (const unsigned int q : this -> quadrature_point_indices()){
        const Point<dimension>&       reference_Q                = this -> quadrature_point_reference(q);
        const Point<space_dimension>& is_physical_Q_push_forward = this -> quadrature_point(q);
              
        std::string ref_q_pnt_string = "At reference quadrature point " 
                                        + std::to_string(q) 
                                        + " with value (";
        for (unsigned int d = 0; d < dimension; d++)
          ref_q_pnt_string += std::to_string(reference_Q(d)) + " ";
        ref_q_pnt_string += ")\n";

        fputs( ref_q_pnt_string.c_str(), forward_result);

        for (unsigned int d = 0; d < dimension; d++)
          parametric_Q(d) = A(d) * reference_Q(d) + b(d);

        double norm = 0;
        for (unsigned int d = 0; d < space_dimension; d++) {
          sh_physical_Q_push_forward(d) = phi.value(parametric_Q, d);
          e_push_forward(d) = std::fabs(is_physical_Q_push_forward(d) - sh_physical_Q_push_forward(d)); 
          if (e_push_forward(d) > TOL){
            // Input error to forward results as summary of which entries are false
            const std::string e_str = "e" + std::to_string(d) + ": ";
            fputs( e_str.c_str(), forward_result);
            fprintf(forward_result, "%+1.8e\n", e_push_forward(d));
            norm = std::max(e_push_forward(d), norm);
          }
        }
        
        if (norm > TOL){
          success = false;
          fputs( "IS:     ", forward_result) ;
          for (unsigned int d = 0; d < space_dimension; d++)
            fprintf(forward_result, "%+1.8e ", is_physical_Q_push_forward(d));
          fputs( "\nSHOULD: ", forward_result) ;
          for (unsigned int d = 0; d < space_dimension; d++)
            fprintf(forward_result, "%+1.8e ", sh_physical_Q_push_forward(d));
          fputs("\n", forward_result);
        } else {
           fputs("No errors. \n", forward_result);
        }
        fputs("\n", forward_result);
      } // for ( q )
    } // for( cell ) // summary of geometry values

    // 2. Write a summary of geometry first derivatives
    fputs( " ============================================================= \n", forward_result);
    fputs( "           Results on geometry derivatives\n", forward_result);
    fputs( " ============================================================= \n", forward_result);
    if (test_inverse) {
      fputs( " ============================================================= \n", inverse_result);
      fputs( "           Results on geometry derivatives\n", inverse_result);
      fputs( " ============================================================= \n", inverse_result);
    }
    for ( const auto& cell : this->tria->active_cell_iterators() ){
      const std::string cell_str = std::to_string(cell -> level()) + "." + std::to_string(cell -> index());
      const std::string ResultOnCell = " >>>>>>>>>>>>>>> Results on cell "  + cell_str + " <<<<<<<<<<<<<<<\n\n";
      fputs( ResultOnCell.c_str(), forward_result);
      if (test_inverse)
        fputs( ResultOnCell.c_str(), inverse_result);

      this -> reinit(cell);
                  
      const Point<dimension>& b = cell -> vertex(0);
      const Point<dimension>& A = cell -> vertex(nvc - 1) + -1. * cell -> vertex(0);

      Point<space_dimension>  sh_physical_Q_push_forward;
      Point<dimension>        parametric_Q;
      Tensor<1, space_dimension, Tensor<1, dimension>> sh_push_forward_gradient;
      Tensor<1, dimension, Tensor<1, space_dimension>> sh_pull_back_gradient;
      Tensor<1, space_dimension, Tensor<1, dimension>> e_push_forward; 
      Tensor<1, dimension, Tensor<1, space_dimension>> e_pull_back;
      DerivativeForm<1, space_dimension, dimension>    is_pull_back_gradient;
      for (const unsigned int q : this -> quadrature_point_indices()){
        // Need points for evaluation
        const Point<dimension>&       reference_Q                = this -> quadrature_point_reference(q);
        for (unsigned int d = 0; d < dimension; d++)
          parametric_Q(d) = A(d) * reference_Q(d) + b(d);

        // Get values to check for correctness
        const DerivativeForm<1, dim, spacedim>& is_push_forward_gradient = this -> jacobian(q);
        if (test_inverse) 
          is_pull_back_gradient = this -> inverse_jacobian(q);
              

        std::string ref_q_pnt_string = "At reference quadrature point " 
                                        + std::to_string(q) 
                                        + " with value (";
        for (unsigned int d = 0; d < dimension; d++)
          ref_q_pnt_string += std::to_string(reference_Q(d)) + " ";
        ref_q_pnt_string += ")\n";
        fputs( ref_q_pnt_string.c_str(), forward_result);
        if (test_inverse) {
          fputs( ref_q_pnt_string.c_str(), inverse_result);
        }


        double norm = 0;
        for (unsigned int d = 0; d < space_dimension; d++) {
          // We need that point for the inverse
          if (test_inverse)
            sh_physical_Q_push_forward(d) = phi.value(parametric_Q, d);

          // This is to be checked in this section
          sh_push_forward_gradient[d]   = phi.gradient(parametric_Q, d);
          for (unsigned int dh = 0; dh < dimension; dh++){
            e_push_forward[d][dh] = std::fabs(sh_push_forward_gradient[d][dh] - is_push_forward_gradient[d][dh]); 
            norm = std::max(e_push_forward[d][dh], norm);
            if (e_push_forward[d][dh] > TOL){
              // Input error to forward results as summary of which entries are false
              const std::string e_str = "e" + std::to_string(d) + std::to_string(dh) + ": ";
              fputs(e_str.c_str(), forward_result);
              fprintf(forward_result, "%+1.8e\n", e_push_forward[d][dh]);
            }
          }
        }
        
        if (norm > TOL){
          success = false;
          fputs( "IS:     ", forward_result);
          for (unsigned int d = 0; d < space_dimension; d++) {
            for (unsigned int dd = 0; dd < dimension; dd++)
              fprintf(forward_result, "%+1.8e ", is_push_forward_gradient[d][dd]);
            fputs("\n        ", forward_result);
          }
          fputs( "\nSHOULD: ", forward_result);
          for (unsigned int d = 0; d < space_dimension; d++) {
            for (unsigned int dd = 0; dd < dimension; dd++)
              fprintf(forward_result, "%+1.8e ", sh_push_forward_gradient[d][dd]);
            fputs("\n        ", forward_result);
          }
          fputs("\n", forward_result);
        } else {
          fputs( "No errors. \n", forward_result);
        }

        if (test_inverse){
          norm = 0;
          for(unsigned int d = 0; d < dimension; d++) {
            sh_pull_back_gradient[d] = inv_phi.gradient(sh_physical_Q_push_forward, d);
            
            for (unsigned int dh = 0; dh < space_dimension; dh++){
              e_pull_back[d][dh] = std::fabs(sh_pull_back_gradient[d][dh] - is_pull_back_gradient[d][dh]);
              norm = std::max(e_pull_back[d][dh], norm);
              if (e_pull_back[d][dh] > TOL){
                fputs(("e" + std::to_string(d) + std::to_string(dh) + ": ").c_str(), inverse_result);
              fprintf(inverse_result, "%+1.8e\n", e_pull_back[d][dh]);
                                
              }
            }
          }
          if (norm > TOL){
            success = false;
            fputs( "IS:     ", inverse_result);
            for (unsigned int d = 0; d < space_dimension; d++) {
              for (unsigned int dd = 0; dd < dimension; dd++)
                fprintf(inverse_result, "%+1.8e ", is_pull_back_gradient[d][dd]);
              fputs("\n        ", inverse_result);
            }
            fputs( "\nSHOULD: ", inverse_result);
            for (unsigned int d = 0; d < space_dimension; d++) {
              for (unsigned int dd = 0; dd < dimension; dd++)
                fprintf(inverse_result, "%+1.8e ", sh_pull_back_gradient[d][dd]);
              fputs("\n        ", inverse_result);
            }
            fputs("\n", inverse_result);
          } else {
            fputs( "No errors. \n" , inverse_result);
          }
        } // if ( test_inverse )
        fputs("\n", forward_result);
        if (test_inverse)
          fputs("\n", inverse_result);

      } // for ( q )
      fputs("\n", forward_result);
      if (test_inverse)
        fputs("\n", inverse_result);
    } // for ( cell ) // test gradients

    // 3. Write a summary of geometry second derivatives
    if (test_second_derivative) { 
      fputs( " ============================================================= \n", forward_result);
      fputs( "           Results on geometry second derivatives\n", forward_result);
      fputs( " ============================================================= \n", forward_result);
      if (test_inverse) {
        fputs( " ============================================================= \n", inverse_result);
        fputs( "           Results on geometry second derivatives\n", inverse_result);
        fputs( " ============================================================= \n", inverse_result);
      }
      for (const auto& cell : this->tria->active_cell_iterators()) {
        const std::string cell_str = std::to_string(cell -> level()) + "." + std::to_string(cell -> index());
        const std::string ResultOnCell = " >>>>>>>>>>>>>>> Results on cell "  + cell_str + " <<<<<<<<<<<<<<<\n\n";
        fputs( ResultOnCell.c_str(), forward_result);
        if (test_inverse)
          fputs( ResultOnCell.c_str(), inverse_result);

        this -> reinit(cell);
                    
        const Point<dimension>& b = cell -> vertex(0);
        const Point<dimension>& A = cell -> vertex(nvc - 1) + -1. * cell -> vertex(0);

        Point<dimension>        parametric_Q;
        Point<space_dimension>  sh_physical_Q_push_forward;
        Tensor<1, space_dimension, Tensor<2, dimension>> sh_push_forward_grad_gradient;
        Tensor<1, dimension, Tensor<2, space_dimension>> sh_pull_back_grad_gradient;
        Tensor<1, space_dimension, Tensor<2, dimension>> e_push_forward; 
        Tensor<1, dimension, Tensor<2, space_dimension>> e_pull_back;
        Tensor<1, dimension, Tensor<2, space_dimension>> is_pull_back_grad_gradient;
        for (const unsigned int q : this->quadrature_point_indices()){
          // Need points for evaluation
          const Point<dimension>&       reference_Q                = this->quadrature_point_reference(q);
          for (unsigned int d = 0; d < dimension; d++)
            parametric_Q(d) = A(d) * reference_Q(d) + b(d);

          // Get values to check for correctness
          const Tensor<1, space_dimension, Tensor<2, dimension>>& 
                  is_push_forward_grad_gradient = this->jacobian_grad(q);
            
          if (test_inverse)
            is_pull_back_grad_gradient = this->inverse_jacobian_grad(q);
                
          std::string ref_q_pnt_string = "At reference quadrature point " 
                                          + std::to_string(q) 
                                          + " with value (";
          for (unsigned int d = 0; d < dimension; d++)
            ref_q_pnt_string += std::to_string(reference_Q(d)) + " ";
          ref_q_pnt_string += ")\n";
          fputs( ref_q_pnt_string.c_str(), forward_result);
          if (test_inverse) {
            fputs( ref_q_pnt_string.c_str(), inverse_result);
          }

          double norm = 0;
          for (unsigned int d = 0; d < space_dimension; d++) {
            // We need that point for the inverse
            if (test_inverse)
              sh_physical_Q_push_forward(d) = phi.value(parametric_Q, d);

            // This is to be checked in this section
            const SymmetricTensor<2, dimension>& tmp = phi.hessian(parametric_Q, d);
            for (unsigned int i = 0; i < dimension; i++) {
              for (unsigned int j = 0; j < i+1; j++){
                sh_push_forward_grad_gradient[d][i][j] = tmp[i][j];
                sh_push_forward_grad_gradient[d][j][i] = tmp[j][i];
              }
            }
            for (unsigned int dh1 = 0; dh1 < dimension; dh1++) {
              for (unsigned int dh2 = 0; dh2 <= dh1; dh2++) {
                e_push_forward[d][dh1][dh2] = 
                        std::fabs(sh_push_forward_grad_gradient[d][dh1][dh2] 
                                        - is_push_forward_grad_gradient[d][dh1][dh2]); 
                norm = std::max(e_push_forward[d][dh1][dh2], norm);
                if (e_push_forward[d][dh1][dh2] > TOL){
                  // Input error to forward results as summary of which entries are false
                  const std::string e_str = ("e" + std::to_string(d) + std::to_string(dh1) + std::to_string(dh2) + ": ");
                  fputs(e_str.c_str(), forward_result);
                  fprintf(forward_result, "%+1.8e\n", e_push_forward[d][dh1][dh2]);
                }
              }
            }
          } // for ( d )
          
          if (norm > TOL){
            success = false;
            std::string space ;
            for (unsigned int d = 0; d < dimension; d++)
              space = space + "         ";
            const std::string IsVShould = ("    IS:    " + space + "vs.   SHOULD: \n" );
            fputs( IsVShould.c_str(), forward_result);
            for (unsigned int d = 0; d < space_dimension; d++){
              fputs(" ----------------------------------------------------------- \n", forward_result);
              fputs( (std::to_string(d) + "|   ").c_str(), forward_result);
              for ( unsigned int dh = 0; dh < dimension-1; dh++){
                for (unsigned int dh1 = 0; dh1 < dimension; dh1++)
                  fprintf(forward_result, "%+1.8e ", is_push_forward_grad_gradient[d][dh][dh1]);
                fputs(( space + std::string("|")).c_str(), forward_result);

                for (unsigned int dh1 = 0; dh1 < dimension; dh1++)
                  fprintf(forward_result, "%+1.8e ", sh_push_forward_grad_gradient[d][dh][dh1]);
                fputs( "\n |   ", forward_result);
              }
              for (unsigned int dh1 = 0; dh1 < dimension; dh1++)
                fprintf(forward_result, "%+1.8e ", sh_push_forward_grad_gradient[d][dimension-1][dh1]);
              fputs( (space + std::string("|")).c_str(), forward_result);

              for (unsigned int dh1 = 0; dh1 < dimension; dh1++)
                fprintf(forward_result, "%+1.8e ", sh_push_forward_grad_gradient[d][dimension-1][dh1]);
              fputs("\n", forward_result);
            }
          } else {
            fputs( "No errors. \n" , forward_result);
          }

          if (test_inverse){
            norm = 0;
            for(unsigned int d = 0; d < dimension; d++) {
              const SymmetricTensor<2, space_dimension>& tmp = inv_phi.hessian(sh_physical_Q_push_forward, d);
              for (unsigned int i = 0; i < space_dimension; i++) {
                for (unsigned int j = 0; j < i+1; j++){
                  sh_pull_back_grad_gradient[d][i][j] = tmp[i][j];
                  sh_pull_back_grad_gradient[d][j][i] = tmp[j][i];
                }
              }
              
              for (unsigned int dh1 = 0; dh1 < space_dimension; dh1++){
                for (unsigned int dh2 = 0; dh2 <= dh1; dh2++) {
                  e_pull_back[d][dh1][dh2] = 
                    std::fabs(sh_pull_back_grad_gradient[d][dh1][dh2] - is_pull_back_grad_gradient[d][dh1][dh2]);
                  norm = std::max(e_pull_back[d][dh1][dh2], norm);
                  if (e_pull_back[d][dh1][dh2] > TOL){
                    const std::string e_str = ("e" + std::to_string(d) + std::to_string(dh1) + std::to_string(dh2) + ": ");
                    fputs(e_str.c_str(), inverse_result);
                    fprintf(inverse_result, "%+1.8e\n", e_pull_back[d][dh1][dh2]);
                                    
                  }
                }
              }
            }
            if (norm > TOL){
              success = false;
              std::string space ;
              for (unsigned int d = 0; d < space_dimension-1; d++)
                space = space + "         ";
              const std::string IsVShould = "    IS:    " + space + "vs.   SHOULD: \n";
              fputs( IsVShould.c_str(), inverse_result);
              for (unsigned int d = 0; d < dimension; d++){
                fputs(" ----------------------------------------------------------- \n", inverse_result);
                fputs( (std::to_string(d) + "|   ").c_str(), inverse_result);
                for ( unsigned int dh = 0; dh < dimension-1; dh++){
                  for (unsigned int dh1 = 0; dh1 < dimension; dh1++)
                    fprintf(inverse_result, "%+1.8e ", is_pull_back_grad_gradient[d][dh][dh1]);
                  fputs( (space + std::string("|")).c_str(), inverse_result);

                  for (unsigned int dh1 = 0; dh1 < dimension; dh1++)
                    fprintf(inverse_result, "%+1.8e ", sh_pull_back_grad_gradient[d][dh][dh1]);
                  fputs("\n |   ", inverse_result);
                }
                for (unsigned int dh1 = 0; dh1 < dimension; dh1++)
                  fprintf(inverse_result, "%+1.8e ", is_pull_back_grad_gradient[d][dimension-1][dh1]);
                fputs( (space + std::string( "|")).c_str(), inverse_result);

                for (unsigned int dh1 = 0; dh1 < dimension; dh1++)
                  fprintf(inverse_result, "%+1.8e ", sh_pull_back_grad_gradient[d][dimension-1][dh1]);
                fputs("\n" , inverse_result);
              }
            } else {
              fputs( "No errors. \n" , inverse_result);
            }
          } // if ( test_inverse )
          fputs("\n", forward_result);
          if (test_inverse)
            fputs("\n", inverse_result);
        } // for ( q )
        fputs("\n", forward_result);
        if (test_inverse)
          fputs("\n", inverse_result);
      } // for ( cell )
    }
    fclose(forward_result);
    fclose(inverse_result);
    
    Assert(success, ExcMessage("The calculated values with TSValues do not match the exact values!\n Checkout the log files for more information."));
#ifndef DEBUG
    if (!success )
      throw ExcMessage("The calculated values with TSValues do not match the exact values!\n Check the log files for more information."); 
#endif
    
    // Revert to the setup before testing
    this -> flags = old_flags;
    this -> init_bernstein_tables();
    this -> set_bernstein_tables();
    this -> init_shape_tables(); 
  } // test_geometry_mapping

  template<int dim, int spacedim, int soldim>
  void TSFaceValues<dim, spacedim, soldim>::test_geometry_mapping(
    const Function<dim>      *phi_ptr    , 
    const Function<spacedim> *inv_phi_ptr, 
    const long double        &TOL,
    const bool                test_second_derivative
  ) { 
          std::cout << "Testing TSFaceValues ... " << std::endl;
    bool success = true; 
    bool test_inverse = inv_phi_ptr != NULL;
    const UpdateFlags old_flags = this -> flags;
    const Functions::ZeroFunction<spacedim> dummy; 
    const Function<dim>         &phi      = *phi_ptr;
    const Function<spacedim>    &inv_phi  = (test_inverse ? *inv_phi_ptr : dummy);
    std::vector< unsigned int >  degrees  = this -> tria->get_degree();
    for (unsigned int& p : degrees) p++;
    const unsigned int nfc = GeometryInfo<dimension>::faces_per_cell;
    const unsigned int nvc = GeometryInfo<dimension>::vertices_per_cell;

    UpdateFlags &flags = this -> flags; 
    flags = update_values | update_gradients | update_quadrature_points |
            update_normal_vectors | update_JxW_values | update_jacobians;
    if (test_inverse)
      flags = flags | update_inverse_jacobians;

    this -> n_gauss_points = degrees; 
    
    // Setup tables for the test
    if (test_second_derivative){
      flags = flags | update_hessians | update_jacobian_grads;
      this -> init_bernstein_tables();
      this -> set_bernstein_tables();
      this -> init_shape_tables(); 
    } else {
      this -> init_bernstein_tables();
      this -> set_bernstein_tables();
      this -> init_shape_tables(); 
    }

    // ToDo: Once clear how exceptions are generated, encase this in a try catch block
    FILE *forward_result, *inverse_result;
    forward_result = fopen("log/test_tsfacevalues_push_forward.txt", "w");
    inverse_result = fopen("log/test_tsfacevalues_pull_back.txt", "w");

    fputs(" ============================================================= \n", forward_result);
    fputs(" ============================================================= \n", forward_result);
    fputs("           Results on push forward mapping\n", forward_result);
    fputs(" ============================================================= \n", forward_result);
    fputs(" ============================================================= \n\n\n", forward_result);
    
    if (test_inverse) {
      fputs(" ============================================================= \n", inverse_result);
      fputs(" ============================================================= \n", inverse_result);
      fputs("           Results on pull back mapping\n", inverse_result);
      fputs(" ============================================================= \n", inverse_result);
      fputs(" ============================================================= \n\n\n", inverse_result);
    }

    // 1. Write a summary of geometry values
    // Note that the inverse values need not be checked at this stage
    fputs( " ============================================================= \n", forward_result);
    fputs( "           Results on geometry values\n", forward_result);
    fputs( " ============================================================= \n", forward_result);
    for ( const auto& cell : this->tria->active_cell_iterators() ){
      for (unsigned int f = 0; f < nfc; f++){
        const std::string cell_str = std::to_string(cell -> level()) + "." + std::to_string(cell -> index()) 
                                      + " at face " + std::to_string(f);
        const std::string ResultOnCell = " >>>>>>>>>>>>>>> Results on cell "  + cell_str + " <<<<<<<<<<<<<<<\n\n";
        fputs( ResultOnCell.c_str(), forward_result);
                    
        this -> reinit(cell, f); 

        const Point<dimension>& b = cell -> vertex(0);
        const Point<dimension>& A = cell -> vertex(nvc - 1) + -1. * cell -> vertex(0);

        Point<dimension>        parametric_Q;
        Point<space_dimension>  sh_physical_Q_push_forward;
        Point<space_dimension>  e_push_forward; 
        for (const unsigned int q : this -> quadrature_point_indices()){
          const Point<dimension>&       reference_Q                = this -> quadrature_point_reference(q);
          const Point<space_dimension>& is_physical_Q_push_forward = this -> quadrature_point(q);
                
          std::string ref_q_pnt_string = "At reference quadrature point " 
                                          + std::to_string(q) 
                                          + " with value (";
          for (unsigned int d = 0; d < dimension; d++)
            ref_q_pnt_string += std::to_string(reference_Q(d)) + " ";
          ref_q_pnt_string += ")\n";

          fputs( ref_q_pnt_string.c_str(), forward_result);

          for (unsigned int d = 0; d < dimension; d++)
            parametric_Q(d) = A(d) * reference_Q(d) + b(d);

          double norm = 0;
          for (unsigned int d = 0; d < space_dimension; d++) {
            sh_physical_Q_push_forward(d) = phi.value(parametric_Q, d);
            e_push_forward(d) = std::fabs(is_physical_Q_push_forward(d) - sh_physical_Q_push_forward(d)); 
            if (e_push_forward(d) > TOL){
              // Input error to forward results as summary of which entries are false
              const std::string e_str = "e" + std::to_string(d) + ": ";
              fputs( e_str.c_str(), forward_result);
              fprintf(forward_result, "%+1.8e\n", e_push_forward(d));
              norm = std::max(e_push_forward(d), norm);
            }
          }
          
          if (norm > TOL){
            success = false;
            fputs( "IS:     ", forward_result) ;
            for (unsigned int d = 0; d < space_dimension; d++)
              fprintf(forward_result, "%+1.8e ", is_physical_Q_push_forward(d));
            fputs( "\nSHOULD: ", forward_result) ;
            for (unsigned int d = 0; d < space_dimension; d++)
              fprintf(forward_result, "%+1.8e ", sh_physical_Q_push_forward(d));
            fputs("\n", forward_result);
          } else {
             fputs("No errors. \n", forward_result);
          }
          fputs("\n", forward_result);
        } // for ( q )
      } // for ( f )
    } // for( cell ) // summary of geometry values

    if (test_inverse) {
      // 1.5 Write a summary of geometry normals
      // The normals are computed from the inverse mapping, but the results 
      // are printed to forward
      fputs( " ============================================================= \n", forward_result);
      fputs( "           Results on normals\n", forward_result);
      fputs( " ============================================================= \n", forward_result);
      for (const auto& cell : this -> tria -> active_cell_iterators()){
        for (unsigned int f = 0; f < nfc; f++){
          if (cell -> face(f) -> at_boundary()){
            const std::string cell_str = std::to_string(cell -> level()) + "." + std::to_string(cell -> index()) 
                                          + " at face " + std::to_string(f);
            const std::string ResultOnCell = " >>>>>>>>>>>>>>> Results on cell "  + cell_str + " <<<<<<<<<<<<<<<\n\n";
            fputs( ResultOnCell.c_str(), forward_result);
                        
            this -> reinit(cell, f); 

            Tensor<1, space_dimension> sh_normal; 
            Tensor<1, space_dimension> e_normal; 
            for (const unsigned int& q : this -> quadrature_point_indices()){
              const Point<dimension>&           reference_Q   = this -> quadrature_point_reference(q);
              const Point<space_dimension>&     physical_Q    = this -> quadrature_point(q);
              const Tensor<1, space_dimension>& is_normal     = this -> normal_vector(q);
                        
              std::string ref_q_pnt_string = "At reference quadrature point " 
                                              + std::to_string(q) 
                                              + " with value (";
              for (unsigned int d = 0; d < dimension; d++)
                ref_q_pnt_string += std::to_string(reference_Q(d)) + " ";
              ref_q_pnt_string += ")\n";
              fputs( ref_q_pnt_string.c_str(), forward_result);
              
              const int s = GeometryInfo<dim>::unit_normal_orientation[f];
              sh_normal = s * inv_phi.gradient(physical_Q, f/2) 
                              / inv_phi.gradient(physical_Q, f/2).norm();

              double norm = 0;
              for (unsigned int d = 0; d < space_dimension; d++){
                e_normal[d] = std::fabs(is_normal[d] - sh_normal[d]);
                norm = std::max(e_normal[d], norm);
                if (e_normal[d] > TOL) {
                  // Input error to forward results as summary of which entries are false
                  const std::string e_str = "e" + std::to_string(d) + ": ";
                  fputs(e_str.c_str(), forward_result);
                  fprintf(forward_result, "%+1.8e\n", e_normal[d]);
                } // if ( TOL )
              } // for ( d )

              if (norm > TOL){
                success = false;
                fputs( "IS:     ", forward_result);
                for (unsigned int d = 0; d < space_dimension; d++) {
                    fprintf(forward_result, "%+1.8e ", is_normal[d]);
                  fputs("\n        ", forward_result);
                }
                fputs( "\nSHOULD: ", forward_result);
                for (unsigned int d = 0; d < space_dimension; d++) {
                    fprintf(forward_result, "%+1.8e ", sh_normal[d]);
                  fputs("\n        ", forward_result);
                }
                fputs("\n", forward_result);
              } else {
                fputs("No errors. \n", forward_result);
              } // if ( norm > TOL )
            } // for ( q )
          } // if ( at_boundary() )
        } // for ( f )
      } // for ( cell )
    } // if (test_inverse)

    // 2. Write a summary of geometry first derivatives
    fputs( " ============================================================= \n", forward_result);
    fputs( "           Results on geometry derivatives\n", forward_result);
    fputs( " ============================================================= \n", forward_result);
    if (test_inverse) {
      fputs( " ============================================================= \n", inverse_result);
      fputs( "           Results on geometry derivatives\n", inverse_result);
      fputs( " ============================================================= \n", inverse_result);
    }
    for ( const auto& cell : this->tria->active_cell_iterators() ){
      for (unsigned int f = 0; f < nfc; f++) {
        const std::string cell_str = std::to_string(cell -> level()) + "." + std::to_string(cell -> index()) 
                                      + " at face " + std::to_string(f);
        const std::string ResultOnCell = " >>>>>>>>>>>>>>> Results on cell "  + cell_str + " <<<<<<<<<<<<<<<\n\n";
        fputs( ResultOnCell.c_str(), forward_result);
        if (test_inverse)
          fputs( ResultOnCell.c_str(), inverse_result);
                    
        this -> reinit(cell, f); 

        const Point<dimension>& b = cell -> vertex(0);
        const Point<dimension>& A = cell -> vertex(nvc - 1) + -1. * cell -> vertex(0);

        Point<space_dimension>  sh_physical_Q_push_forward;
        Point<dimension>        parametric_Q;
        Tensor<1, space_dimension, Tensor<1, dimension>> sh_push_forward_gradient;
        Tensor<1, dimension, Tensor<1, space_dimension>> sh_pull_back_gradient;
        Tensor<1, space_dimension, Tensor<1, dimension>> e_push_forward; 
        Tensor<1, dimension, Tensor<1, space_dimension>> e_pull_back;
        DerivativeForm<1, space_dimension, dimension>    is_pull_back_gradient;
        for (const unsigned int q : this -> quadrature_point_indices()){
          // Need points for evaluation
          const Point<dimension>&       reference_Q                = this -> quadrature_point_reference(q);
          for (unsigned int d = 0; d < dimension; d++)
            parametric_Q(d) = A(d) * reference_Q(d) + b(d);

          // Get values to check for correctness
          const DerivativeForm<1, dim, spacedim>& is_push_forward_gradient = this -> jacobian(q);
          if (test_inverse) 
            is_pull_back_gradient = this -> inverse_jacobian(q);
                

          std::string ref_q_pnt_string = "At reference quadrature point " 
                                          + std::to_string(q) 
                                          + " with value (";
          for (unsigned int d = 0; d < dimension; d++)
            ref_q_pnt_string += std::to_string(reference_Q(d)) + " ";
          ref_q_pnt_string += ")\n";
          fputs( ref_q_pnt_string.c_str(), forward_result);
          if (test_inverse) {
            fputs( ref_q_pnt_string.c_str(), inverse_result);
          }


          double norm = 0;
          for (unsigned int d = 0; d < space_dimension; d++) {
            // We need that point for the inverse
            if (test_inverse)
              sh_physical_Q_push_forward(d) = phi.value(parametric_Q, d);

            // This is to be checked in this section
            sh_push_forward_gradient[d]   = phi.gradient(parametric_Q, d);
            for (unsigned int dh = 0; dh < dimension; dh++){
              e_push_forward[d][dh] = std::fabs(sh_push_forward_gradient[d][dh] - is_push_forward_gradient[d][dh]); 
              norm = std::max(e_push_forward[d][dh], norm);
              if (e_push_forward[d][dh] > TOL){
                // Input error to forward results as summary of which entries are false
                const std::string e_str = "e" + std::to_string(d) + std::to_string(dh) + ": ";
                fputs(e_str.c_str(), forward_result);
                fprintf(forward_result, "%+1.8e\n", e_push_forward[d][dh]);
              }
            }
          }
          
          if (norm > TOL){
            success = false;
            fputs( "IS:     ", forward_result);
            for (unsigned int d = 0; d < space_dimension; d++) {
              for (unsigned int dd = 0; dd < dimension; dd++)
                fprintf(forward_result, "%+1.8e ", is_push_forward_gradient[d][dd]);
              fputs("\n        ", forward_result);
            }
            fputs( "\nSHOULD: ", forward_result);
            for (unsigned int d = 0; d < space_dimension; d++) {
              for (unsigned int dd = 0; dd < dimension; dd++)
                fprintf(forward_result, "%+1.8e ", sh_push_forward_gradient[d][dd]);
              fputs("\n        ", forward_result);
            }
            fputs("\n", forward_result);
          } else {
            fputs( "No errors. \n", forward_result);
          }

          if (test_inverse){
            norm = 0;
            for(unsigned int d = 0; d < dimension; d++) {
              sh_pull_back_gradient[d] = inv_phi.gradient(sh_physical_Q_push_forward, d);
              
              for (unsigned int dh = 0; dh < space_dimension; dh++){
                e_pull_back[d][dh] = std::fabs(sh_pull_back_gradient[d][dh] - is_pull_back_gradient[d][dh]);
                norm = std::max(e_pull_back[d][dh], norm);
                if (e_pull_back[d][dh] > TOL){
                  fputs(("e" + std::to_string(d) + std::to_string(dh) + ": ").c_str(), inverse_result);
                fprintf(inverse_result, "%+1.8e\n", e_pull_back[d][dh]);
                                  
                }
              }
            }
            if (norm > TOL){
              success = false;
              fputs( "IS:     ", inverse_result);
              for (unsigned int d = 0; d < space_dimension; d++) {
                for (unsigned int dd = 0; dd < dimension; dd++)
                  fprintf(inverse_result, "%+1.8e ", is_pull_back_gradient[d][dd]);
                fputs("\n        ", inverse_result);
              }
              fputs( "\nSHOULD: ", inverse_result);
              for (unsigned int d = 0; d < space_dimension; d++) {
                for (unsigned int dd = 0; dd < dimension; dd++)
                  fprintf(inverse_result, "%+1.8e ", sh_pull_back_gradient[d][dd]);
                fputs("\n        ", inverse_result);
              }
              fputs("\n", inverse_result);
            } else {
              fputs( "No errors. \n" , inverse_result);
            }
          } // if ( test_inverse )
          fputs("\n", forward_result);
          if (test_inverse)
            fputs("\n", inverse_result);

        } // for ( q )
        fputs("\n", forward_result);
        if (test_inverse)
          fputs("\n", inverse_result);
      } // for ( f ) 
    } // for ( cell ) // test gradients

    // 3. Write a summary of geometry second derivatives
    if (test_second_derivative) { 
      fputs( " ============================================================= \n", forward_result);
      fputs( "           Results on geometry second derivatives\n", forward_result);
      fputs( " ============================================================= \n", forward_result);
      if (test_inverse) {
        fputs( " ============================================================= \n", inverse_result);
        fputs( "           Results on geometry second derivatives\n", inverse_result);
        fputs( " ============================================================= \n", inverse_result);
      }
      for (const auto& cell : this->tria->active_cell_iterators()) {
        for ( unsigned int f = 0; f < nfc; f++) {
          const std::string cell_str = std::to_string(cell -> level()) + "." + std::to_string(cell -> index());
          const std::string ResultOnCell = " >>>>>>>>>>>>>>> Results on cell "  + cell_str + " <<<<<<<<<<<<<<<\n\n";
          fputs( ResultOnCell.c_str(), forward_result);
          if (test_inverse)
            fputs( ResultOnCell.c_str(), inverse_result);

          this -> reinit(cell, f); 
                      
          const Point<dimension>& b = cell -> vertex(0);
          const Point<dimension>& A = cell -> vertex(nvc - 1) + -1. * cell -> vertex(0);

          Point<dimension>        parametric_Q;
          Point<space_dimension>  sh_physical_Q_push_forward;
          Tensor<1, space_dimension, Tensor<2, dimension>> sh_push_forward_grad_gradient;
          Tensor<1, dimension, Tensor<2, space_dimension>> sh_pull_back_grad_gradient;
          Tensor<1, space_dimension, Tensor<2, dimension>> e_push_forward; 
          Tensor<1, dimension, Tensor<2, space_dimension>> e_pull_back;
          Tensor<1, dimension, Tensor<2, space_dimension>> is_pull_back_grad_gradient;
          for (const unsigned int q : this->quadrature_point_indices()){
            // Need points for evaluation
            const Point<dimension>&       reference_Q                = this->quadrature_point_reference(q);
            for (unsigned int d = 0; d < dimension; d++)
              parametric_Q(d) = A(d) * reference_Q(d) + b(d);

            // Get values to check for correctness
            const Tensor<1, space_dimension, Tensor<2, dimension>>& 
                    is_push_forward_grad_gradient = this->jacobian_grad(q);
              
            if (test_inverse)
              is_pull_back_grad_gradient = this->inverse_jacobian_grad(q);
                  
            std::string ref_q_pnt_string = "At reference quadrature point " 
                                            + std::to_string(q) 
                                            + " with value (";
            for (unsigned int d = 0; d < dimension; d++)
              ref_q_pnt_string += std::to_string(reference_Q(d)) + " ";
            ref_q_pnt_string += ")\n";
            fputs( ref_q_pnt_string.c_str(), forward_result);
            if (test_inverse) {
              fputs( ref_q_pnt_string.c_str(), inverse_result);
            }

            double norm = 0;
            for (unsigned int d = 0; d < space_dimension; d++) {
              // We need that point for the inverse
              if (test_inverse)
                sh_physical_Q_push_forward(d) = phi.value(parametric_Q, d);

              // This is to be checked in this section
              const SymmetricTensor<2, dimension>& tmp = phi.hessian(parametric_Q, d);
              for (unsigned int i = 0; i < dimension; i++) {
                for (unsigned int j = 0; j < i+1; j++){
                  sh_push_forward_grad_gradient[d][i][j] = tmp[i][j];
                  sh_push_forward_grad_gradient[d][j][i] = tmp[j][i];
                }
              }
              for (unsigned int dh1 = 0; dh1 < dimension; dh1++) {
                for (unsigned int dh2 = 0; dh2 <= dh1; dh2++) {
                  e_push_forward[d][dh1][dh2] = 
                          std::fabs(sh_push_forward_grad_gradient[d][dh1][dh2] 
                                          - is_push_forward_grad_gradient[d][dh1][dh2]); 
                  norm = std::max(e_push_forward[d][dh1][dh2], norm);
                  if (e_push_forward[d][dh1][dh2] > TOL){
                    // Input error to forward results as summary of which entries are false
                    const std::string e_str = ("e" + std::to_string(d) + std::to_string(dh1) + std::to_string(dh2) + ": ");
                    fputs(e_str.c_str(), forward_result);
                    fprintf(forward_result, "%+1.8e\n", e_push_forward[d][dh1][dh2]);
                  }
                }
              }
            } // for ( d )
            
            if (norm > TOL){
              success = false;
              std::string space ;
              for (unsigned int d = 0; d < dimension; d++)
                space = space + "         ";
              const std::string IsVShould = ("    IS:    " + space + "vs.   SHOULD: \n" );
              fputs( IsVShould.c_str(), forward_result);
              for (unsigned int d = 0; d < space_dimension; d++){
                fputs(" ----------------------------------------------------------- \n", forward_result);
                fputs( (std::to_string(d) + "|   ").c_str(), forward_result);
                for ( unsigned int dh = 0; dh < dimension-1; dh++){
                  for (unsigned int dh1 = 0; dh1 < dimension; dh1++)
                    fprintf(forward_result, "%+1.8e ", is_push_forward_grad_gradient[d][dh][dh1]);
                  fputs(( space + std::string("|")).c_str(), forward_result);

                  for (unsigned int dh1 = 0; dh1 < dimension; dh1++)
                    fprintf(forward_result, "%+1.8e ", sh_push_forward_grad_gradient[d][dh][dh1]);
                  fputs( "\n |   ", forward_result);
                }
                for (unsigned int dh1 = 0; dh1 < dimension; dh1++)
                  fprintf(forward_result, "%+1.8e ", sh_push_forward_grad_gradient[d][dimension-1][dh1]);
                fputs( (space + std::string("|")).c_str(), forward_result);

                for (unsigned int dh1 = 0; dh1 < dimension; dh1++)
                  fprintf(forward_result, "%+1.8e ", sh_push_forward_grad_gradient[d][dimension-1][dh1]);
                fputs("\n", forward_result);
              }
            } else {
              fputs( "No errors. \n" , forward_result);
            }

            if (test_inverse){
              norm = 0;
              for(unsigned int d = 0; d < dimension; d++) {
                const SymmetricTensor<2, space_dimension>& tmp = inv_phi.hessian(sh_physical_Q_push_forward, d);
                for (unsigned int i = 0; i < space_dimension; i++) {
                  for (unsigned int j = 0; j < i+1; j++){
                    sh_pull_back_grad_gradient[d][i][j] = tmp[i][j];
                    sh_pull_back_grad_gradient[d][j][i] = tmp[j][i];
                  }
                }
                
                for (unsigned int dh1 = 0; dh1 < space_dimension; dh1++){
                  for (unsigned int dh2 = 0; dh2 <= dh1; dh2++) {
                    e_pull_back[d][dh1][dh2] = 
                      std::fabs(sh_pull_back_grad_gradient[d][dh1][dh2] - is_pull_back_grad_gradient[d][dh1][dh2]);
                    norm = std::max(e_pull_back[d][dh1][dh2], norm);
                    if (e_pull_back[d][dh1][dh2] > TOL){
                      const std::string e_str = ("e" + std::to_string(d) + std::to_string(dh1) + std::to_string(dh2) + ": ");
                      fputs(e_str.c_str(), inverse_result);
                      fprintf(inverse_result, "%+1.8e\n", e_pull_back[d][dh1][dh2]);
                                      
                    }
                  }
                }
              }
              if (norm > TOL){
                success = false;
                std::string space ;
                for (unsigned int d = 0; d < space_dimension-1; d++)
                  space = space + "         ";
                const std::string IsVShould = "    IS:    " + space + "vs.   SHOULD: \n";
                fputs( IsVShould.c_str(), inverse_result);
                for (unsigned int d = 0; d < dimension; d++){
                  fputs(" ----------------------------------------------------------- \n", inverse_result);
                  fputs( (std::to_string(d) + "|   ").c_str(), inverse_result);
                  for ( unsigned int dh = 0; dh < dimension-1; dh++){
                    for (unsigned int dh1 = 0; dh1 < dimension; dh1++)
                      fprintf(inverse_result, "%+1.8e ", is_pull_back_grad_gradient[d][dh][dh1]);
                    fputs( (space + std::string("|")).c_str(), inverse_result);

                    for (unsigned int dh1 = 0; dh1 < dimension; dh1++)
                      fprintf(inverse_result, "%+1.8e ", sh_pull_back_grad_gradient[d][dh][dh1]);
                    fputs("\n |   ", inverse_result);
                  }
                  for (unsigned int dh1 = 0; dh1 < dimension; dh1++)
                    fprintf(inverse_result, "%+1.8e ", is_pull_back_grad_gradient[d][dimension-1][dh1]);
                  fputs( (space + std::string( "|")).c_str(), inverse_result);

                  for (unsigned int dh1 = 0; dh1 < dimension; dh1++)
                    fprintf(inverse_result, "%+1.8e ", sh_pull_back_grad_gradient[d][dimension-1][dh1]);
                  fputs("\n" , inverse_result);
                }
              } else {
                fputs( "No errors. \n" , inverse_result);
              }
            } // if ( test_inverse )
            fputs("\n", forward_result);
            if (test_inverse)
              fputs("\n", inverse_result);
          } // for ( q )
          fputs("\n", forward_result);
          if (test_inverse)
            fputs("\n", inverse_result);
        } // for ( f )
      } // for ( cell )
    } // if ( test_second_derivative )
    fclose(forward_result);
    fclose(inverse_result);
    
    Assert(success, ExcMessage("The calculated values with TSValues do not match the exact values!\n Check the log files for more information."));
#ifndef DEBUG
    if (!success )
      throw ExcMessage("The calculated values with TSValues do not match the exact values!\n Check the log files for more information."); 
#endif

    // Revert to the setup before testing
    this -> flags = old_flags;
    this -> init_bernstein_tables();
    this -> set_bernstein_tables();
    this -> init_shape_tables(); 
  } // test_geometry_mapping

#include "ts_values.inst.in"

} // namespace dealt
