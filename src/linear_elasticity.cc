#include <linear_elasticity.h>

namespace Linear_Elasticity {
  template class DataGenerator<2>;
  template class DataGenerator<3>;

  template class Elasticity_SOL<2>;
  template class Elasticity_SOL<3>;

  template class ElasticProblem_SOL<2>;
  template class ElasticProblem_SOL<3>;

  template class ElasticProblem_RHS<2>;
  template class ElasticProblem_RHS<3>;

  template class ElasticProblem_NC_X0<2>;
  template class ElasticProblem_NC_X0<3>;

  template class ElasticProblem_NC_Y0<2>;
  template class ElasticProblem_NC_Y0<3>;

  template class ElasticProblem_NC_Z<2>;
  template class ElasticProblem_NC_Z<3>;

  template class ElasticProblem<2>;
  template class ElasticProblem<3>;


  // =================================================================
  //
  //                  DataGenerator
  //
  // =================================================================
  template<int dim>
  DataGenerator<dim>::DataGenerator(
    const bool /* distortion */
  ) {
    throw ExcNotImplemented();
  } // DataGenerator<dim>::DataGenerator()

  template<>
  DataGenerator<2>::DataGenerator(
    const bool distortion
  ) {
    std::vector< std::vector< double > > kv(2);
    std::vector< Point<3> > cps(9);
    std::vector< unsigned int > deg;
  
    // Define the vector of control points
    const double w = 1.;
    if (distortion) {
      cps[0] = w * Point<2 + 1>(+0.00, +0.00, 1./w);
      cps[1] = w * Point<2 + 1>(+0.50, +0.00, 1./w);
      cps[2] = w * Point<2 + 1>(+1.00, +0.00, 1./w);
      cps[3] = w * Point<2 + 1>(-0.50, +0.50, 1./w);
      cps[4] = w * Point<2 + 1>(+0.00, +0.50, 1./w);
      cps[5] = w * Point<2 + 1>(+0.50, +0.50, 1./w);
      cps[6] = w * Point<2 + 1>(+0.00, +1.00, 1./w);
      cps[7] = w * Point<2 + 1>(+0.50, +1.00, 1./w);
      cps[8] = w * Point<2 + 1>(+1.00, +1.00, 1./w);
    } else { 
      cps[0] = w * Point<2 + 1>(-1.00, -1.00, 1./w);
      cps[1] = w * Point<2 + 1>(+0.00, -1.00, 1./w);
      cps[2] = w * Point<2 + 1>(+1.00, -1.00, 1./w);
      cps[3] = w * Point<2 + 1>(-1.00, +0.00, 1./w);
      cps[4] = w * Point<2 + 1>(+0.00, +0.00, 1./w);
      cps[5] = w * Point<2 + 1>(+1.00, +0.00, 1./w);
      cps[6] = w * Point<2 + 1>(-1.00, +1.00, 1./w);
      cps[7] = w * Point<2 + 1>(+0.00, +1.00, 1./w);
      cps[8] = w * Point<2 + 1>(+1.00, +1.00, 1./w);
    }
  
    // define the knot vectors
    kv[0] = {0, 0, 0, 1, 1, 1};
    kv[1] = {0, 0, 0, 1, 1, 1};
  
    // Define the degree
    deg = {2, 2};
  
    data = IPF_Data<2>(cps, kv, deg);
  } // DataGenerator<2>::DataGenerator()

  template<>
  DataGenerator<3>::DataGenerator(
    const bool distortion
  ) {
    std::vector< std::vector< double > > kv(3);
    std::vector< Point<3> > cps(27);
    std::vector< unsigned int > deg;
  
    // Define the vector of control points
    if ( distortion ) {
      cps[0 ] = Point<3>(+0.0, +0.0, +0.0);
      cps[1 ] = Point<3>(+0.5, +0.0, +0.0);
      cps[2 ] = Point<3>(+1.0, +0.0, +0.0);
      cps[3 ] = Point<3>(-0.5, +0.5, +0.0);
      cps[4 ] = Point<3>(+0.0, +0.5, +0.0);
      cps[5 ] = Point<3>(+0.5, +0.5, +0.0);
      cps[6 ] = Point<3>(+0.0, +1.0, +0.0);
      cps[7 ] = Point<3>(+0.5, +1.0, +0.0);
      cps[8 ] = Point<3>(+1.0, +1.0, +0.0);

      cps[9 ] = Point<3>(+0.0, +0.0, +0.5);
      cps[10] = Point<3>(+0.5, +0.0, +0.5);
      cps[11] = Point<3>(+1.0, +0.0, +0.5);
      cps[12] = Point<3>(-0.5, +0.5, +0.5);
      cps[13] = Point<3>(+0.0, +0.5, +0.5);
      cps[14] = Point<3>(+0.5, +0.5, +0.5);
      cps[15] = Point<3>(+0.0, +1.0, +0.5);
      cps[16] = Point<3>(+0.5, +1.0, +0.5);
      cps[17] = Point<3>(+1.0, +1.0, +0.5);

      cps[18] = Point<3>(+0.0, +0.0, +1.0);
      cps[19] = Point<3>(+0.5, +0.0, +1.0);
      cps[20] = Point<3>(+1.0, +0.0, +1.0);
      cps[21] = Point<3>(-0.5, +0.5, +1.0);
      cps[22] = Point<3>(+0.0, +0.5, +1.0);
      cps[23] = Point<3>(+0.5, +0.5, +1.0);
      cps[24] = Point<3>(+0.0, +1.0, +1.0);
      cps[25] = Point<3>(+0.5, +1.0, +1.0);
      cps[26] = Point<3>(+1.0, +1.0, +1.0);
    } else { 
      cps[0 ] = Point<3>(-1.0, +0.0, +0.0);
      cps[1 ] = Point<3>(+0.0, +0.0, +0.0);
      cps[2 ] = Point<3>(+1.0, +0.0, +0.0);
      cps[3 ] = Point<3>(-1.0, +0.5, +0.0);
      cps[4 ] = Point<3>(+0.0, +0.5, +0.0);
      cps[5 ] = Point<3>(+1.0, +0.5, +0.0);
      cps[6 ] = Point<3>(-1.0, +1.0, +0.0);
      cps[7 ] = Point<3>(+0.0, +1.0, +0.0);
      cps[8 ] = Point<3>(+1.0, +1.0, +0.0);

      cps[9 ] = Point<3>(-1.0, +0.0, +0.5);
      cps[10] = Point<3>(+0.0, +0.0, +0.5);
      cps[11] = Point<3>(+1.0, +0.0, +0.5);
      cps[12] = Point<3>(-1.0, +0.5, +0.5);
      cps[13] = Point<3>(+0.0, +0.5, +0.5);
      cps[14] = Point<3>(+1.0, +0.5, +0.5);
      cps[15] = Point<3>(-1.0, +1.0, +0.5);
      cps[16] = Point<3>(+0.0, +1.0, +0.5);
      cps[17] = Point<3>(+1.0, +1.0, +0.5);

      cps[18] = Point<3>(-1.0, +0.0, +1.0);
      cps[19] = Point<3>(+0.0, +0.0, +1.0);
      cps[20] = Point<3>(+1.0, +0.0, +1.0);
      cps[21] = Point<3>(-1.0, +0.5, +1.0);
      cps[22] = Point<3>(+0.0, +0.5, +1.0);
      cps[23] = Point<3>(+1.0, +0.5, +1.0);
      cps[24] = Point<3>(-1.0, +1.0, +1.0);
      cps[25] = Point<3>(+0.0, +1.0, +1.0);
      cps[26] = Point<3>(+1.0, +1.0, +1.0);
    }

    // define the knot vectors
    kv[0] = {0, 0, 0, 1, 1, 1};
    kv[1] = {0, 0, 0, 1, 1, 1};
    kv[2] = {0, 0, 0, 1, 1, 1};
  
    // Define the degree
    deg = {2, 2, 2};
  
    std::vector< double >     weights(27, 1.);
    data = IPF_Data<3>(cps, weights, kv, deg);
  } // DataGenerator<3>::DataGenerator()


  // =================================================================
  //
  //                  Right-Hand-Side
  //
  // =================================================================
  template<int dim>
  double ElasticProblem_RHS<dim>::value(
    const Point<dim>      &eval,
    const unsigned int     component
  ) const {
    Assert(dim >= 2, ExcNotImplemented());
    
    const double& tmp1 = sol_fcn->gradient_divergence(eval, component);
    const double& tmp2 = sol_fcn->divergence_gradient(eval, component);
    return -1. * lambda.value(eval) * tmp1 +
            -1. * mu.value(eval) * tmp2 +
            -1. * mu.value(eval) * tmp1;
  } // ElasticProblem_RHS<dim>::value()

  template<int dim> 
  void ElasticProblem_RHS<dim>::vector_value(
    const Point<dim>      &eval,
          Vector<double>  &out
  ) const { 
    Assert(dim >= 2, ExcNotImplemented());
    Assert(out.size() == dim, ExcInvalidState());

    Tensor<1, dim> tensor_out;
    Tensor<1, dim> tmp1, tmp2;
    sol_fcn->gradient_divergence_tensor(eval, tmp1);
    sol_fcn->divergence_gradient_tensor(eval, tmp2);
    tensor_out = -1. * lambda.value(eval) * tmp1 +
                 -1. * mu.value(eval) * tmp2 +
                 -1. * mu.value(eval) * tmp1;
    for (unsigned int d = 0; d < dim; d++)
      out(d) = tensor_out[d];
  } // ElasticProblem_RHS<dim>::vector_value()

  template<int dim>
  void ElasticProblem_RHS<dim>::value_list(
    const std::vector< Point<dim> >   &evals,
          std::vector< double >       &out,
    const unsigned int                 component
  ) const {
    Assert(evals.size() == out.size(), ExcDimensionMismatch(evals.size(), out.size()));
    const unsigned int n_pnts = evals.size();
    std::vector< double > tmp1(n_pnts);
    std::vector< double > tmp2(n_pnts);
    
    sol_fcn->gradient_divergence_list(evals, tmp1, component); 
    sol_fcn->divergence_gradient_list(evals, tmp2, component); 
    for (unsigned int n = 0; n < n_pnts; n++){
      out[n] = -1. * lambda.value(evals[n]) * tmp1[n] + 
                 -1. * mu.value(evals[n]) * tmp2[n] +
                 -1. * mu.value(evals[n]) * tmp1[n];
    }
  } // ElasticProblem_RHS<dim>::value_list()
  
  template<int dim>
  void ElasticProblem_RHS<dim>::vector_value_list(
    const std::vector< Point<dim> >     &evals,
          std::vector< Vector<double> > &out
  ) const { 
    Assert(evals.size() == out.size(), ExcDimensionMismatch(evals.size(), out.size()));
    const unsigned int n_pnts = evals.size(); 
    for (unsigned int n = 0; n < n_pnts; n++)
      vector_value(evals[n], out[n]);
  } // ElasticProblem_RHS<dim>::vector_value_list()

  template<int dim> 
  void ElasticProblem_RHS<dim>::tensor_value(
    const Point<dim>      &eval,
          Tensor<1, dim>  &out
  ) const { 
    Assert(dim >= 2, ExcNotImplemented());

    Tensor<1, dim> tmp1, tmp2;
    sol_fcn->gradient_divergence_tensor(eval, tmp1);
    sol_fcn->divergence_gradient_tensor(eval, tmp2);
    out = -1. * lambda.value(eval) * tmp1 +
          -1. * mu.value(eval) * tmp2 +
          -1. * mu.value(eval) * tmp1;
  } // ElasticProblem_RHS<dim>::tensor_value()

  template<int dim>
  void ElasticProblem_RHS<dim>::tensor_value_list(
    const std::vector< Point<dim> >     &evals,
          std::vector< Tensor<1, dim> > &out
  ) const {
    Assert(evals.size() == out.size(), ExcDimensionMismatch(evals.size(), out.size()));
    const unsigned int n_pnts = evals.size(); 
    for (unsigned int n = 0; n < n_pnts; n++)
      tensor_value(evals[n], out[n]);
  } // ElasticProblem_RHS<dim>::tensor_value_list()

  // =================================================================
  //
  //                  Neumann Condition X0
  //
  // =================================================================
  template<int dim>
  double ElasticProblem_NC_X0<dim>::value(
    const Point<dim>      &eval,
    const unsigned int     component
  ) const {
    Point<dim> mapped_eval;
    Tensor<2, dim> geometry_jacobian; 
    for (unsigned int d = 0; d < dim; d++) {
      mapped_eval(d)       = geometry.value(eval, d);
      geometry_jacobian[d] = geometry.gradient(eval, d);
    }
    DerivativeForm<1, dim, dim> tangents(geometry_jacobian);
    tangents = tangents.transpose();


    Tensor<1, dim> normal; 
    switch( dim ) {
      case 2: 
        normal = cross_product_2d(tangents[1] /* x = 0 */);
        break; 
      case 3: 
        normal = cross_product_3d(tangents[1], 
                                  tangents[2]);
        break;
      default:
        throw ExcNotImplemented();
        break;
    }

    // ensure normal points outwards
    if (normal[0] > 0)
      normal *= -1.;
    
    normal /= normal.norm();

             Tensor<2, dim> solution_gradients;
             Tensor<2, dim> stress;
    SymmetricTensor<2, dim> strain; 
    for (unsigned int d = 0; d < dim; d++)
      solution_gradients[d] = solution.gradient(mapped_eval, d);

    strain = symmetrize(solution_gradients);
    for (unsigned int i = 0; i < dim; i++) {
      for (unsigned int j = 0; j < dim; j++) {
        for (unsigned int k = 0; k < dim; k++) {
          for (unsigned int l = 0; l < dim; l++) {
            if (i == j && k == l)
              stress[i][j] += strain[k][l];
            if ( (i == k && j == l) || 
                 (i == l && j == k))
              stress[i][j] += strain[k][l];
          }
        }
      }
    }


    return normal * stress[component];
  } // ElasticProblem_NC_X0<dim>::value()


  template<int dim>
  void ElasticProblem_NC_X0<dim>::tensor_value(
    const Point<dim>      &eval,
          Tensor<1, dim>  &out
  ) const {
    Point<dim> mapped_eval;
    Tensor<2, dim> geometry_jacobian; 
    for (unsigned int d = 0; d < dim; d++) {
      mapped_eval(d)       = geometry.value(eval, d);
      geometry_jacobian[d] = geometry.gradient(eval, d);
    }
    DerivativeForm<1, dim, dim> tangents(geometry_jacobian);
    tangents = tangents.transpose();


    Tensor<1, dim> normal; 
    switch( dim ) {
      case 2: 
        normal = cross_product_2d(tangents[1] /* x = 0 */);
        break; 
      case 3: 
        normal = cross_product_3d(tangents[1], 
                                  tangents[2]);
        break;
      default:
        throw ExcNotImplemented();
        break;
    }

    // ensure normal points outwards
    if (normal[0] > 0)
      normal *= -1.;
    
    normal /= normal.norm();

             Tensor<2, dim> solution_gradients;
             Tensor<2, dim> stress;
    SymmetricTensor<2, dim> strain; 
    for (unsigned int d = 0; d < dim; d++)
      solution_gradients[d] = solution.gradient(mapped_eval, d);

    strain = symmetrize(solution_gradients);
    for (unsigned int i = 0; i < dim; i++) {
      for (unsigned int j = 0; j < dim; j++) {
        for (unsigned int k = 0; k < dim; k++) {
          for (unsigned int l = 0; l < dim; l++) {
            if (i == j && k == l)
              stress[i][j] += strain[k][l];
            if ( (i == k && j == l) || 
                 (i == l && j == k))
              stress[i][j] += strain[k][l];
          }
        }
      }
    }

    for (unsigned int d = 0; d < dim; d++)
      out[d] = normal * stress[d];
  } // ElasticProblem_NC_X0<dim>::tensor_value()

  // =================================================================
  //
  //                  Neumann Condition Y0
  //
  // =================================================================
  template<int dim>
  double ElasticProblem_NC_Y0<dim>::value(
    const Point<dim>      &eval,
    const unsigned int     component
  ) const {
    Point<dim> mapped_eval;
    Tensor<2, dim> geometry_jacobian; 
    for (unsigned int d = 0; d < dim; d++) {
      mapped_eval(d)       = geometry.value(eval, d);
      geometry_jacobian[d] = geometry.gradient(eval, d);
    }
    DerivativeForm<1, dim, dim> tangents(geometry_jacobian);
    tangents = tangents.transpose();


    Tensor<1, dim> normal; 
    switch( dim ) {
      case 2: 
        normal = cross_product_2d(tangents[0] /* y = 0 */);
        break; 
      case 3: 
        normal = cross_product_3d(tangents[0], 
                                  tangents[2]);
        break;
      default:
        throw ExcNotImplemented();
        break;
    }

    // ensure normal points outwards
    if (normal[1] > 0)
      normal *= -1.;

    normal /= normal.norm();

             Tensor<2, dim> solution_gradients;
             Tensor<2, dim> stress;
    SymmetricTensor<2, dim> strain; 
    for (unsigned int d = 0; d < dim; d++)
      solution_gradients[d] = solution.gradient(mapped_eval, d);

    strain = symmetrize(solution_gradients);
    for (unsigned int i = 0; i < dim; i++) {
      for (unsigned int j = 0; j < dim; j++) {
        for (unsigned int k = 0; k < dim; k++) {
          for (unsigned int l = 0; l < dim; l++) {
            if (i == j && k == l)
              stress[i][j] += strain[k][l];
            if ( (i == k && j == l) || 
                 (i == l && j == k))
              stress[i][j] += strain[k][l];
          }
        }
      }
    }

    return normal * stress[component];
  } // ElasticProblem_NC_Y0<dim>::value()

  template<int dim>
  void ElasticProblem_NC_Y0<dim>::tensor_value(
    const Point<dim>      &eval,
          Tensor<1, dim>  &out
  ) const {
    Point<dim> mapped_eval;
    Tensor<2, dim> geometry_jacobian; 
    for (unsigned int d = 0; d < dim; d++) {
      mapped_eval(d)       = geometry.value(eval, d);
      geometry_jacobian[d] = geometry.gradient(eval, d);
    }
    DerivativeForm<1, dim, dim> tangents(geometry_jacobian);
    tangents = tangents.transpose();


    Tensor<1, dim> normal; 
    switch( dim ) {
      case 2: 
        normal = cross_product_2d(tangents[0] /* y = 0 */);
        break; 
      case 3: 
        normal = cross_product_3d(tangents[0], 
                                  tangents[2]);
        break;
      default:
        throw ExcNotImplemented();
        break;
    }

    // ensure normal points outwards
    if (normal[1] > 0)
      normal *= -1.;

    normal /= normal.norm();

             Tensor<2, dim> solution_gradients;
             Tensor<2, dim> stress;
    SymmetricTensor<2, dim> strain; 
    for (unsigned int d = 0; d < dim; d++)
      solution_gradients[d] = solution.gradient(mapped_eval, d);

    strain = symmetrize(solution_gradients);
    for (unsigned int i = 0; i < dim; i++) {
      for (unsigned int j = 0; j < dim; j++) {
        for (unsigned int k = 0; k < dim; k++) {
          for (unsigned int l = 0; l < dim; l++) {
            if (i == j && k == l)
              stress[i][j] += strain[k][l];
            if ( (i == k && j == l) || 
                 (i == l && j == k))
              stress[i][j] += strain[k][l];
          }
        }
      }
    }

    for (unsigned int d = 0; d < dim; d++)
      out[d] = normal * stress[d];
  } // ElasticProblem_NC_Y0<dim>::tensor_value()


  // =================================================================
  //
  //                  Neumann Condition Z 
  //
  // =================================================================
  template<int dim>
  double ElasticProblem_NC_Z<dim>::value(
      const Point<dim>&,
      const unsigned int
  ) const { 
    throw ExcNotImplemented();
  }

  template<int dim>
  void ElasticProblem_NC_Z<dim>::tensor_value(
      const Point<dim>&,
            Tensor<1, dim>&
  ) const { 
    throw ExcNotImplemented();
  }

  template<>
  double ElasticProblem_NC_Z<3>::value(
    const Point<3>        &eval,
    const unsigned int     component
  ) const {
    Point<3> mapped_eval;
    Tensor<2, 3> geometry_jacobian; 
    for (unsigned int d = 0; d < 3; d++) {
      mapped_eval(d)       = geometry.value(eval, d);
      geometry_jacobian[d] = geometry.gradient(eval, d);
    }
    DerivativeForm<1, 3, 3> tangents(geometry_jacobian);
    tangents = tangents.transpose();


    Tensor<1, 3> normal = cross_product_3d(tangents[0], tangents[1]);

    // ensure normal points outwards
    if ((normal[2] > 0 && eval(2) == 0) || 
          normal[2] < 0 && eval(2) == 1)
      normal *= -1.;

    normal /= normal.norm();

             Tensor<2, 3> solution_gradients;
             Tensor<2, 3> stress;
    SymmetricTensor<2, 3> strain; 
    for (unsigned int d = 0; d < 3; d++)
      solution_gradients[d] = solution.gradient(mapped_eval, d);

    strain = symmetrize(solution_gradients);
    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        for (unsigned int k = 0; k < 3; k++) {
          for (unsigned int l = 0; l < 3; l++) {
            if (i == j && k == l)
              strain[i][j] += stress[k][l];
            if ( (i == k && j == l) || 
                 (i == l && j == k))
              strain[i][j] += stress[k][l];
          }
        }
      }
    }

    return normal * stress[component];
  } // ElasticProblem_NC_Z<3>::value()

  template<>
  void ElasticProblem_NC_Z<3>::tensor_value(
    const Point<3>      &eval,
          Tensor<1, 3>  &out
  ) const {
    Point<3> mapped_eval;
    Tensor<2, 3> geometry_jacobian; 
    for (unsigned int d = 0; d < 3; d++) {
      mapped_eval(d)       = geometry.value(eval, d);
      geometry_jacobian[d] = geometry.gradient(eval, d);
    }
    DerivativeForm<1, 3, 3> tangents(geometry_jacobian);
    tangents = tangents.transpose();


    Tensor<1, 3> normal = cross_product_3d(tangents[0], tangents[1]);

    // ensure normal points outwards
    if ((normal[2] > 0 && eval(2) == 0) || 
          normal[2] < 0 && eval(2) == 1)
      normal *= -1.;

    normal /= normal.norm();

             Tensor<2, 3> solution_gradients;
             Tensor<2, 3> stress;
    SymmetricTensor<2, 3> strain; 
    for (unsigned int d = 0; d < 3; d++)
      solution_gradients[d] = solution.gradient(mapped_eval, d);

    strain = symmetrize(solution_gradients);
    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        for (unsigned int k = 0; k < 3; k++) {
          for (unsigned int l = 0; l < 3; l++) {
            if (i == j && k == l)
              stress[i][j] += strain[k][l];
            if ( (i == k && j == l) || 
                 (i == l && j == k))
              stress[i][j] += strain[k][l];
          }
        }
      }
    }

    for (unsigned int d = 0; d < 3; d++)
      out[d] = stress[d] * normal;
  } // ElasticProblem_NC_Z<3>::tensor_value()


  // =================================================================
  // 
  //                  Solution
  //
  // =================================================================
  template<int dim>
  double ElasticProblem_SOL<dim>::value(
    const Point<dim>      &eval,
    const unsigned int     component
  ) const {
    Assert(dim >= 2, ExcNotImplemented());

    const double &x  = eval(0);
    const double &y  = eval(1); 
    const double &pi = numbers::PI;
    const double &frac = 1. / (pi * pi);
    switch (component) {
      case 0: 
        return +1. * frac * std::sin(pi * x) * std::sin(pi * y);
      case 1: 
        return -1. * frac * std::cos(pi/2. * x) * std::cos(pi/2. * y);
      default:
        return 0.;
    } // switch
  } // ElasticProblem_SOL<dim>::value()

  template<int dim> 
  void ElasticProblem_SOL<dim>::vector_value(
    const Point<dim>      &eval,
          Vector<double>  &out
  ) const { 
    Assert(dim >= 2, ExcNotImplemented());
    Assert(out.size() == dim, ExcInvalidState());
    const double &x  = eval(0);
    const double &y  = eval(1); 
    const double &pi = numbers::PI;
    const double &frac = 1. / (pi * pi);

    out(0) = +1. * std::sin(pi * x) * std::sin(pi * y);
    out(1) = -1. * std::cos(pi/2. * x) * std::cos(pi/2. * y);
    if (dim == 3)
      out(2) = 0.;

    out *= frac;
  } // ElasticProblem_SOL<dim>::vector_value()

  template<int dim>
  void ElasticProblem_SOL<dim>::tensor_value(
    const Point<dim>      &eval,
          Tensor<1, dim>  &out
  ) const {
    const double &x  = eval(0);
    const double &y  = eval(1); 
    const double &pi = numbers::PI;
    const double &frac = 1. / (pi * pi);

    out[0] = +1. * std::sin(pi * x) * std::sin(pi * y);
    out[1] = -1. * std::cos(pi/2. * x) * std::cos(pi/2. * y);
    if (dim == 3)
      out[2] = 0.;

    out *= frac;
  } // ElasticProblem_SOL<dim>::tensor_value()

  template<int dim>
  void ElasticProblem_SOL<dim>::value_list(
    const std::vector< Point<dim> >   &evals,
          std::vector< double >       &out,
    const unsigned int                 component
  ) const {
    Assert(evals.size() == out.size(), ExcDimensionMismatch(evals.size(), out.size()));
    const unsigned int n_pnts = evals.size(); 
    for (unsigned int n = 0; n < n_pnts; n++)
      out[n] = value(evals[n], component);
  } // ElasticProblem_SOL<dim>::value_list()
  
  template<int dim>
  void ElasticProblem_SOL<dim>::vector_value_list(
    const std::vector< Point<dim> >     &evals,
          std::vector< Vector<double> > &out
  ) const { 
    Assert(evals.size() == out.size(), ExcDimensionMismatch(evals.size(), out.size()));
    const unsigned int n_pnts = evals.size(); 
    for (unsigned int n = 0; n < n_pnts; n++)
      vector_value(evals[n], out[n]);
  } // ElasticProblem_SOL<dim>::vector_value_list()

  template<int dim>
  Tensor<1, dim> ElasticProblem_SOL<dim>::gradient(
    const Point<dim>    &eval,
    const unsigned int   component
  ) const {
    Assert(dim >= 2, ExcNotImplemented());

    const double         &x  = eval(0);
    const double         &y  = eval(1); 
    const double         &pi = numbers::PI;
    const double         &frac = 1. / pi;
          Tensor<1, dim>  out;
    switch (component) {
      case 0: 
        out[0] = +1. * frac * std::cos(pi * x) * std::sin(pi * y);
        out[1] = +1. * frac * std::sin(pi * x) * std::cos(pi * y);
      break;
      case 1:
        out[0] = +1./2. * frac * std::sin(pi/2. * x) * std::cos(pi/2. * y);
        out[1] = +1./2. * frac * std::cos(pi/2. * x) * std::sin(pi/2. * y);
      break;
      default: 
      break;
    } // switch

    return out;
  } // ElasticProblem_SOL<dim>::gradient()

  template<int dim>
  void ElasticProblem_SOL<dim>::vector_gradient(
    const Point<dim>                    &eval,
          std::vector< Tensor<1, dim> > &out
  ) const {
    Assert(out.size() == dim, ExcDimensionMismatch(out.size(), dim));
    Assert(dim >= 2, ExcNotImplemented());

    const double         &x  = eval(0);
    const double         &y  = eval(1); 
    const double         &pi = numbers::PI;
    out[0][0] = +1. * pi * std::cos(pi * x) * std::sin(pi * y);
    out[0][1] = +1. * pi * std::sin(pi * x) * std::cos(pi * y);
    out[1][0] = +1. * pi/2. * std::sin(pi/2. * x) * std::cos(pi/2. * y);
    out[1][1] = +1. * pi/2. * std::cos(pi/2. * x) * std::sin(pi/2. * y);
    if (dim == 3) {
      out[0][2] = 0.;
      out[1][2] = 0.;
      out[2] = Tensor<1, dim>();
    }
  } // ElasticProblem_SOL<dim>::vector_gradient()

  template<int dim>
  void ElasticProblem_SOL<dim>::gradient_list(
    const std::vector< Point<dim> >         &evals,
          std::vector< Tensor<1, dim> >     &out,
    const unsigned int                       component 
  ) const {
    Assert(evals.size() == out.size(), ExcDimensionMismatch(evals.size(), out.size()));
    const unsigned int n_pnts = evals.size();
    for (unsigned int n = 0; n < n_pnts; n++)
      out[n] = gradient(evals[n], component);
  } // ElasticProblem_SOL<dim>::gradient_list()

  template<int dim>
  void ElasticProblem_SOL<dim>::vector_gradients(
    const std::vector< Point<dim> >                     &evals,
          std::vector< std::vector< Tensor<1, dim> > >  &out
  ) const {
    Assert(evals.size() == out.size(), ExcDimensionMismatch(evals.size(), out.size()));
    const unsigned int n_pnts = evals.size();
    for (unsigned int n = 0; n < n_pnts; n++)
      vector_gradient(evals[n], out[n]); 
  } // ElasticProblem_SOL<dim>::vector_gradients()

  template<int dim>
  double ElasticProblem_SOL<dim>::divergence(
    const Point<dim>    &eval
  ) const {
    return gradient(eval, 0)[0] + gradient(eval, 1)[1]; 
  } // ElasticProblem_SOL<dim>::divergence()

  template<int dim>
  void ElasticProblem_SOL<dim>::divergence_list(
    const std::vector< Point<dim> >     &evals,
          std::vector< double >         &out
  ) const {
    Assert(evals.size() == out.size(), ExcDimensionMismatch(evals.size(), out.size()));
    const unsigned int n_pnts = evals.size();
    for (unsigned int n = 0; n < n_pnts; n++)
      out[n] = divergence(evals[n]);
  } // ElasticProblem_SOL<dim>::divergence_list()

  template<int dim>
  double ElasticProblem_SOL<dim>::gradient_divergence(
    const Point<dim>     &eval,
    const unsigned int    component
  ) const {
    Assert(dim <= 3, ExcNotImplemented());
    const double &x  = eval(0);
    const double &y  = eval(1); 
    const double &pi = numbers::PI;
          double  out;
    switch(component) {
      case 0:
        out =   -1. * std::sin(pi * x) * std::sin(pi * y) 
              + -1./4. * std::sin(pi * x / 2.) * std::sin(pi * y / 2.); 
      break; 
      case 1:
        out =    1. * std::cos(pi * x) * std::cos(pi * y) 
              +  1./4. * std::cos(pi * x / 2.) * std::cos(pi * y / 2.); 
      break;
      default:
        out = 0;
      break;
    } // switch

    return out;
  } // ElasticProblem_SOL<dim>::gradient_divergence()

  template<int dim>
  void ElasticProblem_SOL<dim>::gradient_divergence_list(
    const std::vector< Point<dim> >       &evals,
          std::vector< double >           &out,
    const unsigned int                     component
  ) const {
    Assert(out.size() == evals.size(), ExcDimensionMismatch(out.size(), evals.size()));
    const unsigned int n_pnts = out.size();
    for (unsigned int n = 0; n < n_pnts; n++)
      out[n] = gradient_divergence(evals[n], component);
  } // ElasticProblem_SOL<dim>::gradient_divergence_list()

  template<int dim>
  void ElasticProblem_SOL<dim>::gradient_divergence_tensor(
    const Point< dim >          &eval,
          Tensor<1, dim>        &out
  ) const {
    Assert(dim <= 3, ExcNotImplemented());
    const double &x  = eval(0);
    const double &y  = eval(1); 
    const double &pi = numbers::PI;
    out[0] =     -1. * std::sin(pi * x) * std::sin(pi * y) 
                 -1./4. * std::sin(pi * x / 2.) * std::sin(pi * y / 2.); 
    out[1] =     +1. * std::cos(pi * x) * std::cos(pi * y) 
                 +1./4. * std::cos(pi * x / 2.) * std::cos(pi * y / 2.); 
    if (dim == 3)
      out[2] = 0.;
  } // ElasticProblem_SOL<dim>::gradient_divergence_tensor()

  template<int dim>
  void ElasticProblem_SOL<dim>::gradient_divergence_tensor_list(
    const std::vector< Point<dim> >         &evals,
          std::vector< Tensor<1, dim> >     &out
  ) const {
    Assert(out.size() == evals.size(), ExcDimensionMismatch(out.size(), evals.size()));
    const unsigned int n_pnts = out.size();
    for (unsigned int n = 0; n < n_pnts; n++)
      gradient_divergence_tensor(evals[n], out[n]);
  } // ElasticProblem_SOL<dim>::gradient_divergence_tensor_list()

  template<int dim>
  double ElasticProblem_SOL<dim>::divergence_gradient(
    const Point<dim>     &eval,
    const unsigned int    component
  ) const {
    Assert(dim <= 3, ExcNotImplemented());
    const double &x  = eval(0);
    const double &y  = eval(1); 
    const double &pi = numbers::PI;
          double  out;
    switch(component) {
      case 0:
        out = -2. * std::sin(pi * x) * std::sin(pi * y);
      break; 
      case 1:
        out = +1. / 2. * std::cos(pi / 2. * x) * std::cos(pi / 2. * y);
      break;
      default:
        out = 0;
      break;
    } // switch

    return out;
  } // ElasticProblem_SOL<dim>::divergence_gradient()

  template<int dim>
  void ElasticProblem_SOL<dim>::divergence_gradient_list(
    const std::vector< Point<dim> >       &evals,
          std::vector< double >           &out,
    const unsigned int                     component
  ) const {
    Assert(out.size() == evals.size(), ExcDimensionMismatch(out.size(), evals.size()));
    const unsigned int n_pnts = out.size();
    for (unsigned int n = 0; n < n_pnts; n++)
      out[n] = divergence_gradient(evals[n], component);
  } // ElasticProblem_SOL<dim>::divergence_gradient_list()

  template<int dim>
  void ElasticProblem_SOL<dim>::divergence_gradient_tensor(
    const Point< dim >          &eval,
          Tensor<1, dim>        &out
  ) const {
    Assert(dim <= 3, ExcNotImplemented());
    const double &pi = numbers::PI;
    const double &pi_sq = pi * pi;
    out[0] = -2. * pi_sq * value(eval, 0);
    out[1] = -2. * pi_sq * value(eval, 1);
    if (dim == 3)
      out[2] = 0.;
  } // ElasticProblem_SOL<dim>::divergence_gradient_tensor()

  template<int dim>
  void ElasticProblem_SOL<dim>::divergence_gradient_tensor_list(
    const std::vector< Point<dim> >         &evals,
          std::vector< Tensor<1, dim> >     &out
  ) const {
    Assert(out.size() == evals.size(), ExcDimensionMismatch(out.size(), evals.size()));
    const unsigned int n_pnts = out.size();
    for (unsigned int n = 0; n < n_pnts; n++)
      gradient_divergence_tensor(evals[n], out[n]);
  } // ElasticProblem_SOL<dim>::divergence_gradient_tensor_list()

  template<int dim>
  Tensor<1, dim> ElasticProblem_SOL<dim>::stress(
    const Point<dim>      &eval,
    const unsigned int    &component
  ) const {
             Tensor<2, dim> solution_gradients;
             Tensor<2, dim> stress;
    SymmetricTensor<2, dim> strain; 
    for (unsigned int d = 0; d < dim; d++)
      solution_gradients[d] = gradient(eval, d);

    strain = symmetrize(solution_gradients);
    for (unsigned int i = 0; i < dim; i++) {
      for (unsigned int j = 0; j < dim; j++) {
        for (unsigned int k = 0; k < dim; k++) {
          for (unsigned int l = 0; l < dim; l++) {
            if (i == j && k == l)
              stress[i][j] += strain[k][l];
            if ( (i == k && j == l) || 
                 (i == l && j == k))
              stress[i][j] += strain[k][l];
          }
        }
      }
    }

    return stress[component];
  } // ElasticProblem_SOL<dim>::stress()

  template<int dim>
  Tensor<2, dim> ElasticProblem_SOL<dim>::stress(
    const Point<dim>      &eval
  ) const {
             Tensor<2, dim> solution_gradients;
             Tensor<2, dim> stress;
    SymmetricTensor<2, dim> strain; 
    for (unsigned int d = 0; d < dim; d++)
      solution_gradients[d] = gradient(eval, d);

    strain = symmetrize(solution_gradients);
    for (unsigned int i = 0; i < dim; i++) {
      for (unsigned int j = 0; j < dim; j++) {
        for (unsigned int k = 0; k < dim; k++) {
          for (unsigned int l = 0; l < dim; l++) {
            if (i == j && k == l)
              stress[i][j] += strain[k][l];
            if ( (i == k && j == l) || 
                 (i == l && j == k))
              stress[i][j] += strain[k][l];
          }
        }
      }
    }

    return stress;
  } // ElasticProblem_SOL<dim>::stress()


  // ==============================================================
  //
  //                   Problem
  //
  // ==============================================================
  template<int dim> 
  ElasticProblem<dim>::ElasticProblem(
    const unsigned int ref,
    const unsigned int order_elevate,
    const Strategy     strategy
  ) : data(DataGenerator<dim>().data)
    , ref(ref)
    , order(order_elevate + data.max_degree())
    , tria(data)
    , strategy(strategy)
    , sol_fcn()
    , rhs_fcn(&sol_fcn)
    , neumann_y0(tria.get_IPF())
    , neumann_x0(tria.get_IPF())
    , neumann_z(tria.get_IPF())
    , lambda(1.)
    , mu(1.)
  { 
    Assert(dim == 2 || dim == 3, 
            ExcDimensionMismatch2(dim, 2, 3));
    const std::string name = "linear_elasticity/" + std::to_string(dim) + "d/";
    problem_out = OutputSetup(name, order);

    // Set boundary indicators: 
    // Homogeneous dirichleth boundary at x = 0 or x = 1,
    // In-Homogeneous neumann boundary at y = 0 or y = 1,
    // No boundary conditions for z
    for (const auto& face : tria.active_face_iterators()){
      face -> set_boundary_id(Boundaries::Dirichleth);
      // if ((face -> center()).operator()(0) == 1 || 
      //       (face -> center()).operator()(1) == 1)
      //   face -> set_boundary_id(Boundaries::Dirichleth);
      // else if ((face -> center()).operator()(1) == 0)
      //   face -> set_boundary_id(Boundaries::Neumann_Y0);
      // else if ((face -> center()).operator()(0) == 0)
      //   face -> set_boundary_id(Boundaries::Neumann_X0);
      // else 
      //   face -> set_boundary_id(Boundaries::None);
    } //  for ( face )

    tria.degree_elevate_global(order_elevate);

    tria.prepare_assembly();
  } // ElasticProblem<dim>::ElasticProblem()

  template<int dim>
  void ElasticProblem<dim>::setup_system(
  ) {
    const unsigned int n_dofs = dim * tria.n_active_splines();
    system_rhs.reinit(n_dofs);
    solution.reinit(n_dofs);

    const auto& IEN_array = tria.get_IEN_array(dim);
    Assert(IEN_array.size() > 0, ExcInternalError());

    sparsity_pattern.reinit(
        n_dofs,
        n_dofs,
        n_dofs
    );

    for (const auto& [_, arr] : IEN_array)
      for (unsigned int i : arr)
        for (unsigned int j : arr)
          sparsity_pattern.add(i, j);

    sparsity_pattern.compress();
    system_matrix.reinit(sparsity_pattern);
  } // ElasticProblem<dim>::setup_system()

  template<int dim>
  void ElasticProblem<dim>::assemble_system(
  ) {
    std::vector< unsigned int > degrees = tria.get_degree();
    for (unsigned int d = 0; d < dim; d++)
      degrees[d] += 1;

    TSValues<dim, dim, dim> ts_values(
      &tria, 
      degrees, 
      update_values |
      update_gradients |
      update_quadrature_points |
      update_JxW_values
    );
    TSFaceValues<dim, dim, dim> face_values(
      &tria, 
      degrees, 
      update_values |
      update_gradients |
      update_quadrature_points |
      update_JxW_values | 
      update_normal_vectors
    );

    const unsigned int n_dofs_per_cell = ts_values.n_dofs_per_cell();
    FullMatrix<double> cell_matrix(n_dofs_per_cell, n_dofs_per_cell);
    Vector<double>     cell_rhs(n_dofs_per_cell);

    for (const auto& cell : tria.active_cell_iterators()){
      ts_values.reinit(cell);

      cell_matrix       = 0;
      cell_rhs          = 0;

      // Compute system_matrix part of cell
      for (const unsigned int i : ts_values.dof_indices()){
        const unsigned int comp_i =
                ts_values.system_to_component_index(i).first;
        for (const unsigned int j : ts_values.dof_indices()){
          const unsigned int comp_j =
                  ts_values.system_to_component_index(j).first;
          for (const unsigned int q : ts_values.quadrature_point_indices()){
            cell_matrix(i, j) += (
              ts_values.shape_grad(i, q)[comp_i] *
              ts_values.shape_grad(j, q)[comp_j] * 
              lambda.value(ts_values.quadrature_point(q))
              +
              ts_values.shape_grad(i, q)[comp_j] *
              ts_values.shape_grad(j, q)[comp_i] * 
              mu.value(ts_values.quadrature_point(q))
              +
              (
                (comp_i == comp_j) ? 
                  ts_values.shape_grad(i, q) * 
                  ts_values.shape_grad(j, q) *
                  mu.value(ts_values.quadrature_point(q)) : 
                  0.
              )
            ) * ts_values.JxW(q);
          } // for ( q )
        } // for ( j )
      } // for ( i ) 

      // Compute rhs for this cell
      for (const unsigned int i : ts_values.dof_indices()){
        const unsigned int comp_i =
                ts_values.system_to_component_index(i).first;
        for (const unsigned int q : ts_values.quadrature_point_indices())
          cell_rhs(i) += ts_values.shape_value(i, q) * 
                          rhs_fcn.value(ts_values.quadrature_point(q), comp_i) *
                          ts_values.JxW(q);
      } // for ( i )

      // Check neumann boundary: 
      if (cell -> at_boundary()) {
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; f++) {
          if (!(cell -> face(f) -> at_boundary())
                || cell -> face(f) -> boundary_id() == Boundaries::Dirichleth
                || cell -> face(f) -> boundary_id() == Boundaries::None)
            continue;

          face_values.reinit(cell, f);
          for (const unsigned int q : face_values.quadrature_point_indices()) {
            const Tensor<2, dim>& stress = sol_fcn.stress(
                                            face_values.quadrature_point(q)
                                          );
            for (const unsigned int i : face_values.dof_indices()){
              const unsigned int comp_i =
                face_values.system_to_component_index(i).first;
              cell_rhs(i) += stress[comp_i] *
                              face_values.normal_vector(q) *
                              face_values.shape_value(i, q) *
                              face_values.JxW(q);
            }
          } // for ( i ) 
        } // for ( f )
      } // if ( cell -> at_boundary() )

      // Add values to system
      const std::vector< unsigned int > local_dof_indices =
              tria.get_IEN_array(cell, dim);
      system_matrix.add(local_dof_indices, cell_matrix);
      system_rhs.add(local_dof_indices, cell_rhs);
    } // for ( cell )

    std::map<
        types::global_dof_index, 
        double
      > boundary_values;
    std::map<
        types::boundary_id,
        const Function< dim >*
      > boundary_fcns =
          {{Boundaries::Dirichleth, &sol_fcn}};
    tria.template project_boundary_values<dim>(
      boundary_fcns,
      degrees,
      boundary_values
    );

    MatrixTools::apply_boundary_values(
      boundary_values, 
      system_matrix,
      solution,
      system_rhs
    );

  } // ElasticProblem<dim>::assemble_system()

  template<int dim>
  void ElasticProblem<dim>::solve(
  ) {
    std::cout << "Solving system ... " << std::endl;

    SolverControl solver_control(
                     750 * dim * tria.n_active_splines(), 
                     H1 * 1e-4
                  );
    solver_control.enable_history_data();

    SolverCG<Vector<double>> solver(solver_control);
    // dealii::TrilinosWrappers::SolverCG solver(solver_control);

    PreconditionJacobi<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix);
    
    // TrilinosWrappers::PreconditionIC preconditioner; 
    // preconditioner.initialize(system_matrix);

    solver.solve(system_matrix, solution, system_rhs, preconditioner);

    problem_out.add_values_to_table(
      tria.n_levels() - 1,
      tria.n_active_cells(),
      tria.n_active_splines(),
      solver_control.last_step(),
      *(solver_control.get_history_data().end() - 1)
    );
  
    std::cout << " ... done!" << std::endl;
  } // ElasticProblem<dim>::solve()

  template<int dim>
  void ElasticProblem<dim>::compute_h1_error(
  ) {
    std::vector< unsigned int > degrees = tria.get_degree();
    for (unsigned int& p : degrees)
      p = (p * p + 1);

    TSValues<dim, dim, dim> ts_values(
        &tria,
        degrees,
        update_values |
        update_gradients |
        update_quadrature_points |
        update_hessians |
        update_JxW_values);
  
  
    double h1 = 0;
    double l2 = 0;

    for (const auto& cell : tria.active_cell_iterators()) {
      ts_values.reinit(cell); 
      
      std::vector< unsigned int > local_dof_indices =
        tria.get_IEN_array(cell, dim);

      for (const unsigned int q : ts_values.quadrature_point_indices()) {
        const Point<dim>& mapped_q = ts_values.quadrature_point(q);

        // Get the value differences at quadrature point:
        Tensor<1, dim> u_diff; 
        sol_fcn.tensor_value(mapped_q, u_diff);

        Tensor<2, dim> grad_u_diff;
        for (unsigned int d = 0; d < dim; d++)
          grad_u_diff[d] = sol_fcn.gradient(mapped_q, d);

        // Get the differences for the jacobians
        for (unsigned int i : ts_values.dof_indices()){
          const unsigned int comp_i = ts_values.system_to_component_index(i).first;
          const double       c_i    = solution(local_dof_indices[i]);

          u_diff[comp_i]       -= c_i * ts_values.shape_value(i, q);
          grad_u_diff[comp_i]  -= c_i * ts_values.shape_grad(i, q);          
        } // for( i )

        const double dx = ts_values.JxW(q); 
        for (unsigned int i = 0; i < dim; i++)
          for (unsigned int j = 0; j < dim; j++)
            h1 += grad_u_diff[i][j] * grad_u_diff[i][j] * dx;

        l2 += (      u_diff * u_diff      ) * dx;
      } // for ( q )
    } // for ( cell )

    H1 = std::sqrt(l2 + h1);
    L2 = std::sqrt(l2);
    problem_out.add_values_to_table(L2, H1);
  } // ElasticProblem<dim>::compute_h1_error()

  template<int dim>
  void ElasticProblem<dim>::estimate_mark_refine(
  ) { 
    if (strategy == Strategy::Adaptive) {
      Vector<double> residuals(tria.n_active_cells()) ;
      // std::map< types::boundary_id, Function<dim>* > neumann_data;
      // if ( dim == 2)
      //   neumann_data = {{Boundaries::Neumann_Y0, &neumann_y0}, 
      //                   {Boundaries::Neumann_X0, &neumann_x0}};
      // else {
      //   neumann_data = {{Boundaries::Neumann_Y0, &neumann_y0}, 
      //                   {Boundaries::Neumann_X0, &neumann_x0},
      //                   {Boundaries::Neumann_Z,  &neumann_z}};
      // }
      std::vector< unsigned int > degrees = tria.get_degree();
      for (auto& p : degrees)
        p = p * p + 1;
      tria.linear_elasticity_residual_error_estimate(
          degrees,
          &rhs_fcn, 
          &lambda,
          &mu,
          // neumann_data,
          solution,
          residuals
      );
      tria.refine_fixed_number(residuals, 0.05);
    } else { 
      tria.coarsen_bezier_elements();
      tria.refine_global();
    }

    tria.prepare_assembly();
  } // ElasticProblem<dim>::estimate_mark_refine()

  template<int dim>
  void ElasticProblem<dim>::output_system(
  ) {
    throw ExcNotImplemented();
  } // ElasticProblem<dim>::output_system()

  template<>
  void ElasticProblem<2>::output_system(
  ) {
    std::cout << "Printing system matrix and rhs ... " << std::endl;
  
    const unsigned int level = tria.n_levels() - 1;
    const std::string name = problem_out.degree.string();
    const std::string level_name = name + "l" + std::to_string(level);
  
    std::string matrix = level_name + "_mat.dat" ;
    std::string vector = level_name + "_vec.dat" ;
    std::string soluti = level_name + "_sol.dat" ;
  
    if (tria.n_levels() - 1 < 11) {
      GridOutFlags::Svg svg_flags;
      // svg_flags.label_level_number = true;
      // svg_flags.label_cell_index = true;
      // svg_flags.label_boundary_id = true;
      svg_flags.boundary_line_thickness = 1;
      svg_flags.line_thickness          = 1;
      svg_flags.coloring = GridOutFlags::Svg::Coloring::none;
  
      const std::string name_svg = problem_out.svg.string() + "grid_l" + std::to_string(level);
      
      std::ofstream out(name_svg + ".svg");
      GridOut       grid_out;
      grid_out.set_flags(svg_flags);
      
      grid_out.write_svg(tria, out);

      std::filebuf mat, vec, sol;
      mat.open(matrix.c_str(), std::ios::out);
      vec.open(vector.c_str(), std::ios::out);
      sol.open(soluti.c_str(), std::ios::out);
  
      std::ostream mat_out(&mat);
      std::ostream vec_out(&vec);
      std::ostream sol_out(&sol);
  
      system_matrix.print_formatted(mat_out, 16, true, 1, "0");
      system_rhs.print(vec_out, 16);
      solution.print(sol_out, 16);
  
      mat.close();
      vec.close();
      sol.close();

      const unsigned int N1 = 100; 
      const unsigned int N2 = 100;
      const double xmin = 0; 
      const double ymin = 0; 
      const double xmax = 1; 
      const double ymax = 1; 
      std::vector<Point<2>> evals; 
      FullMatrix< double > E(N1 * N2, 2);
      unsigned int ind = 0;
      for (double j = 0.; j < N2; j++) {
        for (double i = 0.; i < N1; i++) {
          const double x = xmin + (i/(N1-1.)) * (xmax - xmin); 
          const double y = ymin + (j/(N2-1.)) * (ymax - ymin); 
          evals.push_back(Point<2>(x, y));
          E(ind, 0) = x;
          E(ind, 1) = y;
          ind++;
        }
      }

      std::filebuf e_f;
      e_f.open(level_name + "_evals.dat", std::ios::out);
      std::ostream out_e(&e_f);
      E.print_formatted(out_e, 16, true, 1, "0");


      // Print the IPF wireframe:
      // tria.print_grid(name);
      tria.generate_mesh_file<0>(level_name, false, 16);
      tria.generate_mesh_file<0>(level_name, true, 16);
      tria.print_IPF_wireframe(level_name);
      tria.printIPF(2, evals, level_name, 16, true, true);
      tria.coarsen_bezier_elements();
      tria.generate_mesh_file<0>(level_name, false, 16);
      tria.generate_mesh_file<0>(level_name, true, 16);
      tria.print_IPF_wireframe(level_name);
      tria.refine_bezier_elements();
    } 

    // Write the grid to a seperate file: 
    const std::string& name_vtg = problem_out.vtg.string() + "physical_grid_l" + std::to_string(level) + ".vtu";

    // First: Make a copy of the triangulation: 
    Triangulation<2> physical_grid; 
    physical_grid.copy_triangulation(tria);

    // And transform it with the IPF from tria.
    // Note: This will make a linear representation
    // of the boundary and interior nodes.
    const IsoparametricManifold<2> geometry(tria.get_IPF()); 
    GridTools::transform(
      [&geometry](const Point<2>& p) -> Point<2>{
        return geometry.push_forward(p);
      },
      physical_grid 
    );

    // Generate the output object
    DataOut<2> data_out;
    data_out.attach_triangulation(physical_grid); 

    // Estimate the error on each cell
    Vector< double > residuals(tria.n_active_cells());
      // std::map< types::boundary_id, Function<dim>* > neumann_data;
      // if ( dim == 2)
      //   neumann_data = {{Boundaries::Neumann_Y0, &neumann_y0}, 
      //                   {Boundaries::Neumann_X0, &neumann_x0}};
      // else {
      //   neumann_data = {{Boundaries::Neumann_Y0, &neumann_y0}, 
      //                   {Boundaries::Neumann_X0, &neumann_x0},
      //                   {Boundaries::Neumann_Z,  &neumann_z}};
      // }
    std::vector<unsigned int> degrees = tria.get_degree();
    for (unsigned int p : degrees)
      p = p * p + 1;

    tria.linear_elasticity_residual_error_estimate(
        degrees,
        &rhs_fcn, 
        &lambda,
        &mu,
        // neumann_data,
        solution,
        residuals
    );

    data_out.add_data_vector(residuals, "cell_errors");

    // Build patches
    data_out.build_patches(); 

    // Open a file
    std::ofstream vtu_out(name_vtg); 
    data_out.write_vtu(vtu_out);

    // Print the table to a file preemptively
    // output to .txt
    problem_out.write_table_text();
    problem_out.write_table_tex();

    // print the difference pointwise for each cell
    std::cout << " ... output to files done!" << std::endl;
  } // ElasticProblem<dim>::output_system();

  template<int dim>
  void ElasticProblem<dim>::right_hand_side(
    const std::vector< Point<dim> >       &points,
          std::vector< Tensor<1, dim> >   &values
  ) const { 
    Point<dim> point_1, point_2; 
    point_1[0] = +0.5; 
    point_2[0] = -0.5;

    for (unsigned int p = 0; p < points.size(); ++p){
      if ((points[p] - point_1).norm_square() < 0.04 || 
            (points[p] - point_2).norm_square() < 0.04)
        values[p][0] = 1.;
      else 
        values[p][0] = 0.;

      if (points[p].norm_square() < 0.04)
        values[p][1] = 1.;
      else 
        values[p][1] = 0.;
    } // for ( p )
  } // ElasticProblem<dim>::right_hand_side()

  template<int dim>
  void ElasticProblem<dim>::right_hand_side(
    const Point<dim>       &point,
          Tensor<1, dim>   &value
  ) const { 
    Point<dim> point_1, point_2; 
    point_1[0] = +0.5; 
    point_2[0] = -0.5;

    if ((point - point_1).norm_square() < 0.04 || 
          (point - point_2).norm_square() < 0.04)
      value[0] = 1.;
    else 
      value[0] = 0.;

    if (point.norm_square() < 0.04)
      value[1] = 1.;
    else 
      value[1] = 0.;
  } // ElasticProblem<dim>::right_hand_side()

  template<int dim>
  void ElasticProblem<dim>::run(
  ) {
    
    unsigned int cycle = 0;  
    unsigned int old_level = tria.n_levels() - 1; 
    while (tria.n_levels() - 1 < ref) {
      this -> setup_system();
      this -> assemble_system();
      this -> solve();
      this -> compute_h1_error();
      this -> output_system();
      if (cycle < 5) {
        this -> estimate_mark_refine();
      } else { 
        std::cout << "Too many refinements resulted in the same level, enforcing global refinement..." << std::endl;
        tria.coarsen_bezier_elements();
        tria.refine_global(); 
        tria.prepare_assembly();
      }

      if (tria.n_levels() - 1 == old_level) {
        cycle++;
      } else {
        cycle = 0;
        old_level = tria.n_levels() - 1;
      }
    } // while

    // Write the resulting table to line
    problem_out.write_table_text(std::cout);
  } // ElasticProblem<dim>::run()

} // namespace 
