#include <linear_elasticity_inhomogeneous.h>

namespace Linear_Elasticity {
  template class InhomogeneousProblem<2>;
  template class InhomogeneousProblem<3>;

  // ==============================================================
  //
  //                   Solution
  //
  // ==============================================================
  template<int dim>
  double InhomogeneousProblem_SOL<dim>::value(
    const Point<dim>      &eval,
    const unsigned int     component
  ) const {
    Assert(false, ExcNotImplemented());
    return 0;
  } // InhomogeneousProblem_SOL<dim>::value()

  template<>
  double InhomogeneousProblem_SOL<2>::value(
    const Point<2>      &eval,
    const unsigned int     component
  ) const {

    const double &x  = eval(0);
    const double &y  = eval(1); 
    const double &pi = numbers::PI;
    const double &frac = 1. / (pi * pi);
    switch (component) {
      case 0: 
        return frac * (x+0.5) * std::sin(pi * y);
      default:
        return 0.;
    } // switch
  } // InhomogeneousProblem_SOL<2>::value()

  template<>
  double InhomogeneousProblem_SOL<3>::value(
    const Point<3>      &eval,
    const unsigned int     component
  ) const {

    const double &x  = eval(0);
    const double &y  = eval(1); 
    const double &z  = eval(2); 
    const double sqrt_z = std::sqrt( z * ( 1 - z) );
    const double z_term = sqrt_z * (z * (1 - z));
    const double &pi = numbers::PI;
    const double &frac = 1. / (pi * pi);
    switch (component) {
      case 0: 
        return frac * (x+0.5) * z_term * std::sin(pi * y);
      default:
        return 0.;
    } // switch
  } // InhomogeneousProblem_SOL<3>::value()

  template<int dim> 
  void InhomogeneousProblem_SOL<dim>::vector_value(
    const Point<dim>      &eval,
          Vector<double>  &out
  ) const { 
    Assert(dim >= 2, ExcNotImplemented());
    Assert(out.size() == dim, ExcInvalidState());
    out(0) = this -> value(eval, 0);
    out(1) = 0;
    if (dim == 3)
      out(2) = 0;
  } // InhomogeneousProblem_SOL<dim>::vector_value()

  template<int dim>
  void InhomogeneousProblem_SOL<dim>::tensor_value(
    const Point<dim>      &eval,
          Tensor<1, dim>  &out
  ) const {
    out[0] = this -> value(eval, 0);
    out[1] = 0.;
    if (dim == 3)
      out[2] = 0.;

  } // InhomogeneousProblem_SOL<dim>::tensor_value()

  template<int dim>
  void InhomogeneousProblem_SOL<dim>::value_list(
    const std::vector< Point<dim> >   &evals,
          std::vector< double >       &out,
    const unsigned int                 component
  ) const {
    Assert(evals.size() == out.size(), ExcDimensionMismatch(evals.size(), out.size()));
    const unsigned int n_pnts = evals.size(); 
    for (unsigned int n = 0; n < n_pnts; n++)
      out[n] = value(evals[n], component);
  }
  
  template<int dim>
  void InhomogeneousProblem_SOL<dim>::vector_value_list(
    const std::vector< Point<dim> >     &evals,
          std::vector< Vector<double> > &out
  ) const { 
    Assert(evals.size() == out.size(), ExcDimensionMismatch(evals.size(), out.size()));
    const unsigned int n_pnts = evals.size(); 
    for (unsigned int n = 0; n < n_pnts; n++)
      vector_value(evals[n], out[n]);
  } // InhomogeneousProblem_SOL<dim>::vector_value_list()

  template<int dim>
  Tensor<1, dim> InhomogeneousProblem_SOL<dim>::gradient(
    const Point<dim>    &eval,
    const unsigned int   component
  ) const {
    Assert(false, ExcNotImplemented());
  }


  template<>
  Tensor<1, 2> InhomogeneousProblem_SOL<2>::gradient(
    const Point<2>    &eval,
    const unsigned int   component
  ) const {

    const double         &x  = eval(0);
    const double         &y  = eval(1); 
    const double         &pi = numbers::PI;
    const double         &frac = 1. / (pi * pi);
          Tensor<1, 2>  out;
    switch (component) {
      case 0: 
        out[0] = frac * std::sin(pi * y);
        out[1] = (1. / pi) * (x+0.5) * std::cos(pi * y);
      break;
      default: 
      break;
    } // switch

    return out;
  } // InhomogeneousProblem_SOL<2>::gradient()

  template<>
  Tensor<1, 3> InhomogeneousProblem_SOL<3>::gradient(
    const Point<3>    &eval,
    const unsigned int   component
  ) const {
    const double         &x  = eval(0);
    const double         &y  = eval(1); 
    const double         &z  = eval(2); 
    const double      sqrtz  = std::sqrt( z * (1-z) );
    const double      z_term = sqrtz * (z * (1 - z));
    const double         &pi = numbers::PI;
    const double         &frac = 1. / (pi * pi);
          Tensor<1, 3>  out;
    switch (component) {
      case 0: 
        out[0] = frac * std::sin(pi * y) * z_term;
        out[1] = (1. / pi) * (x+0.5) * std::cos(pi * y) * z_term;
        out[2] = frac * 1.5 * (x+0.5) * (1. - 2.*z) * sqrtz * std::sin(pi*y); 
      break;
      default: 
      break;
    } // switch

    return out;
  } // InhomogeneousProblem_SOL<2>::gradient()

  template<int dim>
  void InhomogeneousProblem_SOL<dim>::vector_gradient(
    const Point<dim>                    &eval,
          std::vector< Tensor<1, dim> > &out
  ) const {
    Assert(dim >= 2 && dim <= 3, ExcNotImplemented());

    out[0] = this->gradient(eval, 0);
    out[1] = Tensor<1, dim>();
    if (dim == 3)
      out[2] = Tensor<1, dim>();
  } // InhomogeneousProblem_SOL<dim>::vector_gradient()


  template<int dim>
  void InhomogeneousProblem_SOL<dim>::gradient_list(
    const std::vector< Point<dim> >         &evals,
          std::vector< Tensor<1, dim> >     &out,
    const unsigned int                       component 
  ) const {
    Assert(evals.size() == out.size(), ExcDimensionMismatch(evals.size(), out.size()));
    const unsigned int n_pnts = evals.size();
    for (unsigned int n = 0; n < n_pnts; n++)
      out[n] = gradient(evals[n], component);
  } // InhomogeneousProblem_SOL<dim>::gradient_list()

  template<int dim>
  void InhomogeneousProblem_SOL<dim>::vector_gradients(
    const std::vector< Point<dim> >                     &evals,
          std::vector< std::vector< Tensor<1, dim> > >  &out
  ) const {
    Assert(evals.size() == out.size(), ExcDimensionMismatch(evals.size(), out.size()));
    const unsigned int n_pnts = evals.size();
    for (unsigned int n = 0; n < n_pnts; n++)
      vector_gradient(evals[n], out[n]); 
  } // InhomogeneousProblem_SOL<dim>::vector_gradients()

  template<int dim>
  double InhomogeneousProblem_SOL<dim>::divergence(
    const Point<dim>    &eval
  ) const {
    Assert(dim == 2 || dim == 3, ExcNotImplemented());
    if (dim == 2)
      return gradient(eval, 0)[0] + gradient(eval, 1)[1]; 
    else 
      return gradient(eval, 0)[0] + gradient(eval, 1)[1] + gradient(eval, 2)[2];
  } // InhomogeneousProblem_SOL<dim>::divergence()

  template<int dim>
  void InhomogeneousProblem_SOL<dim>::divergence_list(
    const std::vector< Point<dim> >     &evals,
          std::vector< double >         &out
  ) const {
    Assert(evals.size() == out.size(), ExcDimensionMismatch(evals.size(), out.size()));
    const unsigned int n_pnts = evals.size();
    for (unsigned int n = 0; n < n_pnts; n++)
      out[n] = divergence(evals[n]);
  } // InhomogeneousProblem_SOL<dim>::divergence_list()

  template<int dim>
  double InhomogeneousProblem_SOL<dim>::gradient_divergence(
    const Point<dim>     &eval,
    const unsigned int    component
  ) const {
    Assert(false, ExcInternalError());
    return 0;
  }

  template<>
  double InhomogeneousProblem_SOL<2>::gradient_divergence(
    const Point<2>       &eval,
    const unsigned int    component
  ) const {
    const double &y  = eval(1); 
    const double &pi = numbers::PI;
          double  out;
    switch(component) {
      case 1:
        out = std::cos(pi * y) / pi; 
      break; 
      default:
        out = 0;
      break;
    } // switch

    return out;
  } // InhomogeneousProblem_SOL<2>::gradient_divergence()

  template<>
  double InhomogeneousProblem_SOL<3>::gradient_divergence(
    const Point<3>       &eval,
    const unsigned int    component
  ) const {
    const double &y  = eval(1); 
    const double &z  = eval(2); 
    const double sqrt_z = std::sqrt(z * (1-z));
    const double z_term = sqrt_z * (z*(1-z)); 
    const double &pi = numbers::PI;
          double  out;
    switch(component) {
      case 1:
        out = z_term * std::cos(pi * y) / pi; 
      break; 
      case 2: 
        out = 1.5 * (1.-2.*z) * sqrt_z * std::sin(pi * y) / (pi * pi);
      break;
      default:
        out = 0;
      break;
    } // switch

    return out;
  } // InhomogeneousProblem_SOL<3>::gradient_divergence()

  template<int dim>
  void InhomogeneousProblem_SOL<dim>::gradient_divergence_list(
    const std::vector< Point<dim> >       &evals,
          std::vector< double >           &out,
    const unsigned int                     component
  ) const {
    Assert(out.size() == evals.size(), ExcDimensionMismatch(out.size(), evals.size()));
    const unsigned int n_pnts = out.size();
    for (unsigned int n = 0; n < n_pnts; n++)
      out[n] = gradient_divergence(evals[n], component);
  } // InhomogeneousProblem_SOL<dim>::gradient_divergence_list()


  template<int dim>
  void InhomogeneousProblem_SOL<dim>::gradient_divergence_tensor(
    const Point< dim >          &eval,
          Tensor<1, dim>        &out
  ) const {
    out[0] = 0;
    out[1] = this -> gradient_divergence(eval, 1); 
    if (dim == 3)
      out[2] = this -> gradient_divergence(eval, 2); 
  } // InhomogeneousProblem_SOL<dim>::gradient_divergence_tensor()


  template<int dim>
  void InhomogeneousProblem_SOL<dim>::gradient_divergence_tensor_list(
    const std::vector< Point<dim> >         &evals,
          std::vector< Tensor<1, dim> >     &out
  ) const {
    Assert(out.size() == evals.size(), ExcDimensionMismatch(out.size(), evals.size()));
    const unsigned int n_pnts = out.size();
    for (unsigned int n = 0; n < n_pnts; n++)
      gradient_divergence_tensor(evals[n], out[n]);
  } // InhomogeneousProblem_SOL<dim>::gradient_divergence_tensor_list()

  template<int dim>
  double InhomogeneousProblem_SOL<dim>::divergence_gradient(
    const Point<dim>     &eval,
    const unsigned int    component
  ) const {
    Assert(false, ExcNotImplemented());
    return 0;
  } // InhomogeneousProblem_SOL<dim>::divergence_gradient()

  template<>
  double InhomogeneousProblem_SOL<2>::divergence_gradient(
    const Point<2>     &eval,
    const unsigned int    component
  ) const {
    if (component != 0)
      return 0;
    const double &x  = eval(0); 
    const double &y  = eval(1); 
    const double &pi = numbers::PI;
    return -1. * (x + 0.5) * std::sin(pi * y);
  } // InhomogeneousProblem_SOL<2>::divergence_gradient()

  template<>
  double InhomogeneousProblem_SOL<3>::divergence_gradient(
    const Point<3>     &eval,
    const unsigned int    component
  ) const {
    if (component != 0)
      return 0; 

    const double &x  = eval(0); 
    const double &y  = eval(1); 
    const double &z  = eval(2); 
    const double  sqrt_z = std::sqrt(z * (1-z));
    const double  z_term = sqrt_z * (z * (1-z));
    const double &pi = numbers::PI;
    // const double out = 
    //         (6 * (x+0.5) * (z *z - z + 0.125) * std::sin(pi * y)) / sqrt_z 
    //         - pi * pi * (x+0.5) * z_term * std::sin(pi * y);
    // return out;
    const double out = 
            -1. * z_term * pi * pi - 3. * sqrt_z + 0.75 * (1-2*z) * (1-2*z) / sqrt_z;
    return (x+0.5) * out * std::sin(pi * y) / (pi * pi);
  } // InhomogeneousProblem_SOL<3>::divergence_gradient()

  template<int dim>
  void InhomogeneousProblem_SOL<dim>::divergence_gradient_list(
    const std::vector< Point<dim> >       &evals,
          std::vector< double >           &out,
    const unsigned int                     component
  ) const {
    Assert(out.size() == evals.size(), ExcDimensionMismatch(out.size(), evals.size()));
    const unsigned int n_pnts = out.size();
    for (unsigned int n = 0; n < n_pnts; n++)
      out[n] = divergence_gradient(evals[n], component);
  } // InhomogeneousProblem_SOL<dim>::divergence_gradient_list()

  template<int dim>
  void InhomogeneousProblem_SOL<dim>::divergence_gradient_tensor(
    const Point< dim >          &eval,
          Tensor<1, dim>        &out
  ) const {
    Assert(dim <= 3 && dim >= 2, ExcNotImplemented());
    out[0] = this->divergence_gradient(eval, 0);
    out[1] = 0.;
    if (dim == 3)
      out[2] = 0.;
  } // InhomogeneousProblem_SOL<dim>::divergence_gradient_tensor()

  template<int dim>
  void InhomogeneousProblem_SOL<dim>::divergence_gradient_tensor_list(
    const std::vector< Point<dim> >         &evals,
          std::vector< Tensor<1, dim> >     &out
  ) const {
    Assert(out.size() == evals.size(), ExcDimensionMismatch(out.size(), evals.size()));
    const unsigned int n_pnts = out.size();
    for (unsigned int n = 0; n < n_pnts; n++)
      gradient_divergence_tensor(evals[n], out[n]);
  } // InhomogeneousProblem_SOL<dim>::divergence_gradient_tensor_list()

  template<int dim>
  Tensor<1, dim> InhomogeneousProblem_SOL<dim>::stress(
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
  } // InhomogeneousProblem_SOL<dim>::stress()

  template<int dim>
  Tensor<2, dim> InhomogeneousProblem_SOL<dim>::stress(
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
  } // InhomogeneousProblem_SOL<dim>::stress()

  // ==============================================================
  //
  //                   Problem
  //
  // ==============================================================
  
  template<int dim>
  InhomogeneousProblem<dim>::AssemblyScratchData::AssemblyScratchData(
          TS_Triangulation<dim, dim>  *tria,
    const std::vector<unsigned int>   &degrees
  ) : tria(tria) 
    , ts_values(this->tria, 
                degrees, 
                update_values |
                update_gradients |
                update_JxW_values |
                update_quadrature_points) {}

  template<int dim>
  InhomogeneousProblem<dim>::AssemblyScratchData::AssemblyScratchData(
    const AssemblyScratchData &scratch
  ) : tria(scratch.tria)
    , ts_values(scratch.ts_values) {}


  template<int dim> 
  InhomogeneousProblem<dim>::InhomogeneousProblem(
    const unsigned int ref,
    const unsigned int order_elevate,
    const Strategy     strategy
  ) : data(DataGenerator<dim>(true).data)
    , ref(ref)
    , order(order_elevate + data.max_degree())
    , tria(data)
    , strategy(strategy)
    , sol_fcn()
    , rhs_fcn(&sol_fcn)
    , lambda(1.)
    , mu(1.)
  { 
    Assert(dim == 2 || dim == 3, 
            ExcDimensionMismatch2(dim, 2, 3));
    dealt::OutputSetup::set_dim(3);
    if (strategy == Linear_Elasticity::Strategy::Adaptive){
      const std::string name = "linear_elasticity_inhomogeneous/" + std::to_string(dim) + "d/";
      problem_out = OutputSetup(name, order);
    } else {
      const std::string name = "linear_elasticity_inhomogeneous/" 
                             + std::to_string(dim) 
                             + "d_uniform/";
      problem_out = OutputSetup(name, order);
    }

    tria.degree_elevate_global(order_elevate);

    this->offset = 0;
    tria.refine_global(offset);
          
    // Set boundary indicators: 
    if (dim == 2) {
      for (const auto& face : tria.active_face_iterators())
        if (face -> at_boundary())
          if ((face -> center()).operator()(1) == 0 || 
              (face -> center()).operator()(1) == 1)
            face -> set_boundary_id(Boundaries::Dirichleth0);
          else 
            face ->set_boundary_id(Boundaries::Dirichleth);
    } else {
      for (const auto& face : tria.active_face_iterators())
        if (face -> at_boundary())
          if ((face -> center()).operator()(1) == 0 || 
              (face -> center()).operator()(1) == 1 || 
              (face -> center()).operator()(2) == 0 ||
              (face -> center()).operator()(2) == 1)
            face -> set_boundary_id(Boundaries::Dirichleth0);
          else 
            face ->set_boundary_id(Boundaries::Dirichleth);
    }


    tria.prepare_assembly();
  } // InhomogeneousProblem<dim>::InhomogeneousProblem()

  template<int dim>
  void InhomogeneousProblem<dim>::local_assemble_system(
    const typename Triangulation<dim, dim>::active_cell_iterator &cell,
    AssemblyScratchData                                          &scratch,
    AssemblyCopyData                                             &copy_data
  ) {
    // std::cout << "    Initializing TSValues..." << std::endl;
    scratch.ts_values.reinit(cell);
    // std::cout << "    ... done!" << std::endl;

    const unsigned int n_dofs_per_cell = scratch.ts_values.n_dofs_per_cell();
    copy_data.cell_matrix.reinit(n_dofs_per_cell, n_dofs_per_cell);
    copy_data.cell_rhs.reinit(n_dofs_per_cell);
    copy_data.local_dof_indices.resize(n_dofs_per_cell);

    // Compute system_matrix part of cell
    // std::cout << "    Running assembly routine..." << std::endl;
    for (const unsigned int i : scratch.ts_values.dof_indices()){
      const unsigned int comp_i =
              scratch.ts_values.system_to_component_index(i).first;
      for (const unsigned int j : scratch.ts_values.dof_indices()){
        const unsigned int comp_j =
                scratch.ts_values.system_to_component_index(j).first;
        for (const unsigned int q : scratch.ts_values.quadrature_point_indices()){
          copy_data.cell_matrix(i, j) += (
            scratch.ts_values.shape_grad(i, q)[comp_i] *
            scratch.ts_values.shape_grad(j, q)[comp_j] *
            lambda.value(scratch.ts_values.quadrature_point(q))
            +
            scratch.ts_values.shape_grad(i, q)[comp_j] *
            scratch.ts_values.shape_grad(j, q)[comp_i] * 
            mu.value(scratch.ts_values.quadrature_point(q))
            +
            (
              (comp_i == comp_j) ? 
                scratch.ts_values.shape_grad(i, q) * 
                scratch.ts_values.shape_grad(j, q) * 
                mu.value(scratch.ts_values.quadrature_point(q)) : 
                0.
            )
          ) * scratch.ts_values.JxW(q);
        } // for ( q )
      } // for ( j )

      // Compute rhs for this cell
      for (const unsigned int q : scratch.ts_values.quadrature_point_indices())
        copy_data.cell_rhs(i) += scratch.ts_values.shape_value(i, q) * 
                        rhs_fcn.value(scratch.ts_values.quadrature_point(q), comp_i) *
                        scratch.ts_values.JxW(q);
    } // for ( i ) 
    // std::cout << "    ... done!" << std::endl;

    // std::cout << "    Copying data..." << std::endl;
    copy_data.local_dof_indices = scratch.tria -> get_IEN_array(cell, dim);
    // std::cout << "    ... done!" << std::endl;
  } // InhomogeneousProblem<dim>::local_assemble_system()

  template<int dim>
  void InhomogeneousProblem<dim>::copy_local_to_global(
    const AssemblyCopyData &copy_data
  ) {
    system_matrix.add(copy_data.local_dof_indices, copy_data.cell_matrix);
    system_rhs.add(copy_data.local_dof_indices, copy_data.cell_rhs);
  } // InhomogeneousProblem<dim>::copy_local_to_global();

  template<int dim>
  void InhomogeneousProblem<dim>::assemble_system(
  ) {
    std::vector< unsigned int > degrees = tria.get_degree();
    for (unsigned int d = 0; d < dim; d++)
      degrees[d] += 1;

    AssemblyCopyData copy_data;
    AssemblyScratchData scratch(&tria, degrees);
    for (const auto& cell : tria.active_cell_iterators()) {
      local_assemble_system(
        cell, 
        scratch,
        copy_data
      );
      copy_local_to_global(copy_data);
    }
    // WorkStream::run(
    //   tria.begin_active(),
    //   tria.end(),
    //   *this,
    //   &InhomogeneousProblem::local_assemble_system,
    //   &InhomogeneousProblem::copy_local_to_global,
    //   AssemblyScratchData(&tria, degrees),
    //   AssemblyCopyData(),
    //   1
    // );

    std::cout << "Applying Boundary values ..." << std::endl;
    std::map<
        types::global_dof_index, 
        double
      > boundary_values;
    std::map<
        types::boundary_id,
        const Function< dim >*
      > boundary_fcns =
          {{Boundaries::Dirichleth, &sol_fcn}};
    std::cout << "    Projecting boundary values..." << std::endl;
    tria.template project_boundary_values<dim>(
      boundary_fcns,
      degrees,
      boundary_values
    );
    std::cout << "    ... done! Applying inhomogeneous boundary values ..." << std::endl;

    for (const auto& [dof, _] : boundary_values) 
      if (dof > dim * tria.n_active_splines())
        throw ExcInternalError();

    if (system_matrix.n() != system_matrix.m())
      throw ExcInternalError();

    if (system_matrix.n() != solution.size())
      throw ExcInternalError();


    // output_prior_apply_boundary_values(boundary_values);

    MatrixTools::apply_boundary_values(
      boundary_values, 
      system_matrix,
      solution,
      system_rhs
    );
    std::cout << "    ... done! Applying homogeneous boundary values ..." << std::endl;
    
    // Apply homogeneous boundary conditions
    const auto& boundary_dofs = tria.get_boundary_dofs(dim);
    for (const auto& dof : boundary_dofs.at(Boundaries::Dirichleth0)){
      for (unsigned int j = 0; j < dim * tria.n_active_splines(); j++) {
        system_matrix.set(dof, j, 0);
        system_matrix.set(j, dof, 0);
      }
      system_matrix.set(dof, dof, 1);
      system_rhs(dof) = 0;
    }
    std::cout << "    Homogeneous Data applied" << std::endl;
    std::cout << "... done!" << std::endl;

  } // InhomogeneousProblem<dim>::assemble_system()

  template<int dim>
  void InhomogeneousProblem<dim>::setup_system(
  ) {
    const unsigned int n_dofs = dim * tria.n_active_splines();
    system_rhs.reinit(n_dofs);
    solution.reinit(n_dofs);

    const auto& IEN_array = tria.get_IEN_array(dim);
    Assert(IEN_array.size() > 0, ExcInternalError());

    sparsity_pattern.reinit(
        n_dofs,
        n_dofs,
        tria.get_max_entries_per_row(dim)
    );

    for (const auto& [_, arr] : IEN_array)
      for (unsigned int i : arr)
        for (unsigned int j : arr)
          sparsity_pattern.add(i, j);

    sparsity_pattern.compress();

    system_matrix.reinit(sparsity_pattern);
  } // InhomogeneousProblem<dim>::setup_system()

  template<int dim>
  void InhomogeneousProblem<dim>::solve(
  ) {
    std::cout << "Solving system ... " << std::endl;

    SolverControl solver_control(
                     750 * dim * tria.n_active_splines(), 
                     1e-10
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
      tria.n_levels() - 1 - offset,
      tria.n_active_cells(),
      tria.n_active_splines() * dim,
      solver_control.last_step(),
      *(solver_control.get_history_data().end() - 1)
    );
  
    std::cout << " ... done!" << std::endl;
  } // InhomogeneousProblem<dim>::solve()

  template<int dim>
  void InhomogeneousProblem<dim>::compute_h1_error(
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
  } // InhomogeneousProblem<dim>::compute_h1_error()

  template<int dim>
  void InhomogeneousProblem<dim>::estimate_mark_refine(
  ) { 
    if (strategy == Strategy::Adaptive) {
      Vector<double> residuals(tria.n_active_cells()) ;
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
      tria.refine_fixed_number(residuals, 0.10);
    } else { 
      tria.coarsen_bezier_elements();
      tria.refine_global();
    }

    tria.prepare_assembly();
  } // InhomogeneousProblem<dim>::estimate_mark_refine()

  template<int dim>
  void InhomogeneousProblem<dim>::run(
  ) {
    
    unsigned int cycle = 0;  
    unsigned int old_level = tria.n_levels() - 1; 
    unsigned int r = 0;
    while (tria.n_levels() - 1 < ref) {
      std::cout << "Logging iteration: " << std::to_string(r++) << std::endl;
      std::cout << "Setting up system..." << std::endl;
      this -> setup_system();
      std::cout << "... done! Assembling system..." << std::endl;
      this -> assemble_system();
      std::cout << "... done! Solving system..." << std::endl;
      this -> solve();
      std::cout << "... done! Computing errors..." << std::endl;
      this -> compute_h1_error();
      std::cout << "... done! Printing to files..." << std::endl;
      this -> output_system();
      std::cout << "... done! Refining..." << std::endl;

      if (H1 < 1e-10 )
        break;

      if (cycle < 5) {
        this -> estimate_mark_refine();
      } else { 
        std::cout << "Too many refinements resulted in the same level, enforcing global refinement..." << std::endl;
        tria.coarsen_bezier_elements();
        tria.refine_global(); 
        tria.prepare_assembly();
      }
      std::cout << "... done!" << std::endl;

      if (tria.n_levels() - 1 == old_level) {
        cycle++;
      } else {
        cycle = 0;
        old_level = tria.n_levels() - 1;
      }

    } // while

    // Write the resulting table to line
    problem_out.write_table_text(std::cout);
  } // InhomogeneousProblem<dim>::run()

  template<int dim>
  void InhomogeneousProblem<dim>::output_system(
  ) {
    Assert(dim == 2 || dim == 3, ExcDimensionMismatch2(dim, 2, 3));
  }

  template<int dim>
  void InhomogeneousProblem<dim>::output_prior_apply_boundary_values(
    const std::map<types::global_dof_index, double>& boundary_values
  ) {
    const unsigned int level = tria.n_levels() - 1 - offset;

    const std::string name = problem_out.degree.string();
    const std::string level_name = name + "l" + std::to_string(level);
  
    // std::string matrix_dat = level_name + "_mat.dat" ;
    std::string matrix_blo = level_name + "_mat.block" ;
    // std::string bv_dat     = level_name + "_bv.dat" ;
    std::string bv_blo     = level_name + "_bv.block" ;

    // std::filebuf mat_dat_buf, bv_dat_buf;
    std::filebuf mat_blo_buf, bv_blo_buf;
    // mat_dat_buf.open(matrix_dat.c_str(), std::ios::out);
    mat_blo_buf.open(matrix_blo.c_str(), std::ios::out);
    // bv_dat_buf.open(bv_dat.c_str(), std::ios::out);
    bv_blo_buf.open(bv_blo.c_str(), std::ios::out);
  
    // std::ostream mat_dat_out(&mat_dat_buf), bv_dat_out(&bv_dat_buf);
    std::ostream mat_blo_out(&mat_blo_buf), bv_blo_out(&bv_blo_buf);
  
    // system_matrix.print_formatted(mat_dat_out, 16, true, 1, "0");
    system_matrix.block_write(mat_blo_out); 

    SparsityPattern sparsity_pattern(boundary_values.size(), 2, 2);
    for (unsigned int i = 0; i < boundary_values.size(); i++){
      sparsity_pattern.add(i, 0);
      sparsity_pattern.add(i, 1);
    }
    sparsity_pattern.compress();

    SparseMatrix<double> B(sparsity_pattern);
    auto it = boundary_values.begin(); 
    std::cout << "Translating bv to sparse matrix" << std::endl;
    for (unsigned int i = 0; i < boundary_values.size(); i++){
      B.set(i, 0, it -> first );
      B.set(i, 1, it -> second);
      it++;
    }
    std::cout << "done!" << std::endl;

    // B.print_formatted(bv_dat_out, 16, true, 1, "0");
    B.block_write(bv_blo_out);
  }

  template<>
  void InhomogeneousProblem<3>::output_system(
  ) {
    const unsigned int level = tria.n_levels() - 1 - offset;


    // Write the grid to a seperate file: 
    const std::string& name_vtg = problem_out.vtg.string() + "physical_grid_l" + std::to_string(level) + ".vtk";

    // First: Make a copy of the triangulation: 
    Triangulation<3> physical_grid; 
    physical_grid.copy_triangulation(tria);


    // GridTools::transform(...) does not work with anisotropically refined meshes in 3D
    // Thus, we perform the transformation manually ... *sigh*
    const IsoparametricManifold<3> geometry(tria.get_IPF()); 

    std::vector<bool> treated_vertices(tria.n_vertices(), false);
    for (const auto& cell : physical_grid.active_cell_iterators())
      for (const unsigned int v : cell->vertex_indices())
        if (treated_vertices[cell->vertex_index(v)] == false) {
          cell->vertex(v) = geometry.push_forward(cell->vertex(v));
          treated_vertices[cell->vertex_index(v)] = true;
        }


    // Generate the output object
    DataOut<3> data_out;
    data_out.attach_triangulation(physical_grid); 

    // Estimate the error on each cell
    Vector< double > residuals(tria.n_active_cells());
    std::vector<unsigned int> degrees = tria.get_degree();
    for (unsigned int p : degrees)
      p = p * p + 1;

    tria.linear_elasticity_residual_error_estimate(
        degrees,
        &rhs_fcn, 
        &lambda,
        &mu,
        solution,
        residuals
    );

    data_out.add_data_vector(residuals, "cell_errors");

    // Build patches
    data_out.build_patches(); 

    // Open a file
    std::ofstream vtk_out(name_vtg); 
    data_out.write_vtk(vtk_out);

    problem_out.write_table_text();
    problem_out.write_table_tex();


    for (double z = 0.1; z <= 1; z += 0.3) {
      GridPartialOut::get_cross_section(&tria, 2, z, true, 
                        problem_out.cross_sections_z.string() 
                        + "l" 
                        + std::to_string(level) 
                        + "_physical");
      GridPartialOut::get_cross_section(&tria, 2, z, false, 
                        problem_out.cross_sections_z.string() 
                        + "l" 
                        + std::to_string(level) 
                        + "_parametric");
    }
    for (double y = 0.1; y <= 1; y += 0.3) {
      GridPartialOut::get_cross_section(&tria, 1, y, true, 
                        problem_out.cross_sections_y.string() 
                        + "l" 
                        + std::to_string(level) 
                        + "_physical");
      GridPartialOut::get_cross_section(&tria, 1, y, false, 
                        problem_out.cross_sections_y.string() 
                        + "l" 
                        + std::to_string(level) 
                        + "_parametric");
    }
    for (double x = 0.1; x <= 1; x += 0.3) {
      GridPartialOut::get_cross_section(&tria, 0, x, true, 
                        problem_out.cross_sections_x.string() 
                        + "l" 
                        + std::to_string(level) 
                        + "_physical");
      GridPartialOut::get_cross_section(&tria, 0, x, false, 
                        problem_out.cross_sections_x.string() 
                        + "l" 
                        + std::to_string(level) 
                        + "_parametric");
    }

    if (level > 15)
      return;


    const std::string name = problem_out.degree.string();
    const std::string level_name = name + "l" + std::to_string(level);
  
    std::string matrix = level_name + "_mat.dat" ;
    std::string vector = level_name + "_vec.dat" ;
    std::string soluti = level_name + "_sol.dat" ;

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

    const unsigned int N1 = 50; 
    const unsigned int N2 = 25;
    const unsigned int N3 = 10;
    const double xmin = 0; 
    const double ymin = 0; 
    const double zmin = 0; 
    const double xmax = 1; 
    const double ymax = 1; 
    const double zmax = 1; 
    std::vector<Point<3>> evals; 
    FullMatrix< double > E(N1 * N2 * N3, 3);
    unsigned int ind = 0;
    for (double k = 0.; k < N3; k++){
      for (double j = 0.; j < N2; j++) {
        for (double i = 0.; i < N1; i++) {
          const double x = xmin + (i/(N1-1.)) * (xmax - xmin); 
          const double y = ymin + (j/(N2-1.)) * (ymax - ymin); 
          const double z = zmin + (k/(N3-1.)) * (zmax - zmin); 
          evals.push_back(Point<3>(x, y, z));
          E(ind, 0) = x;
          E(ind, 1) = y;
          E(ind, 2) = z;
          ind++;
        }
      }
    }

    std::filebuf e_f;
    e_f.open(level_name + "_evals.dat", std::ios::out);
    std::ostream out_e(&e_f);
    E.print_formatted(out_e, 16, true, 1, "0");

    tria.printIPF(3, evals, level_name, 16, true, true);
  } // InhomogeneousProblem<3>::output_system()

  template<>
  void InhomogeneousProblem<2>::output_system(
  ) {
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
      [&geometry](const Point<2>& p){
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
  } // InhomogeneousProblem<2>::output_system();

  template<int dim>
  void InhomogeneousProblem<dim>::apply_boundary_values(
    const std::map<types::global_dof_index, double> &boundary_values,
    SparseMatrix<double>                            &matrix,
    Vector<double>                                  &solution,
    Vector<double>                                  &right_hand_side,
    const bool                                       eliminate_columns
  ) {
    Assert(matrix.n() == right_hand_side.size(),
           ExcDimensionMismatch(matrix.n(), right_hand_side.size()));
    Assert(matrix.n() == solution.size(),
           ExcDimensionMismatch(matrix.n(), solution.size()));
    Assert(matrix.n() == matrix.m(),
           ExcDimensionMismatch(matrix.n(), matrix.m()));

    // if no boundary values are to be applied
    // simply return
    if (boundary_values.empty())
      return;


    const types::global_dof_index n_dofs = matrix.m();

    // if a diagonal entry is zero
    // later, then we use another
    // number instead. take it to be
    // the first nonzero diagonal
    // element of the matrix, or 1 if
    // there is no such thing
    double first_nonzero_diagonal_entry = 1;
    std::cout << "      Finding first nonzero diagonal entry..." << std::endl;
    for (unsigned int i = 0; i < n_dofs; ++i)
      if (matrix.diag_element(i) != 0.)
        {
          first_nonzero_diagonal_entry = matrix.diag_element(i);
          break;
        }
    if (first_nonzero_diagonal_entry != matrix.diag_element(0)){
      std::cout << "      ... failed! Log: " << std::endl;
      std::cout << "      Should be " 
                << std::to_string(matrix.diag_element(0))
                << " is "
                << std::to_string(first_nonzero_diagonal_entry) 
                << std::endl;
      throw ExcInternalError();
    }

    std::cout << "      ... done! Applying values ... " << std::endl;
    std::map<types::global_dof_index, double>::const_iterator
      dof  = boundary_values.begin(),
      endd = boundary_values.end();
    for (; dof != endd; ++dof)
      {
        Assert(dof->first < n_dofs, ExcInternalError());
        if (dof->first >= n_dofs){
          std::cout << "      ... failed! Log: " << std::endl;
          std::cout << "      DoF "
                    << std::to_string(dof->first) 
                    << " is out of bounds. Upper bound: "
                    << std::to_string(n_dofs)
                    << std::endl;
          throw ExcInternalError();
        }

        const types::global_dof_index dof_number = dof->first;
        // for each boundary dof:

        // set entries of this line to zero except for the diagonal
        // entry
        for (SparseMatrix<double>::iterator p =
               matrix.begin(dof_number);
             p != matrix.end(dof_number);
             ++p)
          if (p->column() != dof_number)
            p->value() = 0.;

        // set right hand side to
        // wanted value: if main diagonal
        // entry nonzero, don't touch it
        // and scale rhs accordingly. If
        // zero, take the first main
        // diagonal entry we can find, or
        // one if no nonzero main diagonal
        // element exists. Normally, however,
        // the main diagonal entry should
        // not be zero.
        //
        // store the new rhs entry to make
        // the gauss step more efficient
        double new_rhs;
        if (matrix.diag_element(dof_number) != 0.)
          {
            new_rhs = dof->second * matrix.diag_element(dof_number);
            right_hand_side(dof_number) = new_rhs;
          }
        else
          { // Shouldn't happen with us
            std::cout << "      ... failed! Log: " << std::endl;
            std::cout << "      Diagonal of DoF " 
                      << std::to_string(dof_number) 
                      << " is zero!" 
                      << std::endl;
            throw ExcInternalError();

            matrix.set(dof_number, dof_number, first_nonzero_diagonal_entry);
            new_rhs = dof->second * first_nonzero_diagonal_entry;
            right_hand_side(dof_number) = new_rhs;
          }


        // if the user wants to have
        // the symmetry of the matrix
        // preserved, and if the
        // sparsity pattern is
        // symmetric, then do a Gauss
        // elimination step with the
        // present row
        if (eliminate_columns)
          {
            std::cout << "    Eliminating columns..." << std::endl;
            // store the only nonzero entry
            // of this line for the Gauss
            // elimination step
            const double diagonal_entry = matrix.diag_element(dof_number);

            // we have to loop over all rows of the matrix which have
            // a nonzero entry in the column which we work in
            // presently. if the sparsity pattern is symmetric, then
            // we can get the positions of these rows cheaply by
            // looking at the nonzero column numbers of the present
            // row. we need not look at the first entry of each row,
            // since that is the diagonal element and thus the present
            // row
            for (SparseMatrix<double>::iterator q =
                   matrix.begin(dof_number) + 1;
                 q != matrix.end(dof_number);
                 ++q)
              {
                std::cout << "      Entering loop" << std::endl;
                const types::global_dof_index row = q->column();
                std::cout << "        row found" << std::endl;
                if (row >= n_dofs){ 
                  std::cout << "      ...failed! Log: " << std::endl;
                  std::cout << "      Accessing row "
                            << std::to_string(row)
                            << " is out of bounds: "
                            << std::to_string(n_dofs) 
                            << std::endl;
                  throw ExcInternalError();
                }


                // find the position of
                // element
                // (row,dof_number)
                // bool (*comp)(
                //   const SparseMatrix<double>::iterator::value_type p,
                //   const unsigned int column) =
                //   &column_less_than<SparseMatrix<double>::iterator>;
                std::cout << "        Searching column" << std::endl;
                const SparseMatrix<double>::iterator p =
                  Utilities::lower_bound(
                      matrix.begin(row) + 1,
                      matrix.end(row),
                      dof_number,
                      [](const SparseMatrix<double>::iterator::value_type p, 
                         const unsigned int                               column
                        ) {
                          return p.column() < column;
                        }
                  );
                std::cout << "        column found" << std::endl;

                // check whether this line has an entry in the
                // regarding column (check for ==dof_number and !=
                // next_row, since if row==dof_number-1, *p is a
                // past-the-end pointer but points to dof_number
                // anyway...)
                //
                // there should be such an entry! we know this because
                // we have assumed that the sparsity pattern is
                // symmetric and we only walk over those rows for
                // which the current row has a column entry
                Assert((p != matrix.end(row)) && (p->column() == dof_number),
                       ExcMessage(
                         "This function is trying to access an element of the "
                         "matrix that doesn't seem to exist. Are you using a "
                         "nonsymmetric sparsity pattern? If so, you are not "
                         "allowed to set the eliminate_column argument of this "
                         "function, see the documentation."));

                if ( !((p != matrix.end(row)) && (p->column() == dof_number)) ) {
                  std::cout << "      ...failed! Log:" << std::endl;
                  std::cout << "This function is trying to access an element of the "
                            << "matrix that doesn't seem to exist. Are you using a "
                            << "nonsymmetric sparsity pattern? If so, you are not "
                            << "allowed to set the eliminate_column argument of this "
                            << "function, see the documentation.";
                    throw ExcMessage(
                         "This function is trying to access an element of the "
                         "matrix that doesn't seem to exist. Are you using a "
                         "nonsymmetric sparsity pattern? If so, you are not "
                         "allowed to set the eliminate_column argument of this "
                         "function, see the documentation.");
                }


                // correct right hand side
                right_hand_side(row) -=
                  static_cast<double>(p->value()) / diagonal_entry * new_rhs;

                std::cout << "        rhs updated" << std::endl;

                // set matrix entry to zero
                p->value() = 0.;

                std::cout << "        value set to 0" << std::endl;

                std::cout << "      end of loop" << std::endl;
              }
          }

        // preset solution vector
        std::cout << "    Setting solution..." << std::endl;
        solution(dof_number) = dof->second;
        std::cout << "    ...successfull!" << std::endl;
      }
  }

} // namespace
