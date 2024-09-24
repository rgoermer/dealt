#include <poisson.h>


namespace Poisson_Benchmark {
  double Geometry_PF::value(
    const Point<2>&       p,
    const unsigned int    component 
  ) const {
    if (component == 1)
      return p(1); 

    const double& x = p(0);
    const double& y = p(1);

    return x*x*y * (-2. + 2. * y) + x * (1. + 2. * y - 2. * y * y);
  } // Geometry_PF::value()

  Tensor<1, 2> Geometry_PF::gradient(
    const Point<2>&       p, 
    const unsigned int    component
  ) const { 
    Tensor<1, 2> out; 
    if (component == 1) {
      out[1] = 1;
      return out; 
    }

    const double& x = p(0);
    const double& y = p(1);

    out[0] = 1. + (2. - 4. * x) * y + (-2. + 4. * x) * y * y;
    out[1] = x * (2. - 4. * y + x * (-2. + 4. * y));

    return out;
  } // Geometry_PF::gradient()

  SymmetricTensor<2, 2> Geometry_PF::hessian(
    const Point<2>&       p,
    const unsigned int    component
  ) const { 
    SymmetricTensor<2, 2> out; 
    if (component == 1)
      return out;

    const double& x = p(0);
    const double& y = p(1); 

    out[0][0] = -4. * y + 4. * y * y; 
    out[0][1] = 2. - 4. * x - 4. * y + 8. * x * y;
    out[1][1] = -4. * x + 4. * x * x;

    return out;
  } // Geometry_PF::hessian()

  double Geometry_PB::value(
    const Point<2>&       p, 
    const unsigned int    component
  ) const { 
    if (component == 1)
      return p(1); 

    const double& x = p(0); 
    const double& y = p(1); 
    const double& ysq = p(1) * p(1);

    const double sqrt = std::sqrt(1. + 4. * y - 8. * x * y + 8. * x * ysq - 8. * ysq * y + 4. * ysq * ysq);

    const double nom = 0.25 * (-1. - 2. * y + 2. * ysq + sqrt);
    if (std::fabs(nom) < 1e-16)
      return 0; 
    else 
      return nom / (-1. * y + ysq);
  } // Geometry_PB::value()

  Tensor< 1, 2> Geometry_PB::gradient(
    const Point<2>&       p, 
    const unsigned int    component
  ) const { 
    Tensor<1, 2> out; 
    if (component == 1) {
      out[1] = 1; 
      return out;
    }

    const double& x = p(0); 
    const double& y = p(1); 
    const double& ysq = p(1) * p(1);

    const double sqrt = std::sqrt(1. + 4. * y - 8. * x * y + 8. * x * ysq - 8. * ysq * y + 4. * ysq * ysq);

    out[0] = 1. / sqrt;
    out[1] = (0.25 - 1.5 * ysq + y * ysq + x * y * (-1. + 3. * y - 2. * ysq) - 0.25 * sqrt + 0.5 * y * sqrt) 
              / 
            ((-1. + y) * (-1. + y) * ysq * sqrt);

    return out;
  } // Geometry_PB::gradient()

  SymmetricTensor<2, 2> Geometry_PB::hessian(
    const Point<2>&       p, 
    const unsigned int    component
  ) const {
    SymmetricTensor< 2, 2> out;
    if (component == 1)
      return out;

    const double& x = p(0); 
    const double& y = p(1); 
    const double& y1 = (-1. + y);
    const double& y1sq = y1 * y1; 
    const double& y1cu = y1sq * y1;
    const double& ysq = p(1) * p(1);
    const double& ycu = p(1) * p(1) * p(1);
    const double& yft = p(1) * p(1) * p(1) * p(1);

    const double root = (1. + 4. * y - 8. * x * y + 8. * x * ysq - 8. * ycu + 4. * yft);
    const double sqrt = std::sqrt(root);
    const double sqrt3 = root * sqrt; 

    out[0][0] = -(4. * (1. - y) * (1. - y) * y) / ((-1. + y) * sqrt3);
    out[0][1] = (-0.25 + x * (0.5 - y) + 1.5 * ysq - ycu)
                    /
                ((0.125 + 0.5 * y - x * y + x * ysq - ycu + 0.5 * yft) * sqrt);

    const double inner1 = (-1. -2. * y + 2. * ysq + sqrt);
    const double inner2 = (-2. + 4. * y + (2. - 12. * ysq + 8. * ycu + x * (-4. + 8. * y)) / sqrt);
    const double inner3 = (4. - 24. * ysq + 16. * ycu + x * (-8. + 16. * y));
    const double inner4 = (4. - inner3 * inner3 / (4. * sqrt3) + (8. * x + 24. * y * y1)/sqrt);
    const double s1 = 0.5 * y1sq   * inner1;
    const double s2 = 0.5 * y1 * y * inner1;
    const double s3 = 0.5 * ysq    * inner1; 
    const double s4 = 0.5 * y1sq * y * inner2;
    const double s5 = 0.5 * y1 * ysq * inner2;
    const double s6 = 0.25 * y1sq * ysq * inner4;

    out[1][1] = (s1 + s2 + s3 - s4 - s5 + s6)/(y1cu * ycu);

    return out;
  } // Geometry_PB::hessian()

  Poisson_Benchmark::Poisson_Benchmark(int ref, int order) :
          ref(ref),                     // Number of refinements to be done
          order(order),
          data(this->get_IPF_data()),
          tria(data),                   // Initialize the grid
          rhs_fcn(),
          sol_fcn()
  {
    IPF_Data<2, 2> data = get_IPF_data();
    std::vector<unsigned int> degrees = data.deg; 
    tria.degree_elevate_global( order );
  
    problem_out = OutputSetup("poisson_benchmark/", data.max_degree() + order); 
  
    offset = 0; 
    tria.refine_global(offset);
  
    tria.set_boundary_dofs();
  
    // l2_error.reinit(ref + 1);
    // h1_error.reinit(ref + 1);
  
    // mesh_size.reinit(ref + 1);
    // dofs_per_level.reinit(ref + 1);
  
    parametric_mapping = tria.get_IPF();
  }
  
  
  
  
  void Poisson_Benchmark::refine_grid()
  {
    tria.execute_coarsening_and_refinement();
  
    tria.set_boundary_dofs();
  
    // We want to perform calculations on this cell, hence we cut along 
    // the TJunction extensions without producing new splines.
    // This sets is_bezier_mesh = true, and allows further computations 
    // along the way. Note, however, that for adaptive refinement strategies,
    // it is needed to coarse the bezier elements again, since there are no
    // new TSplines i defined on these cells. 
    tria.refine_bezier_elements();
  
    // For the assembly of the system matrix, the bezier extraction operators are
    // crucial, hence we compute them preemptively:
    tria.compute_extraction_operators();
  }
  
  
  void Poisson_Benchmark::setup_system(){
    system_rhs.reinit(tria.n_active_splines());
    solution.reinit(tria.n_active_splines());
    
    // To generate the raw structure of the sparsity matrix, we construct our own SparsityPattern
    // This will be done using the IEN_array from the TS_Triangulation. It already stores the 
    // information which spline has support on which spline. 
    //
    // Currently, the maximum number of entries per row is unknown, hence we initialize it
    // as a FullMatrix and compress it afterwards.
    sparsity_pattern.reinit(tria.n_active_splines(), 
                              tria.n_active_splines() ,
                              tria.n_active_splines()
                           ); 
  
    const auto& IEN_array = tria.get_IEN_array();
    Assert(IEN_array.size() > 0, ExcInternalError());
  
    for (const auto& [_, arr] : IEN_array)
      for (unsigned int i : arr)
        for (unsigned int j : arr)
          sparsity_pattern.add(i, j); 
  
    // Free superfluous space
    sparsity_pattern.compress();
    system_matrix.reinit(sparsity_pattern);
  
  }
  
  void Poisson_Benchmark::assemble_system(){
  
    std::cout << "Assembling system matrix ... " << std::endl;
  
    // Setup initial tables that store the bernstein values / grads / hessians.
    std::vector< unsigned int > degrees = tria.get_degree();
    TSValues<2> ts_values( 
            &tria,
            {degrees[0] + 1, degrees[1] + 1},
            update_values |
            update_gradients |
            update_jacobians |
            update_inverse_jacobians |
            update_quadrature_points |
            update_JxW_values);
  
    // Geometry_PF phi;
    // Geometry_PB inv;
    // ts_values.test_geometry_mapping(
    //   &phi, 
    //   &inv, 
    //   1e-12, 
    //   true
    // );
  
    const unsigned int dofs_per_cell = ts_values.n_dofs_per_cell();
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<double>     cell_rhs(dofs_per_cell);
  
  
    // Get a map that corresponds to boundary values
    const auto& boundary_dofs = 
            tria.get_boundary_dofs();

    for (const auto& cell : tria.active_cell_iterators()){
      // Get the Bernstein values on the cell
      ts_values.reinit(cell); 
  
      // Reset the cell matrix
      cell_matrix       = 0; 
      cell_rhs          = 0;
  
      // Quadrature sum: 
      for (const unsigned int q_index : ts_values.quadrature_point_indices()){
  
        // Build the cell matrix: 
        double dx = ts_values.JxW(q_index); 
        for (const unsigned int i : ts_values.dof_indices())
          for(const unsigned int j : ts_values.dof_indices())
            cell_matrix(i,j) += ( ts_values.shape_grad(i, q_index) 
                                  * ts_values.shape_grad(j, q_index)
                                  * dx ); 
  
        // Map the quadrature point from real cell to parametric cell
        const Point<2>& mapped_q = ts_values.quadrature_point(q_index); 
        const double rhs = rhs_fcn.value(mapped_q); 
        for (const unsigned int i : ts_values.dof_indices())
          cell_rhs(i) += ( ts_values.shape_value(i, q_index) 
                            * rhs
                            * dx ); 
  
      } // for ( q_index )
  
    
      // Add the values to the system
      std::vector< unsigned int > local_dof_indices = 
              tria.get_IEN_array(cell);
  
      system_matrix.add(local_dof_indices, cell_matrix);
      system_rhs.add(local_dof_indices, cell_rhs);
  
    } // for ( cell )
  
    // Use the boundary dofs to eliminate rows and columns with boundary indicators
    unsigned int n_global_dofs = tria.n_active_splines();
    for (const auto& dof : boundary_dofs.at(0)){
      for (unsigned int j = 0; j < n_global_dofs; j++){
        system_matrix.set(dof, j, 0);
        system_matrix.set(j, dof, 0);
      }
  
      system_rhs(dof) = 0; 
      system_matrix.set(dof, dof, 1); 
    }
  
    std::cout << " ... done!" << std::endl; 
  }
  
  void Poisson_Benchmark::output_system()
  {
    std::cout << "Printing system matrix and rhs ... " << std::endl;
  
    const unsigned int level = tria.n_levels() - offset - 1; 
    const std::string l = std::to_string(level);
    const std::string name       = problem_out.degree.string();
    const std::string svg_name   = problem_out.svg.string();
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
    tria.printIPF(evals, level_name, 16, true, true);
    tria.coarsen_bezier_elements();
    tria.generate_mesh_file<0>(level_name, false, 16);
    tria.generate_mesh_file<0>(level_name, true, 16);
    tria.print_IPF_wireframe(level_name);
    tria.refine_bezier_elements();
  
  
    // Print the IPF wireframe:
    tria.print_grid(level_name);
    tria.print_IPF_wireframe(level_name ); 
    
    this -> print_grid(svg_name);
  
  
    // Print the table to a file preemptively
    // output to .txt
    problem_out.write_table_text();
    problem_out.write_table_tex();
  
    std::cout << " ... output to files done!" << std::endl; 
  } // output_system
  
  void Poisson_Benchmark::solve(){
    std::cout << "Solving system ... " << std::endl;
  
    SolverControl            solver_control(750 * tria.n_active_splines(), H1 * 1e-4);
    SolverCG<Vector<double>> solver(solver_control);
  
    PreconditionJacobi<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix);
  
    solver.solve(system_matrix, solution, system_rhs, preconditioner);
  
    problem_out.add_values_to_table(
      tria.n_levels() - offset - 1,
      tria.n_active_cells(),
      tria.n_active_splines(),
      solver_control.last_step(),
      solver_control.tolerance()
    );
  
    std::cout << " ... done!" << std::endl; 
  }
  
  
  void Poisson_Benchmark::compute_h1_error()
  {
    std::cout << "Computing H1-error at Level " 
              << tria.n_levels() - offset - 1 
              << " with degree " 
              << data.max_degree() + order
              << " ... " << std::endl;
  
    std::vector<unsigned int> degrees = tria.get_degree(); 
    degrees[0] = degrees[0] * degrees[0] + 1; 
    degrees[1] = degrees[1] * degrees[1] + 1; 
    TSValues<2> ts_values( 
            &tria,
            degrees,
            update_values |
            update_gradients |
            update_jacobians |
            update_inverse_jacobians |
            update_jacobian_grads |
            update_quadrature_points |
            update_JxW_values);
  
    double h1 = 0;
    double l2 = 0;
    // dofs_per_level(r) = tria.n_active_splines();
    for (const auto& cell : tria.active_cell_iterators()){
      // Get the Bernstein values on the cell
      ts_values.reinit(cell); 
  
      std::vector<unsigned int> local_dof_indices = tria.get_IEN_array(cell);
  
      // Quadrature sum: 
      for (const unsigned int q_index : ts_values.quadrature_point_indices()){
        // Map the quadrature point from real cell to parametric cell
        const Point<2> mapped_q = ts_values.quadrature_point(q_index); 
  
        // Get the value of the approximation
        double u_diff = (-1.) * sol_fcn.value(mapped_q); 
        for (const unsigned int i : ts_values.dof_indices())
          u_diff += (solution(local_dof_indices[i]) *
                      ts_values.shape_value(i, q_index)); 
  
        // Build the gradient value at quadrature point: 
        Tensor<1, 2> grad_u_diff = sol_fcn.gradient(mapped_q); 
        for (const unsigned int i : ts_values.dof_indices())
          grad_u_diff -= ( solution(local_dof_indices[i]) *
                           ts_values.shape_grad(i, q_index));
  
        double dx = ts_values.JxW(q_index); 
        h1 += (( grad_u_diff * grad_u_diff ) * dx) ;
        l2 += ((      u_diff * u_diff      ) * dx); 
      } // for ( q_index )
    } // for ( cell )
    H1 = std::sqrt(l2 + h1);
    L2 = std::sqrt(l2);
    
    problem_out.add_values_to_table(L2, H1);
  
    std::cout << " ... done!" << std::endl; 
  } // compute_h1_error
  
  void Poisson_Benchmark::print_grid(
      const std::string& name
  ) const {
    if (tria.n_levels() - offset - 1 > 16)
      return; 
  
    GridOutFlags::Svg svg_flags;
  //    svg_flags.label_level_number = true;
  //    svg_flags.label_cell_index = true;
  //    svg_flags.label_boundary_id = true;
    svg_flags.boundary_line_thickness = 1;
    svg_flags.line_thickness          = 1;
    svg_flags.coloring = GridOutFlags::Svg::Coloring::none;
  
    std::cout << name << std::endl;
    
    std::ofstream out(name + ".svg");
    GridOut       grid_out;
    grid_out.set_flags(svg_flags);
    
    grid_out.write_svg(tria, out);
  }
  
  void Poisson_Benchmark::compute_cell_errors(){}
  
  
  void Poisson_Benchmark::estimate_and_mark(){
    // special case, if there is only one cell
    if (tria.n_active_cells() == 1) {
      tria.coarsen_bezier_elements(); // switch back to parametric mesh, eventhough nothing really is coarsened
      auto cell = tria.begin_active();
      cell -> set_refine_flag(RefinementCase<2>::cut_axis(cell -> level() % 2 )); 
    } else {
      std::map<TriaIterator<::CellAccessor<2, 2>>, double> local_residuals;
      const std::vector<unsigned int>& degrees = tria.get_degree(); 
      tria.poisson_residual_error_estimate(
                      {(degrees[0] == 1 ? 4 : degrees[0]*degrees[0] + 1), 
                       (degrees[1] == 1 ? 4 : degrees[1]*degrees[1] + 1)}, 
                       &rhs_fcn, 
                       solution,
                       local_residuals
                       );
  
  
  
      std::vector<TriaIterator<::CellAccessor<2, 2>> > mark;
      const auto& bezier = tria.get_bezier_elements();
      double global_residual = 0; 
      for (const auto& [_, r] : local_residuals)
        global_residual += r; 
  
      double weight = 1. / ( tria.n_active_cells() - bezier.size() ); 
      for (const auto& cell : tria.active_cell_iterators()){
        const auto& parent = std::find(bezier.begin(), bezier.end(), cell->parent());
        if (cell -> level() == 0 || bezier.size() == 0){
          if (local_residuals[cell] > 0.75 * weight * global_residual)
            mark.push_back(cell);
        } else if (parent != bezier.end() && *parent == cell->parent()) {
          const auto& parent = cell -> parent(); 
          const double residual = local_residuals[parent -> child(0)] +
                                  local_residuals[parent -> child(1)];
          if (residual > 0.75 * weight * global_residual)
            mark.push_back(parent); 
        } else {
          if (local_residuals[cell] > 0.75 * weight * global_residual)
            mark.push_back(cell);
        }
      }
  
      // Finally prepare the next step, by coarsening the bezier elements
      tria.coarsen_bezier_elements(); 
  
      // Set appropriate refine flags
      tria.set_refine_flags(mark);
    }
    this -> refine_grid();
  }
  
  void Poisson_Benchmark::run(){
  
    tria.refine_bezier_elements();
    tria.compute_extraction_operators();
  
    this -> setup_system(); 
    this -> assemble_system();
    this -> solve();
    this -> compute_h1_error();
    this -> output_system();
    this -> estimate_and_mark();
  
    while (tria.n_levels() - 1 - offset < ref ) {
      this -> setup_system(); 
      this -> assemble_system();
      this -> solve();
      this -> compute_h1_error();
      this -> output_system();
      this -> estimate_and_mark();
    }
  
  
      // Write the resulting table to line
      problem_out.write_table_text(std::cout);
  
  }
  
  IPF_Data<2, 2> Poisson_Benchmark::get_IPF_data()
  {
    std::vector< std::vector< double > > kv; 
    std::vector< Point<2 + 1> > cps; 
    std::vector< unsigned int > deg;
  
    get_IPF_data(kv, cps, deg); 
    IPF_Data<2, 2> out(cps, kv, deg); 
  
    return out; 
  }
  
  void Poisson_Benchmark::get_IPF_data(
                           std::vector< std::vector< double > >& kv,
                           std::vector< Point<2 + 1> >& cps, 
                           std::vector<unsigned int>& deg)
  {
    // Define the vector of control points
    cps = std::vector< Point<3> >(9); 
    const double w = 1.;
    cps[0] = w * Point<3>(0,0,1/w);
    cps[1] = w * Point<3>(0.5,0,1/w);
    cps[2] = w * Point<3>(1.,0,1/w);
    cps[3] = w * Point<3>(0,0.5,1/w);
    cps[4] = w * Point<3>(1.0, 0.5, 1/w);
  //  cps[4] = w * Point<3>(0.5, 0.5,1/w);
    cps[5] = w * Point<3>(1.,0.5,1/w);
    cps[6] = w * Point<3>(0,1.,1/w);
    cps[7] = w * Point<3>(0.5,1.,1/w);
    cps[8] = w * Point<3>(1.,1.,1/w);
  
    // define the knot vectors
    kv = std::vector< std::vector< double > >(2);
    kv[0] = {0, 0, 0, 1, 1, 1}; 
    kv[1] = {0, 0, 0, 1, 1, 1}; 
  
    // Define the degree
    deg = {2, 2}; 
  }
  
  
  double Poisson_RHS::value(
    const Point<2>&     p, 
    const unsigned int /* component */
  ) const {
  //  const double pi = 3.141592653589793;
    const double pi = std::acos(-1); 
    double out = 2. * std::sin( 2 * pi * p(0)) 
                      * std::sin( 2 * pi * p(1)); 
    return out; 
  }
  
  double Poisson_SOL::value(
    const Point<2>&     p, 
    const unsigned int  /* component */
  ) const {
    const double pi = std::acos(-1); 
    double out = 1./(pi * pi * 2. * 2.);
    out *= std::sin(2. * pi * p(0));
    out *= std::sin(2. * pi * p(1));
    return out;
  }
  
  Tensor<1, 2> Poisson_SOL::gradient(
    const Point<2>&     p, 
    const unsigned int /* component */
  ) const {
    Tensor<1, 2> out;
    
    const double pi = std::acos(-1); 
    const double frac = 1. / (2. * pi); 
    out[0] = frac * std::cos(2 * pi * p(0)) * 
                    std::sin(2 * pi * p(1)); 
    out[1] = frac * std::sin(2 * pi * p(0)) * 
                    std::cos(2 * pi * p(1)); 
  
    return out; 
  }
}
