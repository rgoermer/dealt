#include <poisson_neumann.h>

namespace Poisson_Neumann {
  Poisson_Benchmark::Poisson_Benchmark(
    int ref, 
    int order,
    const RefinementStrategy strategy) 
  : ref(ref)
  , order(order)
  , data(this->get_IPF_data())
  , tria(data )
  , rhs_fcn()
  , sol_fcn()
  , neumann_bc()
  , strategy(strategy)
  {
    tria.degree_elevate_global(order);
    tria.refine_bezier_elements();
    const auto& cell = tria.begin_active(); 
    const auto& CPs = tria.get_control_points(cell);
    // for (const auto& cp : CPs){
    //   std::cout << cp << std::endl;
    // }
    tria.coarsen_bezier_elements();

    // tria.degree_elevate(0, 1);
    
    std::string name = "poisson_neumann";
    if (strategy == RefinementStrategy::Adaptive)
      name += "_adaptive/";
    else 
      name += "_uniform/";

    problem_out = OutputSetup(name, data.max_degree() + order); 
  
    // Set boundary indicators
    for (auto& face : tria.active_face_iterators()){
    //   const Point<2>& c = face -> center();
    //   if (face -> at_boundary() &&
    //         (c(1) == 0 || c(1) == 1 ) )
        face -> set_boundary_id(Boundary::Dirichlet);
    }
  
    tria.set_boundary_dofs();
  
    tria.refine_bezier_elements();
    tria.compute_extraction_operators();

  }
  
  
  void Poisson_Benchmark::setup_system(){
    system_rhs.reinit(tria.n_active_splines());
    solution.reinit(tria.n_active_splines());
  
    const auto& IEN_array = tria.get_IEN_array();
    Assert(IEN_array.size() > 0, ExcInternalError());
  
    // The TSplines are our finite elements. However, by bezier decomposition, these
    // are decomposed into a linear combination of Bernstein Polynomials over each cell.
    // Thus, we use the FE_Bernstein<1> class, that represents Bernstein polynomials
    // over [-1, 1]. We can not use FE_Berstein<dim>, because this represents a Bernstein
    // basis with polynomial degree p in every direction. In some cases, this is fine,
    // however, we wish to stress the fact, that TSplines may differ in its polynomial degree
    // along two directions
  
    // To generate the raw structure of the sparsity matrix, we construct our own SparsityPattern
    // This will be done using the IEN_array from the TS_Triangulation. It already stores the
    // information which spline has support on which spline.
    //
    // Currently, the maximum number of entries per row is unknown, hence we initialize it
    // as a FullMatrix and compress it afterwards.
    sparsity_pattern.reinit(
        tria.n_active_splines(),
        tria.n_active_splines(),
        tria.n_active_splines() );
  
    for (const auto& [_, arr] : IEN_array)
      for (unsigned int i : arr)
        for (unsigned int j : arr)
          sparsity_pattern.add(i, j);
  
    // Free superfluous space
    sparsity_pattern.compress();
    system_matrix.reinit(sparsity_pattern);
  
  }
  
  void Poisson_Benchmark::assemble_system(){
    if (tria.n_levels() - 1 < 16) {
      std::string name = problem_out.svg.string() + "step0_bezier_grid_l"
                            + std::to_string(tria.n_levels() - 1)
                            + ".svg";  

      GridOutFlags::Svg svg_flags;
      svg_flags.coloring = GridOutFlags::Svg::Coloring::none;
      // svg_flags.label_level_number  = true;
      // svg_flags.label_cell_index    = true;
      // svg_flags.label_boundary_id   = true;

      std::ofstream out(name);
      GridOut       grid_out;
      grid_out.set_flags(svg_flags);
  
      grid_out.write_svg(tria, out);
    }
  
    std::cout << "Assembling system matrix ... " << std::endl;
  
    // Setup initial tables that store the bernstein values / grads / hessians.
  
    std::vector< unsigned int > degrees = tria.get_degree();
    degrees[0] = degrees[0]  + 1;
    degrees[1] = degrees[1]  + 1;

    TSValues<2> ts_values(
        &tria,
        degrees,
        update_values |
        update_gradients |
        update_quadrature_points |
        update_JxW_values);

    TSFaceValues<2> face_values(
        &tria,
        degrees,
        update_values |
        update_gradients |
        update_quadrature_points |
        update_normal_vectors |
        update_JxW_values);

    // const IPF_Inverse inv_phi; 
    // const IPF_Mapping phi;
    // ts_values.test_geometry_mapping(&phi, &inv_phi, 1e-10, true);
    // face_values.test_geometry_mapping(&phi, &inv_phi, 1e-10, false);
  
    const unsigned int dofs_per_cell = ts_values.n_dofs_per_cell();
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);
  

    for (const auto& cell : tria.active_cell_iterators()){
      // Get the Bernstein values on the cell
      ts_values.reinit(cell);
  
      // Reset the cell matrix
      cell_matrix       = 0;
      cell_rhs          = 0;
  
      // Quadrature sum:
      for (const unsigned int q : ts_values.quadrature_point_indices()){
  
        // Build the cell matrix:
        double dx = ts_values.JxW(q);
        for (const unsigned int i : ts_values.dof_indices())
          for(const unsigned int j : ts_values.dof_indices())
            cell_matrix(i,j) += ( ts_values.shape_grad(i, q)
                                  * ts_values.shape_grad(j, q)
                                  * dx );
  
        // Map the quadrature point from real cell to parametric cell
        const Point<2>& mapped_q = ts_values.quadrature_point(q);
        const double rhs = rhs_fcn.value(mapped_q);
        for (const unsigned int i : ts_values.dof_indices())
          cell_rhs(i) += ( ts_values.shape_value(i, q)
                            * rhs
                            * dx );
  
  
      } // for ( q )
  
      // Check for neumann conditions
      if (cell -> at_boundary()){
        for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; f++){
          if (cell -> face(f) -> at_boundary() &&
                cell -> face(f) -> boundary_id() == Boundary::Neumann){
            face_values.reinit(cell, f);
            for (const unsigned int q : face_values.quadrature_point_indices()){
              const long double g = sol_fcn.gradient(face_values.quadrature_point(q)) *
                                      face_values.normal_vector(q);
              for (const unsigned int i : face_values.dof_indices())
                cell_rhs(i) += g * face_values.shape_value(i, q)
                                 * face_values.JxW(q);
            }
          } // if ( neumann_bc )
        } // for ( f )
      } // if
  
      // Add the values to the system
      std::vector< unsigned int > local_dof_indices =
              tria.get_IEN_array(cell);
  
      system_matrix.add(local_dof_indices, cell_matrix);
      system_rhs.add(local_dof_indices, cell_rhs);
  
    } // for ( cell )
    

    std::map<
        types::global_dof_index, 
        double
      > boundary_values;
    std::map<
        types::boundary_id,
        const Function< 2 >*
      > boundary_fcns =
          {{Boundary::Dirichlet, &sol_fcn}};
    tria.project_boundary_values(
      boundary_fcns,
      degrees,
      boundary_values
    );

    std::cout << "Boundary dofs by face: " << std::endl;
    const auto& boundary_dofs = tria.get_boundary_dofs();
    const auto& splines = tria.get_splines(); 
    for (const auto& dof : boundary_dofs.at(Boundary::Dirichlet)){
      const auto& ts = splines.at(dof); 
      const auto& anchor = ts -> get_anchor();
      std::cout << dof << ": " << 0.5 * anchor.first + 0.5*anchor.second << ", val = " << boundary_values.at(dof) << std::endl;
    } // for ( dof )

    MatrixTools::apply_boundary_values(
      boundary_values, 
      system_matrix,
      solution,
      system_rhs
    );
  
  
    std::cout << " ... done!" << std::endl;
  }
  
  void Poisson_Benchmark::output_system()
  {
    std::cout << "Printing system matrix and rhs ... " << std::endl;
  
    const unsigned int level = tria.n_levels() - 1;
    const std::string name = problem_out.degree.string();
    const std::string level_name = name + "l" + std::to_string(level);
  
    std::string matrix = level_name + "_mat.dat" ;
    std::string vector = level_name + "_vec.dat" ;
    std::string soluti = level_name + "_sol.dat" ;

    // Find some interesting spline index
    if (tria.n_levels() - 1 == 9) {
      const double target_x = 0.625;
      const double target_y_min = 0.875;
      const double target_y_max = 1.;

      const auto& splines = tria.get_splines(); 

      int index = -1;
      for (const auto& t : splines) {
        const auto& A0 = t -> get_anchor(0);
        const auto& A1 = t -> get_anchor(1);
        
        if (A0.first == target_x && 
              A1.first == target_y_min && 
              A1.second == target_y_max)
          index = t -> get_level_index();
      }

      std::cout << "index = " << index << std::endl;
    }
  
    if (tria.n_levels() - 1 < 11) {
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
    Vector< double > cell_errors(tria.n_active_cells());
    std::map<types::boundary_id,
             const Function<2> *>
        neumann_data = {{Boundary::Neumann, &neumann_bc}};
    std::vector<unsigned int> degrees = tria.get_degree();
    degrees[0] = degrees[0] * degrees[0] + 1;
    degrees[1] = degrees[1] * degrees[1] + 1;
   
    ResidualEstimators::Poisson<2>::estimate(
      &tria, 
      degrees,
      solution,
      cell_errors, 
      &rhs_fcn,
      neumann_data
    );

    //for (unsigned int n = 0; n < tria.n_active_cells(); n++)
    //  Assert(std::fabs(new_cell_errors(n) - cell_errors(n)) < 1e-15, 
    //          ExcInternalError());

    data_out.add_data_vector(cell_errors, "cell_errors");
    //data_out.add_data_vector(new_cell_errors, "new_cell_errors");

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
  } // output_system
  
  void Poisson_Benchmark::solve(){
    std::cout << "Solving system ... " << std::endl;

    SolverControl            solver_control(750 * tria.n_active_splines(), H1 * 1e-4);
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
      solver_control.tolerance()
    );
  
    std::cout << " ... done!" << std::endl;
  }
  
  
  void Poisson_Benchmark::compute_h1_error()
  {
    std::cout << "Computing H1-error at Level "
              << tria.n_levels() - 1
              << " with degree "
              << data.max_degree() + order
              << " ... " << std::endl;
  
    std::vector< unsigned int > degrees = tria.get_degree();
    degrees[0] = (degrees[0] * degrees[0]) + 1;
    degrees[1] = (degrees[1] * degrees[1]) + 1;

    TSValues<2> ts_values(
        &tria,
        degrees,
        update_values |
        update_gradients |
        update_quadrature_points |
        update_hessians |
        update_JxW_values);
  
  
    double h1 = 0;
    double l2 = 0;
  
    for (const auto& cell : tria.active_cell_iterators()){
      // Get the Bernstein values on the cell
      ts_values.reinit(cell);
  
      std::vector<unsigned int> local_dof_indices = tria.get_IEN_array(cell);
  
      // Quadrature sum:
      for (const unsigned int q_index : ts_values.quadrature_point_indices()){
        // Map the quadrature point from real cell to parametric cell
        const Point<2>& mapped_q = ts_values.quadrature_point(q_index);
  
        // Get the value of the approximation
        double u_diff = sol_fcn.value(mapped_q);
        for (const unsigned int i : ts_values.dof_indices())
          u_diff -= (solution(local_dof_indices[i]) *
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
  
  
  
  void Poisson_Benchmark::estimate_and_mark(){
    // declare a container to store cells to be marked
    std::vector< TriaIterator<CellAccessor<2, 2>> > mark;
  
    if (strategy == RefinementStrategy::Adaptive) {
      std::map< types::boundary_id,
                const Function<2>* >       
          neumann_data = {{1, &neumann_bc}};
      std::vector<unsigned int> degrees = tria.get_degree();
      degrees[0] = degrees[0] * degrees[0] + 1;
      degrees[1] = degrees[1] * degrees[1] + 1;
      Vector<double> local_residuals(tria.n_active_cells());

      ResidualEstimators::Poisson<2>::estimate(
          &tria,
          degrees,
          solution,
          local_residuals,
          &rhs_fcn,
          neumann_data
      );
       
      tria.refine_fixed_number(local_residuals, 0.10);
    } else { 
      tria.coarsen_bezier_elements();
      tria.refine_global();   
    } // if ( strategy )
  
    if (tria.n_levels() < 15) {
      GridOutFlags::Svg svg_flags;
      svg_flags.coloring = GridOutFlags::Svg::Coloring::none;
      // svg_flags.label_level_number  = true;
      // svg_flags.label_cell_index    = true;
      // svg_flags.label_boundary_id   = true;
  
      std::string name = problem_out.svg.string() + "step0_grid_l"
                            + std::to_string(tria.n_levels() - 1)
                            + ".svg";  
      std::ofstream out(name);
      GridOut       grid_out;
      grid_out.set_flags(svg_flags);
  
      grid_out.write_svg(tria, out);
    }
  
    tria.prepare_assembly(); 
  }
  
  void Poisson_Benchmark::run(){
    unsigned int level_count = 0; 
    cycle = 0;
    while (tria.n_levels() - 1 < ref + 1 && 
            H1 > 1e-13){
      cycle++;
      this -> setup_system();
      this -> assemble_system();
      this -> solve();
      this -> print_error();
      this -> compute_h1_error();
      this -> output_system();
      if (level_count < 5) {
        this -> estimate_and_mark();
        level_count++;
      } else {
        tria.coarsen_bezier_elements();
        tria.refine_global();
        tria.set_boundary_dofs();
        tria.refine_bezier_elements();
        tria.compute_extraction_operators();
        level_count = 0; 
      }
    }

    // Write the resulting table to line
    problem_out.write_table_text(std::cout);
  }
  
  void Poisson_Benchmark::print_error(
  ) {
    std::vector< unsigned int > degrees = tria.get_degree();
    degrees[0] = (degrees[0] * degrees[0]) + 1;
    degrees[1] = (degrees[1] * degrees[1]) + 1;

    TSValues<2> ts_values(
        &tria,
        degrees,
        update_values |
        update_gradients |
        update_quadrature_points |
        update_hessians |
        update_JxW_values);
    const unsigned int nq = ts_values.n_quadrature_points_per_cell();
    const unsigned int nc = tria.n_active_cells();
    FullMatrix<double> B(nq, nc);
    FullMatrix<double> evals(nq*nc, 2);
    unsigned int c = 0, ind = 0;
    for (const auto& cell : tria.active_cell_iterators()){
      // Get the TSpline values on this cell
      ts_values.reinit(cell);
  
      // Get the corresponding global dof indices
      std::vector<unsigned int> local_dof_indices = tria.get_IEN_array(cell);
  
      // For each quadrature point ...
      for (const unsigned int q : ts_values.quadrature_point_indices()){
        // ... compute the value of the solution at the
        // quadrature point and ...
        double diff = -1. * sol_fcn.value(ts_values.quadrature_point(q));
        // double diff = 0.;
  
        // ... build the difference with the numeric solution, to ...
        for (const unsigned int i : ts_values.dof_indices())
          diff += ts_values.shape_value(i, q) * solution(local_dof_indices[i]);
  
        // ... store it in the output matrix:
        B(q, c) = std::fabs(diff);
        // B(q, c) = diff;
        evals(ind, 0) = ts_values.quadrature_point(q)(0);
        evals(ind, 1) = ts_values.quadrature_point(q)(1);
        ind++;
      } // for ( q )
      c++;
    } // for ( cell )
  
    // Finally, print the matrix to a readable .dat file:
  
    // Firstly declare a name in the most complicated but retraceable way:
    const std::string l = std::to_string(tria.n_levels() - 1);
    const std::string name = problem_out.degree.string();
    const std::string error_out_name = name + "l" + l + "_cell_error.dat";
    const std::string evals_out_name = name + "l" + l + "_cell_error_evals.dat";

    Vector<double> data(4);
    data(0) = tria.n_active_cells();
    data(1) = degrees[0];
    data(2) = degrees[1];
    data(3) = tria.n_active_splines();

    std::string data_name = name + "l" + l + "_data.dat";
    std::filebuf dat;
    dat.open(data_name.c_str(), std::ios::out);

    std::ostream data_out(&dat);

    data.print(data_out, 5);
    dat.close();
  
    // Then open a file to write contents into:
    std::filebuf out_file0, out_file1;
    out_file0.open(error_out_name.c_str(), std::ios::out);
    out_file1.open(evals_out_name.c_str(), std::ios::out);
  
    // Then tell C++ to write in that file:
    std::ostream out0(&out_file0);
    std::ostream out1(&out_file1);
  
    // And finally print the contents to that file
    B.print_formatted(out0, 16, true, 1, "0");
    evals.print_formatted(out1, 16, true, 1, "0");
  
    // Free the output file for the system
    out_file0.close();
    out_file1.close();
  }
  
  IPF_Data<2, 2> Poisson_Benchmark::get_IPF_data()
  {
    std::vector< std::vector< double > > kv;
    std::vector< Point<2 + 1> > cps;
    std::vector< unsigned int > deg;
  
    // Define the vector of control points
    cps = std::vector< Point<3> >(9);
    const double w = 1.;
    cps[0] = w * Point<3>(+0.0, +0.0, 1./w);
    // cps[1] = w * Point<3>(+0.2, -0.5, 1./w);
    cps[1] = w * Point<3>(+0.2, +0.5, 1./w);
    // cps[1] = w * Point<3>(+0.5, +0.0, 1./w);
    cps[2] = w * Point<3>(+1.0, +0.0, 1./w);
    cps[3] = w * Point<3>(+0.0, +0.5, 1./w);
    cps[4] = w * Point<3>(+0.2, +0.5, 1./w);
    // cps[4] = w * Point<3>(+0.0, +0.5, 1./w);
    cps[5] = w * Point<3>(+1.0, +0.5, 1./w);
    cps[6] = w * Point<3>(+0.0, +1.0, 1./w);
    // cps[7] = w * Point<3>(+0.2, +1.5, 1./w);
    cps[7] = w * Point<3>(+0.2, +0.5, 1./w);
    // cps[7] = w * Point<3>(+0.5, +1.0, 1./w);
    cps[8] = w * Point<3>(+1.0, +1.0, 1./w);
  
    // define the knot vectors
    kv = std::vector< std::vector< double > >(2);
    kv[0] = {0, 0, 0, 1, 1, 1};
    kv[1] = {0, 0, 0, 1, 1, 1};
  
    // Define the degree
    deg = {2, 2};
  
    IPF_Data<2, 2> out(cps, kv, deg);
  
    return out;
  }
  
  
  double Poisson_RHS::value(
    const Point<2>&     p,
    const unsigned int /* component */
  ) const {
    const double pi = numbers::PI;
    double out = std::sin( 2 * pi * p(0)) *
                 std::cos(     pi * p(1));
    return out;
  }
  
  double Poisson_SOL::value(
    const Point<2>&     p,
    const unsigned int  /* component */
  ) const {
    const double pi = std::acos(-1);
    double out = 1./(pi * pi * 5.) *
                  std::sin(2. * pi * p(0)) *
                  std::cos(     pi * p(1));
    return out;
  }
  
  Tensor<1, 2> Poisson_SOL::gradient(
    const Point<2>&     p,
    const unsigned int /* component */
  ) const {
    Tensor<1, 2> out;
  
    const double pi = numbers::PI;
    const double frac = 1. / (5. * pi);
    out[0] = (+2.) * frac * std::cos(2 * pi * p(0)) *
                            std::cos(    pi * p(1));
    out[1] = (-1.) * frac * std::sin(2 * pi * p(0)) *
                            std::sin(    pi * p(1));
  
    return out;
  }
  
  double Poisson_NC::value(
      const Point<2>&     p,
      const unsigned int  /* component */
  ) const {
    const Poisson_SOL sol_fcn;
    const IPF_Mapping geometry; 
    const Tensor<1, 2> df = sol_fcn.gradient(p);

    Tensor<1, 2> n; 
    const int s = (std::fabs(p(0) - 1./3.) < 1e-15 ? -1 : +1);
    if (p(1) < 0.5) {
      n[0] = +1. * s * (1. - 2. * p(0)); 
      n[1] = -1. * s * (0.4 + 1.2 * p(0));
    } else {
      n[0] = +1. * s * (1. - 2. * p(0)); 
      n[1] = +1. * s * (0.4 + 1.2 * p(0));
    }
  
    return df * (n / n.norm());
  }

  double IPF_Mapping::value(
    const Point<2>& p,
    const unsigned int component
  ) const {
    if (component == 0)
      return (0.4 * p(0) + 0.6 * p(0)*p(0)); 
    else if (component == 1)
      return (p(1) + p(0) * (1. - p(0)) * (1. - 2. * p(1)));
    else 
      return nan("1");
  }

  Tensor<1, 2> IPF_Mapping::gradient(
    const Point<2>& p,
    const unsigned int component
  ) const {
    Tensor<1, 2> out;
    if (component == 0){
      out[0] = 0.4 + 1.2 * p(0);
    } else if (component == 1) {
      out[0] = (1. - 2.*p(1)) * (1. - 2.*p(0));
      out[1] = 1. - (2.*p(0) * (1. - p(0))); 
    } else {
      out[0] = nan("1");
      out[1] = nan("1");
    }

    return out;
  }

  SymmetricTensor<2, 2> IPF_Mapping::hessian(
    const Point<2>& p,
    const unsigned int component
  ) const {
    SymmetricTensor<2, 2> out; 
    if (component == 0) {
      out[0][0] = 1.2;
    } else if (component == 1) {
      out[0][0] = -2. * (1. - 2.*p(1));
      out[0][1] = -2. * (1. - 2.*p(0));
      out[1][1] = 0.;
    }
    return out;
  }

  double IPF_Inverse::value(
    const Point<2>& p,
    const unsigned int component
  ) const {
    const double root = std::sqrt(1. + 15. * p(0));
    double out;
    if (component == 0) {
      out = (root-1.) / 3.;
    } else {
      const double D = (30. * p(0) - 10. * root + 19.);
      out = (9. * p(1) - 5. * root + 15. * p(0) + 5.) / D;
    }
    return out;
  }

  Tensor<1, 2> IPF_Inverse::gradient(
    const Point<2>& p,
    const unsigned int component
  ) const {
    Tensor<1, 2> out; 
    const double root = std::sqrt(1. + 15. * p(0));
    if (component == 0) {
      out[0] = 5. / (2. * root);  
    } else if (component == 1) {
      const double D9 = (30. * p(0) - 10. * root + 19.);
      const double D = D9 / 9.;
      out[0] = -1. * 135. * (2. * root - 5.) * (2. * p(1) - 1.) / (2. * D9 * D9 * root);  
      out[1] = 1. / D;
    }
    return out;
  }

  SymmetricTensor<2, 2> IPF_Inverse::hessian(
    const Point<2>& p,
    const unsigned int component
  ) const {
    SymmetricTensor<2, 2> out; 
    const double root = std::sqrt(1. + 15. * p(0));
    const double root3 = std::pow(1. + 15. * p(0), 1.5);
    if (component == 0) {
      out[0][0] = -75. / (4. * root3);
    } else if (component == 1) {
      const double D9 = (30. * p(0) - 10. * root + 19.);
      const double D93 = D9 * D9 * D9;
      out[0][0] = (2025. * (16. * root3 + 150. * root - 1350. * p(0) -175.) * (2.*p(1) -1.)) 
                  / (4. * root3 * D93);
      out[0][1] = 90. * (-3. + 15./(2. * root))
                  / (D9 * D9);
      out[1][1] = 0.;
    }

    return out;

  }
}
