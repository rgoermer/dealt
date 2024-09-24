//First Grid





#include <poisson_neumann_standard.h>

namespace Poisson_Neumann { 
  double Poisson_RHS2::value(
    const Point<2>&     p,
    const unsigned int /* component */
  ) const {
    const double& x = p(0);
    const double& y = p(1);
    double out = 6 * (x + 2. * y * y);
                 
    return -1. * out;
  
  }

  double Poisson_SOL2::value(
    const Point<2>&     p,
    const unsigned int  /* component */
  ) const {
    const double& x = p(0);
    const double& y = p(1);
    double out = x * x * x 
                  + y * y * y * y;
                 
                 
    return out;
  }
  
  Tensor<1, 2> Poisson_SOL2::gradient(
    const Point<2>&     p,
    const unsigned int /* component */
  ) const {
    const double& x = p(0);
    const double& y = p(1);
  
    Tensor<1, 2> out;
    out[0] = 3. * x * x;
    out[1] = 4. * y * y * y;
  
    return out;
  }

  // Define functions for SquishedGeometry class
  Point<2> SquishedGeometry::push_forward(
    const Point<2> &chart_point
  ) const {
  	const double x = chart_point[0];
  	const double y = chart_point[1];
  	Point<2> phi_xy;
  
  	phi_xy[0] = 0.4*x + 0.6*x*x;
  	phi_xy[1] = x * (1. - x) * (1. - 2.*y) + y;
    // phi_xy[1] = y + x * (1 - x) * (1. - 2.*y + y*y - y*y*y);
  
  	return phi_xy;
  }

  DerivativeForm<1, 2, 2> SquishedGeometry::push_forward_gradient (
    const Point<2> &chart_point
  ) const {
  	const double x = chart_point[0];
  	const double y = chart_point[1];
  	Tensor<2, 2> Jacobi;
  	
  	Jacobi[0][0] = 1.2*x + 0.4;
  	Jacobi[0][1] = 0.;
  	Jacobi[1][0] = (1. - 2.*x) * (1. - 2.*y);
  	Jacobi[1][1] = 1. - 2.*x + 2.*x*x;
  	
  	DerivativeForm<1, 2, 2> phi_xy = Jacobi;
  	
  	return phi_xy;
  }
  
  Point<2> SquishedGeometry::pull_back(
    const Point<2> &space_point
  ) const {
  	const double x = space_point[0];
  	const double y = space_point[1];
  	Point<2> xy;
  
  	const double root = std::sqrt(15.*x +1.);
  	const double denominator = (30.*x - 10.*root + 19.);
  
  	xy[0] = (root - 1) / 3.;
  	xy[1] = (9.*y - 5.*(root - 3.*x - 1.)) / denominator;
    // const double temp = 
    //           ( std::pow(
    //                 1.98903 + 60.4938*x + 280.556 * x*x + 50.9259 *x*x*x 
    //                 - 1.98903 * root  - 45.5761 * x * root - 50.9259 * x * x * root 
    //                 + 16.6667 * y + 175. * x * y + 75. * x*x*y - 16.6667 * root * y
    //                 - 50. * x * root * y + 
    //                 std::sqrt(
    //                   -3175.33 * std::pow(
    //                         -0.513333 - 4.04 * x - 1.5 * x * x 
    //                         + 0.513333 * root + x * root
    //                       ,3) 
    //                   + std::pow(
    //                         1.98903 + 50.9259 * x*x*x - 1.98903 * root 
    //                         + (16.6667 - 16.6667*root) * root 
    //                         + x * x * (280.556 - 50.9259 * root + 75. *y)
    //                         + x * (60.4938 - 45.576 * root + ( 175. - 50. * root) * y)
    //                       ,2)
    //                 )
    //               , 1./3.));
  
    //xy[1] = 1./3. + (3.59311 + 10.4993 * x * x - 3.59311 * root + x * (28.2782 - 6.99956 * root)) /
    //          ((-1. -3. * x + root ) * temp) 
    //          - 0.47622 * temp / (-1. -3. * x + root);

  	return xy;
  }
  
  
  std::unique_ptr<Manifold<2, 2>> SquishedGeometry::clone(
  ) const {
  	return std::make_unique<SquishedGeometry>();
  }
  
  
  
  
  
  
  
  
  // dependency on order
  Poisson_Benchmark_Standard::Poisson_Benchmark_Standard(
    const unsigned int& p
  ) : fe(p) 
  	, dof_handler(triangulation)
    , map( )
    , sol_fcn()
    , rhs_fcn()
    , nc_fcn()
  {
    problem_out = OutputSetup("poisson_neumann_standard_no_nc/", fe.degree); 
  }
  
  void Poisson_Benchmark_Standard::make_grid(
  ) {
  	SquishedGeometry geometry;
  	
  	const std::vector<unsigned int> repetition = { 1, 1 };
  	GridGenerator::subdivided_hyper_rectangle(triangulation,
  											  repetition,
  											  Point<2> (0, +0.0),
  											  Point<2> (1, +1.0));
  
  	for (const auto& face : triangulation.active_face_iterators()){
  		if (face -> at_boundary()){
  			if ((face -> center()).operator()(0) == 0 || 
  					(face -> center()).operator()(0) == 1 )
  				face -> set_boundary_id(Boundary::Dirichleth); 
  			else 
  				face -> set_boundary_id(Boundary::Dirichleth);
  		}
  	}
  	
  	// GridTools::transform(
  	//   [&geometry](const Point<2> &chart_point) {
  	//   	return geometry.push_forward(chart_point);
  	//   },
  	//   triangulation
    // );
  		
    std::string name = "out/poisson_neumann_standard/00svg/o"
                          + std::to_string(fe.degree)
                          + "/step0_standard_grid_l0"
                          + ".svg";  
  	std::ofstream out(name);
    
    GridOutFlags::Svg svg_flags;
    svg_flags.coloring = GridOutFlags::Svg::Coloring::none;
    // svg_flags.label_level_number  = true;
    // svg_flags.label_cell_index    = true;
    // svg_flags.label_boundary_id   = true;

  	GridOut 	grid_out;
    grid_out.set_flags(svg_flags);
  	grid_out.write_svg(triangulation, out);
  	
  	// explain the triangulation to use our geometry whenever a new point is needed for refining
    triangulation.set_all_manifold_ids(1);
  	triangulation.set_manifold(1, geometry);
  }
 
  void Poisson_Benchmark_Standard::
    test_face_quadrature_points(
  ) const { 

    QGauss<1> face_quadrature_formula(2 * fe.degree + 1);
  	FEFaceValues<2> fe_face_values_cell(
      fe,
      face_quadrature_formula,
      update_values | update_quadrature_points |
  		update_normal_vectors |
  		update_JxW_values
    );
  	FEFaceValues<2> fe_face_values_neighbor(
      fe,
      face_quadrature_formula,
      update_values | update_quadrature_points |
  		update_normal_vectors |
  		update_JxW_values
    );

    unsigned int f = 0;
    auto cell = triangulation.begin_active();
    auto neighbor = cell -> neighbor( f );

    while ((neighbor.state() == IteratorState::invalid || 
            neighbor.state() == IteratorState::past_the_end)
            && f < 4) {
      f++;
      neighbor = cell -> neighbor(f); 
    }

    Assert(f < 4, ExcInternalError());
    Assert(neighbor.state() == IteratorState::valid, ExcInternalError());
    const unsigned int fb = cell -> neighbor_face_no(f);

    std::cout << " Checking quadrature points on face " << cell -> face(f) -> index()
              << " at " << cell -> face(f) -> vertex(0) << " x " << cell -> face(f) -> vertex(1)
              << std::endl
              << "from cell " << cell -> level() << "." << cell -> index() 
              << " with neighbor " << neighbor -> level() << "." << neighbor -> index() 
              << std::endl;

    fe_face_values_cell.reinit(cell, f); 
    fe_face_values_neighbor.reinit(neighbor, fb); 

    for (const auto& q : fe_face_values_cell.quadrature_point_indices()){
      std::cout << "q = " << q << ": " << std::endl;
      std::cout << "c: " << fe_face_values_cell.quadrature_point(q) << std::endl;
      std::cout << "n: " << fe_face_values_neighbor.quadrature_point(q) << std::endl;
    }
  } // test_face_quadrature_points 
  
  
  
  // step-6: Ordering differs from before because of injection constraints to other functions
  void Poisson_Benchmark_Standard::setup_system(
  ) {
  	// First: dofs
  	dof_handler.distribute_dofs(fe);
  
  	solution.reinit(dof_handler.n_dofs());
  	system_rhs.reinit(dof_handler.n_dofs());
  	
  	// Second: constraints
  	constraints.clear();
  	DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  	
  	VectorTools::interpolate_boundary_values(
      map,
      dof_handler,
  		Boundary::Dirichleth,
  		sol_fcn,
  		constraints
    );
  	constraints.close();
  	
  	// Third: sparse System Matrix
  	DynamicSparsityPattern dsp(dof_handler.n_dofs());
  
  	DoFTools::make_sparsity_pattern(dof_handler,
  									dsp,
  									constraints,
  									/*keep_constrained_dofs = */ false);
  									
  	
  
  	sparsity_pattern.copy_from(dsp);
  	system_matrix.reinit(sparsity_pattern);


    const std::string sp_svg_name = problem_out.svg.string() + "step0_sp_l"
                          + std::to_string(triangulation.n_levels() - 1)
                          + ".svg";
    const std::string sp_name = problem_out.degree.string() + "l"
                          + std::to_string(triangulation.n_levels() - 1)
                          + "_sp.dat";

    std::ofstream out_sp_svg(sp_svg_name), 
                  out_sp(sp_name);
    sparsity_pattern.print_svg(out_sp_svg);
    sparsity_pattern.print(out_sp);
  }
  
  
  
  void Poisson_Benchmark_Standard::assemble_system(
  ) {
  	const QGauss<2> quadrature_formula(2 * fe.degree + 1);
  	const QGauss<1> face_quadrature_formula(2 * fe.degree + 1);
  	
    FEValues<2> fe_values(
      map,
      fe,
      quadrature_formula,
      update_values | update_gradients |
      update_quadrature_points | update_JxW_values
    );
                            
  
  	FEFaceValues<2> fe_face_values(
      map,
      fe,
      face_quadrature_formula,
      update_values | update_quadrature_points |
  		update_normal_vectors |
  		update_JxW_values
    );

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
   
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);
      
  	// contribution with local numbering, need to transfer to global numbering
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    
    
    
    // Calculation for all cells
    for (const auto &cell : dof_handler.active_cell_iterators()) {
  	  // Calculation on an List of dofs first, then insert to the right entries
    	cell_matrix = 0;
    	cell_rhs    = 0;
    
      fe_values.reinit(cell);
  
      // Quadrature
      for (const unsigned int q_index : fe_values.quadrature_point_indices()) {
        // Local contribution to system matrix
        for (const unsigned int i : fe_values.dof_indices()) {
          for (const unsigned int j : fe_values.dof_indices())
  			    cell_matrix(i, j) +=
                   ( fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                      fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                      fe_values.JxW(q_index));           // dx
  
          // Local contribution for RHS
          const auto &x_q = fe_values.quadrature_point(q_index);
          cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q)
                          rhs_fcn.value(x_q) *    	  // f(x_q)			
                          fe_values.JxW(q_index));        // dx
        } // for ( i )
      } // for ( q_index ) 
      
      if (cell -> at_boundary()){ 
  	    for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; f++) {
          const auto& face = cell -> face(f); 
  	    	if (face->at_boundary() && (face->boundary_id() == Boundary::Neumann)) {
  	    		fe_face_values.reinit(cell, face);
  	    		for (const unsigned int& q_point : fe_face_values.quadrature_point_indices()) {
  	    			const double neumann_value = sol_fcn.gradient(fe_face_values.quadrature_point(q_point)) *
  	    			                              fe_face_values.normal_vector(q_point);
   
  	    			for (unsigned int i = 0; i < dofs_per_cell; ++i)
  	    				cell_rhs(i) +=
  	    					(fe_face_values.shape_value(i, q_point) * // phi_i(x_q)
  	    					 neumann_value *                          // g(x_q)
  	    					 fe_face_values.JxW(q_point));            // dx
            } // for ( q_pnt )
          } // if ( neumann_bc )
        }  // for ( f )
      } // if 
      // cell->get_dof_indices(local_dof_indices);
          /*for (unsigned int i = 0; i < dofs_per_cell; ++i)
  	  {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                system_matrix.add(local_dof_indices[i],
                                  local_dof_indices[j],
                                  cell_matrix(i, j));
   
              system_rhs(local_dof_indices[i]) += cell_rhs(i);
  	  }
          //*/
          
          
  	  
  	  // step-6: assembly of the system matrix with constraints
  	  // constraint matrix takes care of BC and hanging nodes
  	  cell->get_dof_indices(local_dof_indices);
  
  	  constraints.distribute_local_to_global(
  	  			cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
  	  
  	
    }	// end loop over cells
    
    
    
    /*
  	constraints.condense(system_matrix);
  	constraints.condense(system_rhs);
  
  	std::map<types::global_dof_index, double> boundary_values;
  	VectorTools::interpolate_boundary_values(dof_handler,
  											 Boundary::Dirichleth,
  											 Functions::ZeroFunction<2>(),
  											 boundary_values);
  											 
  	MatrixTools::apply_boundary_values(boundary_values,
  									   system_matrix,
  									   solution,
  									   system_rhs);
    */
  }
  
  
  
  void Poisson_Benchmark_Standard::solve(
  ) {
  	SolverControl            solver_control(100000, H1 * 1e-4);
    solver_control.enable_history_data();
  	SolverCG<Vector<double>> solver(solver_control);
  	
  	
  	PreconditionSSOR<SparseMatrix<double>> preconditioner;
  	preconditioner.initialize(system_matrix, 1.2);

  	
  	solver.solve(system_matrix, solution, system_rhs, preconditioner);

    problem_out.add_values_to_table(
      triangulation.n_levels() - 1,
      triangulation.n_active_cells(),
      dof_handler.n_dofs(),
      solver_control.last_step(),
      *(solver_control.get_history_data().end()-1)
      // solver_control.tolerance()
    );
  
  	constraints.distribute(solution);                          
  }
  
  
  void Poisson_Benchmark_Standard::refine_grid(
  ) {
    if (triangulation.n_levels() != 1){ 
  	  Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
  	  const Poisson_NC fcn;
  	  const Function<2>* fcn_ptr = &fcn;
  	  std::map<types::boundary_id,  const Function<2>*> neumann_bc = 
          {{Boundary::None,  fcn_ptr}};
  	  KellyErrorEstimator<2>::estimate(
                       map,
                       dof_handler,
  	   								 QGauss<2 - 1>(2 * fe.degree + 1),
  	  								 neumann_bc,
  	  								 solution,
  	  								 estimated_error_per_cell);
  	  
  	  // refine top 30 percent of errors of cells
  	  // for (unsigned int i = 0; i < triangulation.n_active_cells(); i++)
  	  //   std::cout <<  estimated_error_per_cell(i) <<  std::endl;
  	  GridRefinement::refine_and_coarsen_fixed_number(triangulation,
  	  												estimated_error_per_cell,
  	  												0.3,
  	  												0.0);
  
      /* 
      for (auto cell : triangulation.active_cell_iterators())
        cell -> set_refine_flag();
      */
        
    } else { 
      for (auto& cell : triangulation.active_cell_iterators())
        cell -> set_refine_flag();
    }

    SolutionTransfer<2> solution_trans(dof_handler); 
    const Vector<double> old_solution = solution;

    triangulation.prepare_coarsening_and_refinement();
    solution_trans.prepare_for_coarsening_and_refinement(old_solution);

  	triangulation.execute_coarsening_and_refinement();

    setup_system();

    solution_trans.interpolate(old_solution, solution);
    constraints.distribute(solution);
  
  }
  
  
  void Poisson_Benchmark_Standard::print_numerical_solution(
  ) const { 
    const unsigned int N1 = 100;
    const unsigned int N2 = 100; 
    std::vector< Point< 2 > > evals; 

    SquishedGeometry geometry; 

    for (unsigned int j = 0; j < N1; j++)
      for (unsigned int i = 0; i < N2; i++)
        evals.push_back(
          geometry.push_forward(
            Point<2>(
              0. + i / (N1 - 1.),
              0. + j / (N2 - 1.)
            )
          )
        );   

    const TensorProductPolynomials<2>& base(
      Polynomials::generate_complete_Lagrange_basis(
        QGaussLobatto<1>(fe.degree + 1).get_points()
      )
    );
    const unsigned int& n_poly = base.n();

    Vector<double> sol(N1 * N2); 
    // MappingQ<2, 2> map(1); // linear mapping
    std::vector<double> values(n_poly);
    std::vector< Tensor<1, 2> > grads;
    std::vector< Tensor<2, 2> > grad_grads;
    std::vector< Tensor<3, 2> > third_derivatives;
    std::vector< Tensor<4, 2> > fourth_derivatives;
    std::vector< types::global_dof_index > global_dof_indices(n_poly);
    for (const auto& cell : dof_handler.active_cell_iterators()){
      cell -> get_dof_indices(global_dof_indices);
      for (unsigned int n = 0; n < N1* N2; n++){
        if (cell -> point_inside(evals[n])){
          const Point<2>& unit_point = map.transform_real_to_unit_cell(cell, evals[n]);
          base.evaluate(unit_point, values, grads, grad_grads, third_derivatives, fourth_derivatives);        
          for (unsigned int i = 0; i < n_poly; i++)
            sol(n) += values[i] * solution(global_dof_indices[i]);
        } // if ( )
      } // for ( n )
    } // for ( cell )

    // Output: 
    const unsigned int level = triangulation.n_levels() - 1; 
    const std::string name = problem_out.degree.string() 
                              + "l" + std::to_string(level)
                              + "_numerical_solution.dat";

    std::filebuf f;
    f.open(name.c_str(), std::ios::out);
    std::ostream out(&f);

    sol.print(out, 16);
  } // print_numerical_solution()
  
  void Poisson_Benchmark_Standard::print_mesh_file(
  ) const {
    FullMatrix<double> cell_list(triangulation.n_active_cells(), 4);
    FullMatrix<double> vert_list(triangulation.n_used_vertices(), 2);
    unsigned int ind = 0;
    for (const auto& cell : triangulation.active_cell_iterators()){
      for (unsigned int v = 0; v < 4; v++){
        cell_list(ind, v) = cell -> vertex_index(v);
        const Point<2>& vertex = cell -> vertex(v);
        vert_list(cell->vertex_index(v), 0) = vertex(0); 
        vert_list(cell->vertex_index(v), 1) = vertex(1); 
      }
      ind++;
    }

    const std::string mesh = problem_out.degree.string() 
                              + "l"
                              + std::to_string(triangulation.n_levels() - 1)
                              + "_mesh.txt";
    std::ofstream mesh_out(mesh, std::ios::out | std::ios::trunc);
    mesh_out << 4 
             << " " 
             << triangulation.n_active_cells() 
             << " " 
             << triangulation.n_used_vertices() 
             << " " 
             << 0 
             << " " 
             << 2 
             << std::endl;

    for (unsigned int i = 0; i < triangulation.n_active_cells(); i++){
      for (unsigned int d = 0; d < 4; d++)
        mesh_out << cell_list(i, d) << " ";

      mesh_out << std::endl;
    }

    for (unsigned int i = 0; i < triangulation.n_used_vertices(); i++){
      for (unsigned int d = 0; d < 2; d++)
        mesh_out << vert_list(i, d) << " ";
      
      mesh_out << std::endl;
    } 

    mesh_out.close();

  } // print_mesh_file

  void Poisson_Benchmark_Standard::process_results(
    const unsigned int cycle
  ) {
  	//.vtu for visualization
    std::cout << "        Printing to files" << std::endl;
  	Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
  	const Poisson_NC fcn;
  	const Function<2>* fcn_ptr = &fcn;
  	std::map<types::boundary_id,  const Function<2>*> neumann_bc = {{1,  fcn_ptr}};
  	KellyErrorEstimator<2>::estimate(
                     // map,
                     dof_handler,
  	 								 QGauss<2 - 1>(2 * fe.degree + 1),
  									 neumann_bc,
  									 solution,
  									 estimated_error_per_cell);

  	DataOut<2> data_out;
  	data_out.attach_dof_handler(dof_handler);
  	data_out.add_data_vector(solution, "solution");
    data_out.add_data_vector(estimated_error_per_cell, "estimated_error_no_map");

  	KellyErrorEstimator<2>::estimate(
                     map,
                     dof_handler,
  	 								 QGauss<2 - 1>(2 * fe.degree + 1),
  									 neumann_bc,
  									 solution,
  									 estimated_error_per_cell);
    data_out.add_data_vector(estimated_error_per_cell, "estimated_error");

  	data_out.build_patches(map, 0);

    const std::string folder_vtg = problem_out.vtg.string();
    const std::string folder_svg = problem_out.svg.string();
    const std::string level      = std::to_string(triangulation.n_levels() - 1);

    std::string vtu_name = folder_vtg
                          + "/step0_standard_grid_l"
                          + level
                          + ".vtu";
    std::string svg_name = folder_svg
                          + "/step0_standard_grid_l"
                          + level 
                          + ".svg";


    GridOutFlags::Svg svg_flags;
    svg_flags.coloring = GridOutFlags::Svg::Coloring::none;
    // svg_flags.label_level_number  = true;
    // svg_flags.label_cell_index    = true;
    // svg_flags.label_boundary_id   = true;
  
    std::ofstream vtu_out(vtu_name);

  	data_out.write_vtu(vtu_out);
    GridOut       grid_out;
    grid_out.set_flags(svg_flags);
  
    if (triangulation.n_levels()-1 < 11) {
  	  std::ofstream svg_out(svg_name);
      grid_out.write_svg(triangulation, svg_out);
  	  print_numerical_solution();
      print_mesh_file();
    }
  	
    // std::cout << "        Integrating differences per cell..." << std::endl;	
  	// Table for Comparison
   	Vector<float> difference_per_cell(dof_handler.get_triangulation().n_active_cells());
   	VectorTools::integrate_difference(
                      map, 
                      dof_handler,
   									  solution,
   									  sol_fcn,
   									  difference_per_cell,
   									  QGauss<2>(fe.degree + 1),
   									  VectorTools::L2_norm);
   
  
    std::cout << "        Computing global error..." << std::endl;
  	L2 =  VectorTools::compute_global_error(triangulation,
  							difference_per_cell,
  							VectorTools::L2_norm);
  
  	VectorTools::integrate_difference(
      map,
  		dof_handler,
  		solution,
  		sol_fcn,
  		difference_per_cell,
  		QGauss<2>(fe.degree + 1),
  		VectorTools::H1_norm);

    H1 =  VectorTools::compute_global_error(
                              triangulation, 
                              difference_per_cell,
                              VectorTools::H1_norm);
  
  
    problem_out.add_values_to_table(L2, H1);
    
  	// const double H1_error = VectorTools::compute_global_error(triangulation,
  	// 						difference_per_cell,
  	// 						VectorTools::H1_seminorm);
    std::cout << "        ... done" << std::endl;
  
  
  	const unsigned int n_active_cells = triangulation.n_active_cells();
  	const unsigned int n_dofs         = dof_handler.n_dofs();
  
  	// output in console
  	std::cout << "Refinement " << cycle << ':' << std::endl
  	<< "   Number of active cells:       " << n_active_cells << std::endl
  	<< "   Number of degrees of freedom: " << n_dofs << std::endl;
  
  	
  	
  	
  	// output in .text 
    problem_out.write_table_text();
    problem_out.write_table_tex();
  }
  
  
  
  void Poisson_Benchmark_Standard::compute_h1_error(
  ) { 
    const QGauss<2> quadrature_formula(fe.degree + 1);

    FEValues<2> fe_values(
      map,
      fe,
      quadrature_formula,
      update_values | update_gradients |
      update_quadrature_points | update_JxW_values
    );

    double h1 = 0;
    double l2 = 0;
    for (const auto& cell : dof_handler.active_cell_iterators()){
      fe_values.reinit(cell); 
      std::vector<types::global_dof_index> local_dof_indices(fe.n_dofs_per_cell());
      cell -> get_dof_indices(local_dof_indices);
      for (const unsigned int q : fe_values.quadrature_point_indices()) {
        const double dx                = fe_values.JxW(q);
              double u_diff            = sol_fcn.value(fe_values.quadrature_point(q)); 
              Tensor<1, 2> grad_u_diff = sol_fcn.gradient(fe_values.quadrature_point(q));
        for (const unsigned int i : fe_values.dof_indices()){
          u_diff       -= solution(local_dof_indices[i]) * fe_values.shape_value(i, q);
          grad_u_diff  -= solution(local_dof_indices[i]) * fe_values.shape_grad(i, q);
        } // for ( i )
        
        h1 += (( grad_u_diff * grad_u_diff ) * dx);
        l2 += ((      u_diff * u_diff      ) * dx);
      } // for ( q ) 
    } // for ( cell )
    
    H1 = std::sqrt(l2 + h1);
    L2 = std::sqrt(l2);

    problem_out.add_values_to_table(L2, H1);
  }
  
  
  
  
  void Poisson_Benchmark_Standard::run(
    const unsigned int ref
  ) {
  	this -> make_grid();

    std::cout << "    Setting system... " << std::endl;
  	setup_system();

  	//refinement cycle
    unsigned int cycle = 0;
    while (triangulation.n_levels() < ref + 1 && 
                    H1 > 1e-14){
  		std::cout << "Cycle " << cycle << ':' << std::endl;


      std::cout << "    Assembling system ..." << std::endl;
  	  assemble_system();

      std::cout << "    Solving system..." << std::endl;
  	  solve();

      std::cout << "    Processing results..." << std::endl;
  	  process_results(cycle++);

      std::cout << "    Refining grid..." << std::endl;
  	  refine_grid();
  	}

    problem_out.write_table_text(std::cout);
	}
  	//outputting in tables

} // namespace Poisson_Neumann



