
#include <lshape_standard.h>

namespace LShape { 
  
  void LShapePostprocessor::evaluate_scalar_field(
    const DataPostprocessorInputs::Scalar<2>  &inputs, 
          std::vector< Vector<double> >       &outputs
  ) const {
    AssertDimension(outputs.size(), inputs.solution_values.size());
    const unsigned int &n_outputs = outputs.size(); 
    for (unsigned int p = 0; p < n_outputs; p++){
      const auto& eval = inputs.evaluation_points[p];
      outputs[p](0) = 
        std::fabs(inputs.solution_values[p] - sol_fcn.value(eval));
      outputs[p](1) = 
        (inputs.solution_gradients[p] - sol_fcn.gradient(eval)).norm_square();
    }


    //for (unsigned int p = 0; p < n_outputs; p++){
    //  printf("Input:  %+1.4e \nOutput: %+1.4e\n", inputs.solution_values[p], outputs[p](0));
    //  printf("sol_fcn: %+1.4e\n", sol_fcn.value(inputs.evaluation_points[p]));
    //  printf("eval:   [%+1.4e, %+1.4e]\n", inputs.evaluation_points[p](0), inputs.evaluation_points[p](1));
    //}
  }

  std::vector< std::string > LShapePostprocessor::get_names(
  ) const {
    std::vector< std::string > out;
    out.emplace_back("pointwise_difference");
    out.emplace_back("pw_diff_grads");
    return out;
  }

  UpdateFlags LShapePostprocessor::get_needed_update_flags(
  ) const {
    return update_values | update_gradients | update_quadrature_points;
  }

  // dependency on order
  LShape_Benchmark_Standard::LShape_Benchmark_Standard(
    const unsigned int order   , 
    const Specifier    type    , 
    const Strategy     strategy 
  ) : type(type)
    , strategy(strategy)
    , fe(order)
    , dof_handler(triangulation)
    , rhs_fcn(type)
    , sol_fcn(type)
    , neumann_bc1(type)
    , neumann_bc2(type)
    , neumann_bc3(type)
    , neumann_bc4(type)
  {
    const FiniteElementData<2>::Conformity conf = fe.conforming_space;
    switch (conf){
      case (FiniteElementData<2>::Conformity::unknown): 
        std::cout << "conformity: unknown" << std::endl; 
        break;
      case (FiniteElementData<2>::Conformity::L2): 
        std::cout << "conformity: L2" << std::endl; 
        break;
      case (FiniteElementData<2>::Conformity::Hcurl): 
        std::cout << "conformity: Hcurl" << std::endl; 
        break;
      case (FiniteElementData<2>::Conformity::Hdiv): 
        std::cout << "conformity: Hdiv" << std::endl; 
        break;
      case (FiniteElementData<2>::Conformity::H1): 
        std::cout << "conformity: H1" << std::endl; 
        break;
      case (FiniteElementData<2>::Conformity::H2):
        std::cout << "conformity: H2" << std::endl; 
        break;
      default:
        Assert(false, ExcInternalError());
        break;
    }
    Assert(fe.conforms(FiniteElementData<2>::Conformity::H1), ExcInternalError());

    // Setup filesystem for currrent run:
    const std::string t_str = (type == Specifier::classic ? "c" : "b"); 
    const std::string s_str = (strategy == Strategy::adaptive ? "a" : "u"); 
    const std::string out_init = "lshape_standard_no_nc/lshape_" + s_str + t_str + "/";

    // create filesystem for output
    problem_out = OutputSetup(out_init, fe.degree);
  }
  
  void LShape_Benchmark_Standard::make_grid(
  ) {
    const std::vector<unsigned int> repetition = { 2, 2 };
    const std::vector<int> remove = { 1, 1 };
  
    GridGenerator::subdivided_hyper_L(
      triangulation,
      repetition,
      Point<2> (-1, -1),
      Point<2> (1, 1), 
      remove
    );
  


    for (const auto& face : triangulation.active_face_iterators()){
      if (face -> at_boundary()){
        face -> set_boundary_id(Boundary::dirichleth);
        // const auto& c = face -> center();
        // if ( std::fabs(c(0) + 1.) < 1e-15 )
        //   face -> set_boundary_id(Boundary::neumann_1);
        // else if ( std::fabs(c(1) - 1.) < 1e-15 )
        //   face -> set_boundary_id(Boundary::neumann_2);
        // else if ( std::fabs(c(0) - 1.) < 1e-15 )
        //   face -> set_boundary_id(Boundary::neumann_3); 
        // else if ( std::fabs(c(1) + 1.) < 1e-15 )
        //   face -> set_boundary_id(Boundary::neumann_4);
        // else if (std::fabs(c(0) + 0.) < 1e-15)
        //   face ->set_boundary_id(Boundary::dirichleth);
        // else if (std::fabs(c(1) + 0.) < 1e-15)
        //   face ->set_boundary_id(Boundary::dirichleth);
        // else 
        //   face -> set_boundary_id(Boundary::none);
      }
    }
      
    std::string name = problem_out.svg.string() + "grid_l0.svg";
    std::ofstream out(name);
    
    GridOutFlags::Svg svg_flags;
    svg_flags.coloring = GridOutFlags::Svg::Coloring::none;
    // svg_flags.label_level_number  = true;
    // svg_flags.label_cell_index    = true;
    // svg_flags.label_boundary_id   = true;

    GridOut   grid_out;
    grid_out.set_flags(svg_flags);
    grid_out.write_svg(triangulation, out);
    
  }
  
  
  
  
  // step-6: Ordering differs from before because of injection constraints to other functions
  void LShape_Benchmark_Standard::setup_system(
  ) {
    // First: dofs
    dof_handler.distribute_dofs(fe);
  
    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
    
    // Second: constraints
    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    
    VectorTools::interpolate_boundary_values(
      dof_handler,
      Boundary::dirichleth,
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
  }
  
  
  
  void LShape_Benchmark_Standard::assemble_system(
  ) {
    const QGauss<2> quadrature_formula(fe.degree + 1);
    const QGauss<1> face_quadrature_formula(fe.degree + 1);
    
    FEValues<2> fe_values(
      fe,
      quadrature_formula,
      update_values | update_gradients |
      update_quadrature_points | update_JxW_values
    );
                            
  
    FEFaceValues<2> fe_face_values(
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
                   ( fe_values.shape_grad(i, q_index) *   // grad phi_i(x_q)
                      fe_values.shape_grad(j, q_index) *  // grad phi_j(x_q)
                      fe_values.JxW(q_index));            // dx
  
          // Local contribution for RHS
          const auto &x_q = fe_values.quadrature_point(q_index);
          cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q)
                          rhs_fcn.value(x_q) *                // f(x_q)     
                          fe_values.JxW(q_index));            // dx
        } // for ( i )
      } // for ( q_index ) 
    
    for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; f++) 
      {
      const auto& face = cell -> face(f); 
        if (face->at_boundary() && (face->boundary_id() == Boundary::neumann_1 ||
                                    face->boundary_id() == Boundary::neumann_2 ||
                                    face->boundary_id() == Boundary::neumann_3 ||
                                    face->boundary_id() == Boundary::neumann_4)) {
          // std::cout << "At cell " << cell->level() << "." << cell -> index() << " and face " << f << std::endl;
          fe_face_values.reinit(cell, face);
   
          // std::cout << "Normals given by: " << std::endl;
          for (const unsigned int& q_point : fe_face_values.quadrature_point_indices()) {
            const double neumann_value = sol_fcn.gradient(fe_face_values.quadrature_point(q_point)) *
                                          fe_face_values.normal_vector(q_point);
            // const auto& normal = fe_face_values.normal_vector(q_point);
            // const double neumann_value = g.value(fe_face_values.quadrature_point(q_point));
            // printf("%2u: [%+1.8e %+1.8e]\n", q_point, normal[0], normal[1]);
   
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              cell_rhs(i) +=
                (fe_face_values.shape_value(i, q_point) * // phi_i(x_q)
                 neumann_value *                          // g(x_q)
                 fe_face_values.JxW(q_point));            // dx
          }
        }
      } 
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
      
    
    } // end loop over cells
    
  }
  
  
  
  void LShape_Benchmark_Standard::solve(
  ) {
    SolverControl            solver_control(100000, H1 * 1e-4);
    SolverCG<Vector<double>> solver(solver_control);
    
    
    PreconditionJacobi<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    
    solver.solve(system_matrix, solution, system_rhs, preconditioner);

    problem_out.add_values_to_table(
      triangulation.n_levels() - 1,
      triangulation.n_active_cells(),
      dof_handler.n_dofs(),
      solver_control.last_step(),
      solver_control.tolerance()
    );
  
    constraints.distribute(solution);                          
  }
  
  
  void LShape_Benchmark_Standard::refine_grid(
  ) {
    if (strategy == Strategy::adaptive  
            /* && (triangulation.n_levels() - 1) % 4 != 0 */ ){
      Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
      
      std::map<types::boundary_id,  const Function<2>*> neumann_bc =
        {{7,  &neumann_bc1}, 
         {8,  &neumann_bc2}, 
         {9,  &neumann_bc3},
         {10,  &neumann_bc4}
        };

      KellyErrorEstimator<2>::estimate(dof_handler,
                       QGauss<2 - 1>(fe.degree * fe.degree + 1),
                       neumann_bc,
                       solution,
                       estimated_error_per_cell, 
                       {},
                       nullptr,
                       numbers::invalid_unsigned_int,
                       numbers::invalid_subdomain_id,
                       numbers::invalid_material_id,
                       KellyErrorEstimator<2>::Strategy::face_diameter_over_twice_max_degree);
      
      // refine top 30 percent of errors of cells
      // for (unsigned int i = 0; i < triangulation.n_active_cells(); i++)
      //   std::cout <<  estimated_error_per_cell(i) <<  std::endl;
      GridRefinement::refine_and_coarsen_fixed_fraction(triangulation,
                              estimated_error_per_cell,
                              1./3.,
                              0.0);
    } else {
      // perform global refinement
      for (auto cell : triangulation.active_cell_iterators())
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
  
  
  
  void LShape_Benchmark_Standard::process_results(
  ) {
    DataOut<2> data_out;
    LShapePostprocessor postprocessor(type);
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "solution");
    data_out.add_data_vector(solution, postprocessor);
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
    
    std::map<types::boundary_id,  const Function<2>*> neumann_bc =
      {{Boundary::neumann_1,  &neumann_bc1}, 
       {Boundary::neumann_2,  &neumann_bc2}, 
       {Boundary::neumann_3,  &neumann_bc3},
       {Boundary::neumann_4,  &neumann_bc4}
      };

    KellyErrorEstimator<2>::estimate(dof_handler,
                     QGauss<1>(fe.degree * fe.degree + 1),
                     neumann_bc,
                     solution,
                     estimated_error_per_cell, 
                     {},
                     nullptr,
                     numbers::invalid_unsigned_int,
                     numbers::invalid_subdomain_id,
                     numbers::invalid_material_id,
                     KellyErrorEstimator<2>::Strategy::face_diameter_over_twice_max_degree);
    
    data_out.add_data_vector(estimated_error_per_cell, "estimated_error");

    std::cout << "        Integrating differences per cell..." << std::endl;  
    // Table for Comparison
    Vector<float> difference_per_cell(triangulation.n_active_cells());
    VectorTools::integrate_difference(
                      dof_handler,
                      solution,
                      sol_fcn,
                      difference_per_cell,
                      QGauss<2>(fe.degree + 1),
                      VectorTools::L2_norm);
  
  
    std::cout << "        Computing global error..." << std::endl;
    // const double L2_error = VectorTools::compute_global_error(triangulation,
    //             difference_per_cell,
    //             VectorTools::L2_norm);
    
    data_out.add_data_vector(difference_per_cell, "Cellwise_L2");
  
    difference_per_cell = 0;
    VectorTools::integrate_difference(
      dof_handler,
      solution,
      sol_fcn,
      difference_per_cell,
      QGauss<2>(fe.degree + 1),
      VectorTools::H1_norm);
  
    this -> compute_h1_error();
  
    // H1 = VectorTools::compute_global_error(triangulation,
    //             difference_per_cell,
    //             VectorTools::H1_norm);

    data_out.add_data_vector(difference_per_cell, "Cellwise_H1");
    std::cout << "        ... done" << std::endl;

    //.vtu for visualization
    std::cout << "        Printing to files" << std::endl;

    const std::string vtu_name = problem_out.vtg.string();
    const std::string svg_name = problem_out.svg.string();
    std::ofstream vtu_out(vtu_name);

    data_out.build_patches();
    data_out.write_vtu(vtu_out);
                          
    if (triangulation.n_levels() - 1 < 11 && strategy == Strategy::adaptive) {
      std::string svg_name = problem_out.svg.string() + "grid_l" + std::to_string(triangulation.n_levels() - 1) + ".svg";
      GridOutFlags::Svg svg_flags;
      svg_flags.coloring = GridOutFlags::Svg::Coloring::none;
      // svg_flags.label_level_number  = true;
      // svg_flags.label_cell_index    = true;
      // svg_flags.label_boundary_id   = true;
  
      std::ofstream svg_out(svg_name);

      GridOut       grid_out;
      grid_out.set_flags(svg_flags);
  
      grid_out.write_svg(triangulation, svg_out);
    } 
    
    
  
  
    const unsigned int n_active_cells = triangulation.n_active_cells();
    const unsigned int n_dofs         = dof_handler.n_dofs();
  
    // output in console
    std::cout << "Cycle " << cycle++ << ":" << std::endl
    << "   Number of active cells:       " << n_active_cells << std::endl
    << "   Number of degrees of freedom: " << n_dofs << std::endl;
  
    // output in .text 
    problem_out.write_table_text();
    problem_out.write_table_tex();
  }
  
  
  
  void LShape_Benchmark_Standard::compute_h1_error(
  ) { 
    const QGauss<2> quadrature_formula(fe.degree + 1);

    FEValues<2> fe_values(
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

    problem_out.add_values_to_table(H1, L2);
  }
  
  
  
  void LShape_Benchmark_Standard::run(
    const unsigned int ref
  ) {
    this -> make_grid();

    std::cout << "    Setting system... " << std::endl;
    setup_system();

    //refinement cycle
    while (triangulation.n_levels() < ref + 1 ){
      std::cout << "    Assembling system ..." << std::endl;
      assemble_system();

      std::cout << "    Solving system..." << std::endl;
      solve();

      std::cout << "    Processing results..." << std::endl;
      process_results();

      std::cout << "    Refining grid..." << std::endl;
      refine_grid();
    }

    std::cout << std::endl;
    problem_out.write_table_text(std::cout);
  }
    //outputting in tables

} // namespace LShape



