#include <beam.h>

namespace Classical {
  
  // ==============================================================
  //
  //                   RHS / Neumann
  //
  // ==============================================================
  double ClassicalProblemRHS::value(
    const Point<3>        &/* eval */,
    const unsigned int     component
  ) const { 
    return 0;
    // return component == 2 ? -1. : 0;
  }

  double ClassicalProblemNeumann::value(
    const Point<3>        &/* eval */,
    const unsigned int     component
  ) const { 
    return component == 2 ? -1. : 0;
    // return 0.;
  }

  // ==============================================================
  //
  //                   Problem
  //
  // ==============================================================
  
  ClassicalProblem::AssemblyScratchData::AssemblyScratchData(
          TS_Triangulation<3>         *tria,
    const std::vector<unsigned int>   &degrees
  ) : tria(tria) 
    , ts_values(tria, 
                degrees, 
                update_values |
                update_gradients |
                update_JxW_values |
                update_quadrature_points) 
    , face_values(tria, 
                  degrees,
                  update_values |
                  update_JxW_values |
                  update_quadrature_points){}

  ClassicalProblem::AssemblyScratchData::AssemblyScratchData(
    const AssemblyScratchData &scratch
  ) : tria(scratch.tria)
    , ts_values(scratch.ts_values) 
    , face_values(scratch.face_values){}


  ClassicalProblem::ClassicalProblem(
    const unsigned int ref,
    const unsigned int order_elevate,
    const Type         problem
  ) : data(DataGenerator(problem).data)
    , ref(ref)
    , order(order_elevate + data.max_degree())
    , tria(data)
    , rhs_fcn()
    , nc_fcn()
    , lambda(1.)
    , mu(1.)
  { 
    dealt::OutputSetup::set_dim(3);
    if (problem == Type::Beam){
      tria.degree_elevate(1, 1);
      tria.degree_elevate(2, 1);
      problem_out = OutputSetup("beam/", order);
    } else {
      tria.degree_elevate_global(1);
      problem_out = OutputSetup("slab/", order + 1);
    }
    tria.degree_elevate_global(order_elevate);


    this->offset = 0;
    tria.refine_global(offset);
          
    if (problem == Type::Beam)
      for (const auto& face : tria.active_face_iterators()) {
        if (face -> at_boundary())
          if ((face -> center()).operator()(0) == 0 )
            face -> set_boundary_id(Boundaries::Dirichlet);
          else if ((face -> center()).operator()(2) == 1 )
            face -> set_boundary_id(Boundaries::Neumann);
          else 
            face -> set_boundary_id(Boundaries::HomogeneousNeumann);
      }
    else 
      for (const auto& face: tria.active_face_iterators()) {
        if (face->at_boundary())
          if ((face->center()).operator()(2) == 0)
            face -> set_boundary_id(Boundaries::Dirichlet);
          else if ((face->center()).operator()(2) == 1)
            face -> set_boundary_id(Boundaries::Neumann);
          else 
            face -> set_boundary_id(Boundaries::HomogeneousNeumann);
      }

    run_times.set_auto_fill_mode(true);

    run_times.declare_column("Level");        // Done and done

    run_times.declare_column("Lin. Solve");   // Done and done
    run_times.declare_column("Assembly");     // Done and done
    run_times.add_column_to_supercolumn("Lin. Solve", "Solve");
    run_times.add_column_to_supercolumn("Assembly",   "Solve");
    

    run_times.declare_column("Estimate");             // Done and done

    run_times.declare_column("Mark and Refine");      // Done and done
    run_times.declare_column("Set Boundary DoFs");    // Done and done
    run_times.declare_column("Refine Bezier");        // Done and done
    run_times.declare_column("Extraction Operators"); // Done and done

    run_times.add_column_to_supercolumn("Set Boundary DoFs",    "Refine");
    run_times.add_column_to_supercolumn("Refine Bezier",        "Refine");
    run_times.add_column_to_supercolumn("Extraction Operators", "Refine");
    run_times.add_column_to_supercolumn("Mark and Refine",      "Refine");

    prepare_assembly_and_measure_time();

    // tria.prepare_assembly();
  } // ClassicalProblem::ClassicalProblem()

  void ClassicalProblem::prepare_assembly_and_measure_time(
  ) {
    auto t_set_boundary_dofs_start = high_resolution_clock::now();
    tria.set_boundary_dofs();
    auto t_set_boundary_dofs_end   = high_resolution_clock::now();

    auto t_refine_bezier_elements_start = high_resolution_clock::now();
    tria.refine_bezier_elements();
    auto t_refine_bezier_elements_end   = high_resolution_clock::now();

    auto t_compute_extraction_operators_start = high_resolution_clock ::now();
    tria.compute_extraction_operators();
    auto t_compute_extraction_operators_end   = high_resolution_clock::now();

    duration<double> t_boundary_dofs = 
            t_set_boundary_dofs_end - t_set_boundary_dofs_start;
    duration<double> t_refine_bezier_elements = 
            t_refine_bezier_elements_end - t_refine_bezier_elements_start;
    duration<double> t_compute_extraction_operators = 
            t_compute_extraction_operators_end - t_compute_extraction_operators_start;

    run_times.add_value(
        "Set Boundary DoFs",
        convert_time(t_boundary_dofs)
    );
    run_times.add_value(
        "Refine Bezier",
        convert_time(t_refine_bezier_elements)
    );
    run_times.add_value(
        "Extraction Operators",
        convert_time(t_compute_extraction_operators)
    );
    run_times.add_value("Level", tria.n_levels()-1);
  } // ClassicalProblem::prepare_assembly_and_measure_time();

  const std::string ClassicalProblem::convert_time(
    const std::chrono::duration<double>& time
  ) const {
    std::string out;
    // time is given in seconds! 
    if (time < std::chrono::seconds(1)) {
      out = std::to_string(duration_cast<std::chrono::milliseconds>(time).count()) + "ms";
    } else if (time < std::chrono::minutes(1)) {
      out = std::to_string(time.count()) + "s";
    } else if (time < std::chrono::hours(1)){
      out = std::to_string(duration_cast<std::chrono::minutes>(time).count()) + "m";
    } else {
      out = std::to_string(duration_cast<std::chrono::hours>(time).count()) + "h";
    }

    return out; 
  }

  void ClassicalProblem::local_assemble_system(
    const typename Triangulation<3, 3>::active_cell_iterator     &cell,
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

    // Apply Neumann Boundary:
    if (cell->at_boundary()) {
      for (unsigned int face = 0; face < GeometryInfo<3>::faces_per_cell; face++){
        if (cell->face(face)->at_boundary() && 
            cell->face(face)->boundary_id() == Boundaries::Neumann) {
          scratch.face_values.reinit(cell, face);
          for (const unsigned int i : scratch.face_values.dof_indices()){
            const unsigned int comp_i =
                    scratch.face_values.system_to_component_index(i).first;
            for (const unsigned int q : scratch.face_values.quadrature_point_indices())
              copy_data.cell_rhs(i) += 
                  scratch.face_values.shape_value(i, q) *
                  nc_fcn.value(scratch.face_values.quadrature_point(q), comp_i) *
                  scratch.ts_values.JxW(q);
          }
        }
      }
    }
    // std::cout << "    ... done!" << std::endl;

    // std::cout << "    Copying data..." << std::endl;
    copy_data.local_dof_indices = scratch.tria -> get_IEN_array(cell, 3);
    // std::cout << "    ... done!" << std::endl;
  } // ClassicalProblem::local_assemble_system()

  void ClassicalProblem::copy_local_to_global(
    const AssemblyCopyData &copy_data
  ) {
    system_matrix.add(copy_data.local_dof_indices, copy_data.cell_matrix);
    system_rhs.add(copy_data.local_dof_indices, copy_data.cell_rhs);
  } // ClassicalProblem::copy_local_to_global();

  void ClassicalProblem::assemble_system(
  ) {
    auto t_assemble_system_start = high_resolution_clock::now();

    std::vector< unsigned int > degrees = tria.get_degree();
    degrees[0] += 1; degrees[1] += 1; degrees[2] += 1;

    AssemblyCopyData    copy_data;
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


    std::cout << "    ... done! Applying homogeneous boundary values ..." << std::endl;
    
    // Apply homogeneous boundary conditions
    const auto& boundary_dofs = tria.get_boundary_dofs(3);
    for (const auto& dof : boundary_dofs.at(Boundaries::Dirichlet)){
      for (unsigned int j = 0; j < 3 * tria.n_active_splines(); j++) {
        system_matrix.set(dof, j, 0);
        system_matrix.set(j, dof, 0);
      }
      system_matrix.set(dof, dof, 1);
      system_rhs(dof) = 0;
    }
    // std::cout << "    Homogeneous Data applied" << std::endl;
    // std::cout << "... done!" << std::endl;

    auto t_assemble_system_end = high_resolution_clock::now();
    duration<double> t_assemble_system = 
            t_assemble_system_end - t_assemble_system_start;

    run_times.add_value(
        "Assembly",
        convert_time(t_assemble_system)
    );

  } // ClassicalProblem::assemble_system()

  void ClassicalProblem::setup_system(
  ) {
    const unsigned int n_dofs = 3 * tria.n_active_splines();
    system_rhs.reinit(n_dofs);
    solution.reinit(n_dofs);

    const auto& IEN_array = tria.get_IEN_array(3);
    Assert(IEN_array.size() > 0, ExcInternalError());

    sparsity_pattern.reinit(
        n_dofs,
        n_dofs,
        tria.get_max_entries_per_row(3)
    );

    for (const auto& [_, arr] : IEN_array)
      for (unsigned int i : arr)
        for (unsigned int j : arr)
          sparsity_pattern.add(i, j);

    sparsity_pattern.compress();

    system_matrix.reinit(sparsity_pattern);
  } // ClassicalProblem::setup_system()

  void ClassicalProblem::solve(
  ) {
    std::cout << "Solving system ... " << std::endl;

    auto t_solve_start = high_resolution_clock::now();
    SolverControl solver_control(
                     750 * 3 * tria.n_active_splines(), 
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

    auto t_solve_end = high_resolution_clock::now();

    duration<double> t_solve = 
            t_solve_end - t_solve_start; 
    run_times.add_value(
       "Lin. Solve",
       convert_time(t_solve)
    );

    problem_out.add_values_to_table(
      tria.n_levels() - 1 - offset,
      tria.n_active_cells(),
      tria.n_active_splines() * 3,
      solver_control.last_step(),
      *(solver_control.get_history_data().end() - 1)
    );
  
    std::cout << " ... done!" << std::endl;
  } // ClassicalProblem::solve()

  void ClassicalProblem::compute_residual(
  ) { 
    std::vector< unsigned int > degrees = tria.get_degree();
    degrees[0] += 1; degrees[1] += 1; degrees[2] += 1;
    const IsoparametricManifold<3> geometry(tria.get_IPF()); 

    TSValues<3, 3, 3> ts_values(
      &tria, 
      degrees, 
      update_values |
      update_gradients |
      update_hessians  |
      update_JxW_values |
      update_quadrature_points
    );

    residual = 0.;
    for (const auto& cell : tria.active_cell_iterators()) {
      double cell_residual = 0;
      const std::vector< unsigned int > local_IEN_array
        = tria.get_IEN_array(cell, 3);
      ts_values.reinit(cell);
      
      for (const unsigned int q : ts_values.quadrature_point_indices()) {
        Tensor<1, 3>    integrand;
        const Point<3>& Q        = ts_values.quadrature_point(q);
        const double    lambda_q = lambda.value(Q);
        const double    mu_q     = mu.value(Q);
        for (const unsigned int i : ts_values.dof_indices()) {
          const double c_i = solution(local_IEN_array[i]);
          const Tensor<2, 3>& hess = 
            ts_values.shape_hessian(i, q);
          const unsigned int comp_i = 
            ts_values.system_to_component_index(i).first;

          for (unsigned int d = 0; d < 3 ; d++) {
            integrand[d]      += c_i * (lambda_q + mu_q)  *
                                        hess[comp_i][d];  // \lambda \nabla (\nabla \cdot u)
            integrand[comp_i] += c_i * mu_q * hess[d][d]; // \mu (\nabla \cdot \nabla u)
          }
        } // for ( i )

        for (unsigned int d = 0; d < 3; d++)
          integrand[d] += rhs_fcn.value(Q, d);

        cell_residual += integrand * integrand * ts_values.JxW(q);
      } // for ( q )
      const double cell_width = (geometry.push_forward(cell->vertex(0))).distance(
                                 geometry.push_forward(cell->vertex(7)));
      cell_residual *= cell_width * cell_width / 24.;

      residual += cell_residual;
    } // cell
    residual = std::sqrt(residual);
  } // ClassicalProblem::compute_residual

  void ClassicalProblem::estimate_mark_refine(
  ) { 
    Vector<double> residuals(tria.n_active_cells()) ;
    std::vector< unsigned int > degrees = tria.get_degree();
    for (auto& p : degrees)
      p = p + 2;
    const Functions::ConstantFunction<3> homogeneous_neumann(0., 3);
    std::map< types::boundary_id,
              const Function<3>* > neumann_data = 
              {{Boundaries::Neumann, &nc_fcn},
               {Boundaries::HomogeneousNeumann, &homogeneous_neumann}};


    auto t_estimate_start = high_resolution_clock::now();
    tria.linear_elasticity_residual_error_estimate(
        degrees,
        &rhs_fcn, 
        &lambda,
        &mu,
        neumann_data,
        solution,
        residuals
    );
    auto t_estimate_end = high_resolution_clock::now();

    compute_residual();
    max_residual = residuals.linfty_norm();
    problem_out.add_values_to_table(max_residual, residual);


    
    auto t_mark_and_refine_start = high_resolution_clock::now();
    tria.refine_fixed_number(residuals, 0.10);
    auto t_mark_and_refine_end = high_resolution_clock::now();

    duration<double> t_estimate = 
            t_estimate_end - t_estimate_start;
    duration<double> t_mark_and_refine = 
            t_mark_and_refine_end - t_mark_and_refine_start;

    run_times.add_value(
        "Estimate",
        convert_time(t_estimate)
    );
    run_times.add_value(
        "Mark and Refine",
        convert_time(t_mark_and_refine)
    );

    // tria.prepare_assembly();
  } // ClassicalProblem::estimate_mark_refine()

  void ClassicalProblem::run(
  ) {
    
    unsigned int cycle = 0;  
    unsigned int old_level = tria.n_levels() - 1; 
    unsigned int r = 0;

    std::string run_times_tex_out_name(problem_out.degree.string() + "run_times.tex");
    std::string run_times_text_out_name(problem_out.degree.string() + "run_times.txt");


    while (tria.n_levels() - 1 < ref) {
      this -> setup_system();
      this -> assemble_system();
      this -> solve();
      // this -> output_system();

      if (residual < 1e-10 )
        break;

      if (cycle < 5) {
        this -> estimate_mark_refine();
      } else { 
        std::cout << "Too many refinements resulted in the same level, enforcing global refinement..." << std::endl;
        tria.refine_global(); 
      }
      std::cout << "... done!" << std::endl;

      if (tria.n_levels() - 1 == old_level) {
        cycle++;
      } else {
        cycle = 0;
        old_level = tria.n_levels() - 1;
      }

      std::filebuf file_tex_out;
      std::filebuf file_text_out;

      file_tex_out.open(run_times_tex_out_name.c_str(), std::ios::out);
      file_text_out.open(run_times_text_out_name.c_str(), std::ios::out);

      std::ostream run_times_tex_out(&file_tex_out);
      std::ostream run_times_text_out(&file_text_out);

      run_times.write_tex(run_times_tex_out, false);
      run_times.write_text(run_times_text_out, dealii::TableHandler::TextOutputFormat::org_mode_table);

      file_tex_out.close();
      file_text_out.close();

      problem_out.write_table_text();
      problem_out.write_table_tex();

      run_times.start_new_row();
      this -> prepare_assembly_and_measure_time();
    } // while

    // Write the resulting table to line
    problem_out.write_table_text(std::cout);
  } // ClassicalProblem::run()


  void ClassicalProblem::output_system(
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
    GridTools::transform(
      [&geometry](const Point<3>& p){ return geometry.push_forward(p);},
      physical_grid
    );


    // Generate the output object
    DataOut<3> data_out;
    data_out.attach_triangulation(physical_grid); 

    // Estimate the error on each cell
    Vector< double > residuals(tria.n_active_cells());
    std::vector<unsigned int> degrees = tria.get_degree();
    for (unsigned int p : degrees)
      p = p * p + 1;

    const Functions::ConstantFunction<3> homogeneous_neumann(0.);
    std::map< types::boundary_id,
              const Function<3>* > neumann_data = 
              {{Boundaries::Neumann, &nc_fcn},
               {Boundaries::HomogeneousNeumann, &homogeneous_neumann}};

    tria.linear_elasticity_residual_error_estimate(
        degrees,
        &rhs_fcn, 
        &lambda,
        &mu,
        neumann_data,
        solution,
        residuals
    );

    Vector<double> levels (tria.n_active_cells());
    auto cell = tria.begin_active();
    for (unsigned int n = 0; n < tria.n_active_cells(); n++){
      levels(n) = cell->level();
      cell++;
    }

    data_out.add_data_vector(residuals, "cell_errors");
    data_out.add_data_vector(levels, "levels");

    // Build patches
    data_out.build_patches(); 

    // Open a file
    std::ofstream vtk_out(name_vtg); 
    data_out.write_vtk(vtk_out);

    // residual = residuals.l2_norm() ; 
    // max_residual = residuals.linfty_norm();
    // problem_out.add_values_to_table(max_residual, residual);

    // problem_out.write_table_text();
    // problem_out.write_table_tex();

    if (level > 1)
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

    const unsigned int N1 = 100;
    const unsigned int N2 = 25;
    const unsigned int N3 = 25;
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
    tria.coarsen_bezier_elements();
    tria.print_IPF_wireframe(level_name);
    tria.refine_bezier_elements();
  } // ClassicalProblem::output_system()

  DataGenerator::DataGenerator(
    const Type problem
  ) { 
    if (problem == Type::Slab)
      generate_slab_data();
    else 
      generate_beam_data();


  } // DataGenerator

  void DataGenerator::generate_beam_data(
  ) { 
    std::vector< std::vector< double > > kv(3); 
    kv[1] = {0., 0., 1., 1.};
    kv[2] = kv[1];
    kv[0] = {0., 0., 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 10., 10.};

    std::vector<unsigned int> degrees = {2, 1, 1}; 

    std::vector< Point<3> > controls(48, Point<3>());
    for (unsigned int i = 0; i < 12; i++)
      controls[i] = Point<3>(0 + 10. * i / 11., 0, 0);

    for (unsigned int i = 12; i < 24; i++){
      controls[i] = controls[i-12];
      controls[i](1) = 1;
    }

    for (unsigned int i = 24; i < 48; i++){
      controls[i] = controls[i-24];
      controls[i](2) = 1; 
    }

    std::vector<double> weights(48, 1.);
    this -> data = IPF_Data<3>(controls, weights, kv, degrees);
  } // DataGenerator::generate_beam_data

  void DataGenerator::generate_slab_data(
  ) {
    std::vector< std::vector<double> > kv(3, {0., 0., 1., 1.});

    std::vector<unsigned int> degrees = {1, 1, 1};

    std::vector< Point<3> > controls(8, Point<3>());
    const double height = 0.25;
    controls[0] = Point<3>(0, 0, 0);
    controls[1] = Point<3>(1, 0, 0);

    controls[2] = Point<3>(0, 1, 0);
    controls[3] = Point<3>(1, 1, 0);
    
    controls[4] = Point<3>(0, 0, height);
    controls[5] = Point<3>(1, 0, height);

    controls[6] = Point<3>(0, 1, height);
    controls[7] = Point<3>(1, 1, height);

    std::vector<double> weights(8, 1.);
    this -> data = IPF_Data<3>(controls, weights, kv, degrees);
  }

} // namespace
