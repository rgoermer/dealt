#include <poisson_benchmark.h>
#include <residual_error_estimators.h>

namespace Poisson_Benchmark {
        Poisson_Benchmark_3D::Poisson_Benchmark_3D(
    int ref, 
    int order
  ) : ref(ref)
    , order(order)
    , data(this->get_data(0, 1))
    , tria(data)
    , rhs_fcn(1)
    , nc_fcn(0)
    , dirichlet_fcn(0)
    , eps_fcn(1e-3) 
  {
    tria.degree_elevate(2, 1); 
    tria.degree_elevate_global(order);

    for (const auto& face : tria.active_face_iterators())
      if (face -> at_boundary())
        if ((face->center()).operator()(0) == 1 ||
            (face->center()).operator()(1) == 1 || 
            (face->center()).operator()(2) == 1)
          face -> set_boundary_id(Boundary_IDs::Neumann);
        else 
          face -> set_boundary_id(Boundary_IDs::Dirichlet);

    this->offset = 0;
    tria.refine_global(offset);

    const std::string name = "poisson_benchmark_3d/";
    dealt::OutputSetup::set_dim(3);
    problem_out = OutputSetup(name, 2+order);

    tria.set_boundary_dofs(); 
    tria.prepare_assembly();
  } // Poisson_Benchmark_3D

  void Poisson_Benchmark_3D::run(
  ) {
    unsigned int cycle = 0;  
    unsigned int old_level = tria.n_levels() - 1; 
    while (tria.n_levels() - 1 < ref && 
            residual > 1e-10) {
      this -> setup_system();
      this -> assemble_system();
      this -> solve();
      this -> output_system();

      if (cycle < 5) {
        this -> estimate_mark_refine();
      } else { 
        std::cout << "Too many refinements resulted in the same level, enforcing global refinement..." << std::endl;
        tria.coarsen_bezier_elements();
        tria.refine_global(); 
      }
      tria.prepare_assembly();

      if (tria.n_levels() - 1 == old_level) {
        cycle++;
      } else {
        cycle = 0;
        old_level = tria.n_levels() - 1;
      }
    } // while

    // Write the resulting table to line
    problem_out.write_table_text(std::cout);
  } // Poisson_Benchmark_3D::run()

  void Poisson_Benchmark_3D::setup_system(
  ) { 
    const unsigned int n_dofs = tria.n_active_splines();
    system_rhs.reinit(n_dofs);
    solution.reinit(n_dofs);

    const auto& IEN_array = tria.get_IEN_array();
    Assert(IEN_array.size() > 0, ExcInternalError());

    sparsity_pattern.reinit(
        n_dofs,
        n_dofs,
        tria.get_max_entries_per_row() 
    );

    for (const auto& [_, arr] : IEN_array)
      for (unsigned int i : arr)
        for (unsigned int j : arr)
          sparsity_pattern.add(i, j);

    sparsity_pattern.compress();
    system_matrix.reinit(sparsity_pattern);
  } // Poisson_Benchmark_3D::setup_system()

  void Poisson_Benchmark_3D::assemble_system(
  ) {
    std::vector< unsigned int > degrees = tria.get_degree();
    degrees[0] += 1; degrees[1] += 1; degrees[2] += 1;

    WorkStream::run(
      tria.begin_active(),
      tria.end(),
      *this,
      &Poisson_Benchmark_3D::local_assemble_system,
      &Poisson_Benchmark_3D::copy_local_to_global,
      AssemblyScratchData(&tria, degrees),
      AssemblyCopyData()
    );

    std::map<
        types::global_dof_index, 
        double
      > boundary_values;
    std::map<
        types::boundary_id,
        const Function< 3 >*
      > boundary_fcns =
          {{Boundary_IDs::Dirichlet, &dirichlet_fcn}};
    tria.project_boundary_values(
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
  } // Poisson_Benchmark_3D::assemble_system()

  void Poisson_Benchmark_3D::local_assemble_system(
    const typename Triangulation<3, 3>::active_cell_iterator        &cell,
    AssemblyScratchData                                             &scratch,
    AssemblyCopyData                                                &copy_data
  ) { 

    scratch.ts_values.reinit(cell);

    const unsigned int n_dofs_per_cell = scratch.ts_values.n_dofs_per_cell();
    copy_data.cell_matrix.reinit(n_dofs_per_cell, n_dofs_per_cell);
    copy_data.cell_rhs.reinit(n_dofs_per_cell);
    copy_data.local_dof_indices = scratch.tria -> get_IEN_array(cell);

    // Compute system_matrix part of cell
    const double eps = eps_fcn.value(Point<3>());
    const double rhs = rhs_fcn.value(Point<3>());
    for (const unsigned int i : scratch.ts_values.dof_indices()){
      for (const unsigned int j : scratch.ts_values.dof_indices()){
        for (const unsigned int q : scratch.ts_values.quadrature_point_indices()){
          copy_data.cell_matrix(i, j) += (
            eps * 
            scratch.ts_values.shape_grad(i, q) *
            scratch.ts_values.shape_grad(j, q) 
            +
            scratch.ts_values.shape_value(i, q) *
            scratch.ts_values.shape_value(j, q) * 
            1.
          ) * scratch.ts_values.JxW(q);
        } // for ( q )
      } // for ( j )

      // Compute rhs for this cell
      for (const unsigned int q : scratch.ts_values.quadrature_point_indices())
        copy_data.cell_rhs(i) += 
                        rhs *
                        scratch.ts_values.shape_value(i, q) * 
                        scratch.ts_values.JxW(q);
    } // for ( i ) 

    // if (!cell -> at_boundary())
    //   return; 

    // const double nc_value = nc_fcn.value(Point<3>());
    // for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; f++)
    //   if (cell->face(f)->at_boundary())
    //     if (cell->face(f)->boundary_id() == Boundary_IDs::Neumann) {
    //       scratch.face_values.reinit(cell, f); 
    //       for (const unsigned int i : scratch.face_values.dof_indices())
    //         for (const unsigned int q : scratch.face_values.quadrature_point_indices())
    //           copy_data.cell_rhs(i) += 
    //                   nc_value * 
    //                   scratch.face_values.shape_value(i, q) *
    //                   scratch.face_values.JxW(q);
    //     } // if (boundary_id)

  } // Poisson_Benchmark_3D::local_asemble_system()

  void Poisson_Benchmark_3D::copy_local_to_global(
    const AssemblyCopyData &copy_data
  ) { 
    system_matrix.add(copy_data.local_dof_indices, copy_data.cell_matrix);
    system_rhs.add(copy_data.local_dof_indices, copy_data.cell_rhs);
  } // Poisson_Benchmark::copy_local_to_global()

  void Poisson_Benchmark_3D::solve(
  ) {
    std::cout << "Solving system ... " << std::endl;

    SolverControl solver_control(
                     750 *  tria.n_active_splines(), 
                     1e-10
                  );
    solver_control.enable_history_data();

    SolverCG<Vector<double>> solver(solver_control);

    PreconditionJacobi<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix);
    
    solver.solve(system_matrix, solution, system_rhs, preconditioner);

    problem_out.add_values_to_table(
      tria.n_levels() - 1 - offset,
      tria.n_active_cells(),
      tria.n_active_splines() ,
      solver_control.last_step(),
      *(solver_control.get_history_data().end() - 1)
    );
  
    std::cout << " ... done!" << std::endl;
  } // Poisson_Benchmark_3D::solve()

  void Poisson_Benchmark_3D::estimate_mark_refine(
  ) { 
    Vector<double> residuals(tria.n_active_cells()) ;
    std::vector< unsigned int > degrees = tria.get_degree();
    degrees[0] = degrees[0] * degrees[0] + 1;
    degrees[1] = degrees[1] * degrees[1] + 1;
    degrees[2] = degrees[2] * degrees[2] + 1;

    std::map< 
              types::boundary_id, 
        const Function<3>* 
      > neumann_bc = {{Boundary_IDs::Neumann, &nc_fcn}};
    Functions::ConstantFunction<3> a(1.);
    ResidualEstimators::Poisson<3>::estimate(
            &tria,
            degrees,
            solution,
            residuals,
            &rhs_fcn,
            neumann_bc,
            &eps_fcn,
            &a
        );


    residual = residuals.l2_norm(); 
    max_residual = residuals.linfty_norm();
    tria.refine_fixed_number(residuals, 0.10);
  } // void Poisson_Benchmark_3D::estimate_mark_refine()

  void Poisson_Benchmark_3D::output_system(
  ) {
    const unsigned int level = tria.n_levels() - 1;

    // Write the grid to a seperate file: 
    const std::string& name_vtg = problem_out.vtg.string() + "physical_grid_l" + std::to_string(level) + ".vtk";

    // First: Make a copy of the triangulation: 
    Triangulation<3> physical_grid; 
    physical_grid.copy_triangulation(tria);


    // GridTools::transform(...) does not work with anisotropically refined meshes in 3D
    // Thus, we perform the transformation manually ... *sigh*
    const IsoparametricManifold<3> geometry(tria.get_IPF()); 
    GridTools::transform(
      [&geometry](const Point<3> p){
        return geometry.push_forward(p);
      }, 
      physical_grid
    );

    // Generate the output object
    DataOut<3> data_out;
    data_out.attach_triangulation(physical_grid); 

    // Build patches
    data_out.build_patches(); 

    // Open a file
    std::ofstream vtk_out(name_vtg); 
    data_out.write_vtk(vtk_out);


    problem_out.add_values_to_table(max_residual, residual);

    problem_out.write_table_text();
    problem_out.write_table_tex();

    for (double z = 0.1; z <= 1; z += 0.3) {
      GridPartialOut::get_cross_section(&tria, 2, z, true, 
                        problem_out.cross_sections_z.string() 
                        + "l" 
                        + std::to_string(level) 
                        + "_physical");
    }
    for (double y = 0.1; y <= 1; y += 0.3) {
      GridPartialOut::get_cross_section(&tria, 1, y, true, 
                        problem_out.cross_sections_y.string() 
                        + "l" 
                        + std::to_string(level) 
                        + "_physical");
    }
    for (double x = 0.1; x <= 1; x += 0.3) {
      GridPartialOut::get_cross_section(&tria, 0, x, true, 
                        problem_out.cross_sections_x.string() 
                        + "l" 
                        + std::to_string(level) 
                        + "_physical");
    }

    if (level > 15)
      return;

    const std::string name = problem_out.degree.string();
    const std::string level_name = name + "l" + std::to_string(level-offset);
  
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

    tria.printIPF(evals, level_name, 16, true, true);
  } // Poisson_Benchmark_3D::output_system()


  IPF_Data<3, 3> Poisson_Benchmark_3D::get_data(
    const double r,
    const double R
  ) {
    std::vector< std::vector< double > > kv(3);
    std::vector< Point<3> >              points(18);
    std::vector< double >                weights(18, 1.);
    std::vector< unsigned int >          deg = {2, 2, 1};
    const double sqrt2 = 1. / std::sqrt(2.);

    // Define vector of control points
    // points[0] = r * Point<3>(1, 0, 0); 
    // points[1] = r * Point<3>(1, 0, 1); weights[1] = sqrt2;
    // points[2] = r * Point<3>(0, 0, 1);
    // points[3] = r * Point<3>(1, 1, 0); weights[3] = sqrt2;
    // points[4] = r * Point<3>(1, 1, 1); weights[4] = 0.5;
    // points[5] = r * Point<3>(0, 0, 1); weights[5] = sqrt2;
    // points[6] = r * Point<3>(0, 1, 0);
    // points[7] = r * Point<3>(0, 1, 1); weights[7] = sqrt2;
    // points[8] = r * Point<3>(0, 0, 1);

    // points[ 9] = R * Point<3>(1, 0, 0);
    // points[10] = R * Point<3>(1, 0, 1); weights[10] = sqrt2;
    // points[11] = R * Point<3>(0, 0, 1);
    // points[12] = R * Point<3>(1, 1, 0); weights[12] = sqrt2;
    // points[13] = R * Point<3>(1, 1, 1); weights[13] = 0.5; 
    // points[14] = R * Point<3>(0, 0, 1); weights[14] = sqrt2;
    // points[15] = R * Point<3>(0, 1, 0); 
    // points[16] = R * Point<3>(0, 1, 1); weights[16] = sqrt2;
    // points[17] = R * Point<3>(0, 0, 1);
    points[0 ] = Point<3>(0.0, 0.0, 0.0);
    points[1 ] = Point<3>(0.5, 0.0, 0.0);
    points[2 ] = Point<3>(1.0, 0.0, 0.0);
    points[3 ] = Point<3>(0.0, 0.5, 0.0);
    points[4 ] = Point<3>(0.5, 0.5, 0.0);
    points[5 ] = Point<3>(1.0, 0.5, 0.0);
    points[6 ] = Point<3>(0.0, 1.0, 0.0);
    points[7 ] = Point<3>(0.5, 1.0, 0.0);
    points[8 ] = Point<3>(1.0, 1.0, 0.0);

    points[9 ] = Point<3>(0.0, 0.0, 1.0);
    points[10] = Point<3>(0.5, 0.0, 1.0);
    points[11] = Point<3>(1.0, 0.0, 1.0);
    points[12] = Point<3>(0.0, 0.5, 1.0);
    points[13] = Point<3>(0.5, 0.5, 1.0);
    points[14] = Point<3>(1.0, 0.5, 1.0);
    points[15] = Point<3>(0.0, 1.0, 1.0);
    points[16] = Point<3>(0.5, 1.0, 1.0);
    points[17] = Point<3>(1.0, 1.0, 1.0);


    kv[0] = {0, 0, 0, 1, 1, 1};
    kv[1] = {0, 0, 0, 1, 1, 1};
    kv[2] = {0, 0, 1, 1};

    return IPF_Data<3, 3>(points, weights, kv, deg);
  } // Poisson_Benchmark_3D::get_data()
} // namespace Poisson_Benchmark
