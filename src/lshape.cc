#include <lshape.h>
#include <residual_error_estimators.h>

namespace LShape {
  template class LShape_Problem<2>;
  template class LShape_Problem<3>;

  template class LShape_NC3<2>;
  template class LShape_NC4<2>;
  
  
  
  // =========================================================
  // =========================================================
  //              Instantiations
  // =========================================================
  // =========================================================
  //
  // =========================================================
  // =========================================================
  //              DataGenerator implementation 
  // =========================================================
  // =========================================================
  template<int dim, int spacedim>
  DataGenerator<dim, spacedim>::DataGenerator(
      const unsigned int /* n_elements */, 
      const bool         /* no_c0_edges */
  ) { 
    Assert(dim > 1, ExcNotImplemented());
    Assert(dim < 4, ExcNotImplemented());
  }
  

  // =========================================================
  // =========================================================
  //             LShape_RHS Implementation 
  // =========================================================
  // =========================================================
  template<int dim>
  double LShape_RHS_benchmark<dim>::value(
    const Point<dim>&       /* p */,
    const unsigned int    /* component */
  ) const { 
    Assert(dim > 1, ExcNotImplemented());
    Assert(dim < 4, ExcNotImplemented());
    return 0.;
  }
  
  
  template<int dim>
  double LShape_RHS_classic<dim>::value(
    const Point<dim>&       /* p */,
    const unsigned int    /* component */
  ) const  {
    Assert(dim > 1, ExcNotImplemented());
    Assert(dim < 4, ExcNotImplemented());
    return 0.;
  }
  
  template<int dim>
  double LShape_RHS<dim>::value(
      const Point<dim>&     p,
      const unsigned int  /* component */
  ) const {
    double out;
    switch (type) {
      case classic: {
        //LShape_RHS_classic sol;
        //out = sol.value(p);
        out = 0;
        break;
      }
      case benchmark: {
        LShape_RHS_benchmark<dim> sol;
        out = sol.value(p);
        break;
      }
      default: {
        out = nan("1");
      }
    } // switch
    return out;
  }
  
  // =========================================================
  // =========================================================
  //             LShape_SOL Implementation 
  // =========================================================
  // =========================================================
  template<int dim>
  double LShape_SOL_benchmark<dim>::value(
    const Point<dim>&     p,
    const unsigned int  /* component */
  ) const {
    Assert(dim > 1, ExcNotImplemented());
    Assert(dim < 4, ExcNotImplemented());
    return 0.;
  }
  
  template<int dim>
  Tensor<1, dim> LShape_SOL_benchmark<dim>::gradient(
    const Point<dim>&     p,
    const unsigned int  /* component */
  ) const {
    Assert(dim > 1, ExcNotImplemented());
    Assert(dim < 4, ExcNotImplemented());
    Tensor<1, dim> out;
    return out;
  }
  
  template<int dim>
  double LShape_SOL_classic<dim>::value(
    const Point<dim>&     p,
    const unsigned int  /* component */
  ) const { 
    Assert(dim > 1, ExcNotImplemented());
    Assert(dim < 4, ExcNotImplemented());
    return 0.;
  }
  
  template<int dim>
  Tensor<1, dim> LShape_SOL_classic<dim>::gradient(
    const Point<dim>&     p,
    const unsigned int  /* component */
  ) const {
    Assert(dim > 1, ExcNotImplemented());
    Assert(dim < 4, ExcNotImplemented());
    Tensor<1, dim> out;
    return out;
  }
   
  template<int dim>
  double LShape_SOL<dim>::value(
    const Point<dim>&    p,
    const unsigned int /* component */
  ) const {
    double out;
    switch (type) {
      case classic: {
        LShape_SOL_classic<dim> sol;
        out = sol.value(p);
        break;
      }
      case benchmark: {
        LShape_SOL_benchmark<dim> sol;
        out = sol.value(p);
        break;
      }
      default: {
        out = nan("1");
      }
    } // switch
    return out;
  }
  
  template<int dim>
  Tensor<1, dim> LShape_SOL<dim>::gradient(
    const Point<dim>&    p,
    const unsigned int /* component */
  ) const {
    Tensor<1, dim> out;
    switch (type) {
      case classic: {
        LShape_SOL_classic<dim> sol;
        out = sol.gradient(p);
        break;
      }
      case benchmark: {
        LShape_SOL_benchmark<dim> sol;
        out = sol.gradient(p);
        break;
      }
      default: {
        std::cout << "Unspecified type! " << std::endl;
        out[0] = nan("1");
        out[1] = nan("1");
        (dim > 2 ) ? out[2] = nan("1") : out[1] = nan("1");
      }
    } // switch
    AssertIsFinite(out[0]);
    AssertIsFinite(out[1]);
    return out;
  }
  
  // =========================================================
  // =========================================================
  //             LShape_Problem Implementation 
  // =========================================================
  // =========================================================
  template<int dim>
  void LShape_Problem<dim>::set_boundary_ids(
  ){
    Assert(false, ExcNotImplemented());
  }

  template<int dim>
  void LShape_Problem<dim>::print_grid(
      const std::string& name
  ){} // does nothing

  template<>
  void LShape_Problem<2>::print_grid(
      const std::string& name
  ){
    // Don't print grid if level too high.
    if (tria.n_levels() - offset - 1 > 15)
      return; 

    GridOutFlags::Svg svg_flags;
//    svg_flags.label_level_number = true;
//    svg_flags.label_cell_index = true;
//    svg_flags.label_boundary_id = true;
    svg_flags.boundary_line_thickness = 1;
    svg_flags.line_thickness          = 1;
    svg_flags.coloring = GridOutFlags::Svg::Coloring::none;
  
    std::ofstream out(name + ".svg");
    GridOut       grid_out;
    grid_out.set_flags(svg_flags);
  
    grid_out.write_svg(tria, out);

    tria.generate_mesh_file<0>(name, true, 16); 
    tria.generate_mesh_file<0>(name, false, 16);
  }

  template<>
  void LShape_Problem<3>::print_grid(
      const std::string& name
  ){
    GridOutFlags::Eps<3> eps_flags;
    eps_flags.azimut_angle = 45;
    eps_flags.turn_angle = 45;
  
    std::ofstream out(name + ".eps");
    GridOut       grid_out;
    grid_out.set_flags(eps_flags);
  
    grid_out.write_eps(tria, out);
    tria.print_grid(name + ".dat", 10);
  }
  
  template<>
  void LShape_Problem<2>::set_boundary_ids(
  ) {
    for (const auto& face : 
        tria.active_face_iterators()){
      if (face -> at_boundary()) {
        const auto& c = face -> center();
        if ( std::fabs(c(0) + 0.) < 1e-15 )
          face -> set_boundary_id(Boundary::neumann_1);
        else if ( std::fabs(c(0) - 1.) < 1e-15 )
          face -> set_boundary_id(Boundary::dirichleth);
        else if ( std::fabs(c(1) + 0.) < 1e-15 )
          face -> set_boundary_id(Boundary::dirichleth); 
        else if ( std::fabs(c(1) - 1.) < 1e-15 )
          face -> set_boundary_id(Boundary::neumann_2);
        else 
          face -> set_boundary_id(Boundary::none);
      }
    }
  }
  
  template<>
  void LShape_Problem<3>::set_boundary_ids(
  ) {
          // Do nothing
  }

  template<int dim>
  LShape_Problem<dim>::LShape_Problem(
    const unsigned int ref, 
    const unsigned int order,
    const Specifier    type,
    const Strategy     strategy
  ) : type(type),
      strategy(strategy),
      ref(ref),
      order(order),
      data(DataGenerator<dim, dim>(1, true).data),
      tria(data),
      rhs_fcn(type),
      sol_fcn(type),
      neumann_bc1(type),
      neumann_bc2(type),
      neumann_bc_diag(type)
  {
    // initialize std::strings for output files:
    const std::string t_str = (type == Specifier::benchmark) ? "b" : "c";
    const std::string s_str = (strategy == Strategy::adaptive) ? "a" : "u";
    const std::string d_suffix = (dim == 2 ? "" : "_3d"); 
  
    const std::string out_init = "lshape" 
                + d_suffix 
                + "/lshape_" 
                + s_str + t_str + "/";

    problem_out = OutputSetup(out_init, data.max_degree() + order);
  
    // set the boundary id for the corresponding dimension
    set_boundary_ids();
  
    tria.degree_elevate_global( order );
  
    // Assure the assumption is fulfilled, i.e. that
    // there are at least 2 by 2 cells in the active region
    offset = 3 * (this->data.max_degree() + order) / 2;
    tria.refine_global(offset);
  

    const std::string name = problem_out.svg.string() + "l0";
    print_grid(name);
  
    // Initialize information about boundary dofs
    tria.set_boundary_dofs();
  
    // This call does nothing at the beginning, but it switches the mesh to
    // a bezier mesh, i.e. TS_Triangulation::is_bezier = true. 
    // We are only allowed to compute on the bezier mesh
    tria.refine_bezier_elements();
    tria.compute_extraction_operators();
  }
  
  
  // =========================================================
  // =========================================================
  //             LShape_NC Implementation 
  // =========================================================
  // =========================================================
  template<int dim>
  double LShape_NC1<dim>::value(
    const Point<dim>&     p,
    const unsigned int  /* component */
  ) const {
    Tensor<1, dim> df;
    if (type == classic){
      const LShape_SOL_classic<dim> f;
      df = f.gradient(p);
    } else if (type == benchmark){
      const LShape_SOL_benchmark<dim> f;
      df = f.gradient(p);
    }
  
    if (dim == 2){
      return -1. * df[0];
    }
  
    return nan("1");
  }

  template<int dim>
  double LShape_NC2<dim>::value(
    const Point<dim>&     p,
    const unsigned int  /* component */
  ) const {
    Tensor<1, dim> df;
    if (type == classic){
      const LShape_SOL_classic<dim> f;
      df = f.gradient(p);
    } else if (type == benchmark){
      const LShape_SOL_benchmark<dim> f;
      df = f.gradient(p);
    }
  
    if (dim == 2){
      return +1. * df[1];
    }
  
    return nan("1");
  }

  template<int dim>
  double LShape_NC3<dim>::value(
    const Point<dim>&     p,
    const unsigned int  /* component */
  ) const {
    Tensor<1, dim> df;
    if (type == classic){
      const LShape_SOL_classic<dim> f;
      df = f.gradient(p);
    } else if (type == benchmark){
      const LShape_SOL_benchmark<dim> f;
      df = f.gradient(p);
    } else {
      Assert(false, ExcInternalError());
    }
  
    if (dim == 2){
      return +1. * df[0];
    }
  
    return nan("1");
  }

  template<int dim>
  double LShape_NC4<dim>::value(
    const Point<dim>&     p,
    const unsigned int  /* component */
  ) const {
    Tensor<1, dim> df;
    if (type == classic){
      const LShape_SOL_classic<dim> f;
      df = f.gradient(p);
    } else if (type == benchmark){
      const LShape_SOL_benchmark<dim> f;
      df = f.gradient(p);
    }
  
    if (dim == 2){
      return -1. * df[1];
    }
  
    return nan("1");
  }
  template<int dim>
  double LShape_NC_Diag<dim>::value(
    const Point<dim>&     p,
    const unsigned int  /* component */
  ) const {
    Tensor<1, dim> df;
    if (type == classic){
      const LShape_SOL_classic<dim> f;
      df = f.gradient(p);
    } else if (type == benchmark){
      const LShape_SOL_benchmark<dim> f;
      df = f.gradient(p);
    }
  
    if (dim == 2){
      return (+1. * df[0] + -1. * df[1]) / std::sqrt(2);
    }
  
    return nan("1");
  }
  
  // =========================================================
  // =========================================================
  //             LShape_Problem Implementation
  // =========================================================
  // =========================================================
  template<int dim>
  void LShape_Problem<dim>::print_solution(
  ) {
    Assert(dim > 1, ExcNotImplemented());
    Assert(dim < 4, ExcNotImplemented());
  }
  
  // Print the numeric solution at the boundary
  template<int dim>
  void LShape_Problem<dim>::print_numeric_solution(
  ) {
    Assert(dim > 1, ExcNotImplemented());
    Assert(dim < 4, ExcNotImplemented());
  }
  
  template<int dim>
  void LShape_Problem<dim>::print_numeric_solution(
      const double /* x */
  ) {
    Assert(dim > 1, ExcNotImplemented());
    Assert(dim < 3, ExcNotImplemented());
  }
  
  template<>
  void LShape_Problem<2>::print_solution(
  ) {
    const unsigned int N = 100;
  
    FullMatrix<double> B_x0(N, N), B_y0(N, N);
    FullMatrix<double> B_x1(N, N), B_y1(N, N);
    FullMatrix<double> B_x2(N, N), B_y2(N, N);
  
  
    // Declare output matrices for the corresponding patches
    FullMatrix<double> B0(N, N), B1(N, N), B2(N, N);
  
    for (unsigned int j = 0; j < N; j++){
      for (unsigned int i = 0; i < N; i++){
  
        B0(j, i) = sol_fcn.value( Point<2>( -1. + i/(N-1.), +0. + j/(N-1.) ));
        B1(j, i) = sol_fcn.value( Point<2>( +0. + i/(N-1.), +0. + j/(N-1.) ));
        B2(j, i) = sol_fcn.value( Point<2>( +0. + i/(N-1.), -1. + j/(N-1.) ));
  
  
        B_x0(i, j) = -1. + j/(N-1.); B_y0(i, j) = +0. + i/(N-1.);
        B_x1(i, j) = +0. + j/(N-1.); B_y1(i, j) = +0. + i/(N-1.);
        B_x2(i, j) = +0. + j/(N-1.); B_y2(i, j) = -1. + i/(N-1.);
  
      }
    }
    const std::string name0 = problem_out.degree.string() + "solution.0.dat";
    const std::string name1 = problem_out.degree.string() + "solution.1.dat";
    const std::string name2 = problem_out.degree.string() + "solution.2.dat";
  
  
    const std::string name_x0 = problem_out.degree.string() + "solution.x0.dat";
    const std::string name_x1 = problem_out.degree.string() + "solution.x1.dat";
    const std::string name_x2 = problem_out.degree.string() + "solution.x2.dat";
  
    const std::string name_y0 = problem_out.degree.string() + "solution.y0.dat";
    const std::string name_y1 = problem_out.degree.string() + "solution.y1.dat";
    const std::string name_y2 = problem_out.degree.string() + "solution.y2.dat";
  
  
    std::filebuf f0, f1, f2;
    std::filebuf f_x0, f_x1, f_x2, f_y0, f_y1, f_y2;
    f0.open(name0.c_str(), std::ios::out);
    f1.open(name1.c_str(), std::ios::out);
    f2.open(name2.c_str(), std::ios::out);
  
  
    f_x0.open(name_x0.c_str(), std::ios::out);
    f_x1.open(name_x1.c_str(), std::ios::out);
    f_x2.open(name_x2.c_str(), std::ios::out);
  
    f_y0.open(name_y0.c_str(), std::ios::out);
    f_y1.open(name_y1.c_str(), std::ios::out);
    f_y2.open(name_y2.c_str(), std::ios::out);
  
  
    std::ostream out0(&f0), out1(&f1), out2(&f2);
    std::ostream out_x0(&f_x0), out_x1(&f_x1), out_x2(&f_x2);
    std::ostream out_y0(&f_y0), out_y1(&f_y1), out_y2(&f_y2);
  
    B0.print_formatted(out0, 16, true, 1, "0");
    B1.print_formatted(out1, 16, true, 1, "0");
    B2.print_formatted(out2, 16, true, 1, "0");
  
  
    B_x0.print_formatted(out_x0, 16, true, 1, "0");
    B_x1.print_formatted(out_x1, 16, true, 1, "0");
    B_x2.print_formatted(out_x2, 16, true, 1, "0");
  
    B_y0.print_formatted(out_y0, 16, true, 1, "0");
    B_y1.print_formatted(out_y1, 16, true, 1, "0");
    B_y2.print_formatted(out_y2, 16, true, 1, "0");
  
  
    f0.close(); f1.close(); f2.close();
  
    f_x0.close(); f_x1.close(); f_x2.close();
    f_y0.close(); f_y1.close(); f_y2.close();
  }
  
  
  template<>
  void LShape_Problem<2>::print_numeric_solution(
  ) {
    // Seperately, print the gradient on the boundary
    const unsigned int N = 100;
  
    // Declare a container to store the values on the boundary
    FullMatrix<double> B_num(4, N);
    FullMatrix<double> B_exa(4, N);
  
    const auto& splines = tria.get_splines();
    const auto& mapping = tria.get_IPF();
    const auto& kv      = data.kv;
    const unsigned int n0 = kv[0].size() - 1;
    const unsigned int n1 = kv[1].size() - 1;
    for (unsigned int j = 0; j < N; j++){
      const Point<2> P0(kv[0][0], kv[1][0] + (kv[1][n1] -  kv[1][0]) * j / (N-1.));
      const Point<2> P1(kv[0][n0], kv[1][0] + (kv[1][n1] -  kv[1][0]) * j / (N-1.));
      const Point<2> P2(kv[0][0] + (kv[0][n0] -  kv[0][0]) * j / (N-1.),kv[1][0]);
      const Point<2> P3(kv[0][0] + (kv[0][n0] -  kv[0][0]) * j / (N-1.),kv[1][n1]);
  
      for (unsigned int i = 0; i < tria.n_active_splines(); i++){
        B_num(0, j) += solution[i] * splines[i] -> value(P0);
        B_num(1, j) += solution[i] * splines[i] -> value(P1);
        B_num(2, j) += solution[i] * splines[i] -> value(P2);
        B_num(3, j) += solution[i] * splines[i] -> value(P3);
      }
  
      B_exa(0, j) = sol_fcn.value(mapping.point_value(P0));
      B_exa(1, j) = sol_fcn.value(mapping.point_value(P1));
      B_exa(2, j) = sol_fcn.value(mapping.point_value(P2));
      B_exa(3, j) = sol_fcn.value(mapping.point_value(P3));
    }
  
    std::string name_num = problem_out.degree.string() 
                            + "solution.num.bdry.l" 
                            + std::to_string(tria.n_levels() - offset - 1) 
                            + ".dat";
    std::string name_exa = problem_out.degree.string() 
                            + "solution.exa.bdry.l" 
                            + std::to_string(tria.n_levels() - offset - 1) 
                            + ".dat";
    std::filebuf f0, f1;
    f0.open(name_num.c_str(), std::ios::out);
    f1.open(name_exa.c_str(), std::ios::out);
  
    std::ostream out0(&f0), out1(&f1);
    B_num.print_formatted(out0, 16, true, 1, "0");
    B_exa.print_formatted(out1, 16, true, 1, "0");
  
    f0.close(), f1.close();
  }
  
  
  template<>
  void LShape_Problem<2>::print_numeric_solution(
      const double x
  ) {
    const unsigned int N = 100;
    const unsigned int n_dofs = tria.n_active_splines();
    const auto& splines = tria.get_splines();
    const auto& mapping = tria.get_IPF();
    Vector<double> num_pi4(N);
    Vector<double> sol_pi4(N);
    for (unsigned int i = 0; i < N; i++){
      Point<2> P(x, 0. + i/(N-1.));
      Vector<double> Phat;
      mapping.vector_value(P, Phat);
      sol_pi4(i) = sol_fcn.value(Point<2>(Phat(0), Phat(1)));
      for (unsigned int j = 0; j < n_dofs; j++){
         num_pi4(i) += solution(j) * splines[j] -> value(P);
      }
    }
  
    std::string pi40 = problem_out.degree.string() + "sol_pi4.dat";
    std::string pi41 = problem_out.degree.string() + "num_pi4.dat";
  
    std::filebuf f0, f1;
    f0.open(pi40.c_str(), std::ios::out);
    f1.open(pi41.c_str(), std::ios::out);
  
    std::ostream pi40_out(&f0);
    std::ostream pi41_out(&f1);
  
    sol_pi4.print(pi40_out, 16);
    num_pi4.print(pi41_out, 16);
  }
  
  
  
  
  template<int dim>
  void LShape_Problem<dim>::run(){
  
    // For debugging purposes: Print out the solution
    print_solution();
    unsigned int nlr = 0;
    unsigned int old_level = 0; 
    cycle = 0;
    while ((tria.n_levels() - offset - 1) < (ref + 1)){
      nlr++; 
      this -> setup_system(); 
      this -> assemble_system();
      this -> solve();
      this -> compute_h1_error();
      this -> print_error();
      this -> output_system();
      this -> estimate_and_refine();
      if (tria.n_levels() - offset - 1 != old_level)
        nlr = 0;
      else if (nlr == 5){
        std::cout << "Too many refinements resulted in the same level, global refinement is enforced!" << std::endl;
        // enforce a global refinement step: 
        tria.coarsen_bezier_elements();
        tria.refine_global();
        tria.set_boundary_dofs();
        tria.refine_bezier_elements();
        tria.compute_extraction_operators();
        nlr = 0; 
      }
      cycle++;
    }
    problem_out.write_table_text(std::cout);
  }
  
  template<int dim>
  void LShape_Problem<dim>::solve(){
    std::cout << "Solving system ... " << std::endl;
  
    SolverControl            solver_control(750 * tria.n_active_splines(), H1 * 1e-4);
    SolverCG<Vector<double>> solver(solver_control);

    PreconditionJacobi<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix);

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
  
  template<int dim>
  void LShape_Problem<dim>::assemble_system(){
    const std::string name = problem_out.svg.string() 
                              + "l"
                              + std::to_string(tria.n_levels() - offset - 1);
    print_grid(name);

    std::cout << "Assembling system matrix ... " << std::endl;
  
    // Setup initial tables that store the bernstein values / grads / hessians.
    std::vector< unsigned int > degrees = tria.get_degree();
    for (unsigned int d = 0; d < dim; d++)
      degrees[d] = degrees[d]*degrees[d] + 1;
  
    TSValues<dim> ts_values( 
          &tria,
          degrees, 
          update_values |
          update_gradients |
          update_quadrature_points |
          update_JxW_values);
  
    TSFaceValues<dim> face_values( 
          &tria,
          degrees, 
          update_values |
          update_gradients |
          update_quadrature_points |
          update_normal_vectors |
          update_JxW_values);
  
    const unsigned int  dofs_per_cell = ts_values.n_dofs_per_cell();
    FullMatrix<double>  cell_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<double>      cell_rhs(dofs_per_cell);
  
    for (const auto& cell : tria.active_cell_iterators()){
      // Get the Bernstein values on the cell
      ts_values.reinit(cell); 
  
      // Reset the cell matrix
      cell_matrix       = 0; 
      cell_rhs          = 0;
  
      // Quadrature sum: 
      for (const unsigned int q_index : ts_values.quadrature_point_indices()){
        // Build the cell matrix: 
        const double dx = ts_values.JxW(q_index);
        const double f  = rhs_fcn.value(ts_values.quadrature_point(q_index));
        for (const unsigned int i : ts_values.dof_indices()){
          for(const unsigned int j : ts_values.dof_indices())
            cell_matrix(i,j) += ( ts_values.shape_grad(i, q_index) 
                                  * ts_values.shape_grad(j, q_index)
                                  * dx ); 
  
  
          cell_rhs(i) += (ts_values.shape_value(i, q_index) *
                            f * dx);
        } // for ( i )
      } // for ( q_index )
  
  
      if (cell -> at_boundary()){
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; f++){
          if (cell -> face(f) -> at_boundary()
                  && ( cell -> face(f) -> boundary_id() == Boundary::neumann_1 
                        || cell -> face(f) -> boundary_id() == Boundary::neumann_2
                        || cell -> face(f) -> boundary_id() == Boundary::neumann_diag)){
            face_values.reinit(cell, f);
            for (const unsigned int q : face_values.quadrature_point_indices()){
              for (const unsigned int i : face_values.dof_indices()){
                cell_rhs(i) += (sol_fcn.gradient(face_values.quadrature_point(q)) *
                                  face_values.normal_vector(q) *
                                  face_values.shape_value(i, q) *
                                  face_values.JxW(q));
              }
            }
          } // if (cell -> face(f) -> at_boundary())
        } // for ( f )
      } // if (cell -> at_boundary())
    
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
        const Function< dim >*
      > boundary_fcns =
          {{Boundary::dirichleth, &sol_fcn}};
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
  
  
    std::cout << " ... done!" << std::endl; 
  } // assemble system
  
  template<int dim>
  void LShape_Problem<dim>::setup_system(
  ) {
    system_rhs.reinit(tria.n_active_splines());
    solution.reinit(tria.n_active_splines());
  
    // To generate the raw structure of the sparsity matrix, we construct our own SparsityPattern
    // This will be done using the IEN_array from the TS_Triangulation. It already stores the 
    // information which spline has support on which cell. 
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
  
  
  template<int dim>
  void LShape_Problem<dim>::compute_h1_error()
  {
    std::cout << "Computing H1-error at Level " 
              << tria.n_levels() - offset - 1 
              << " with degree " 
              << data.max_degree() + order
              << " ... " << std::endl;
  
    std::vector<unsigned int> degrees = tria.get_degree(); 
    for (unsigned int d = 0; d < dim; d++)
      degrees[d] = degrees[d]*degrees[d] + 1;
  
    TSValues<dim> ts_values(
                &tria,
                degrees, 
                update_values |
                update_gradients |
                update_quadrature_points |
                update_JxW_values);
  
    // In case we have computed the error on this level alreday
    double h1 = 0;
    double l2 = 0;
    for (const auto& cell :
            tria.active_cell_iterators()){

      // Get the TSpline values on this cell
      ts_values.reinit(cell); 
  
      // Get the corresponding global dof indices
      std::vector<unsigned int> local_dof_indices =
              tria.get_IEN_array(cell);
  
      // Quadrature sum: 
      for (const unsigned int q_index :
             ts_values.quadrature_point_indices()){
        // Map the quadrature point from real cell to parametric cell
        const Point<dim> mapped_q = 
                ts_values.quadrature_point(q_index); 
  
        // Get the value of the approximation
        double u_diff = sol_fcn.value(mapped_q);
        for (const unsigned int i : ts_values.dof_indices())
          u_diff -= (solution(local_dof_indices[i]) *
                      ts_values.shape_value(i, q_index)); 
  
        // Build the gradient value at quadrature point: 
        Tensor<1, dim> grad_u_diff = sol_fcn.gradient(mapped_q);
        for (const unsigned int i : ts_values.dof_indices())
          grad_u_diff -= ( solution(local_dof_indices[i]) *
                           ts_values.shape_grad(i, q_index));
  
        const double dx = ts_values.JxW(q_index);
        h1 += (( grad_u_diff * grad_u_diff ) * dx) ;
        l2 += ((      u_diff * u_diff      ) * dx);
      } // for ( q_index )
    } // for ( cell )
    H1 = std::sqrt(l2 + h1);
    L2 = std::sqrt(l2);

    problem_out.add_values_to_table(L2, H1);
  
    std::cout << " ... done!" << std::endl; 
  } // compute_h1_error

  template<int dim>
  void LShape_Problem<dim>::estimate_and_refine(){
          // does nothing
  } 

  template<>
  void LShape_Problem<2>::estimate_and_refine(){
    constexpr unsigned int dim = 2;
    if (strategy == Strategy::uniform){
      // For uniform refinement no estimator needs to be invoked,
      // this saves some computations
      tria.coarsen_bezier_elements();
      tria.refine_global();
    } else if (strategy == Strategy::adaptive) {
      // declare a container to store cells to be marked
      std::vector< TriaIterator<CellAccessor<dim, dim>> > mark;
  
      // special case, if there is only one cell
      if (tria.n_active_cells() == 1) {
        mark.push_back(tria.begin_active());
      } else {
        Vector<double> local_residuals(tria.n_active_splines());
        std::map< 
            types::boundary_id,
            const Function<dim>*
            > neumann_data = 
            {{Boundary::neumann_1, &neumann_bc1}, 
              {Boundary::neumann_2, &neumann_bc2}, 
              {Boundary::neumann_diag, &neumann_bc_diag}};
        std::vector<unsigned int> degrees = tria.get_degree();
        for(unsigned int d = 0; d < dim; d++)
          degrees[d] = degrees[d]*degrees[d] + 1;

        ResidualEstimators::Poisson<dim>::estimate(
            &tria,
            degrees,
            solution,
            local_residuals,
            &rhs_fcn,
            neumann_data
        );
  
        tria.refine_fixed_number(local_residuals, 0.10);

      } // if ( special case )
  
      // Finally prepare the next step, by coarsening the bezier elements
      tria.coarsen_bezier_elements();
  
      // And then mark all cells in mark
      tria.set_refine_flags(mark);
  
      // Finally, perform refinement with the following call
      tria.execute_coarsening_and_refinement();
    }
  
    std::string name = problem_out.svg.string() + "l"
                        + std::to_string(tria.n_levels() - offset - 1);
    print_grid(name);
  
  
    // Before we switch to the Bezier grid, calculate information about
    // boundary dofs
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

  template<>
  void LShape_Problem<3>::estimate_and_refine(){
    if (tria.n_active_cells() == 1){
      tria.coarsen_bezier_elements();
      tria.refine_global();
    } else {
      std::vector< TriaIterator<CellAccessor<3,3>> > mark;
      for (const auto& cell : tria.active_cell_iterators()){
        const Point<2> on_line(-1, 0);
        for (unsigned int v = 0; v < 8; v++){
          const Point<3> vert = cell->vertex(v);
          if ((on_line + Point<2>(vert(0), vert(1))).norm() < 1e-15)  
            mark.push_back(cell);
        }
      }

      tria.coarsen_bezier_elements();
      tria.set_refine_flags(mark); 
      tria.execute_coarsening_and_refinement();
    }

  
    // Before we switch to the Bezier grid, calculate information about
    // boundary dofs
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
  
  
  
  template<int dim>
  void LShape_Problem<dim>::output_system()
  {
    std::cout << "Printing system matrix and rhs ... " << std::endl;
  
    const unsigned int level = tria.n_levels() - offset - 1; 
    const std::string l = std::to_string(level);
    const std::string name = problem_out.degree.string(); 
  
    // std::string matrix = name + "_mat.dat" ;
    // std::string vector = name + "_vec.dat" ;
    // std::string soluti = name + "_sol.dat" ;
    // 
    // std::filebuf mat, vec, sol;
    // mat.open(matrix.c_str(), std::ios::out); 
    // vec.open(vector.c_str(), std::ios::out); 
    // sol.open(soluti.c_str(), std::ios::out); 
  
    // std::ostream mat_out(&mat);
    // std::ostream vec_out(&vec);
    // std::ostream sol_out(&sol);
  
    // system_matrix.print_formatted(mat_out, 16, true, 1, "0");
    // system_rhs.print(vec_out, 16);
    // solution.print(sol_out, 16);
  
    // mat.close();
    // vec.close();
    // sol.close();
  
  
    // tria.printIPF(name, 8, true, true);
  
    // Print the IPF wireframe:
    // tria.print_grid(name, 8);
    tria.print_IPF_wireframe(name + "l" + l, 8, 2);
    tria.coarsen_bezier_elements();
    tria.print_IPF_wireframe(name + "l" + l, 8, 2);
    tria.refine_bezier_elements();
    
  
    // Print the function at degree phi = pi / 4, i.e. x = 1
    // in case of one parametric element
    print_numeric_solution( 1. );
  
    // Print the function at the boundary
    print_numeric_solution( );

    problem_out.write_table_text(); 
    problem_out.write_table_tex();
  
    std::cout << " ... output to files done!" << std::endl; 
  } // output_system
  
  
  
  
  
  template<int dim>
  void LShape_Problem<dim>::print_error(
  ) {
    std::vector<unsigned int> degrees = tria.get_degree(); 
    for (unsigned int d = 0; d < dim; d++)
      degrees[d] = degrees[d]*degrees[d] + 1;
    TSValues<dim> ts_values(
                &tria,
                degrees, 
                update_values |
                update_gradients |
                update_quadrature_points |
                update_JxW_values);

    const unsigned int nq = ts_values.n_quadrature_points_per_cell();
    const unsigned int nc = tria.n_active_cells();
    FullMatrix<double> B(nq, nc);
    unsigned int c = 0;
    for (const auto& cell : 
            tria.active_cell_iterators()){
      // Get the TSpline values on this cell
      ts_values.reinit(cell);
  
      // Get the corresponding global dof indices
      std::vector<unsigned int> local_dof_indices =
            tria.get_IEN_array(cell);
  
      // For each quadrature point ...
      for (const unsigned int q :
              ts_values.quadrature_point_indices()){
        // ... compute the value of the solution at the
        // quadrature point and ...
        double diff = sol_fcn.value(ts_values.quadrature_point(q));
  
        // ... build the difference with the numeric solution, to ...
        for (const unsigned int i :
                ts_values.dof_indices())
          diff -= ts_values.shape_value(i, q) 
                   * solution(local_dof_indices[i]);
  
        // ... store it in the output matrix:
        B(q, c) = std::fabs(diff);
      } // for ( q )
      c++;
    } // for ( cell )
  
    // Finally, print the matrix to a readable .dat file:
  
    // Firstly declare a name in the most complicated but retraceable way:
    const std::string l = 
            std::to_string(tria.n_levels() - offset - 1);
    const std::string name = problem_out.degree.string() 
                              + "lshape_o" 
                              + std::to_string(data.max_degree() + order) 
                              + "_l" + l
                              + "_cell_error.dat";
  
    // Then open a file to write contents into:
    std::filebuf out_file;
    out_file.open(name.c_str(), std::ios::out);
  
    // Then tell C++ to write in that file:
    std::ostream out(&out_file);
  
    // And finally print the contents to that file
    B.print_formatted(out, 16, true, 1, "0");
  
    // Free the output file for the system
    out_file.close();
  }
  
  
  // =========================================================
  // =========================================================
  //              Specializations
  // =========================================================
  //
  // =========================================================
  // =========================================================
  //              DataGenerator implementation 
  // =========================================================
  // =========================================================
  
  // dim = 2 {
  template<>
  DataGenerator<2, 2>::DataGenerator(
      const unsigned int n_elements,
      const bool no_c0_edges
   ){
    if (n_elements == 1) {
      // In case of 1 element, the solution is mirrored along the x = y axis,
      // and since for 1 parametric element we automatically have no c0 edges,
      // it is not needed in this setup.
  
      // Return 2D Data needed for half the LShape
      std::vector< Point< 2 + 1> >        cps =
        {
          Point<2 + 1>(-1.0, +0.0, 1.),
          Point<2 + 1>(+0.0, +0.0, 1.),
          Point<2 + 1>(-1.0, +1.0, 1.),
          Point<2 + 1>(+1.0, +1.0, 1.)
        };
      std::vector< std::vector<double> >  kvs =
           { {0., 0., 1., 1.}, {0., 0., 1., 1. } };
      std::vector< unsigned int >         deg = {1, 1};
  
      data = IPF_Data<2, 2>(cps, kvs, deg);
    } else if (n_elements == 2 && no_c0_edges){
      // Return 2D Data needed for LShape
      std::vector< Point< 2 + 1> >        cps =
        {
          Point<2 + 1>(-1.0, +0.0, 1.),
          Point<2 + 1>(+0.0, +0.0, 1.),
          Point<2 + 1>(+0.0, +0.0, 1.),
          Point<2 + 1>(+0.0, -1.0, 1.),
          Point<2 + 1>(-1.0, +1.0, 1.),
          Point<2 + 1>(+1.0, +1.0, 1.),
          Point<2 + 1>(+1.0, +1.0, 1.),
          Point<2 + 1>(+1.0, -1.0, 1.)
        };
      std::vector< std::vector<double> >  kvs =
        { {0., 0., 0., 1., 2., 2., 2.}, {0., 0., 1., 1. } };
      std::vector< unsigned int >         deg = {2, 1};
  
      data = IPF_Data<2, 2>(cps, kvs, deg);
    } else if (n_elements == 2 && !(no_c0_edges)) {
      // Return 2D Data needed for LShape
      std::vector< Point< 2 + 1> >        cps =
        {
          Point<2 + 1>(-1.0, +0.0, 1.),
          Point<2 + 1>(+0.0, +0.0, 1.),
          Point<2 + 1>(+0.0, -1.0, 1.),
          Point<2 + 1>(-1.0, +1.0, 1.),
          Point<2 + 1>(+1.0, +1.0, 1.),
          Point<2 + 1>(+1.0, -1.0, 1.)
        };
      std::vector< std::vector<double> >  kvs =
        { {0., 0., 1., 2., 2.}, {0., 0., 1., 1. } };
      std::vector< unsigned int >         deg = {1, 1};
  
      data = IPF_Data<2, 2>(cps, kvs, deg);
    } else if (n_elements == 3 && no_c0_edges) {
      // Return 2D Data needed for LShape
      std::vector< Point< 2 + 1> >        cps =
        {
          Point<2 + 1>(-1.0, +0.0, 1.), // P0
          Point<2 + 1>(+0.0, +0.0, 1.), // P1
          Point<2 + 1>(+0.0, +0.0, 1.), // P2
          Point<2 + 1>(+0.0, +0.0, 1.), // P3
          Point<2 + 1>(+0.0, -1.0, 1.), // P4
          Point<2 + 1>(-1.0, +1.0, 1.), // P5
          Point<2 + 1>(+0.0, +1.0, 1.), // P6
          Point<2 + 1>(+1.0, +1.0, 1.), // P7
          Point<2 + 1>(+0.0, +1.0, 1.), // P8
          Point<2 + 1>(+1.0, -1.0, 1.)  // P9
        };
      std::vector< std::vector<double> >  kvs =
        { {0., 0., 0., 1., 2., 3., 3., 3.}, {0., 0., 1., 1. } };
      std::vector< unsigned int >         deg = {2, 1};
  
      data = IPF_Data<2, 2>(cps, kvs, deg);
    } else if (n_elements == 4 && !no_c0_edges) {
      // Return 2D Data needed for LShape
      std::vector< Point< 2 + 1> >        cps =
        {
          Point<2 + 1>(-1.0, +0.0, 1.), // P0
          Point<2 + 1>(+0.0, +0.0, 1.), // P1
          Point<2 + 1>(+0.0, +0.0, 1.), // P2
          Point<2 + 1>(+0.0, +0.0, 1.), // P3
          Point<2 + 1>(+0.0, -1.0, 1.), // P4
          Point<2 + 1>(-1.0, +1.0, 1.), // P5
          Point<2 + 1>(+0.0, +1.0, 1.), // P6
          Point<2 + 1>(+1.0, +1.0, 1.), // P7
          Point<2 + 1>(+0.0, +1.0, 1.), // P8
          Point<2 + 1>(+1.0, -1.0, 1.)  // P9
        };
      std::vector< std::vector<double> >  kvs =
        { {0., 0., 1., 2., 3., 4., 4.}, {0., 0., 1., 1. } };
      std::vector< unsigned int >         deg = {1, 1};
  
      data = IPF_Data<2, 2>(cps, kvs, deg);
    }
  }
  // } dim = 2, DataGenerator
  
  // dim = 3 {
  template<>
  DataGenerator<3, 3>::DataGenerator(
      const unsigned int /* n_elements */,
      const bool         /* no_c0_edges */
  ){
      // In case of 1 element, the solution is mirrored along the x = y axis,
      // and since for 1 parametric element we automatically have no c0 edges,
      // it is not needed in this setup.
  
      // Return 2D Data needed for half the LShape
      std::vector< Point< 3 > >        cps =
        {
          Point<3>(-1.0, +0.0, -1.),
          Point<3>(+0.0, +0.0, -1.),
          Point<3>(-1.0, +1.0, -1.),
          Point<3>(+1.0, +1.0, -1.),
          Point<3>(-1.0, +0.0, +1.),
          Point<3>(+0.0, +0.0, +1.),
          Point<3>(-1.0, +1.0, +1.),
          Point<3>(+1.0, +1.0, +1.)
        };
      std::vector< double > w = {1., 1., 1., 1., 1., 1., 1., 1.};
  
      std::vector< std::vector<double> >  kvs =
           { {0., 0., 1., 1.}, {0., 0., 1., 1. }, {0., 0., 1., 1. } };
      std::vector< unsigned int >         deg = {1, 1, 1};
  
      data = IPF_Data<3, 3>(cps, w, kvs, deg);
  }
  // } dim = 3, Datagenerator
  
  
  // =========================================================
  // =========================================================
  //             LShape_RHS Implementation 
  // =========================================================
  // =========================================================
  
  // dim = 2 {
  template<>
  double LShape_RHS_benchmark<2>::value(
    const Point<2>&       p,
    const unsigned int    /* component */
  ) const  {
    const long double pi = 1.; //std::acos(-1);
    const long double sx = std::sin(pi * p(0));
    const long double sy = std::sin(pi * p(1));
    return 2 * pi * pi * sx * sy;
  }
  // } dim = 2
  
  // dim = 3 {
  template<>
  double LShape_RHS_benchmark<3>::value(
    const Point<3>&       p,
    const unsigned int    /* component */
  ) const  {
    // ToDo
    return 0 * p[0] * p[1];
  }
  // } dim = 3
  
  // dim = 2 {
  template<>
  double LShape_RHS_classic<2>::value(
    const Point<2>&       /* p */,
    const unsigned int    /* component */
  ) const  {
    return 0.;
  }
  // } dim = 2
  
  // dim = 3 {
  template<>
  double LShape_RHS_classic<3>::value(
    const Point<3>&       /* p */,
    const unsigned int    /* component */
  ) const  {
    return 0.;
  }
  // } dim = 3
  
  
  // =========================================================
  // =========================================================
  //             LShape_SOL Implementation 
  // =========================================================
  // =========================================================
  
  // dim = 2 {
  template<>
  double LShape_SOL_benchmark<2>::value(
    const Point<2>&     p,
    const unsigned int  /* component */
  ) const {
    const long double pi = 1.; //std::acos(-1);
    const long double sx = std::sin(pi * p(0));
    const long double sy = std::sin(pi * p(1));
  
    return (sx * sy);
  //  if (p.norm() == 0)
  //    return 0;
  //
  //  const long double pi = std::acos(-1);
  //  const long double degree  = std::atan2(p(1), p(0));
  //  const long double length  = std::pow(p.norm(), 2./3.);
  //
  //  return length * std::sin((2. * degree + pi) / 3.);
  }
  // } dim = 2
  
  // dim = 3 {
  template<>
  double LShape_SOL_benchmark<3>::value(
    const Point<3>&     p,
    const unsigned int  /* component */
  ) const {
    // ToDo
    return 0. * p[0] * p[1];
  }
  // } dim = 3
  
  // dim = 2 {
  template<>
  Tensor<1, 2> LShape_SOL_benchmark<2>::gradient(
    const Point<2>&     p,
    const unsigned int  /* component */
  ) const {
    const long double pi = 1.; //std::acos(-1);
    const long double sx = std::sin(pi * p(0));
    const long double sy = std::sin(pi * p(1));
    const long double cx = std::cos(pi * p(0));
    const long double cy = std::cos(pi * p(1));
  
    Tensor<1, 2> out;
    out[0] = pi * cx * sy;
    out[1] = pi * sx * cy;
  
    return out;
  }
  // } dim = 2
  
  // dim = 3 {
  template<>
  Tensor<1, 3> LShape_SOL_benchmark<3>::gradient(
    const Point<3>&     p,
    const unsigned int  /* component */
  ) const {
    // ToDo
    Tensor<1, 3> out;
    return out * p[0];
  }
  // } dim = 3
  
  // dim = 2 {
  template<>
  double LShape_SOL_classic<2>::value(
    const Point<2>&     p,
    const unsigned int  /* component */
  ) const { 
    if (p.norm() == 0)
      return 0;
  
    const long double pi = std::acos(-1);
    const long double degree  = std::atan2(p(1), p(0));
    const long double length  = std::pow(p.norm(), 2./3.);
    AssertIsFinite(pi);
    AssertIsFinite(degree);
    AssertIsFinite(length);
  
    return length * std::sin((2. * degree + pi) / 3.);
  }
  // } dim = 2
  
  // dim = 3 {
  template<>
  double LShape_SOL_classic<3>::value(
    const Point<3>&     p,
    const unsigned int  /* component */
  ) const { 
    if (p.norm() == 0)
      return 0;
  
    return 0.;
  }
  // } dim = 3
  
  // dim = 2 {
  template<>
  Tensor<1, 2> LShape_SOL_classic<2>::gradient(
    const Point<2>&     p,
    const unsigned int  /* component */
  ) const {
    const long double pi = std::acos(-1);
    const double degree  = std::atan2(p(1), p(0));
    const double length  = p(1) * p(1) + p(0) * p(0);
  
    const long double S = std::sin((2. * degree + pi) / 3.);
    const long double C = std::cos((2. * degree + pi) / 3.);
  
    Tensor<1, 2> out;
    if (length < 1e-14)
      return out;

    out[0] = 2. * (
                p(0) * S - p(1) * C
              ) / (3. * std::pow(length, 2./3.));
    out[1] = 2. * (
                p(1) * S + p(0) * C
              ) / (3. * std::pow(length, 2./3.));
    return out;
  }
  // } dim = 2
  
  // dim = 3 {
  template<>
  Tensor<1, 3> LShape_SOL_classic<3>::gradient(
    const Point<3>&     p,
    const unsigned int  /* component */
  ) const {
    Tensor<1, 3> out;
    return out * p[0];
  }
  // } dim = 3
  
}
