#include <ts_laplace.h>

template class TS_Laplace<2>;
template class TS_Laplace<3>;


template<int spacedim>
TS_Laplace<spacedim>::TS_Laplace(int ref, int order) : 
        ref(ref),                     // Number of refinements to be done
        tria(this->get_IPF_data()),   // Initialize the grid
        dof_handler(tria),            // Assign the dof handler -- This needs adjustment
        fe(this->get_IPF_data().max_degree()+order) // we use Bernstein Polynomials of degree 2
{
  IPF_Data<2, spacedim> data = get_IPF_data();
  unsigned int max_deg = data.max_degree();
  std::vector<unsigned int> degrees = data.deg; 
  for (unsigned int d = 0; d < 2 /*dim*/; d++)
    tria.degree_elevate(d, max_deg - degrees[d] + order);
}


template<int spacedim>
IPF_Data<2, spacedim> TS_Laplace<spacedim>::get_IPF_data()
{
  IPF_Data<2, spacedim> out; 

  Assert(false, ExcImpossibleInDim(spacedim));

  return out; 
}


template<int spacedim>
void TS_Laplace<spacedim>::get_IPF_data(
                         std::vector< std::vector< double > >& kv,
                         std::vector< Point<spacedim + 1> >& cps, 
                         std::vector<unsigned int>& deg)
{
  Assert(false, ExcImpossibleInDim(spacedim));
}

template<int spacedim>
void TS_Laplace<spacedim>::refine_grid()
{
  // coarsen the bezier elements for next refinement step
  tria.coarsen_bezier_elements();
  if (tria.n_levels() > 2){
    const Point<2> p(-0.5, -1.);
    for (const auto& cell : tria.active_cell_iterators()){
      RefinementCase<2> ref = RefinementCase<2>::cut_axis(cell -> level() % 2); 
      for (unsigned int n = 0; n < GeometryInfo<2>::vertices_per_cell; n++){
        if ( (p + cell -> vertex(n)).norm() < std::fabs(1e-10)){
          cell -> set_refine_flag(ref); 
          break; // The other vertices won't satisfy this condition
        } // if ( ... )
      } // for ( n < n_vertices )
    } // for ( cell )

    tria.execute_coarsening_and_refinement();
  } else {
    tria.refine_global(1); 
  }
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

  // retrieve the isoparametric function from the triangulation
  // This might be superfluous ... 
  parametric_mapping = tria.get_IPF();
}

template<int spacedim>
void TS_Laplace<spacedim>::setup_system(){
  // set dimensions for system:
  dof_handler.distribute_dofs(fe);

  system_rhs.reinit(tria.n_active_splines());
  solution.reinit(tria.n_active_splines());

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
  sparsity_pattern.reinit(tria.n_active_splines(), 
                            tria.n_active_splines(), 
                            tria.get_max_nnz_supports()
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

  constraints.clear(); 

  DoFTools::make_hanging_node_constraints(dof_handler, constraints); 

  Functions::ZeroFunction<2>   fcn(1);
  std::map<types::boundary_id, const Function<2>* > 
      boundary_functions{{0, &fcn}};
  VectorTools::project_boundary_values(dof_handler, 
                                        boundary_functions, 
                                        QGauss<1>(fe.degree + 1), 
                                        constraints 
                                        );
  
  constraints.close();
}

template<int spacedim>
void TS_Laplace<spacedim>::assemble_system(){
  std::cout << "Assembling system matrix ... " << std::endl;
  QGauss<2>   quadrature_formula(fe.degree + 1);
  FEValues<2> fe_values(fe,
                        quadrature_formula,
                        update_values | update_gradients | update_quadrature_points |
                        update_covariant_transformation |
                        update_JxW_values);


  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
  FullMatrix<double> bernstein_matrix(dofs_per_cell, dofs_per_cell); 
  Vector<double>     cell_rhs(dofs_per_cell);
  Vector<double>     bernstein_rhs(dofs_per_cell);


  const auto&  boundary_dofs = 
          tria.get_boundary_dofs();

  for (const auto& cell : dof_handler.active_cell_iterators()){
    // Get the Bernstein values on the cell
    fe_values.reinit(cell); 
    const std::vector< unsigned int > local_dof_indices =
            tria.get_IEN_array( cell );

    // Reset the cell matrix
    cell_matrix = 0; 
    bernstein_matrix = 0;
    cell_rhs    = 0;

    // Quadrature sum: 
    for (const unsigned int q_index : fe_values.quadrature_point_indices()){
      // Build the cell matrix: 
      double dx = fe_values.JxW(q_index); 
      for (const unsigned int i : fe_values.dof_indices())
        for(const unsigned int j : fe_values.dof_indices())
          cell_matrix(i,j) += fe_values.shape_grad(i, q_index) * 
                                fe_values.shape_grad(j, q_index) * 
                                dx; 

      for (const unsigned int i : fe_values.dof_indices())
        cell_rhs(i) += fe_values.shape_value(i, q_index) 
                        * 1. 
                        * dx; 

    } // for ( q_index )

  
    // Add the values to the system
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

template<int spacedim>
void TS_Laplace<spacedim>::output_system()
{
  std::cout << "Printing system matrix and rhs ... " << std::endl;

  int level = tria.n_levels() - 1; 
  std::string matrix = "out/ts_laplace_mat_" + std::to_string(level) + ".dat";
  std::string vector = "out/ts_laplace_rhs_" + std::to_string(level) + ".dat";
  std::string soluti = "out/ts_laplace_sol_" + std::to_string(level) + ".dat";
  
  std::filebuf mat, vec, sol;
  mat.open(matrix.c_str(), std::ios::out); 
  vec.open(vector.c_str(), std::ios::out); 
  sol.open(soluti.c_str(), std::ios::out); 

  std::ostream mat_out(&mat);
  std::ostream vec_out(&vec);
  std::ostream sol_out(&sol);

  system_matrix.print_formatted(mat_out, 5, true, 1, "0"); 
  system_rhs.print(vec_out, 5); 
  solution.print(sol_out, 5); 

  mat.close();
  vec.close();
  sol.close();

  // Print the Spline values
  unsigned int N1 = 201;
  unsigned int N2 = 201;

  std::vector< Point<2> > Points(N1*N2);
  IPF_Data<2, spacedim> data = get_IPF_data();
  const auto& Xi = data.kv; 
  double x0 = Xi[0][0];
  double x1 = Xi[0][Xi[0].size() - 1]; 
  double y0 = Xi[1][0];
  double y1 = Xi[1][Xi[1].size() - 1]; 

  unsigned int pos = 0;
  for (unsigned int i = 0; i < N1; i ++)
    for (unsigned int j = 0; j < N2; j++)
      Points[pos++] = Point<2>(x0 + i/(N1-1.) * (x1 - x0), 
                                y0 + j/(N2-1.) * (y1 - y0));

  std::string name; 
  if (spacedim == 2)
    name = "disc_";
  else 
    name = "hem2d_";

  tria.printIPF(Points, name, true, true);

  // Print the IPF wireframe:
  tria.print_IPF_wireframe("_" + name + std::to_string(tria.n_levels() - 1)); 

  std::cout << " ... done!" << std::endl; 
} // output_system

template<int spacedim>
void TS_Laplace<spacedim>::solve(){
  std::cout << "Solving system ... " << std::endl;

  SolverControl            solver_control(1000, 1e-6 * system_rhs.l2_norm());
  SolverCG<Vector<double>> solver(solver_control);
  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());

  std::cout << " ... done!" << std::endl; 
}

template<int spacedim>
void TS_Laplace<spacedim>::output_results(){
  DataOut<2> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");
  data_out.build_patches();

  std::ofstream output("out/ts_laplace_hem2d_solution.vtk");
  data_out.write_vtk(output);
}

template<int spacedim>
void TS_Laplace<spacedim>::run(){
  tria.refine_bezier_elements();
  tria.compute_extraction_operators();
  parametric_mapping = tria.get_IPF();

  this -> setup_system(); 
  this -> assemble_system();
  this -> solve();
  this -> output_system();
  for (int r = 0; r < ref; r++){
    this -> refine_grid(); 
    this -> setup_system(); 
    this -> assemble_system();
    this -> solve();
    this -> output_system();
  }
//  this -> output_results();
}


// Specializations
// Data for disc:
template<>
void TS_Laplace<2>::get_IPF_data(
                         std::vector< std::vector< double > >& kv,
                         std::vector< Point<2 + 1> >& cps, 
                         std::vector<unsigned int>& deg)
{
  // Define the vector of control points
  cps = std::vector< Point<3> >(9); 
  cps[0] = Point<3>(0,0,1);
  cps[1] = Point<3>(1,0,1);
  cps[2] = Point<3>(2,0,1);
  cps[3] = Point<3>(0,1,1);
  cps[4] = Point<3>(1,1,1);
//  cps[4] = Point<3>(1.5,1.7,0.7);
  cps[5] = Point<3>(2,1,1);
  cps[6] = Point<3>(0,2,1);
  cps[7] = Point<3>(1,2,1);
  cps[8] = Point<3>(2,2,1);

  // define the knot vectors
  kv = std::vector< std::vector< double > >(2);
  kv[0] = {0, 0, 0, 1, 1, 1}; 
  kv[1] = {0, 0, 0, 1, 1, 1}; 

  // Define the degree
  deg = {2, 2}; 

//  cps = std::vector< Point<3> >(18);
//  cps[0] = Point<3>(1.,0.,1.);
//  cps[1] = Point<3>(1./std::sqrt(2),-1./std::sqrt(2),1./std::sqrt(2));
//  cps[2] = Point<3>(0.,-1.,1.);
//  cps[3] = Point<3>(-1./std::sqrt(2),-1./std::sqrt(2),1./std::sqrt(2));
//  cps[4] = Point<3>(-1.,0.,1.);
//  cps[5] = Point<3>(-1./std::sqrt(2),1./std::sqrt(2),1./std::sqrt(2));
//  cps[6] = Point<3>(0.,1.,1.);
//  cps[7] = Point<3>(1./std::sqrt(2),1./std::sqrt(2),1./std::sqrt(2));
//  cps[8] = Point<3>(1.,0.,1.);
//  cps[9] = Point<3>(2.,0.,1.);
//  cps[10] = Point<3>(2./std::sqrt(2),-2./std::sqrt(2),1./std::sqrt(2));
//  cps[11] = Point<3>(0.,-2.,1.);
//  cps[12] = Point<3>(-2./std::sqrt(2),-2./std::sqrt(2),1./std::sqrt(2));
//  cps[13] = Point<3>(-2.,0.,1.);
//  cps[14] = Point<3>(-2./std::sqrt(2),2./std::sqrt(2),1./std::sqrt(2));
//  cps[15] = Point<3>(0.,2.,1.);
//  cps[16] = Point<3>(2./std::sqrt(2),2./std::sqrt(2),1./std::sqrt(2));
//  cps[17] = Point<3>(2.,0.,1.);
//
//  
//  // Define the knot vectors: 
//  kv = std::vector< std::vector<double>>(2); 
//  kv[1] = std::vector<double>({0,0,1,1});
//  kv[0] = std::vector<double>({0,0,0,1,1,2,2,3,3,4,4,4});
//
//  // Laslty, define the degree
//  deg = {2, 1};
}

template<>
void TS_Laplace<3>::get_IPF_data(
                         std::vector< std::vector< double > >& kv,
                         std::vector< Point<3 + 1> >& cps, 
                         std::vector<unsigned int>& deg)
{
  // Define the vector of control points: 
  cps = std::vector<Point<4>>(9); 
  std::vector<std::vector<double>> o_cp(9);

  double w = 1./std::sqrt(2.);
  o_cp[0] = {9.98, 0., 0., 1.};
  o_cp[1] = {w*9.98, w*9.98, 0., w};
  o_cp[2] = {0., 9.98, 0., 1.};
  o_cp[3] = {w*9.98, 0., w*9.98, w};
  o_cp[4] = {9.98/2., 9.98/2., 9.98/2., 1./2.};
  o_cp[5] = {0., w*9.98, w*9.98, w};
  o_cp[6] = {0., 0., 9.98, 1.};
  o_cp[7] = {0., 0., w*9.98, w};
  o_cp[8] = {0., 0., 9.98, 1.};

  for (int i = 0; i < 9; i++){
    for (int j = 0; j < 4; j++){
      cps[i](j) = o_cp[i][j];
    }
  }

  // Define the knot vectors: 
  kv = std::vector< std::vector<double>>(2); 
  std::vector<double> o_kv = {0, 0, 0, 1, 1, 1};

  kv[0] = o_kv;
  kv[1] = o_kv;

  // Laslty, define the degree
  deg = {2, 2};
}

template<>
IPF_Data<2, 3> TS_Laplace<3>::get_IPF_data()
{
  std::vector< Point<3 + 1> > cps;
  std::vector< std::vector< double > > kv; 
  std::vector< unsigned int > deg; 
  get_IPF_data(kv, cps, deg);

  IPF_Data<2, 3> out(cps, kv, deg); 
  return out;
}

template<>
IPF_Data<2, 2> TS_Laplace<2>::get_IPF_data()
{
  std::vector< Point<2 + 1> > cps;
  std::vector< std::vector< double > > kv; 
  std::vector< unsigned int > deg; 
  get_IPF_data(kv, cps,  deg);

  IPF_Data<2, 2> out(cps, kv, deg); 
  return out;
}


