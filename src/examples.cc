#include <examples.h>

std::vector< std::vector<double> > kv_smooth_lshape(){
  std::vector< std::vector<double> > out(2); 
  std::vector<double> kv0 = {0., 0., 0., .5, 1., 1., 1.};
  std::vector<double> kv1 = {0., 0., 0., 1., 1., 1.};

  out[0] = kv0;
  out[1] = kv1;

  return out;
}

std::vector< Point<2 + 1> > cp_smooth_lshape(){
  std::vector< Point<2 + 1> > out(12); 

  out[0] = Point<2 + 1>(0, 0, 1); 
  out[4] = Point<2 + 1>(-1, 0, 1); 
  out[8] = Point<2 + 1>(-2, 0, 1); 
  
  out[1] = Point<2 + 1>(0, 1, 1); 
  out[5] = Point<2 + 1>(-1, 2, 1); 
  out[9] = Point<2 + 1>(-2, 2, 1); 
  
  out[2] = Point<2 + 1>(1, 1.5, 1); 
  out[6] = Point<2 + 1>(1, 4, 1); 
  out[10] = Point<2 + 1>(1, 5, 1); 

  out[3] = Point<2 + 1>(3, 1.5, 1); 
  out[7] = Point<2 + 1>(3, 4, 1); 
  out[11] = Point<2 + 1>(3, 5, 1); 

  return out; 
}

std::vector< std::vector<double> > kv_circle(){
  std::vector< std::vector<double> > out(1);
  std::vector< double > kv2(12);

  kv2[0] = 0;
  kv2[1] = 0;
  kv2[2] = 0;
  kv2[3] = 1;
  kv2[4] = 1;
  kv2[5] = 2;
  kv2[6] = 2;
  kv2[7] = 3;
  kv2[8] = 3;
  kv2[9] = 4;
  kv2[10] = 4;
  kv2[11] = 4;

  out[0] = kv2;

  return out;
}

std::vector< std::vector<double> > kv_disk(){
  std::vector< std::vector<double> > out(2);
  std::vector< double > kv1(4);
  std::vector< double > kv2(12);

  kv1[0] = 0;
  kv1[1] = 0;
  kv1[2] = 1;
  kv1[3] = 1;

  kv2[0] = 0;
  kv2[1] = 0;
  kv2[2] = 0;
  kv2[3] = 1;
  kv2[4] = 1;
  kv2[5] = 2;
  kv2[6] = 2;
  kv2[7] = 3;
  kv2[8] = 3;
  kv2[9] = 4;
  kv2[10] = 4;
  kv2[11] = 4;

  out[1] = kv1;
  out[0] = kv2;

  return out;
}

std::vector< std::vector<double> > kv_cylinder(){
  std::vector< std::vector<double> > out(3);
  std::vector< double > kv1(4);
  std::vector< double > kv2(12);
  std::vector< double > kv3(6);
  kv1[0] = 0;
  kv1[1] = 0;
  kv1[2] = 1;
  kv1[3] = 1;

  kv2[0] = 0;
  kv2[1] = 0;
  kv2[2] = 0;
  kv2[3] = 1;
  kv2[4] = 1;
  kv2[5] = 2;
  kv2[6] = 2;
  kv2[7] = 3;
  kv2[8] = 3;
  kv2[9] = 4;
  kv2[10] = 4;
  kv2[11] = 4;

  kv3[0] = 0;
  kv3[1] = 0;
  kv3[2] = 0;
  kv3[3] = 1;
  kv3[4] = 1;
  kv3[5] = 1;

  out[1] = kv1;
  out[0] = kv2;
  out[2] = kv3;

  return out;
}

std::vector< Point<3> > cp_circle(){
  std::vector< Point<3> > out(9);

  out[0] = Point<3>(0.,1.,1.);
  out[1] = Point<3>(1./std::sqrt(2),1./std::sqrt(2),1./std::sqrt(2));
  out[2] = Point<3>(1.,0.,1.);
  out[3] = Point<3>(1./std::sqrt(2),-1./std::sqrt(2),1./std::sqrt(2));
  out[4] = Point<3>(0.,-1.,1.);
  out[5] = Point<3>(-1./std::sqrt(2),-1./std::sqrt(2),1./std::sqrt(2));
  out[6] = Point<3>(-1.,0.,1.);
  out[7] = Point<3>(-1./std::sqrt(2),1./std::sqrt(2),1./std::sqrt(2));
  out[8] = Point<3>(0.,1.,1.);


  return out;
}


std::vector< Point<3> > cp_disk(){
  std::vector< Point<3> > out(18);

  out[0] = Point<3>(1.,0.,1.);
  out[1] = Point<3>(1./std::sqrt(2),-1./std::sqrt(2),1./std::sqrt(2));
  out[2] = Point<3>(0.,-1.,1.);
  out[3] = Point<3>(-1./std::sqrt(2),-1./std::sqrt(2),1./std::sqrt(2));
  out[4] = Point<3>(-1.,0.,1.);
  out[5] = Point<3>(-1./std::sqrt(2),1./std::sqrt(2),1./std::sqrt(2));
  out[6] = Point<3>(0.,1.,1.);
  out[7] = Point<3>(1./std::sqrt(2),1./std::sqrt(2),1./std::sqrt(2));
  out[8] = Point<3>(1.,0.,1.);
  out[9] = Point<3>(2.,0.,1.);
  out[10] = Point<3>(2./std::sqrt(2),-2./std::sqrt(2),1./std::sqrt(2));
  out[11] = Point<3>(0.,-2.,1.);
  out[12] = Point<3>(-2./std::sqrt(2),-2./std::sqrt(2),1./std::sqrt(2));
  out[13] = Point<3>(-2.,0.,1.);
  out[14] = Point<3>(-2./std::sqrt(2),2./std::sqrt(2),1./std::sqrt(2));
  out[15] = Point<3>(0.,2.,1.);
  out[16] = Point<3>(2./std::sqrt(2),2./std::sqrt(2),1./std::sqrt(2));
  out[17] = Point<3>(2.,0.,1.);


  return out;
}

std::vector< Point<4> > cp_cylinder(){
  std::vector< std::vector<double> > out(54);

  double h1 = 4;
  double h2 = 2.5;
  double h3 = 1;
  double w = 1./std::sqrt(2);

  out[0] = {0.,1., h1, 1.};
  out[1] = {w , w, h1*w, w};
  out[2] = {1., 0., h1, 1.};
  out[3] = {w , -w, h1*w, w};
  out[4] = {0.,-1., h1, 1.};
  out[5] = {-w, -w , h1*w, w};
  out[6] = {-1., 0., h1, 1.};
  out[7] = {-w, w, h1*w, w};
  out[8] = {0., 1., h1, 1.};

  out[9] = {0.,2., h1, 1.};
  out[10] = {2*w , 2*w, h1*w, w};
  out[11] = {2., 0., h1, 1.};
  out[12] = {2*w , -2*w, h1*w, w};
  out[13] = {0.,-2., h1, 1.};
  out[14] = {-2*w, -2*w , h1*w, w};
  out[15] = {-2., 0., h1, 1.};
  out[16] = {-2*w, 2*w, h1*w, w};
  out[17] = {0., 2., h1, 1.};


  out[18] = {0.,1., h2, 1.};
  out[19] = {w , w, h2*w, w};
  out[20] = {1., 0., h2, 1.};
  out[21] = {w , -w, h2*w, w};
  out[22] = {0.,-1., h2, 1.};
  out[23] = {-w, -w , h2*w, w};
  out[24] = {-1., 0., h2, 1.};
  out[25] = {-w, w, h2*w, w};
  out[26] = {0., 1., h2, 1.};

  out[27] = {0.,2., h2, 1.};
  out[28] = {2*w , 2*w, h2*w, w};
  out[29] = {2., 0., h2, 1.};
  out[30] = {2*w , -2*w, h2*w, w};
  out[31] = {0.,-2., h2, 1.};
  out[32] = {-2*w, -2*w , h2*w, w};
  out[33] = {-2., 0., h2, 1.};
  out[34] = {-2*w, 2*w, h2*w, w};
  out[35] = {0., 2., h2, 1.};


  out[36] = {0.,1., h3, 1.};
  out[37] = {w , w, h3*w, w};
  out[38] = {1., 0., h3, 1.};
  out[39] = {w , -w, h3*w, w};
  out[40] = {0.,-1., h3, 1.};
  out[41] = {-w, -w , h3*w, w};
  out[42] = {-1., 0., h3, 1.};
  out[43] = {-w, w, h3*w, w};
  out[44] = {0., 1., h3, 1.};

  out[45] = {0.,2., h3, 1.};
  out[46] = {2*w , 2*w, h3*w, w};
  out[47] = {2., 0., h3, 1.};
  out[48] = {2*w , -2*w, h3*w, w};
  out[49] = {0.,-2., h3, 1.};
  out[50] = {-2*w, -2*w , h3*w, w};
  out[51] = {-2., 0., h3, 1.};
  out[52] = {-2*w, 2*w, h3*w, w};
  out[53] = {0., 2., h3, 1.};


  std::vector< Point<4> > o(54);

  for (unsigned int i = 0; i < 54; i++){
    for (unsigned int j = 0; j < 4; j++){
      o[i](j) = out[i][j];
    }
  }

  return o;
}

std::vector< std::vector<double>> kv_hemisphere()
{
  std::vector< std::vector<double>> out(2); 
  std::vector<double> o = {0, 0, 0, 1, 1, 1};

  out[0] = o;
  out[1] = o;

  return out;
}

std::vector< Point<4, double> > cp_hemisphere()
{
  std::vector<Point<4>> out(9); 
  std::vector<std::vector<double>> o(9);

  double w = 1./std::sqrt(2.);
  o[0] = {9.98, 0., 0., 1.};
  o[1] = {w*9.98, w*9.98, 0., w};
  o[2] = {0., 9.98, 0., 1.};
  o[3] = {w*9.98, 0., w*9.98, w};
  o[4] = {9.98/2., 9.98/2., 9.98/2., 1./2.};
  o[5] = {0., w*9.98, w*9.98, w};
  o[6] = {0., 0., 9.98, 1.};
  o[7] = {0., 0., w*9.98, w};
  o[8] = {0., 0., 9.98, 1.};

  for (int i = 0; i < 9; i++){
    for (int j = 0; j < 4; j++){
      out[i](j) = o[i][j];
    }
  }

  return out;
}

void print_example(std::string ex, int ref, int order){

  if (ex == "disc") {
      std::cout << " Running TSpline demo for disk\n ";
      const unsigned int dim = 2;
      const unsigned int spacedim = 2;

      std::vector< Point<spacedim+1,double> > cp;
      std::vector< std::vector<double> > kv;

      unsigned int N1 = 101;
      unsigned int N2 = 51;

      cp = cp_disk();
      kv = kv_disk();

      std::vector< unsigned int > deg = {2, 1};

      TS_Triangulation<dim,spacedim> T(kv,
                                   deg,
                                   cp);

      std::vector< Point<dim> > Points(N1*N2);

      unsigned int pos = 0;
      for (unsigned int i = 0; i < N1; i ++){
        for (unsigned int j = 0; j < N2; j++){
          Points[pos++] = Point<dim>(i/(N1-1.) * 4., j/(N2-1.) );
        }
      }

      const Point<dim> p0(-2., 0.);
      const Point<dim> p1(-2., -1.);
      T.degree_elevate(1, order); 
      T.degree_elevate(0, order + 1);
//      T.degree_elevate_global(order);

      T.printIPF(Points, "", false);
      T.print_grid("_disc_0");
      T.print_IPF_wireframe("_disc_0");
      for (int r = 0; r < ref; r++){
        // condition for refinement...
        for (auto& cell : T.active_cell_iterators())
          for (unsigned int n = 0; n < GeometryInfo<dim>::vertices_per_cell; n++ )
            if ((cell -> vertex(n) + p1).norm() < std::fabs(1e-16))
              cell -> set_refine_flag(RefinementCase<dim>::cut_axis(cell->level()%dim));
          
        T.execute_coarsening_and_refinement();
        T.print_grid("_disc_" + std::to_string(r+1));
        T.print_IPF_wireframe("_disc_" + std::to_string(r+1));
      }
    }
  else if (ex == "hem2d_glob") {
      std::cout << " Running TSpline demo for hemisphere\n ";
      const unsigned int dim = 2;
      const unsigned int spacedim = 3;

      std::vector< Point<spacedim+1,double> > cp;
      std::vector< std::vector<double> > kv;

      unsigned int N1 = 101;
      unsigned int N2 = 51; 

      cp = cp_hemisphere();
      kv = kv_hemisphere();

      std::vector< unsigned int > deg = {2, 2};

      TS_Triangulation<dim,spacedim> T(kv,
                                   deg,
                                   cp);


      std::vector< Point<dim> > Points(N1*N2);

      unsigned int pos = 0;
      for (unsigned int i = 0; i < N1; i ++){
        for (unsigned int j = 0; j < N2; j++){
          Points[pos++] = Point<dim>(i/(N1-1.) * 1., j/(N2-1.) );
        }
      }

      GridOutFlags::Svg svg_flags;
      svg_flags.coloring = GridOutFlags::Svg::Coloring::none;
//      svg_flags.label_level_number = true;
//      svg_flags.label_cell_index = true;

      std::string name = "out/grid_hem2d_glob_" + std::to_string(ref) + ".svg";
      std::ofstream out(name);
      GridOut       grid_out; 
      grid_out.set_flags(svg_flags);


      // Refinement: 
      T.degree_elevate_global(order);
      T.refine_global(ref);


      T.print_IPF_wireframe(ex);
      T.print_grid(ex);

      grid_out.write_svg(T, out);

      T.printIPF(Points, "", false);
    }
  else if (ex == "hem2d_fsnglrty") {
      std::cout << " Running TSpline demo for hemisphere\n ";
      const unsigned int dim = 2;
      const unsigned int spacedim = 3;

      std::vector< Point<spacedim+1,double> > cp;
      std::vector< std::vector<double> > kv;

      unsigned int N1 = 101;
      unsigned int N2 = 51; 

      cp = cp_hemisphere();
      kv = kv_hemisphere();

      std::vector< unsigned int > deg = {2, 2};

      TS_Triangulation<dim,spacedim> T(kv,
                                   deg,
                                   cp);

      T.degree_elevate_global(order); 

      T.print_IPF_wireframe("_level_0");
      T.print_grid("_level_0");

      std::vector< Point<dim> > Points(N1*N2);

      unsigned int pos = 0;
      for (unsigned int i = 0; i < N1; i ++){
        for (unsigned int j = 0; j < N2; j++){
          Points[pos++] = Point<dim>(i/(N1-1.) * 1., j/(N2-1.) );
        }
      }

      GridOutFlags::Svg svg_flags;
      svg_flags.coloring = GridOutFlags::Svg::Coloring::none;
      svg_flags.label_level_number = true;
      svg_flags.label_cell_index = true;

      std::string name = "out/grid_hem2d_fsnglrty_" + std::to_string(0) + ".svg";
      std::ofstream out(name);
      GridOut       grid_out; 
      grid_out.set_flags(svg_flags);

      grid_out.write_svg(T, out);

//      const Point<dim> p(-0.5, -1.);
      unsigned int n_vertices = GeometryInfo<dim>::vertices_per_cell;

      T.begin_active(0) -> set_refine_flag(RefinementCase<dim>::cut_x);
      T.execute_coarsening_and_refinement();
      T.print_IPF_wireframe("_level_1");
      T.print_grid("_level_1"); 
      name = "out/grid_hem2d_fsnglrty_" + std::to_string(1) + ".svg";
      std::ofstream out1(name);
      grid_out.write_svg(T, out1);
      // Refinement: 
      for (int i = 1; i < ref; i ++){
        for (const auto& cell : T.Triangulation<dim>::active_cell_iterators()){
        RefinementCase<dim> ref = RefinementCase<dim>::cut_axis(cell -> level() % dim); 
          for (unsigned int n = 0; n < n_vertices; n++){
            const auto& p = cell -> vertex(n); 
            if ( std::abs(p(dim-1) - 1) < 1e-16){
              cell -> set_refine_flag(ref); 
              break;
            }
          }
        }

        T.execute_coarsening_and_refinement();
        T.print_IPF_wireframe("_level_" + std::to_string(i + 1));
        T.print_grid("_level_" + std::to_string(i+1));
        T.print_bezier_grid("_level_" + std::to_string(i+1));

        name = "out/grid_hem2d_fsnglrty_" + std::to_string(i+1) + ".svg";
        std::ofstream out(name);

        grid_out.write_svg(T, out);

        T.refine_bezier_elements(); 
        name = "out/grid_hem2d_fsnglrty_bezier_" + std::to_string(i+1) + ".svg";
        
        std::ofstream out2(name);

        grid_out.write_svg(T, out2);
        T.coarsen_bezier_elements(); 
      }

      T.printIPF(Points, "", false);
    }
  else if (ex == "hem2d_csnglrty"){
      std::cout << " Running TSpline demo for hemisphere\n ";
      const unsigned int dim = 2;
      const unsigned int spacedim = 3;

      std::vector< Point<spacedim+1,double> > cp;
      std::vector< std::vector<double> > kv;

      unsigned int N1 = 101;
      unsigned int N2 = 51; 

      cp = cp_hemisphere();
      kv = kv_hemisphere();

      std::vector< unsigned int > deg = {2, 2};
      std::vector< unsigned int > sparsity_max_rows;
      std::vector< unsigned int > splines_on_level;
      std::vector< unsigned int > cells_on_level;

      TS_Triangulation<dim,spacedim> T(kv,
                                   deg,
                                   cp);

      T.degree_elevate_global(order); 

      std::vector< Point<dim> > Points(N1*N2);
      unsigned int pos = 0;
      for (unsigned int i = 0; i < N1; i ++){
        for (unsigned int j = 0; j < N2; j++){
          Points[pos++] = Point<dim>(i/(N1-1.) * 1., j/(N2-1.) );
        }
      }

      T.print_IPF_wireframe("_level_0");
      T.print_grid("_level_0");
      T.printIPF(Points, "hem2d_csnglrty_", true, false);

      GridOutFlags::Svg svg_flags;
      svg_flags.coloring = GridOutFlags::Svg::Coloring::none;
      svg_flags.label_level_number = true;
      svg_flags.label_cell_index = true;

      std::string name = "out/grid_hem2d_csnglrty_" + std::to_string(0) + ".svg";
      std::ofstream out(name);
      GridOut       grid_out; 
      grid_out.set_flags(svg_flags);

      grid_out.write_svg(T, out);

      const Point<dim> p(-0.5, -1.);
      unsigned int n_vertices = GeometryInfo<dim>::vertices_per_cell;

      T.begin_active(0) -> set_refine_flag(RefinementCase<dim>::cut_x);
      T.execute_coarsening_and_refinement();
      T.print_IPF_wireframe("_level_1");
      T.print_grid("_level_1"); 
      T.printIPF(Points, "hem2d_csnglrty_", true, false);

      T.refine_bezier_elements(); 
      SparsityPattern sp(T.n_active_splines(), T.n_active_splines(), T.n_active_splines()); 
      const auto& IEN = T.get_IEN_array();

      for (const auto& [_, arr] : IEN)
        for (unsigned int i : arr)
          for (unsigned int j : arr)
            sp.add(i, j); 

      sp.compress();
      sparsity_max_rows.push_back(sp.max_entries_per_row());
      splines_on_level.push_back(T.n_active_splines());
      cells_on_level.push_back(T.n_active_cells());

      std::string sn = "out/hem2d_csnglrty_sparsity_pattern_" + std::to_string(1) + ".svg";
      std::ofstream so(sn);

      sp.print_svg(so);
      T.coarsen_bezier_elements();

      name = "out/grid_hem2d_csnglrty_" + std::to_string(1) + ".svg";
      std::ofstream out1(name);
      grid_out.write_svg(T, out1);
      // Refinement: 
      for (int i = 1; i < ref; i ++){
        for (const auto& cell : T.Triangulation<dim>::active_cell_iterators()){
        RefinementCase<dim> ref = RefinementCase<dim>::cut_axis(cell -> level() % dim); 
          for (unsigned int n = 0; n < n_vertices; n++){
            if ( (p + cell -> vertex(n)).norm() < 1e-16){
              cell -> set_refine_flag(ref); 
              break;
            }
          }
        }

        T.execute_coarsening_and_refinement();
        T.print_IPF_wireframe("_level_" + std::to_string(i + 1));
        T.print_grid("hem2d_csnglrty_");
        T.print_bezier_grid("hem2d_csnglrty_");
        T.printIPF(Points, "hem2d_csnglrty_", true, false);
        

        name = "out/grid_hem2d_csnglrty_" + std::to_string(i+1) + ".svg";
        std::ofstream out(name);

        grid_out.write_svg(T, out);

        T.refine_bezier_elements(); 

        T.compute_extraction_operators();
        name = "out/grid_hem2d_csnglrty_bezier_" + std::to_string(i+1) + ".svg";
        
        std::ofstream out2(name);

        grid_out.write_svg(T, out2);

        SparsityPattern sparsity_pattern(T.n_active_splines(), T.n_active_splines(), T.n_active_splines()); 
        const auto& IEN_array = T.get_IEN_array();

        for (const auto& [_, arr] : IEN_array)
          for (unsigned int i : arr)
            for (unsigned int j : arr)
              sparsity_pattern.add(i, j); 

        sparsity_pattern.compress();
        sparsity_max_rows.push_back(sparsity_pattern.max_entries_per_row());
        splines_on_level.push_back(T.n_active_splines());
        cells_on_level.push_back(T.n_active_cells());

        std::string sparsity_name = "out/hem2d_csnglrty_sparsity_pattern_" + std::to_string(i+1) + ".svg";
        std::ofstream sparsity_out = std::ofstream(sparsity_name);
        sparsity_pattern.print_svg(sparsity_out);

        T.coarsen_bezier_elements(); 

        T.test_bezier();
      }
      T.printIPF(Points, "", false);
      std::cout << "lev : nnz | n_ts | n_cels " << std::endl;
      unsigned int N = sparsity_max_rows.size(); 
      for (unsigned int i = 0; i < N; i++){
        std::cout << i+1 << " : " << sparsity_max_rows[i] 
                         << " | " << splines_on_level[i] 
                         << " | " << cells_on_level[i] 
                         << std::endl;
      }
    }
  else if (ex == "hem3d_basic") {
      std::cout << " Running TSpline demo for 3D hemisphere\n ";
      const unsigned int dim = 3;
      const unsigned int spacedim = 3;

      std::vector<Point<spacedim+1, double>> cp; 
      std::vector<std::vector<double>> kv; 

      unsigned int N1 = 51;
      unsigned int N2 = 26;
      unsigned int N3 = 11; 

      cp = cp_hemisphere(); 
      kv = kv_hemisphere();

      Point<spacedim+1> P0 = 0.04*Point<spacedim+1, double>::unit_vector(0);
      Point<spacedim+1> P1 = 0.04*Point<spacedim+1, double>::unit_vector(1);
      Point<spacedim+1> P2 = 0.04*Point<spacedim+1, double>::unit_vector(2);
      
      double w = 1./std::sqrt(2.); 
      cp.push_back(cp[0] + P0); 
      cp.push_back(cp[1] + w*P0 + w*P2); 
      cp.push_back(cp[2] + P2); 
      cp.push_back(cp[3] + w*P0 + w*P1); 
      cp.push_back(cp[4] + 1./2.*(P0 + P1 + P2)); 
      cp.push_back(cp[5] + w*P2); 
      cp.push_back(cp[6] + P1); 
      cp.push_back(cp[7] + w*P1 + w*P2); 
      cp.push_back(cp[8] + P2); 

      kv.push_back({0,0,1,1});
      
      std::vector< unsigned int > deg = {2, 2, 1};

      TS_Triangulation<dim,spacedim> T(kv,
                                   deg,
                                   cp);

      std::vector< Point<dim> > Points(N1*N2*N3);

      unsigned int pos = 0;
      for (unsigned int x = 0; x < N1; x++){
        for (unsigned int y = 0; y < N2; y++){
          for (unsigned int z = 0; z < N3; z++) {
            Points[pos++] = Point<dim>(x/(N1-1.) * 1., y/(N2-1.) * 1., z/(N3-1.) * 1.);
          }
        }
      }
      // T.printIPF(Points);

      // Try some refinement: 
      for (int i = 0; i < ref; i ++){
        std::cout << "Refining grid in iteration " << i << "... " << std::endl;
        RefinementCase<dim> ref = RefinementCase<dim>::cut_axis(i % dim); 
        T.Triangulation<dim>::begin_active(i) -> set_refine_flag(ref); 
        T.execute_coarsening_and_refinement();

        std::cout << "Saving intermediate results, i.e. grid in knot domain and parametric grid ... " << std::endl;
        T.print_IPF_vertices("_level_" + std::to_string(i + 1));
        T.print_grid("_level_" + std::to_string(i+1));
        T.print_bezier_grid("_level_" + std::to_string(i+1));
        T.test_bezier();
      }
      T.printIPF(Points, "", false);

      GridOut       grid_out; 
      T.print_grid(ex);
      std::ofstream test("test.vtk");
      grid_out.write_vtk(T, test);
      std::cout << "Grid written to grid.*" << std::endl;
    }
  else if (ex == "smooth_lshape") {
    const std::vector< std::vector<double> >& kv = kv_smooth_lshape();
    const std::vector< Point<2 + 1, double> >& cp = cp_smooth_lshape();
    std::vector< unsigned int > deg = {2, 2}; 
    
    const unsigned int dim = 2; 
    const unsigned int spacedim = 2; 

    TS_Triangulation<dim, spacedim> T(kv, deg, cp);

    unsigned int N1 = 101; 
    unsigned int N2 = 51;

    std::vector< Point<dim> > Points(N1*N2);

    unsigned int pos = 0;
    for (unsigned int i = 0; i < N1; i ++){
      for (unsigned int j = 0; j < N2; j++){
        Points[pos++] = Point<dim>(i/(N1-1.) * 1., j/(N2-1.) );
      }
    }
   
    GridOutFlags::Svg svg_flags;
    svg_flags.coloring = GridOutFlags::Svg::Coloring::none;
    svg_flags.label_level_number = false;
    svg_flags.label_cell_index = false;

    std::string name = "out/grid_smooth_lshape_" + std::to_string(0) + ".svg";
    std::ofstream out(name);
    GridOut       grid_out; 
    grid_out.set_flags(svg_flags);

    grid_out.write_svg(T, out);

    T.degree_elevate_global(order); 
    T.printIPF(Points, "_smooth_lshape_0", false);
    T.print_IPF_wireframe("_smooth_lshape_0");
    T.print_grid("_smooth_lshape_0");

    Point<spacedim> p(0.5, 0.); 
    for (int i = 0; i < ref; i++){
      for (auto& cell : T.active_cell_iterators())
        if (cell -> point_inside(p))
          cell -> set_refine_flag(RefinementCase<dim>::cut_axis(cell -> level() % dim ));
    
      T.execute_coarsening_and_refinement();
      T.print_IPF_wireframe("_smooth_lshape_" + std::to_string(i + 1));
      T.print_grid("_smooth_lshape_" + std::to_string(i+1));

      name = "out/grid_smooth_lshape_" + std::to_string(i+1) + ".svg";
      std::ofstream out(name);

      grid_out.write_svg(T, out);

      T.refine_bezier_elements(); 
      name = "out/grid_smooth_lshape_bezier_" + std::to_string(i+1) + ".svg";
      
      std::ofstream out2(name);

      grid_out.write_svg(T, out2);
      T.coarsen_bezier_elements(); 
    }
    std::cout << "Example done..." << std::endl;
  } else if (ex == "ts_laplace"){
    TS_Laplace<2> ts_laplace(ref, order);
    ts_laplace.run(); 
  } else if ( ex == "benchmark"){
      Poisson_Benchmark::Poisson_Benchmark benchmark_test(ref,  order);
      benchmark_test.run();
  } else if ( ex == "lshape_ac"){
    LShape::LShape_Problem<2> l(
        ref, order,
        LShape::Specifier::classic,
        LShape::Strategy::adaptive);
    l.run();
  } else if ( ex == "lshape_uc"){
    LShape::LShape_Problem<2> l(
        ref, order,
        LShape::Specifier::classic,
        LShape::Strategy::uniform);
    l.run();
  } else if ( ex == "lshape_ab"){
    LShape::LShape_Problem<2> l(
        ref, order,
        LShape::Specifier::benchmark,
        LShape::Strategy::adaptive);
    l.run();
  } else if ( ex == "lshape_ub"){
    LShape::LShape_Problem<2> l(
        ref, order,
        LShape::Specifier::benchmark,
        LShape::Strategy::uniform);
    l.run();
  } else if (ex == "lshape_standard_ac"){
    LShape::LShape_Benchmark_Standard test(
      order, 
      LShape::Specifier::classic,
      LShape::Strategy::adaptive
    );
    test.run(ref);
  } else if (ex == "lshape_standard_uc"){
    LShape::LShape_Benchmark_Standard test(
      order, 
      LShape::Specifier::classic,
      LShape::Strategy::uniform
    );
    test.run(ref);
  } else if (ex == "lshape_standard_ab"){
    LShape::LShape_Benchmark_Standard test(
      order, 
      LShape::Specifier::benchmark,
      LShape::Strategy::adaptive
    );
    test.run(ref);
  } else if (ex == "lshape_standard_ub"){
    LShape::LShape_Benchmark_Standard test(
      order, 
      LShape::Specifier::benchmark,
      LShape::Strategy::uniform
    );
    test.run(ref);
  } else if (ex == "poisson_nc"){
    Poisson_Neumann::Poisson_Benchmark test(ref, order);
    test.run();
  } else if (ex == "poisson_nc_uniform") {
    Poisson_Neumann::Poisson_Benchmark test(
      ref, order, 
      Poisson_Neumann::RefinementStrategy::Uniform
    );
    test.run();
  } else if (ex == "poisson_nc_adaptive") {
    Poisson_Neumann::Poisson_Benchmark test(
      ref, order, 
      Poisson_Neumann::RefinementStrategy::Adaptive
    );
    test.run();
  } else if (ex == "poisson_nc_standard") {
    Poisson_Neumann::Poisson_Benchmark_Standard test(order);
    test.run(ref); 
  } else if (ex == "lshape3d_ac"){
    LShape::LShape_Problem<3> l(
        ref, order,
        LShape::Specifier::classic,
        LShape::Strategy::adaptive);
    l.run();
  } else if (ex == "poisson_3d") {
    Poisson_Benchmark::Poisson_Benchmark_3D p(
        ref, order
    );
    p.run();
  } else if (ex == "linear_elasticity_2d") {
    Linear_Elasticity::ElasticProblem<2> 
      problem(ref, order);
    problem.run();
  } else if (ex == "linear_elasticity_2d_uniform") {
    Linear_Elasticity::ElasticProblem<2> 
      problem(ref, order, Linear_Elasticity::Strategy::Uniform);
    problem.run();
  } else if (ex == "elasticity_inhomogeneous_3d") {
    Linear_Elasticity::InhomogeneousProblem<3>
      problem(ref, order, Linear_Elasticity::Strategy::Adaptive);
    problem.run();
  } else if (ex == "elasticity_inhomogeneous_3d_uniform") {
    Linear_Elasticity::InhomogeneousProblem<3>
      problem(ref, order, Linear_Elasticity::Strategy::Uniform);
    problem.run();
  } else if (ex == "elasticity_inhomogeneous_2d") {
    Linear_Elasticity::InhomogeneousProblem<2>
      problem(ref, order, Linear_Elasticity::Strategy::Adaptive);
    problem.run();
  } else if (ex == "linear_elasticity_3d") {
    Linear_Elasticity::ElasticProblem<3> 
      problem(ref, order);
    problem.run();
  } else if (ex == "linear_elasticity_3d_uniform") {
    Linear_Elasticity::ElasticProblem<3> 
      problem(ref, order, Linear_Elasticity::Strategy::Uniform);
    problem.run();
  } 
  else 
    Assert(false, ExcMessage("The chosen case is not implemented"));
}


