/*
 * Utilities.cc
 *
 *  Created on: Mar 26, 2020
 *      Author: goermer
 */

#include <utilities.h>

#ifdef DEBUG
void debug_info_start(){
  std::cout << "##############  DEBUG INFO  #####################" << std::endl;
  std::cout << "\n\n";
}
#endif

#ifdef DEBUG
void debug_info_end(){
  std::cout << std::endl;
  std::cout << "############  END DEBUG INFO  ###################" << std::endl;
  std::cout << std::endl << std::endl;
}
#endif

std::string indent(int n){
  int N = 4*n;
  std::string out;
  for (int i = 0; i < N; i++)
    out += " ";

  return out;
}

namespace dealt {
  OutputSetup::OutputSetup(
    const std::string& problem,
    const unsigned int degree,
    const std::vector<std::string>& columns,
    const std::vector<std::string>& tex_captions,
    const std::vector<bool>&        scientific, 
    const std::vector<unsigned int> precision,
    const std::vector<std::string>& super_column_names ,
    const std::vector<std::vector<std::string>> super_columns
  ) { 
    std::filesystem::path cp = std::filesystem::current_path();
    std::filesystem::path out = cp / std::filesystem::path("out/");
    this -> problem           = out / std::filesystem::path(problem);
    this -> degree            = this -> problem / std::filesystem::path("o" + std::to_string(degree) + "/");
    this -> svg               = this -> degree / std::filesystem::path("00svg/");
    this -> vtg               = this -> degree / std::filesystem::path("01vtg/");
    if (dim == 3) {
      this -> cross_sections    = this -> vtg / std::filesystem::path("cross_sections/");
      this -> cross_sections_x  = this -> cross_sections / std::filesystem::path("x/");
      this -> cross_sections_y  = this -> cross_sections / std::filesystem::path("y/");
      this -> cross_sections_z  = this -> cross_sections / std::filesystem::path("z/");
    }

    init_filesystem(cp, out, problem);
    init_table(
      columns, 
      tex_captions, 
      scientific, 
      precision,
      super_column_names,
      super_columns
    );
  }

  void OutputSetup::init_table(
    const std::vector<std::string>& columns,
    const std::vector<std::string>& tex_captions,
    const std::vector<bool>&        scientific, 
    const std::vector<unsigned int> precision,
    const std::vector<std::string>& super_column_names ,
    const std::vector<std::vector<std::string>> super_columns
  ) {
    Assert(super_column_names.size() == super_columns.size(), 
              dealii::StandardExceptions
              ::ExcMessage("Please provide an equal amount of super columns!"));
    Assert(columns.size() == tex_captions.size() || tex_captions.size() == 0, 
              dealii::StandardExceptions
              ::ExcMessage("Please make sure that either no tex caption is provided"
                         " or that each column is provided with a caption."));
    Assert(columns.size() == scientific.size() || scientific.size() == 0,
              dealii::StandardExceptions
              ::ExcMessage("Please make sure that either no column is specified"
                         " for scientific output or that each column is "
                         "specified with scientific output."));
    Assert(columns.size() == precision.size() || precision.size() == 0,
              dealii::StandardExceptions
              ::ExcMessage("Please make sure that either no column is specified"
                         " with some precision or that each column is "
                         "specified with some precision."));

    const unsigned int n_cols = columns.size();
    const unsigned int n_caps = tex_captions.size();
    const unsigned int n_scie = scientific.size(); 
    const unsigned int n_prec = precision.size();
    const unsigned int n_supe = super_columns.size();
    for (unsigned int c = 0; c < n_cols; c++)
      table.declare_column(columns[c]);

    for (unsigned int c = 0; c < n_caps; c++)
      table.set_tex_caption(columns[c], tex_captions[c]);

    for (unsigned int c = 0; c < n_scie; c++)
      table.set_scientific(columns[c], scientific[c]);

    for (unsigned int c = 0; c < n_prec; c++)
      table.set_precision(columns[c], precision[c]);
   
    for (unsigned int sc = 0; sc < n_supe; sc++)
      for (const auto& c : super_columns[sc])
        table.add_column_to_supercolumn(c, super_column_names[sc]);
  }
  
  OutputSetup::OutputSetup(
    const OutputSetup& other
  ) {
    this -> problem           = other.problem; 
    this -> degree            = other.degree;
    this -> svg               = other.svg;
    this -> vtg               = other.vtg;
    this -> table             = other.table;
    this -> cross_sections    = other.cross_sections;
    this -> cross_sections_x  = other.cross_sections_x;
    this -> cross_sections_y  = other.cross_sections_y;
    this -> cross_sections_z  = other.cross_sections_z;
  }
  
  OutputSetup& OutputSetup::operator=(
    const OutputSetup&& other
  ) { 
    this -> problem           = other.problem; 
    this -> degree            = other.degree;
    this -> svg               = other.svg;
    this -> vtg               = other.vtg;
    this -> table             = other.table;
    this -> cross_sections    = other.cross_sections;
    this -> cross_sections_x  = other.cross_sections_x;
    this -> cross_sections_y  = other.cross_sections_y;
    this -> cross_sections_z  = other.cross_sections_z;

    return *this;
  }

  void OutputSetup::init_filesystem(
    const std::filesystem::path& cp,
    const std::filesystem::path& out,
    const std::string& problem
  ) { 
    if (std::filesystem::create_directory(out, cp))
      std::cout << "Created output directory " << out.string() << std::endl;
  
  
    // Create each subpath if necessary
    std::stringstream     problem_path(problem);
    std::filesystem::path sub_path = out;
    std::string           sub_path_str; 
    while (std::getline(problem_path, sub_path_str, '/')) {
      sub_path = sub_path / std::filesystem::path(sub_path_str); 
      if (std::filesystem::create_directory(sub_path, cp))
        std::cout << "Created output directory " << sub_path.string() << std::endl;
    }
  
    if (std::filesystem::create_directory(this -> degree, cp))
      std::cout << "Created output directory " << this -> degree.string() << std::endl;
  
    if (std::filesystem::create_directory(this->svg, cp))
      std::cout << "Created output directory " << this->svg.string() << std::endl;
  
    if (std::filesystem::create_directory(this->vtg, cp))
      std::cout << "Created output directory " << this->vtg.string() << std::endl;
  
    if (dim == 3) {
      if (std::filesystem::create_directory(this->cross_sections, cp))
        std::cout << "Created output directory " << this->cross_sections.string() << std::endl;

      if (std::filesystem::create_directory(this->cross_sections_x, cp))
        std::cout << "Created output directory " << this->cross_sections_x.string() << std::endl;

      if (std::filesystem::create_directory(this->cross_sections_y, cp))
        std::cout << "Created output directory " << this->cross_sections_y.string() << std::endl;

      if (std::filesystem::create_directory(this->cross_sections_z, cp))
        std::cout << "Created output directory " << this->cross_sections_z.string() << std::endl;
    }
  }
  
  void OutputSetup::add_values_to_table(
    const unsigned int  &level, 
    const unsigned int  &n_cells,
    const unsigned int  &n_dofs,
    const unsigned int  &k_cg,
    const double        &eps_cg
  ) {
    if (prev_level != 0 && level == prev_level) {
      //remove row from table
      table.clear_current_row(); 
      cycles++;
    } else {
      prev_level = level;
      cycles = 1;
    }
  
  
    // Add values to table
    table.add_value("Level",  level);
    table.add_value("Cells",  n_cells);
    table.add_value("DoFs",   n_dofs);
    table.add_value("Cycles", cycles);
    table.add_value("k_{CG}",   k_cg);
    table.add_value("TOL_{CG}", eps_cg);
  }
  
  void OutputSetup::add_values_to_table(
    const double        &l2_error,
    const double        &h1_error
  ) {
    table.add_value("L2", l2_error);
    table.add_value("H1", h1_error);
  }

  void OutputSetup::add_values_to_table(
    const std::vector<std::string>& columns,
    const std::vector<double>&      values,
    const bool&                     clear
  ) { 
    Assert(columns.size() == values.size(), 
            dealii::StandardExceptions
            ::ExcMessage("Each column in columns needs exactly "
                       "one value from values."));

    if (clear)
      table.clear_current_row();

    const unsigned int n_cols = columns.size();
    for (unsigned int c = 0; c < n_cols; c++)
      table.add_value(columns[c], values[c]);
  }

  void OutputSetup::write_table_text(
    std::ostream& out
  ) const {
    table.write_text(out, dealii::TableHandler::TextOutputFormat::org_mode_table);
  }
  
  void OutputSetup::write_table_tex(
    std::ostream& out
  ) const {
    table.write_tex(out);
  }
  
  void OutputSetup::write_table_text(
    const std::string& name
  ) const {
    std::string out_name;
    if (name == "")
      out_name = degree.string() + "table";
    else 
      out_name = name;
  
    std::ofstream out(out_name + ".txt");
    table.write_text(out, dealii::TableHandler::TextOutputFormat::org_mode_table);
  }
  
  void OutputSetup::write_table_tex(
    const std::string& name 
  ) const {
    std::string out_name;
    if (name == "")
      out_name = degree.string() + "table";
    else 
      out_name = name;
  
    std::ofstream out(out_name + ".tex");
    table.write_tex(out);
  }

  void OutputSetup::set_dim(
    const unsigned int other_dim
  ) {
    dim = other_dim;
  }
  unsigned int OutputSetup::dim = 0;
}
