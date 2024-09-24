/*
 * Utilities.h
 *
 *  Created on: Mar 24, 2020
 *      Author: goermer
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_


#include <deal.II/base/convergence_table.h>

#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <filesystem>


#ifdef DEBUG
void extern debug_info_start();
#endif

#ifdef DEBUG
void extern debug_info_end();
#endif

std::string extern indent(int n = 1);


namespace dealt {
  struct OutputSetup {
  private: 
    static unsigned int dim;

  public: 
    static void set_dim(const unsigned int dim);

    OutputSetup(){}
  
    OutputSetup(
      const std::string& problem,
      const unsigned int degree,
      const std::vector<std::string>& columns =
        {"Cycles", "Level", "Cells", "DoFs", "L2", "H1", "k_{CG}", "TOL_{CG}"},
      const std::vector<std::string>& tex_captions = 
        {"Cycles", "Level", "\\# cells", "\\# dofs", "$L_2$-error", "$H^1$-error", "CG iters", "CG last tol"},
      const std::vector<bool>& scientific =
        {false, false, true, true, true, true, true, true},
      const std::vector<unsigned int> precision = 
        {0, 0, 0, 0, 8, 8, 0, 8},
      const std::vector<std::string>& super_column_names = 
        {"CG-solver", "Errors", "Grid Info"},
      const std::vector<std::vector<std::string>> super_columns =
        {{"k_{CG}", "TOL_{CG}"}, {"L2", "H1"}, {"Cycles", "Level", "Cells", "DoFs"}}
    );
  
    OutputSetup(
      const OutputSetup& other
    );
  
    OutputSetup& operator=(
      const OutputSetup&& other
    );

    void init_filesystem(
      const std::filesystem::path& cp,
      const std::filesystem::path& out,
      const std::string& problem
    );
    void init_table(
      const std::vector<std::string>& columns,
      const std::vector<std::string>& tex_captions,
      const std::vector<bool>&        scientific, 
      const std::vector<unsigned int> precision,
      const std::vector<std::string>& super_column_names ,
      const std::vector<std::vector<std::string>> super_columns
    );
  
    void add_values_to_table(
      const unsigned int  &level, 
      const unsigned int  &n_cells,
      const unsigned int  &n_dofs,
      const unsigned int  &k_cg,
      const double        &eps_cg
    );
  
    void add_values_to_table(
      const double        &l2_error,
      const double        &h1_error
    );

    void add_values_to_table(
      const std::vector<std::string> &columns,
      const std::vector<double>      &values,
      const bool                     &clear = false
    );
  
    void write_table_text(
      std::ostream& out
    ) const;
  
    void write_table_tex(
      std::ostream& out
    ) const;
  
    void write_table_text(
      const std::string& name = ""  
    ) const;
  
    void write_table_tex(
      const std::string& name = ""  
    ) const;
  
  public:                 
    // These are the main paths to be used to output the results
    std::filesystem::path problem; 
    std::filesystem::path degree; 
    std::filesystem::path svg;    
    std::filesystem::path vtg;    
    std::filesystem::path cross_sections;
    std::filesystem::path cross_sections_x;
    std::filesystem::path cross_sections_y;
    std::filesystem::path cross_sections_z;
  
    // This is the convergence table. It removes the most recent
    // row if the input parameter level coincides with prev_level
    dealii::ConvergenceTable    table;
    
    unsigned int        prev_level = 0;
    unsigned int        cycles     = 1;
  
  } ;
}

#endif /* INC_UTILITIES_H_ */
