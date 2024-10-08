/*
 * residual_error_estimatos.h
 * 
 *  Created on: Mar 05, 2024
 *      Author: hoai nguyen
 * 
 */

#ifndef RESIDUAL_ERROR_ESTIMATOR_H_
#define RESIDUAL_ERROR_ESTIMATOR_H_

// STL header
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include <cmath>
#include <memory>
#include <algorithm>
#include <unordered_map>
#include <variant>
#include <vector>

#include <utility>

// Deal header
#include <deal.II/base/config.h>
#include <deal.II/grid/filtered_iterator.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>

#include <deal.II/base/tensor.h>
#include <deal.II/base/function.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/polynomials_bernstein.h>
#include <deal.II/base/quadrature_lib.h>

// Own Header
#include <utilities.h>
#include <tspline_function.h>
#include <ts_triangulation.h>
#include <isoparametric_function.h>

namespace dealt {
  using namespace dealii;
  
  namespace ResidualEstimators {
    
    template<int dim, int spacedim = dim>
    class Poisson{

    protected:
      static constexpr unsigned int dimension       = dim; 
      static constexpr unsigned int space_dimension = spacedim;
      
      using face_iterator           = TriaIterator<dealii::TriaAccessor<dim-1, dim, dim>>;
      using cell_iterator           = TriaIterator<dealii::CellAccessor<dim, dim>>;
      using active_face_iterator    = TriaActiveIterator<dealii::TriaAccessor<dim-1, dim, dim>>;
      using active_cell_iterator    = TriaActiveIterator<dealii::CellAccessor<dim, dim>>;
      
      
    public:
      static void estimate(
        const TS_TriangulationBase<dim, spacedim>*  ts_tria, 
        const std::vector< unsigned int >&          n_gauss_points,
        const Vector<double>&                       solution,
              Vector<double>&                       residuals,
        const Function<spacedim>*                   rhs_fcn,
        const std::map< 
                      types::boundary_id,
                const Function<spacedim>* 
              >&                                    neumann_bc,
        const Function<spacedim>*                   sigma = NULL,
        const Function<spacedim>*                   a = NULL
      );
    
    private: 
      static double estimate_one_cell(
        const TS_TriangulationBase<dim, spacedim>*  ts_tria,
        const active_cell_iterator                  cell,
        const std::vector< unsigned int >&          n_gauss_points,
        const Function<spacedim>*                   rhs_fcn,
        
        const Vector<double>&                       solution,
        const std::map< 
                      types::boundary_id,
                const Function<spacedim>* 
              >&                                    neumann_bc,
        const Function<spacedim>*                   sigma,
        const Function<spacedim>*                   a
      ) ;

      static double estimate_face_residuals(
          const TS_TriangulationBase<dim, spacedim>*  ts_tria,
          const active_cell_iterator                  cell,
          const std::vector<unsigned int>&            n_gauss_points,
          const Vector<double>&                       solution,
          const std::map<
                        types::boundary_id,
                  const Function<spacedim>*
                >&                                    neumann_bc
      ) ;

      static double estimate_cell_residual(
          const TS_TriangulationBase<dim, spacedim>*  ts_tria,
          const active_cell_iterator                  cell,
          const std::vector< unsigned int >&          n_gauss_points,
          const Function<spacedim>*                   rhs_fcn,
          
          const Vector<double>&                       solution,
          const Function<spacedim>*                   sigma,
          const Function<spacedim>*                   a
      ) ;

      static double estimate_cell_residual_sigma_const(
          const TS_TriangulationBase<dim, spacedim>*  ts_tria,
          const active_cell_iterator                  cell,
          const std::vector<unsigned int>&            n_gauss_points,
          const Function<spacedim>*                   rhs_fcn,
          const Vector<double>&                       solution,
          const Function<spacedim>*                   a
      ) ;
    }; // Class Poisson



    template<int dim, int spacedim = dim>
    class Linear_Elasticity{

    protected:
      static constexpr unsigned int dimension       = dim; 
      static constexpr unsigned int space_dimension = spacedim;
      
      using face_iterator           = TriaIterator<dealii::TriaAccessor<dim-1, dim, dim>>;
      using cell_iterator           = TriaIterator<dealii::CellAccessor<dim, dim>>;
      using active_face_iterator    = TriaActiveIterator<dealii::TriaAccessor<dim-1, dim, dim>>;
      using active_cell_iterator    = TriaActiveIterator<dealii::CellAccessor<dim, dim>>;
    
    public:
      static void estimate(
        const TS_TriangulationBase<dim, spacedim>*  ts_tria, 
        const std::vector< unsigned int >&          n_gauss_points,
        const Vector<double>&                       solution,
              Vector<double>&                       residuals,
        const Function<spacedim>*                   rhs_fcn,
        const std::map< 
                      types::boundary_id,
                const Function<spacedim>* 
              >&                                    neumann_bc,
        const Function<spacedim>*                   lambda = NULL,
        const Function<spacedim>*                   mu = NULL
      );

    private:
      static double estimate_one_cell(
        const TS_TriangulationBase<dim, spacedim>*  ts_tria,
        const active_cell_iterator                  cell,
        const std::vector< unsigned int >&          n_gauss_points,
        const Function<spacedim>*                   rhs_fcn,
        
        const Vector<double>&                       solution,
        const std::map< 
                      types::boundary_id,
                const Function<spacedim>* 
              >&                                    neumann_bc,
        const Function<spacedim>*                   lambda,
        const Function<spacedim>*                   mu
      ) ;

      static double estimate_face_residuals(
          const TS_TriangulationBase<dim, spacedim>*  ts_tria,
          const active_cell_iterator                  cell,
          const std::vector<unsigned int>&            n_gauss_points,
          const Vector<double>&                       solution,
          const std::map<
                        types::boundary_id,
                  const Function<spacedim>*
                >&                                    neumann_bc
      ) ;

      static double estimate_cell_residual(
          const TS_TriangulationBase<dim, spacedim>*  ts_tria,
          const active_cell_iterator                  cell,
          const std::vector< unsigned int >&          n_gauss_points,
          const Function<spacedim>*                   rhs_fcn,
          
          const Vector<double>&                       solution,
          const Function<spacedim>*                   lambda,
          const Function<spacedim>*                   mu
      ) ;

      static double estimate_cell_residual_lambda_const(
          const TS_TriangulationBase<dim, spacedim>*  ts_tria,
          const active_cell_iterator                  cell,
          const std::vector<unsigned int>&            n_gauss_points,
          const Function<spacedim>*                   rhs_fcn,
          const Vector<double>&                       solution,
          const Function<spacedim>*                   mu
      ) ;

      static double estimate_cell_residual_mu_const(
          const TS_TriangulationBase<dim, spacedim>*  ts_tria,
          const active_cell_iterator                  cell,
          const std::vector<unsigned int>&            n_gauss_points,
          const Function<spacedim>*                   rhs_fcn,
          const Vector<double>&                       solution,
          const Function<spacedim>*                   lambda
      ) ;
    
        


    }; // Class Linear_Elasticity

  } // namespace ResidualEstimators

} // namespace dealt

#endif