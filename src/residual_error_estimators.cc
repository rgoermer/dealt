/*
 * residual_error_estimators.cc
 * 
 *  Created on: Mar 05, 2024
 *      Author: hoai nguyen
 * 
 */

#include <residual_error_estimators.h> 



template class dealt::ResidualEstimators::Poisson<2, 2>;
template class dealt::ResidualEstimators::Poisson<2, 3>;
template class dealt::ResidualEstimators::Poisson<3, 3>;

namespace dealt {
  using namespace dealii;
  namespace ResidualEstimators {

    template <int dim, int spacedim>
    void Poisson<dim, spacedim>::estimate(
        const TS_TriangulationBase<dim, spacedim>*    ts_tria,
        const std::vector<unsigned int>&              n_gauss_points,
        const Vector<double>&                         solution,
              Vector<double>&                         residuals,
        const Function<spacedim>*                     rhs_fcn,
        const std::map<
                      types::boundary_id,
                const Function<spacedim>*
              >&                                      neumann_bc,
        const Function<spacedim>*                     sigma,
        const Function<spacedim>*                     a
      ) {
      Assert(dim <= spacedim, ExcNotImplemented());
      Assert(dim >= 2 && dim <= 3, ExcNotImplemented());
      Assert(spacedim >= 2 && spacedim <= 3, ExcNotImplemented());

      // Loop over all cells
      unsigned int counter = 0;
      for (const auto &cell : ts_tria->active_cell_iterators())
      {
        residuals(counter++) = estimate_one_cell(ts_tria, cell,
                                                 n_gauss_points,
                                                 rhs_fcn, solution,
                                                 neumann_bc,
                                                 sigma, a);
      }
    } // estimate

      template<int dim, int spacedim>
      double Poisson<dim, spacedim>::estimate_one_cell(
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
      ) {

        double face_residuals = 0.;
        double cell_residual = 0.;
        if (!neumann_bc.empty())
          face_residuals = estimate_face_residuals(ts_tria, cell,
                                              n_gauss_points, 
                                              solution, neumann_bc);
        // check, if sigma is NULL
        bool sigma_const = false;
        try
        {
          if (sigma != NULL)
            sigma->gradient(Point<spacedim>());
        }
        catch (...)
        {
          sigma_const = true;
        }
        if (sigma_const == true) {
          cell_residual = estimate_cell_residual_sigma_const(ts_tria, cell,
                                                      n_gauss_points,
                                                      rhs_fcn, solution,
                                                      a);
        }
        else {
          cell_residual = estimate_cell_residual(ts_tria, cell,
                                          n_gauss_points,
                                          rhs_fcn, solution,
                                          sigma, a);
        }
        return std::sqrt(face_residuals + cell_residual); /*needs adjustment, but what exactly?*/
      } // estimate_one_cell

      template<int dim, int spacedim>
      double Poisson<dim, spacedim>::estimate_face_residuals(
          const TS_TriangulationBase<dim, spacedim>*  ts_tria,
          const active_cell_iterator                  cell,
          const std::vector< unsigned int >&          n_gauss_points,
          const Vector<double>&                       solution,
          const std::map< 
                        types::boundary_id,
                  const Function<spacedim>* 
                >&                                    neumann_bc
      ) {
        TSFaceValues<dimension, space_dimension> face_values(
            ts_tria, n_gauss_points,
            update_values |
            update_gradients |
            update_quadrature_points |
            update_JxW_values |
            update_normal_vectors);

        const std::vector<unsigned int>& p              = ts_tria -> get_degree();
        const std::map<unsigned int, unsigned int>& mof = ts_tria -> get_mof();
        std::vector<unsigned int> local_dof_indices     = ts_tria -> get_IEN_array(cell);
        const unsigned int nvf = GeometryInfo<dimension>::vertices_per_face;
        double face_residuals = 0.;

        for (unsigned int f = 0;
              f < GeometryInfo<dimension>::faces_per_cell;
              f++) {
          // Get orientation of current face
          const auto& face = cell -> face(f);
          const Point<dim>& c = (-1.) *
                                face->vertex(0) +
                                face->vertex(nvf - 1);
          unsigned int orientation = 0;
          for ( ; orientation < dimension && c(orientation) != 0; orientation++);
          
          // Compute the Jumps of faces at C0 continuity
          if (face->has_children())
          {
            // This is merely a special case, as refined faces are considered
            // in-active and hence their index is not represented in mof
            if (mof.at(face->child(0)->index()) == p[orientation])
            {
              // C0 Edge detected:
              face_values.reinit(cell, f);
              //const unsigned int ppf = 2 * face_values.n_quadrature_points_per_face();
              for (const unsigned int q_index : face_values.quadrature_point_indices())
              {
                double g = 0;
                for (const unsigned int i : face_values.dof_indices())
                  g += solution(local_dof_indices[i]) *
                       face_values.shape_grad(i, q_index) *
                       face_values.normal_vector(q_index);
                //const unsigned int ch = q_index > ppf ? 1 : 0;
                face_residuals += 0.25 * g * g * face_values.JxW(q_index);
              } // for ( q_index )
            } // otherwise it is atleast a C1 edge, and cannot be at the boundary
          }
          else if (mof.at(face->index()) == p[orientation])
          {
            // C0 edge detected
            face_values.reinit(cell, f);
            // Since this face is not refined, we can simply compute
            // the jump terms along it and store the values at the
            // corresponding place in face_residuals
            for (const unsigned int q_index : face_values.quadrature_point_indices())
            {
              double g = 0;
              for (const unsigned int i : face_values.dof_indices())
                g += solution(local_dof_indices[i]) *
                     face_values.shape_grad(i, q_index) *
                     face_values.normal_vector(q_index);
              face_residuals += 0.25 * g * g * face_values.JxW(q_index);
            } // for ( q_index )
          }
          else if (face->at_boundary())
          {
            const auto &bc = neumann_bc.find(face->boundary_id());
            if (bc == neumann_bc.end())
              continue;

            face_values.reinit(cell, f);
            // If the face is at the boundary, compute the jump from
            // the face's boundary id
            for (const unsigned int q_index : face_values.quadrature_point_indices())
            {
              double g = neumann_bc.at(face->boundary_id())
                             ->value(face_values.quadrature_point(q_index));
              for (const unsigned int i : face_values.dof_indices())
                g -= solution(local_dof_indices[i]) *
                     face_values.shape_grad(i, q_index) *
                     face_values.normal_vector(q_index);
              face_residuals += g * g * face_values.JxW(q_index);
            } // for ( q_index )
          } // if ( ... )
        } // for (faces)

        return face_residuals;
      } // estimate_face_residual_sigma_const

      template<int dim, int spacedim>
      double Poisson<dim, spacedim>::estimate_cell_residual(
          const TS_TriangulationBase<dim, spacedim>*  ts_tria,
          const active_cell_iterator                  cell,
          const std::vector< unsigned int >&          n_gauss_points,
          const Function<spacedim>*                   rhs_fcn,
          const Vector<double>&                       solution,
          const Function<spacedim>*                   sigma,
          const Function<spacedim>*                   a
      ) {
        TSValues<dimension, space_dimension> ts_values(
            ts_tria, n_gauss_points,
            update_values |
            update_gradients |
            update_quadrature_points |
            update_JxW_values |
            update_hessians);
        const unsigned int nvc = GeometryInfo<dimension>::vertices_per_cell;

        // Compute the residual on the cell
        double cell_residual = 0.;
        const std::vector<unsigned int> &local_dof_indices = ts_tria -> get_IEN_array(cell);
        ts_values.reinit(cell);
        for (const unsigned int q_index : ts_values.quadrature_point_indices())
        {
          double slu  = 0;
          double gsgu = 0;
          double au   = 0;
          const Point<space_dimension> &Q = ts_values.quadrature_point(q_index);
          const Tensor<1, space_dimension> &sigma_grad = 
                                sigma==NULL 
                                ? Tensor<1, space_dimension>() 
                                : sigma->gradient(Q);
          const double sigma_val = (sigma==NULL ? 1 : sigma->value(Q));
          const double a_val = (a==NULL ? 0 : a->value(Q));
          for (const unsigned int i : ts_values.dof_indices())
          {
            double l = 0;
            for (unsigned int d = 0; d < space_dimension; d++)
              l  += ts_values.shape_hessian(i, q_index)[d][d];
            slu  += solution(local_dof_indices[i]) * sigma_val * l;
            au   += solution(local_dof_indices[i]) * 
                    a_val * 
                    ts_values.shape_value(i, q_index);
            gsgu += solution(local_dof_indices[i]) *
                    sigma_grad *
                    ts_values.shape_grad(i, q_index);
          } // for ( i )
          double g = (slu + gsgu - au + rhs_fcn->value(Q));
          cell_residual += g * g * ts_values.JxW(q_index);
        } // for ( q_index )
        const double cell_width = (cell->vertex(0)).distance(cell->vertex(nvc - 1));

        return cell_residual * cell_width * cell_width;
      } // estimate_cell_residual


      template<int dim, int spacedim>
      double Poisson<dim, spacedim>::estimate_cell_residual_sigma_const(
          const TS_TriangulationBase<dim, spacedim>*  ts_tria,
          const active_cell_iterator                  cell,
          const std::vector<unsigned int>&            n_gauss_points,
          const Function<spacedim>*                   rhs_fcn,
          const Vector<double>&                       solution,
          const Function<spacedim>*                   a
      ) {
        TSValues<dimension, space_dimension> ts_values(
            ts_tria, n_gauss_points,
            update_values |
            update_gradients |
            update_quadrature_points |
            update_JxW_values |
            update_hessians);

        // Compute the residual on the cell
        double cell_residual = 0.;
        const unsigned int nvc = GeometryInfo<dimension>::vertices_per_cell;
        const std::vector<unsigned int> &local_dof_indices = ts_tria -> get_IEN_array(cell);
        ts_values.reinit(cell);
        for (const unsigned int q_index : ts_values.quadrature_point_indices())
        {
          double slu  = 0;
          double gsgu = 0;
          double au   = 0;
          const Point<space_dimension> &Q = ts_values.quadrature_point(q_index);
          const Tensor<1, space_dimension> &sigma_grad = Tensor<1, space_dimension>(); 
          const double a_val = (a==NULL ? 0 : a->value(Q));
          for (const unsigned int i : ts_values.dof_indices())
          {
            double l = 0;
            for (unsigned int d = 0; d < space_dimension; d++)
              l  += ts_values.shape_hessian(i, q_index)[d][d];
            slu  += solution(local_dof_indices[i]) /** sigma_val*/ * l;
            au   += solution(local_dof_indices[i]) * 
                    a_val * 
                    ts_values.shape_value(i, q_index);
            gsgu += solution(local_dof_indices[i]) *
                    sigma_grad *        
                    ts_values.shape_grad(i, q_index);
          } // for ( i )
          double g = (slu + gsgu - au + rhs_fcn->value(Q));
          cell_residual += g * g * ts_values.JxW(q_index);
        } // for ( q_index )
        const double cell_width = (cell->vertex(0)).distance(cell->vertex(nvc - 1));

        return cell_residual * cell_width * cell_width;
      } // estimate_cell_residual_sigmal_const


      template<int dim, int spacedim>
      void Linear_Elasticity<dim, spacedim>::estimate(
          const TS_TriangulationBase<dim, spacedim>*    ts_tria,
          const std::vector<unsigned int>&              n_gauss_points,
          const Vector<double>&                         solution,
                Vector<double>&                         residuals,
          const Function<spacedim>*                     rhs_fcn,
          const std::map<
                        types::boundary_id,
                  const Function<spacedim>*
                >&                                      neumann_bc,
          const Function<spacedim>*                     lambda,
          const Function<spacedim>*                     mu
        ) {


        // Loop over all cells
        unsigned int counter = 0;
        for (const auto &cell : ts_tria->active_cell_iterators())
        {
          residuals(counter++) = estimate_one_cell(ts_tria, cell,
                                                   n_gauss_points,
                                                   rhs_fcn, solution,
                                                   neumann_bc,
                                                   sigma, a);
        }
      } // estimate
  } // namespace ResidualEstimators

} // namespace dealt

