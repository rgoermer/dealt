/*
 * TS_Triangulation.h
 *
 *  Created on: Mar 26, 2020
 *      Author: goermer
 */

#ifndef TS_TRIANGULATION_H_
#define TS_TRIANGULATION_H_

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
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>

#include <deal.II/base/tensor.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/polynomials_bernstein.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/numerics/data_out.h>


// Own Header
#include <utilities.h>
#include <tspline_function.h>
#include <isoparametric_function.h>

namespace dealt {
  using namespace dealii;

  template<int dim, int spacedim = dim>
  struct IPF_Data {
    std::vector< Point<spacedim+1> >      cps;
    std::vector< std::vector< double > >  kv;
    std::vector< unsigned int >           deg;

    IPF_Data( const std::vector< Point<spacedim+1> >      &other_cps,
              const std::vector< std::vector< double > >  &other_kv,
              const std::vector< unsigned int >           &other_deg
    ) : cps(other_cps), kv(other_kv), deg(other_deg) {}
    IPF_Data(
        const std::vector< Point<spacedim> >        &other_cps,
        const std::vector< double >                 &other_w,
        const std::vector< std::vector< double > >  &other_kv,
        const std::vector< unsigned int >           &other_deg
    ) : kv(other_kv), deg(other_deg) {
#ifdef DEBUG
      const unsigned int n_cps = other_cps.size();
      const unsigned int n_w   = other_w.size();
      AssertDimension(n_cps, n_w);
#endif
      cps.resize(other_cps.size());
      for (unsigned int i = 0; i < other_cps.size(); i++){
        for (unsigned int d = 0; d < spacedim; d++)
          cps[i](d) = other_w[i] * other_cps[i](d);
        cps[i](spacedim) = other_w[i];
      }
    }

    IPF_Data() = default;

    unsigned int max_degree() const {
      unsigned int max = deg[0];
      for(unsigned int d = 0; d < dim; d++)
        if (max < deg[d])
          max = deg[d];

      return max;
    }
  };

  template<int dim, int spacedim = dim>
  class TS_TriangulationBase : public Triangulation<dim, dim> {
  // =============================================================================
  // =============================================================================
  //                 typename declarations  
  // =============================================================================
  // =============================================================================
  protected: 
    static constexpr unsigned int dimension       = dim; 
    static constexpr unsigned int space_dimension = spacedim;

    //  Define iterator shortcuts to be used throughout this class -- these names are
    //  according to the deal.II library and directly derived from there
    using face_iterator           = TriaIterator<dealii::TriaAccessor<dim-1, dim, dim>>;
    using cell_iterator           = TriaIterator<dealii::CellAccessor<dim, dim>>;
    using active_face_iterator    = TriaActiveIterator<dealii::TriaAccessor<dim-1, dim, dim>>;
    using active_cell_iterator    = TriaActiveIterator<dealii::CellAccessor<dim, dim>>;

    using vertex_iterator         = TriaIterator<dealii::TriaAccessor<0, dimension, dimension>>;
    using active_vertex_iterator  = TriaActiveIterator<dealii::TriaAccessor<0, dimension, dimension>>;
    using IteratorSelector        = dealii::internal::TriangulationImplementation::Iterators<dimension, dimension>;
    using line_iterator           = typename IteratorSelector::line_iterator;
    using active_line_iterator    = typename IteratorSelector::active_line_iterator;
    using quad_iterator           = typename IteratorSelector::quad_iterator;
    using active_quad_iterator    = typename IteratorSelector::active_quad_iterator;
    using hex_iterator            = typename IteratorSelector::hex_iterator;
    using active_hex_iterator     = typename IteratorSelector::active_hex_iterator;

    // A sorthand notation for the anchor bounds as a tupel of points
    using anchor_bounds    = std::pair<Point<dimension>, Point<dimension>>;

    // A shorthand notation for TSplineFunction and a shared pointer thereof
    using TSpline          = TSplineFunction<dimension, space_dimension>;
    using ts_ptr           = std::shared_ptr<TSpline>;
    using ControlPoint     = Point<space_dimension+ 1>;


  // =============================================================================
  // =============================================================================
  //                 member declarations 
  // =============================================================================
  // =============================================================================
  protected:
    // Save the entirety of TSplines in a vector, i.e. the basis
    // ToDo: This vector might be unnecessary, since the parents are
    // stored within a parent pointer from new T-Splines.
    std::vector< ts_ptr > old_splines;

    // Save the active splines in a seperate array
    std::vector< ts_ptr > active_splines;


    // Since the vector of splines contains all splines over all levels,
    // we store the number of active splines in a seperate variable
    unsigned int n_splines;

    // For an easy implementation we save the bezier mesh B(T) together with the T-mesh,
    // where the extra bezier elements generated by extending hanging interfaces, will be
    // stored in an extra array. Note however, that we save the parent elements, i.e.
    // the elements in the parametric mesh T, that will be cut by extending hanging
    // interfaces. Thus, to compute the integrals arising in the bilinear form a(u, v)
    // or in the right hand side f(v), the user is supposed to perform the following
    // operations:
    // 1. Cut all elements in bezier_elements accordingly; no splines are generated
    // 2. Perform computations using standard finite element procedure and expressing
    //    TSplines locally as Bernstein-Polynomials
    // 3. Coarsen the refined elements immediately after refinement, to get back to the
    //    parametric T-mesh.
    //
    // This procedure enables us to use Bezier-Decomposition for T-Splines, because
    // on every cell we have the same amount of T-Splines.
    //
    // Consider the following parametric grid. Note, that the double lines on the
    // right indicate lines with multiplicity 2.
    //
    //    + - - - + - - - + - - - - - - - + +
    //    |       |       |               | |
    //    |       |       |               | |
    //    |       |       |               | |
    //    + - - - + - - - + - - - - - - - + +
    //    |       |       |               | |
    //    |       |       |       E       | |
    //    |       |       |               | |
    //    + - - - + - - - + - - - + - - - + +
    //    |       |       |       |       | |
    //    |       |       |       |       | |
    //    |       |       |       |       | |
    //    + - - - + - - - + - - - + - - - + +
    //
    // With px = py = 2, the T-mesh element E has 10 DoFs, but from the Bezier Extraction
    // Operator we infer 9 = (2 + 1)*(2 + 1) DoFs. So instead we add the T-Junction axtension
    // to obtain the Bezier mesh by step 1
    //
    //    + - - - + - - - + - - - - - - - + +
    //    |       |       |               | |
    //    |       |       |               | |
    //    |       |       |               | |
    //    + - - - + - - - + - - - + - - - + +
    //    |       |       |       :       | |
    //    |       |       |       :       | |
    //    |       |       |       :       | |
    //    + - - - + - - - + - - - + - - - + +
    //    |       |       |       |       | |
    //    |       |       |       |       | |
    //    |       |       |       |       | |
    //    + - - - + - - - + - - - + - - - + +
    //
    // where no new T-Splines are generated. Each Bezier Element now has 9 DoFs and
    // we can compute the integrals for each Bezier element. Note, that only the cut
    // element E changed, that is stored seperately, and each old cell including the
    // children of E are now the bezier elements. After integration [step 2] we then
    // reverse this process by coarsening the bezier elements.
    //
    // We provide functions for refining the bezier elements and coarsening.
    //
    // Further, we will use this array to switch between the parametric T-mesh and
    // the Bezier T-mesh
    std::vector< cell_iterator > bezier_elements;

    // Vector of degrees
    std::vector<unsigned int> p;

    // Save the distinct base knots and its multiplicities. This will be used
    // to raise the degree of the basis functions before refinement. While
    // the base knots will only be initially used for degree elevation, 
    // the vector of multiplicities is used to declare the multiplicity
    // of faces.
    std::vector<std::vector<double>>        base_knots;
    std::vector<std::vector<unsigned int>>  multiplicities;

    // Assign each active cell of the bezier mesh to a set 
    // of globa DoFs. This may then be used for matrix assembly
    std::map< 
      const active_cell_iterator, 
      std::vector< types::global_dof_index > 
    > IEN_array;

    // Assign each boundary id a set of global DoFs. 
    std::map< 
      types::boundary_id,
      std::vector< types::global_dof_index>
    > boundary_dofs;

    // the extraction operators in nD are computed from 1D extraction operators. 
    // To save some storage, we assign each active cell of the bezier mesh a 
    // container where each entry corresponds to the 1D extraction operator rows
    // of this spline. Combined with the IEN_array, 
    //      extraction_operators[cell][i]
    // returns the dim 1D extraction operators on cell for the global DoF 
    // IEN_array[cell][i]. 
    // Using tensor products of vectors, the complete container for a cell
    // will later be transformed to a FullMatrix<double> and passed onto the 
    // TSValues or TSFaceValues class. See compute_extraction_operators(). 
    std::map< 
      active_cell_iterator,
      std::vector< std::array< Vector<double>, dimension> >
    > extraction_operators;

    // This is the same as before, but since faces may be split in two parts,
    // we need new extraction operator rows for faces. 
    std::map< 
      active_cell_iterator,
      std::map<
        unsigned int,
        std::vector< std::array< Vector<double>, dimension> >
      >
    > face_operators;

    // To ease on checks for bezier mesh, we store a variable that is used internally
    // only as a safeguard to ensure that some computations can only be done if
    // this mesh is a bezier mesh.
    bool                                    is_bezier_mesh = false;

    // Save for each face the multiplicity
    std::map< 
      unsigned int /* face index */,
      unsigned int /* mult */
    > mof;

    // Upper and lower bounds of knot vectors
    // ToDo: This may be derived from base_knots. So it may be removed.
    Point<dimension> kv_lower;
    Point<dimension> kv_upper;

    // IPF to be constructed offline for the geometry mapping
    IsoparametricFunction<dimension, space_dimension> IPF;

    // for the sparsity pattern, we sort the tsplines specified by a pattern
    // It is initialized to sort by barycentriccoordiantes of tsplines, however,
    // if renumber_splinese is called with a different argument, the pattern is
    // overwritten in order to sort by the new pattern after each refinement
    // step.
    //
    // Currently the following patterns are allowed:
    // 1. barycentric     - sort tsplines by barycentric coordinates
    // 2. volumetric      - sort tsplines by sice of its support, starting with the smallest.
    //                      If the supports of two splines are of the same size, we sort by
    //                      barycentric coordinates
    // 3. inv_volumetric  - sort tsplines by sice of its support, starting with the largest.
    //                      If the supports of two splines are of the same size, we sort by
    //                      barycentric coordinates
    // ToDo: Deprecated. Remove functionality.
    std::string pattern = "barycentric";

    // In case we start with anisotropic cells, we want to cut the longes edge first:
    unsigned int level_offset;

  // =============================================================================
  // =============================================================================
  //                constructor  
  // =============================================================================
  // =============================================================================
  public:
    // Standard constructor to initialize dealii::Triangulation<dim>
    // other values will be defined by derived classes
    TS_TriangulationBase(
    ) : Triangulation<dimension>(
          Triangulation<dimension>::MeshSmoothing::none, true
        ) {}

    ~TS_TriangulationBase() = default;

    TS_TriangulationBase(
      const TS_TriangulationBase<dim, spacedim>& other
    ) = delete;

    TS_TriangulationBase<dim, spacedim>& operator=(
      TS_TriangulationBase<dim, spacedim>&& other
    ) noexcept = default;

  // =============================================================================
  // =============================================================================
  //              Getter and Setter functions
  // =============================================================================
  // =============================================================================
  public:
    // get the underlying Isoparametric mapping
    const IsoparametricFunction<dim, spacedim>& 
        get_IPF() const ;

    // get the number of active TSplines
    const unsigned int& 
        n_active_splines()const;

    // get the vector of polynomial degrees
    const std::vector< unsigned int >& 
        get_degree() const ;

    // get the degree at a specific direction
    const unsigned int& 
        get_degree(const unsigned int d) const; 

    // get the map of multiplicities
    const std::map< unsigned int, unsigned int >& 
        get_mof() const;

    // get IEN array
    const std::map<
            const active_cell_iterator, 
            std::vector<types::global_dof_index> 
      >& get_IEN_array() const ;

    // get IEN array for a specific cell
    const std::vector< types::global_dof_index >& 
        get_IEN_array(
      const active_cell_iterator& cell
    ) const ;

    // Get IEN array for vector valued systems
    std::map<
        const active_cell_iterator, 
        std::vector<types::global_dof_index> 
      > get_IEN_array(
      const unsigned int n_components
    ) const;

    // get IEN array for vector valued problems on a specific cell
    std::vector< types::global_dof_index > 
        get_IEN_array(
      const active_cell_iterator& cell,
      const unsigned int n_components
    ) const;

    // get active splines on specific cell
    const std::vector< ts_ptr > 
        get_splines(
      const active_cell_iterator& cell
    ) const ;

    // get all active splines
    const std::vector< ts_ptr >& get_splines() const ;

    // return points that define a box around this triangulation
    std::pair< const Point<dim>&, const Point<dim>& > 
        get_bounding_box(
    ) const ;

    // get Bezier coefficients on a specific cell
    virtual const FullMatrix<double> get_bezier_coefficients(
        const active_cell_iterator& cell
    ) const = 0;

    // get Bezier coefficients on a specific cell at a speecific face
    // Note, that the face index corresponds to a face at the boundary of
    // this cell
    virtual const std::vector< FullMatrix<double> >
      get_bezier_coefficients(
        const active_cell_iterator& cell,
        const unsigned int face
    ) const = 0;

    // return boundary dof indices with boundary ids
    const std::map<
        types::boundary_id,
        std::vector< types::global_dof_index >
      >& get_boundary_dofs(
    ) const ;

    // return boundary dof indices with boundary ids
    std::map<
        types::boundary_id,
        std::vector< types::global_dof_index >
      > get_boundary_dofs(
      const unsigned int n_components
    ) const ;

    // return a container to the bezier elements
    const std::vector<cell_iterator>& get_bezier_elements(
    ) const ;

    // Return the control points of the splines of cell specified by its id
    const std::vector< ControlPoint > get_control_points(
        const active_cell_iterator& cell
    ) const ;

    // get a (heuristic) number of non-zero matrix entries for sparse matrices
    unsigned int get_max_nnz_supports() const ;

    std::vector< unsigned int > get_max_entries_per_row(
      const unsigned int n_components = 1
    ) const;

  // =============================================================================
  // =============================================================================
  //                Routines for matrix assembly
  // =============================================================================
  // =============================================================================
  public:
    void prepare_assembly();

    void set_refine_flags(
        const std::vector< cell_iterator >& mark
    ) ;

    template<int soldim = 1>
    void project_boundary_values(
      const std::map< 
                types::boundary_id, 
                const Function<spacedim, double>* 
            >& boundary_functions,
      const std::vector< unsigned int>& n_gauss_points,
            std::map< 
                types::global_dof_index,
                double 
            >& boundary_values
    ) const;


    /// Computes the residual error for poisson-like problems of the form
    ///  \f[
    ///    \left\lbrace
    ///    \f{aligned}{
    ///      -\grad\cdot(\sigma \grad u) &= f     &&\text{in} \Omega 
    ///      \frac{\partial u}{\partial n} &= g_i^N &&\text{on} \Gamma_i^N \subset \partial \Omega
    ///      u &= g^N_j && \text{on} \Gamma^D_j \subset\partial\Omega
    ///    \f}
    ///    \right.
    ///  \f]
    /// where \f$\sigma\colon \mathbb{R}^d\to \mathbb{R}\f$ is a continuous differentiable 
    /// function, \f$\colon\mathbb{R}^d \to \mathbb{R}\f$ is a square integrable function, 
    /// \f$ g_i^N\colon\\mathbb{R}^d \to \mathbb{R} \f$, \f$i=1,\dots,m_N\f$, and \f$ 
    /// g_j^D\colon\mathbb{R}^d \to \mathbb{R}\f$, \f$ j = 1,\dots, m_D\f$ are 
    /// Neumann, resp. Dirichleth conditions. Note that 
    ///  \f[
    ///    \bigcup_{i=1}^m_N \Gamma_i^N \cap \bigcup_{j=1}^m_D \Gamma_j^D = \emptyset
    ///  \f]
    /// and further 
    ///  \f[
    ///    \Gamma_i^a \cap \Gamma_j^a = \emptyset,\, i=1,\dots,m_a, \, a\in{N,D}.
    ///  \f]
    void poisson_residual_error_estimate(
        const std::vector< unsigned int >&          n_gauss_points,
        const Function<space_dimension>*            rhs_fcn,
        const Function<space_dimension>*            sigma,
        const std::map< 
                      types::boundary_id,
                const Function<space_dimension>* 
              >&                                    neumann_bc,
        const Vector<double>&                       solution,
        std::map< cell_iterator,
                  double >&                         residuals
    ) const;


    /// Computes the residual error for poisson-like problems of the form
    ///  \f[
    ///    \left\lbrace
    ///    \f{aligned}{
    ///      -\grad\cdot(\sigma \grad u) + au &= f     &&\text{in} \Omega 
    ///      \frac{\partial u}{\partial n} &= g_i^N &&\text{on} \Gamma_i^N \subset \partial \Omega
    ///      u &= g^N_j && \text{on} \Gamma^D_j \subset\partial\Omega
    ///    \f}
    ///    \right.
    ///  \f]
    /// where \f$\sigma\colon \mathbb{R}^d\to \mathbb{R}\f$ is a continuous differentiable 
    /// function, \f$\colon\mathbb{R}^d \to \mathbb{R}\f$ is a square integrable function, 
    /// \f$ g_i^N\colon\\mathbb{R}^d \to \mathbb{R} \f$, \f$i=1,\dots,m_N\f$, and \f$ 
    /// g_j^D\colon\mathbb{R}^d \to \mathbb{R}\f$, \f$ j = 1,\dots, m_D\f$ are 
    /// Neumann, resp. Dirichleth conditions. Note that 
    ///  \f[
    ///    \bigcup_{i=1}^m_N \Gamma_i^N \cap \bigcup_{j=1}^m_D \Gamma_j^D = \emptyset
    ///  \f]
    /// and further 
    ///  \f[
    ///    \Gamma_i^a \cap \Gamma_j^a = \emptyset,\, i=1,\dots,m_a, \, a\in{N,D}.
    ///  \f]
    void poisson_residual_error_estimate(
        const std::vector< unsigned int >&          n_gauss_points,
        const Function<space_dimension>*            rhs_fcn,
        const Function<space_dimension>*            sigma,
        const Function<space_dimension>*            a,
        const std::map< 
                      types::boundary_id,
                const Function<space_dimension>* 
              >&                                    neumann_bc,
        const Vector<double>&                       solution,
              Vector<double>&                       residuals
    ) const;

    void poisson_residual_error_estimate(
        const std::vector< unsigned int >&          n_gauss_points,
        const Function<space_dimension>*            rhs_fcn,
        const Functions::ConstantFunction<space_dimension>*    sigma,
        const Function<space_dimension>*            a,
        const std::map< 
                      types::boundary_id,
                const Function<space_dimension>* 
              >&                                    neumann_bc,
        const Vector<double>&                       solution,
              Vector<double>&                       residuals
    ) const;

    void poisson_residual_error_estimate(
        const std::vector< unsigned int >&          n_gauss_points,
        const Function<space_dimension>*            rhs_fcn,
        const Function<space_dimension>*            sigma,
        const Functions::ConstantFunction<space_dimension>*    a,
        const std::map< 
                      types::boundary_id,
                const Function<space_dimension>* 
              >&                                    neumann_bc,
        const Vector<double>&                       solution,
              Vector<double>&                       residuals
    ) const;

    void poisson_residual_error_estimate(
        const std::vector< unsigned int >&          n_gauss_points,
        const Function<space_dimension>*            rhs_fcn,
        const Functions::ConstantFunction<space_dimension>*    sigma,
        const Functions::ConstantFunction<space_dimension>*    a,
        const std::map< 
                      types::boundary_id,
                const Function<space_dimension>* 
              >&                                    neumann_bc,
        const Vector<double>&                       solution,
              Vector<double>&                       residuals
    ) const;

    /// Same as Before, but remove const-ness from neumann_bc map
    void poisson_residual_error_estimate(
        const std::vector< unsigned int >&          n_gauss_points,
        const Function<space_dimension>*            rhs_fcn,
        const Function<space_dimension>*            sigma,
        const std::map< 
                      types::boundary_id,
                      Function<space_dimension>* 
              >&                                    neumann_bc,
        const Vector<double>&                       solution,
        std::map< cell_iterator,
                  double >&                         residuals
    ) const;

    // Same as before, but without Neumann boundary conditions, i.e. 
    //  \f[
    //    \gamm_i^N \equiv 0.
    //  \f]
    void poisson_residual_error_estimate(
        const std::vector< unsigned int >&         n_gauss_points,
        const Function<space_dimension>*           rhs_fcn,
        const Function<space_dimension>*           sigma,
        const Vector<double>&                      solution,
        std::map< cell_iterator,
                  double >&                        residuals
    ) const;

    // Same as before, but with the condition that \f$\sigma \equiv 1\f$
    void poisson_residual_error_estimate(
        const std::vector< unsigned int >&         n_gauss_points,
        const Function<space_dimension>*           rhs_fcn,
        const std::map< 
                      types::boundary_id,
                const Function<space_dimension>* 
              >&                                   neumann_bc,
        const Vector<double>&                      solution,
        std::map< cell_iterator,
                  double >&                        residuals
    ) const;

    // Same as before, but with the condition that \f$\sigma \equiv 1\f$
    void poisson_residual_error_estimate(
        const std::vector< unsigned int >&         n_gauss_points,
        const Function<space_dimension>*           rhs_fcn,
        const std::map< 
                      types::boundary_id,
                      Function<space_dimension>* 
              >&                                   neumann_bc,
        const Vector<double>&                      solution,
        std::map< cell_iterator,
                  double >&                        residuals
    ) const;

    // Same as before, but with the condition that \f$\sigma \equiv 1\f$
    // and without Neumann boundary conditions.
    void poisson_residual_error_estimate(
        const std::vector< unsigned int >&         n_gauss_points,
        const Function<space_dimension>*           rhs_fcn,
        const Vector<double>&                      solution,
        std::map< cell_iterator,
                  double >&                        residuals
    ) const;

    // Same as before but now with Vector< type > for output
    void poisson_residual_error_estimate(
        const std::vector< unsigned int >&          n_gauss_points,
        const Function<space_dimension>*            rhs_fcn,
        const Function<space_dimension>*            sigma,
        const std::map< 
                      types::boundary_id,
                const Function<space_dimension>* 
              >&                                    neumann_bc,
        const Vector<double>&                       solution,
              Vector<double>&                       residuals
    ) const;

    /// Same as Before, but remove const-ness from neumann_bc map
    void poisson_residual_error_estimate(
        const std::vector< unsigned int >&          n_gauss_points,
        const Function<space_dimension>*            rhs_fcn,
        const Function<space_dimension>*            sigma,
        const std::map< 
                      types::boundary_id,
                      Function<space_dimension>* 
              >&                                    neumann_bc,
        const Vector<double>&                       solution,
              Vector<double>&                       residuals
    ) const;

    // Same as before, but without Neumann boundary conditions, i.e. 
    //  \f[
    //    \gamm_i^N \equiv 0.
    //  \f]
    void poisson_residual_error_estimate(
        const std::vector< unsigned int >&         n_gauss_points,
        const Function<space_dimension>*           rhs_fcn,
        const Function<space_dimension>*           sigma,
        const Vector<double>&                      solution,
              Vector<double>&                     residuals
    ) const;

    // Same as before, but with the condition that \f$\sigma \equiv 1\f$
    void poisson_residual_error_estimate(
        const std::vector< unsigned int >&         n_gauss_points,
        const Function<space_dimension>*           rhs_fcn,
        const std::map< 
                      types::boundary_id,
                const Function<space_dimension>* 
              >&                                   neumann_bc,
        const Vector<double>&                      solution,
              Vector<double>&                      residuals
    ) const;

    // Same as before, but with the condition that \f$\sigma \equiv 1\f$
    void poisson_residual_error_estimate(
        const std::vector< unsigned int >&         n_gauss_points,
        const Function<space_dimension>*           rhs_fcn,
        const std::map< 
                      types::boundary_id,
                      Function<space_dimension>* 
              >&                                   neumann_bc,
        const Vector<double>&                      solution,
              Vector<double>&                      residuals
    ) const;

    // Same as before, but with the condition that \f$\sigma \equiv 1\f$
    // and without Neumann boundary conditions.
    void poisson_residual_error_estimate(
        const std::vector< unsigned int >&         n_gauss_points,
        const Function<space_dimension>*           rhs_fcn,
        const Vector<double>&                      solution,
              Vector<double>&                      residuals
    ) const;


    // Isotropic materials defined by coefficients /lambda and /mu
    void linear_elasticity_residual_error_estimate(
        const std::vector< unsigned int >&          n_gauss_points,
        const Function<space_dimension>*            rhs_fcn,
        const Function<space_dimension>*            lambda,
        const Function<space_dimension>*            mu,
        const std::map< 
                      types::boundary_id,
                const Function<space_dimension>* 
              >&                                    neumann_bc,
        const Vector<double>&                       solution,
        std::map< cell_iterator,
                  double >&                         residuals
    ) const;

    void linear_elasticity_residual_error_estimate(
        const std::vector< unsigned int >&          n_gauss_points,
        const Function<space_dimension>*            rhs_fcn,
        const Function<space_dimension>*            lambda,
        const Function<space_dimension>*            mu,
        const std::map< 
                      types::boundary_id,
                const Function<space_dimension>* 
              >&                                    neumann_bc,
        const Vector<double>&                       solution,
              Vector<double>&                       residuals
    ) const;

    void linear_elasticity_residual_error_estimate(
        const std::vector< unsigned int >&          n_gauss_points,
        const Function<space_dimension>*            rhs_fcn,
        const Function<space_dimension>*            lambda,
        const Function<space_dimension>*            mu,
        const std::map< 
                      types::boundary_id,
                      Function<space_dimension>* 
              >&                                    neumann_bc,
        const Vector<double>&                       solution,
              Vector<double>&                       residuals
    ) const;

    void linear_elasticity_residual_error_estimate(
        const std::vector< unsigned int >&          n_gauss_points,
        const Function<space_dimension>*            rhs_fcn,
        const Function<space_dimension>*            lambda,
        const Function<space_dimension>*            mu,
        const Vector<double>&                       solution,
              Vector<double>&                       residuals
    ) const;

    // Special cases where lambda and mu are constant functions
    // Isotropic materials defined by coefficients /lambda and /mu
    void linear_elasticity_residual_error_estimate(
        const std::vector< unsigned int >&          n_gauss_points,
        const Function<space_dimension>*            rhs_fcn,
        const Functions::ConstantFunction<space_dimension>*    lambda,
        const Functions::ConstantFunction<space_dimension>*    mu,
        const std::map< 
                      types::boundary_id,
                const Function<space_dimension>* 
              >&                                    neumann_bc,
        const Vector<double>&                       solution,
        std::map< cell_iterator,
                  double >&                         residuals
    ) const;

    void linear_elasticity_residual_error_estimate(
        const std::vector< unsigned int >&          n_gauss_points,
        const Function<space_dimension>*            rhs_fcn,
        const Functions::ConstantFunction<space_dimension>*    lambda,
        const Functions::ConstantFunction<space_dimension>*    mu,
        const std::map< 
                      types::boundary_id,
                const Function<space_dimension>* 
              >&                                    neumann_bc,
        const Vector<double>&                       solution,
              Vector<double>&                       residuals
    ) const;

    void linear_elasticity_residual_error_estimate(
        const std::vector< unsigned int >&          n_gauss_points,
        const Function<space_dimension>*            rhs_fcn,
        const Functions::ConstantFunction<space_dimension>*    lambda,
        const Functions::ConstantFunction<space_dimension>*    mu,
        const std::map< 
                      types::boundary_id,
                      Function<space_dimension>* 
              >&                                    neumann_bc,
        const Vector<double>&                       solution,
              Vector<double>&                       residuals
    ) const;

    void linear_elasticity_residual_error_estimate(
        const std::vector< unsigned int >&          n_gauss_points,
        const Function<space_dimension>*            rhs_fcn,
        const Functions::ConstantFunction<space_dimension>*    lambda,
        const Functions::ConstantFunction<space_dimension>*    mu,
        const Vector<double>&                       solution,
              Vector<double>&                       residuals
    ) const;
  // =============================================================================
  // =============================================================================
  //                  Printer routines
  // =============================================================================
  // =============================================================================
  public:
    // Print the values of the IPF in a matrix. Specify, wether or not to print splines
    // as well as the IPF. Thus, we can either print all Splines or simply the IPF
    void printIPF(
        const std::vector< Point<dim> >& Points,
        const std::string& out_name,
        const unsigned int precision = 4,
        const bool print_splines = true,
        const bool print_IPF = true
    ) const ;

    void printIPF(
        const unsigned int n_components,
        const std::vector< Point<dim> >& Points,
        const std::string& out_name,
        const unsigned int precision = 4,
        const bool print_splines = true,
        const bool print_IPF = true
    ) const ;

    void printIPF(
        const std::string& out_name,
        const unsigned int precision = 4,
        const bool print_splines = true,
        const bool print_IPF = false
    ) const ;

    // Print the wireframe of the parametric T-mesh mapped to the parametric domain,
    // that is, every line in the grid is mapped to a line in the domain.
    // Each line in the grid is seperated in a set of N intermediate points, and each
    // point is mapped to the parametric domain.
    // The result is stored in a matrix, with dimension (n x 2*dim), where the first dim
    // columns correspond to the map of the first vertex from a line, and the last dim columns
    // correspond o the last vertex from the same line. The amount of rows is determined
    // by the amount intermediate points N and the amount of vertices, i.e. n = N * n_vertices.
    void print_IPF_wireframe(
        const std::string& out_name,
        const unsigned int precision = 4,
        const unsigned int n_intermediate_points = 11
    ) const ;

    void print_IPF_surface_wireframe(
        const std::string& out_name,
        const unsigned int precision = 4,
        const unsigned int n_intermediate_points = 11
    ) const ;

    // Same as above, but without intermediate points
    void print_IPF_vertices(
        const std::string& out_name,
        const unsigned int precision = 4
    ) const ;

    // Same as print_IPF_vertices, but without the mapping to the parametric domain.
    void print_grid(
        const std::string& out_name,
        const unsigned int precision = 4
    ) const ;

    // Same as above, but the grid switches to the bezier T-mesh before printing.
    void print_bezier_grid(
        const std::string& out_name,
        const unsigned int precision = 4
    ) const ;

    // Convert a given knot vector to a string for output
    const std::string kv_to_string(
        const std::vector<double>& kx
    ) const;

    template<int soldim = 0>
    void generate_mesh_file(
      const std::string&  out_name, 
      const bool          parametric    = true,
      const unsigned int  precision     = 8,
      const std::map<unsigned int, Point<soldim>>& 
                          vertex_values = {}
    ) const ;

  // =============================================================================
  // =============================================================================
  //                 Manipulating the grid 
  // =============================================================================
  // =============================================================================
  // To manipulate the grid, the user is [per deal.II standard] supposed to call
  // execute_coarsening_and_refinement(), execute_refinement() or execute_coarsening()
  // to manipulate the grid with previously flagged cells.
  //
  // In order to keep analysis suitability, we have to employ the coarse neighborhood
  // introduced by Morgenstern et Al. and refine each cell within as long as there
  // are coarser cells within before refining the initially marked cells.
  //
  // Before refining, however, we need to save the splines, which will be affected,
  // i.e. the splines that have to change their local knot vector. The splines will be
  // updated after refinement is done, together with generating new splines
  public:
    // From the marked cells, calculate the coarse neighborhood, then the affected splines,
    // i.e. the splines that will become inactive and from which we have to generate new
    // splines.
    void refine_fixed_number(
      const Vector< double >& residuals, 
      const double            percentile
    );

    void refine_fixed_number(
      const std::map< 
              active_cell_iterator, 
              double>&        residuals,
      const unsigned int      percentile
    );

    void execute_coarsening_and_refinement(
    ) override;

    // Refine every cell within this grid 'times' time.
    void refine_global(
        const unsigned int times = 1
    );

    // Raise the polynomial degree of this T-mesh, and hence the TSplines by 'times' time
    // in each direction
    void degree_elevate_global(
        const unsigned int times = 1
    );

    // Elevate the polynomial degree of this T-mesh in a specified direction by 'times'
    void degree_elevate(
        const unsigned int d,
        const unsigned int times = 1
    );

    // set the boundary dofs according to the boundary ids
    void set_boundary_dofs();

    // refine all cells, that were identified as bezier cells.
    void refine_bezier_elements();

    // Coarsen all previously refined bezier cells
    void coarsen_bezier_elements();

    // Calculate the bezier extraction operator for each cell
    // This call has to be made from the user. The bool is_bezier_mesh acts as a
    // safeguard to prevent its use if the mesh has not previously refined.
    // However, the user can force its computation, in which case this function
    // calls refine_bezier_elements(). This feature is to be used with caution
    void compute_extraction_operators(
        // const bool& force = false
    );

    // A function to test if we have found all bezier elements. If the
    // mesh is a bezier mesh, then every cell is supported by
    //      \prod_{d=0}^{dim-1}(p[d]+1)
    // splines
    void test_bezier();

  // =============================================================================
  // =============================================================================
  //                Helper section 
  // =============================================================================
  // =============================================================================
  // Other helper functions to help build the member variables, that are purely
  // for internal use only
    // Find all bezier elements, by searching T-Junction extensions
    virtual void find_bezier_elements() = 0; 

    // Compute the IEN array
    // This function is automatically called, whenever refine_bezier_elements() is called.
    // This assures input faults, as the IEN_array and the extraction_operators can only
    // be used for bezier meshes
    void compute_IEN_array();

    // Find the coarse neighborhood of the marked cells
    virtual void get_coarse_neighborhood(
      std::vector<active_cell_iterator>& coarse_neighborhood
    ) = 0;

    // get splines that are to be set in-active and from which we generate new splines
    // 
    // The idea of this function is rather simple: for each cell marked for refinement in the coarse
    // neighborhood, we collect the TSplines around the cell with the following conditions, where we
    // assume a (single) cut in direction d
    //    1. Check if TSpline ts has the line segment in its d-th local knot vector
    //    2. Check if support(ts) intersects with the refined cell
    // In detail, consider the following 2D example for sake of simplicity
    //
    //              + - - - - - + - - - - - + - - - - - +
    //              |           |           |           |
    //              |           + = = + = = +           |
    //              |           |     |     |           |
    //   y{l+2}     + - - - - - + = = + = = + - - - - - +
    //              |           |           |           |
    //              |           |           |           |
    //              |           |           |           |
    //   y{l+1}     + ~ ~ ~ ~ ~ X - - - - - X - - - - - +
    //              |           | marked    |           |
    //              |           | for y-cut |           |
    //              |           |           |           |
    //   y{l}       + - - - - - X - - - - - X - - - - - +
    //              |           |           |           |
    //              |           |           |           |
    //              |           |           |           |
    //   y{l-1}     + - - - - - + - - - - - + - - - - - +
    //              x{k-2}      x{k-1}   x{k}=x{k+1}      x{k+2}
    //
    // In this example we consider the horizontal lines to be the anchors, i.e.
    // px is even and py is odd. Every TSpline that has the segment (x{k-1}, x{k})
    // in the middle of their (first) index vector with support on the cell
    // then needs to be updated, as well as their corresponding ControlPoints, e.g.
    // the splines located on the anchors marked by =
    // Note, that we actually check if the x-segment of the anchor our spline is located on
    // is contained within the x-segment of the cell marked for refinement, i.e. if
    // A = a{1} x ... x a{d} is the anchor of a spline and C = c{1} x ... x{d} is the
    // cell marked for refinement along axis k, then there must hold a{k} \subset c{k}.
    //
    // Note further, that this criterion ensures that splines adjacent to the cell will
    // not be refined, e.g. the splines on the lines marked by ~
    void get_affected_splines(
      const std::vector<active_cell_iterator>& coarse_neighborhood,
      std::vector< ts_ptr >& to_be_updated
    ) const ;

    // from the splines returned by the previous function, get all splines with support
    // on the specified face, and sort them along direction dr, where each sub-array
    // is aligned with less then p[dr] + 1 splines
    std::vector<std::vector< ts_ptr >> get_affected_splines(
      const active_face_iterator&   face,
      const std::vector<ts_ptr>&    to_be_updated,
            unsigned int&           dr
    ) const ;

    // merge kv1 into kv2
    void merge_kv(
      const std::vector<double>& kv1,
            std::vector<double>& kv2,
      const unsigned int d
    ) const ;

    // return a vector of the inner faces from a cell. Here, 'cell' is a cell from
    // the coarse neighborhood of the initially marked cells, and was refined
    // before new TSplines were computed.
    const std::vector< face_iterator > get_inner_faces(
      const cell_iterator& cell
    ) ;

    // Check whether or not a given knot vector\f$\mathbf{x}\f$ with a given degree \f$p\f$ is
    // \f$p\f$-open and append necessary knots at the front or remove superfluous
    // knots from the interior. The condition we enforce on any 1D knot
    // vector is 
    // \f[
    //      x_0 = ... = x_p < x_{p+1} <= ... <= x_n < x_{n+1} = ... = x_{n+p+1}
    // \f]
    void ensure_open_knot_vector(
        std::vector< double >& kv,
        const unsigned int p
    ) const ;

    // setup the multiplicty of face array
    void setup_mof();

    // Get the vector valued distance between a cell and a point
    void vector_dist(
        const cell_iterator& cell,
        const Point<dim>& z,
              Point<dim>& distance
    ) const;

    // Get the vector valued distance between two cells
    void vector_dist(
        const cell_iterator& cell1,
        const cell_iterator& cell2,
              Point<dim>& distance
    ) const;
  }; // TS_TriangulationBase<dim, spacedim>
  
  // This class serves the purpose to define TS_Triangulation with 2 template arguments. 
  // As the actual implementation will be split for 2D and 3D in template specializations,
  // this instantiation will only throw errors whenever the user tries to create 
  // objects of type TS_Triangulation<dim, spacedim> with dim > 3 or dim < 2.
  template<int, int>
  class TS_Triangulation;
  
  // TS_Triangulation specialization with dim = 2
  template<int spacedim>
  class TS_Triangulation<2, spacedim> : public TS_TriangulationBase<2, spacedim> {
  private:
    static constexpr unsigned int dimension       = 2; 
    static constexpr unsigned int space_dimension = spacedim;

    // Import namespaces fro base class
    using face_iterator           = typename TS_TriangulationBase<dimension, space_dimension>::face_iterator;
    using cell_iterator           = typename TS_TriangulationBase<dimension, space_dimension>::cell_iterator;
    using active_face_iterator    = typename TS_TriangulationBase<dimension, space_dimension>::active_face_iterator;
    using active_cell_iterator    = typename TS_TriangulationBase<dimension, space_dimension>::active_cell_iterator;

    using anchor_bounds    = typename TS_Triangulation<dimension, space_dimension>::anchor_bounds;

    using TSpline          = typename TS_Triangulation<dimension, space_dimension>::TSpline;
    using ts_ptr           = typename TS_Triangulation<dimension, space_dimension>::ts_ptr;
    using ControlPoint     = typename TS_Triangulation<dimension, space_dimension>::ControlPoint;

  // =============================================================================
  // =============================================================================
  //                constructor  
  // =============================================================================
  // =============================================================================
  public:
    TS_Triangulation(
      const std::vector< std::vector<double> >&      knot_vectors , 
      const std::vector< unsigned int >&             degree ,
      const std::vector< Point<spacedim+1> >&        wcps 
    );

    TS_Triangulation(
      const IPF_Data<dimension, space_dimension>& data
    );

    TS_Triangulation(
    ) ;

    ~TS_Triangulation() = default;

    TS_Triangulation(
      const TS_Triangulation<dimension, space_dimension>& other
    ) = delete;

    TS_Triangulation<dimension, space_dimension>& operator=(
      TS_Triangulation<dimension, space_dimension>&& other
    ) noexcept = delete;

    void create_triangulation(
      const std::vector< std::vector<double> >& knot_vectors,
      const std::vector< unsigned int >&        degree,
      const std::vector< Point<spacedim+1> >&   wcps
    );

    void create_triangulation(
      const IPF_Data<dimension, space_dimension>& data
    );

  // =============================================================================
  // =============================================================================
  //               Overrides: 
  // =============================================================================
  // =============================================================================
  private:  
    // Find all bezier elements, by searching T-Junction extensions
    void find_bezier_elements() override; 

    // Find the coarse neighborhood of the marked cells
    void get_coarse_neighborhood(
      std::vector<active_cell_iterator>& coarse_neighborhood
    ) override;
 
  public: 
    // get Bezier coefficients on a specific cell
    const FullMatrix<double> get_bezier_coefficients(
        const active_cell_iterator& cell
    ) const override;

    // get Bezier coefficients on a specific cell at a speecific face
    // Note, that the face index corresponds to a face at the boundary of
    // this cell
    const std::vector< FullMatrix<double> >
      get_bezier_coefficients(
        const active_cell_iterator& cell,
        const unsigned int face
    ) const override;
  }; // TS_Triangulation<2, spacedim>

  // TS_Triangulation specialization with dim = 3
  template<int spacedim>
  class TS_Triangulation<3, spacedim> : public TS_TriangulationBase<3, spacedim> {
  private:
    static constexpr unsigned int dimension       = 3; 
    static constexpr unsigned int space_dimension = spacedim;

    using face_iterator           = typename TS_TriangulationBase<dimension, space_dimension>::face_iterator;
    using cell_iterator           = typename TS_TriangulationBase<dimension, space_dimension>::cell_iterator;
    using active_face_iterator    = typename TS_TriangulationBase<dimension, space_dimension>::active_face_iterator;
    using active_cell_iterator    = typename TS_TriangulationBase<dimension, space_dimension>::active_cell_iterator;

    using anchor_bounds    = typename TS_Triangulation<dimension, space_dimension>::anchor_bounds;

    using TSpline          = typename TS_Triangulation<dimension, space_dimension>::TSpline;
    using ts_ptr           = typename TS_Triangulation<dimension, space_dimension>::ts_ptr;
    using ControlPoint     = typename TS_Triangulation<dimension, space_dimension>::ControlPoint;

  // =============================================================================
  // =============================================================================
  //                constructor  
  // =============================================================================
  // =============================================================================
  public:
    TS_Triangulation(
      const std::vector< std::vector<double> >&      knot_vectors , 
      const std::vector< unsigned int >&             degree ,
      const std::vector< Point<spacedim+1> >&        wcps 
    ) ;

    TS_Triangulation(
      const IPF_Data<dimension, space_dimension>& data
    ) ;

    TS_Triangulation() ;

    ~TS_Triangulation() = default;

    TS_Triangulation(
      const TS_Triangulation<dimension, space_dimension>& other
    ) = delete;

    TS_Triangulation<dimension, space_dimension>& operator=(
      TS_Triangulation<dimension, space_dimension>&& other
    ) noexcept = delete;

    void create_triangulation(
      const std::vector< std::vector<double> >&          knot_vectors,
      const std::vector< unsigned int >&                 degree,
      const std::vector< Point<spacedim+1> >&            wcps
    );

    void create_triangulation(
      const IPF_Data<dimension, space_dimension>& data
    );

  // =============================================================================
  // =============================================================================
  //               Overrides: 
  // =============================================================================
  // =============================================================================
  private:  
    // Find all bezier elements, by searching T-Junction extensions
    void find_bezier_elements() override; 

    // Find the coarse neighborhood of the marked cells
    void get_coarse_neighborhood(
      std::vector<active_cell_iterator>& coarse_neighborhood
    ) override;
 
  public: 
    // get Bezier coefficients on a specific cell
    const FullMatrix<double> get_bezier_coefficients(
        const active_cell_iterator& cell
    ) const override;

    // get Bezier coefficients on a specific cell at a speecific face
    // Note, that the face index corresponds to a face at the boundary of
    // this cell
    const std::vector< FullMatrix<double> >
      get_bezier_coefficients(
        const active_cell_iterator& cell,
        const unsigned int face
    ) const override;
  
  }; // TS_Triangulation<3, spacedim>
  

  template<int dim, int spacedim = dim, int soldim = 1>
  class TSValuesBase {
  protected:
    //  Define iterator shortcuts to be used throughout this class -- these names are
    //  according to the deal.II library and directly derived from there
    using face_iterator           = TriaIterator<dealii::TriaAccessor<dim-1, dim, dim>>;
    using cell_iterator           = TriaIterator<dealii::CellAccessor<dim, dim>>;
    using active_face_iterator    = TriaActiveIterator<dealii::TriaAccessor<dim-1, dim, dim>>;
    using active_cell_iterator    = TriaActiveIterator<dealii::CellAccessor<dim, dim>>;

    static constexpr unsigned int dimension          = dim;
    static constexpr unsigned int space_dimension    = spacedim;
    static constexpr unsigned int solution_dimension = soldim;
    int                        face_no             = -1;
    unsigned int               offset              = 0;
    unsigned int               face_quadrature_points; 
    std::array<unsigned int, dim> quadrature_points_per_face;


    Table<2, double>           bernstein_values;        // Define values on reference element
    Table<2, Tensor<1, dim> >  bernstein_grads;         // Define gradients on reference element
    Table<2, Tensor<2, dim> >  bernstein_hessians;      // Define hessians on reference element
    // Table<2, Tensor<3, dim> >  bernstein_3rd;

    unsigned int               dofs_per_cell;                 // Number of splines with support on one cell
    unsigned int               bernstein_quadrature_points;   // number of quadrature points used for bernstein polynomial
    unsigned int               shape_quadrature_points;       // number of quadrature points used for shape functions
    unsigned int               total_quadrature_points;       // number of quadrature points used in total [incl. subfaces]

    Table<1, double>           quadrature_weights;
    Table<1, Point<dim> >      quadrature_points;

    UpdateFlags                flags;

    std::vector< unsigned int >      degrees;
    std::vector< unsigned int >      n_gauss_points;

    Table<2, double>                 shape_values;
    Table<2, Tensor<1, spacedim> >   shape_grads;
    Table<2, Tensor<2, spacedim> >   shape_hessians;
    Table<2, Tensor<3, spacedim> >   third_derivative;

    std::vector< Point<spacedim> >   mapped_quadrature_points;

    std::vector< DerivativeForm< 1, dim, spacedim> >    J; // jacobian on cell
    std::vector< DerivativeForm< 1, spacedim, dim> >    I; // inverse jacobian on cell
    std::vector< Tensor<1, spacedim, Tensor<2, dim> > > H; // Derivative of Jacobian, i.e. Hessian of transformation
    std::vector< Tensor<1, dim, Tensor<2, spacedim> > > HI; // Hessian of inverse transformation

    std::vector< double >                               dx; // For integration; represents JxW from FEValues
  
    std::vector< Tensor<1, spacedim> >                  normals; // For the normals on faces.

    TS_TriangulationBase<dim, spacedim>*    tria;

    Quadrature< dim >                   cell_quadrature;
    std::array<Quadrature<dim>, dim>    face_quadratures;
  public:

    TSValuesBase(
      const TS_TriangulationBase<dim, spacedim>* tria,
      const unsigned int                    n_gauss_points,
      const UpdateFlags                     flags,
      const int&                            face_no = -1
    );

    TSValuesBase(
      const TS_TriangulationBase<dim, spacedim>* tria,
      const std::vector<unsigned int>&      n_gauss_points,
      const UpdateFlags                     flags,
      const int&                            face_no = -1
    );

    // Do we want these constructors? 
    // We might need it for hp capabilities. 
    // But do we need hp-capabilities? Probably not ...
    TSValuesBase() = default;
    TSValuesBase(const TSValuesBase& other) = default;
    TSValuesBase& operator=(const TSValuesBase& other) = default;

    const double& shape_value(
      const unsigned int function_no,
      const unsigned int point_no
    ) const ;

    const Tensor<1, spacedim>& shape_grad(
      const unsigned int function_no,
      const unsigned int point_no
    ) const ;

    const Tensor<2, spacedim>& shape_hessian(
      const unsigned int function_no,
      const unsigned int point_no
    ) const ;

    const Tensor<3, spacedim>& shape_3rd_derivative(
      const unsigned int function_no,
      const unsigned int point_no
    ) const ;

    const DerivativeForm<1, dim, spacedim>& jacobian(
      const unsigned int point_no
    ) const ;

    const DerivativeForm<1, spacedim, dim>& inverse_jacobian(
      const unsigned int point_no
    ) const ;

    const Tensor<1, spacedim, Tensor<2, dim>>& jacobian_grad(
      const unsigned int point_no
    ) const ;

    const Tensor<1, dim, Tensor<2, spacedim>>& inverse_jacobian_grad(
      const unsigned int point_no
    ) const ;

    const double& JxW(
      const unsigned int point_no
    ) const ;

    std_cxx20::ranges::iota_view<unsigned int, unsigned int>
            dof_indices() const;

    std_cxx20::ranges::iota_view<unsigned int, unsigned int>
            quadrature_point_indices() const;

    const Point<spacedim>& quadrature_point(
      const unsigned int q
    ) const ;

    const Point<dim>& quadrature_point_reference(
      const unsigned int q
    ) const ;

    const unsigned int& n_quadrature_points_per_cell() const ;

    const unsigned int& n_quadrature_points_per_face() const;

    const std::vector<unsigned int>& n_q_points_per_cell() const;

    unsigned int n_dofs_per_cell() const;


    virtual const Tensor<1, spacedim>& normal_vector(
      const unsigned int q
    ) const = 0;

    std::pair<unsigned int, unsigned int> system_to_component_index(
      const unsigned int index
    ) const ;

    unsigned int component_to_system_index(
      const unsigned int component, 
      const unsigned int index
    ) const;

    const std::vector< Point<spacedim> >& get_quadrature_points(
    ) const;

    const unsigned int& n_quadrature_points() const;


  // protected member functions to be used for 
  // derived classes, in particular TSFaceValuesBase
  protected:
    // Initialize shape tables
    void init_shape_tables(
      const unsigned int n_subfaces = 1
    );
    
  protected: 
    // private member functions to help define the necessary tables
    
    // Set all tables to the corresponding size, determined by
    // the amount of gauss points, and the degrees
    void init_bernstein_tables();

    
    // Pre-emptively set the bernstein tables to be used for 
    // Bezier decomposition later.
    void set_bernstein_tables();
  }; // TSValuesBase class

  // This is a dummy-class, that enables the proper
  // layout for TSValues. This class is mostly empty,
  // some member functions are defined but will throw
  // an error if the user tries to call it.
  template<int dim, int spacedim = dim, int soldim = 1>
  class TSValues : public TSValuesBase<dim, spacedim, soldim> {
  private:
    // internal usage
    static constexpr unsigned int dimension          = dim; 
    static constexpr unsigned int space_dimension    = spacedim;
    static constexpr unsigned int solution_dimension = soldim;

    // Import typenames from base class
    using face_iterator        
            = typename TSValuesBase<dimension, space_dimension, solution_dimension>::face_iterator; 
    using cell_iterator        
            = typename TSValuesBase<dimension, space_dimension, solution_dimension>::cell_iterator; 
    using active_face_iterator 
            = typename TSValuesBase<dimension, space_dimension, solution_dimension>::active_face_iterator; 
    using active_cell_iterator 
            = typename TSValuesBase<dimension, space_dimension, solution_dimension>::active_cell_iterator; 

  public: 
    // Default standard constructor
    TSValues() = default;

    // Initialize this object using the constructor
    // of the base class for isotropic polynomials
    TSValues(
      const TS_TriangulationBase<dim, spacedim>*  tria,
      const unsigned int                      n_gauss_points,
      const UpdateFlags                       flags
    );

    // Initialize this object using the constructor
    // of the base class for anisotropic polynomials
    TSValues(
      const TS_TriangulationBase<dim, spacedim>*  tria,
      const std::vector<unsigned int>&        n_gauss_points,
      const UpdateFlags                       flags
    );

    const Tensor<1, spacedim>& normal_vector(
        const unsigned int /* q */
    ) const override {
      Assert(false, ExcMessage("Trying to access normals for Cells, which makes no sense."));
      // To silence warnings
      return this->normals[0];
    }

    // Set the values defined by flags on a specific cell
    void reinit(
      const active_cell_iterator& cell
    ) ;

    // perform a test of TSValues on each cell. 
    // Prints a comprehensive summary to a log directory, that 
    // is created if not already given
    // Note, that phi is provided in the library. To perform 
    // an extensive test including the derivatives and values of 
    // the inverse, provide the inverse to this function
    void test_geometry_mapping(
      const Function<dim>       *phi_ptr,
      const Function<spacedim>  *inv_phi_ptr = NULL,
      const long double         &TOL = 1e-15,
      const bool                 test_second_derivative = true
    ) ;

  private:
    // Set values of shape tables on a
    // cell and specific face using 
    // already defined values offset and
    // face_no
    void set_grad_tables(
      const active_cell_iterator& cell
    );

    void set_grad_and_hessian_tables(
      const active_cell_iterator& cell
    );
  }; // TSValues< dim, spacedim>


  // Define the base class for TSFaceValues:
  template<int dim, int spacedim = dim, int soldim = 1>
  class TSFaceValuesBase : public TSValuesBase<dim, spacedim, soldim> {
  protected: 
    // Are we using an anisotropic quadrature?
    bool anisotropic_quadrature = true;

  public:
    TSFaceValuesBase(
      const TS_TriangulationBase<dim, spacedim>* tria,
      const std::vector<unsigned int>&           n_gauss_points,
      const UpdateFlags                          flags
    );

    TSFaceValuesBase(
      const TS_TriangulationBase<dim, spacedim>* tria,
      const unsigned int&                   n_gauss_points,
      const UpdateFlags                     flags
    );

    TSFaceValuesBase() = default; 
    ~TSFaceValuesBase() = default;

    const Tensor<1, spacedim>& normal_vector(
        const unsigned int q
    ) const override;

  }; // TSFaceValuesBase< dim, spacedim>

  template<int dim, int spacedim = dim, int soldim = 1>
  class TSFaceValues : public TSFaceValuesBase<dim, spacedim, soldim> {
  private:
    // internal usage
    static constexpr unsigned int dimension          = dim; 
    static constexpr unsigned int space_dimension    = spacedim;
    static constexpr unsigned int solution_dimension = soldim;

    // Import typenames from base class
    using face_iterator        = 
            typename TSValuesBase<dimension, space_dimension, solution_dimension>::face_iterator; 
    using cell_iterator        = 
            typename TSValuesBase<dimension, space_dimension, solution_dimension>::cell_iterator; 
    using active_face_iterator =
            typename TSValuesBase<dimension, space_dimension, solution_dimension>::active_face_iterator; 
    using active_cell_iterator =
            typename TSValuesBase<dimension, space_dimension, solution_dimension>::active_cell_iterator; 

  public: 
    // Initialize this object using the constructor
    // of the base class for isotropic polynomials
    TSFaceValues(
      const TS_TriangulationBase<dim, spacedim>*  tria,
      const unsigned int                      n_gauss_points,
      const UpdateFlags                       flags
    );

    // Initialize this object using the constructor
    // of the base class for anisotropic polynomials
    TSFaceValues(
      const TS_TriangulationBase<dim, spacedim>*  tria,
      const std::vector<unsigned int>&        n_gauss_points,
      const UpdateFlags                       flags
    );


    // Set the values defined by flags on a specific cell
    void reinit(
      const active_cell_iterator& cell,
      const          int&         face_no
    ) ;

    // perform a test of TSValues on each cell. 
    // Prints a comprehensive summary to a log directory, that 
    // is created if not already given
    // Note, that phi is provided in the library. To perform 
    // an extensive test including the derivatives and values of 
    // the inverse, provide the inverse to this function
    void test_geometry_mapping(
      const Function<dim>       *phi_ptr,
      const Function<spacedim>  *inv_phi_ptr = NULL,
      const long double         &TOL = 1e-15,
      const bool                 test_second_derivative = false
    ) ;
    
  private:
    // Set values of shape tables on a
    // cell and specific face using 
    // already defined values offset and
    // face_no
    void set_grad_tables(
      const active_cell_iterator& cell,
      const unsigned int sub_face
    );

    void set_grad_and_hessian_tables(
      const active_cell_iterator& cell,
      const unsigned int sub_face
    );
  }; // TSFaceValues<dim, spacedim>


} // namespace dealt

#endif /* TSPLINEBASIS_H_ */
