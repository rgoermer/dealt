/*
 *
 *  Created on: Apr 1, 2020
 *      Author: goermer
 */

#ifndef TSPLINEFUNCTION_H_
#define TSPLINEFUNCTION_H_

#include <memory>
#include <variant>

#include <utilities.h>

#include <deal.II/base/function.h>
#include <deal.II/base/auto_derivative_function.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/base/exceptions.h>

#include <utility>


namespace dealt {
  using namespace dealii;

  template<int dim, int spacedim>
  class TSplineFunction : public Function<dim>{

  // Type definitions for readability
  private:
    //  Define iterator shortcuts to be used throughout this class -- these names are
    //  according to the deal.II library and directly derived from there
    using face_iterator           = TriaIterator<dealii::TriaAccessor<dim-1, dim, dim>>;
    using cell_iterator           = TriaIterator<dealii::CellAccessor<dim, dim>>;
    using active_face_iterator    = TriaActiveIterator<dealii::TriaAccessor<dim-1, dim, dim>>;
    using active_cell_iterator    = TriaActiveIterator<dealii::CellAccessor<dim, dim>>;

    using vertex_iterator         = TriaIterator<dealii::TriaAccessor<0, dim, dim>>;
    using active_vertex_iterator  = TriaActiveIterator<dealii::TriaAccessor<0, dim, dim>>;


    // A shorthand notation for TSplineFunction
    using TSpline                 = TSplineFunction<dim, spacedim>;

    // A sorthand notation for the anchor bounds as a tupel of points
    using anchor_bounds = std::pair<Point<dim>, Point<dim>>;

    // Each Spline is associated with a ControlPoint. The entirety of TSplines then gives rise to the
    // Isoparametric Function. A ControlPoint is simply a Point in space extended with its weight, yielding
    // a Point of dimension (spacedim + 1)
    using ControlPoint = Point<spacedim + 1>;
    using uint_vec = std::vector<unsigned int>;
    using uint = unsigned int;

  // private members of this class
  private:

    // The attributes we need for refinement and definition of the current Spline.
    unsigned int  lev;             // level for refinement
    unsigned int  level_index;     // index of spline on a specific level
    unsigned int  ind;             // index of spline over all generated splines
    bool          active = true;   // Is the spline used or not?
    ControlPoint  cp;              // We store the control point within the TSpline data
    anchor_bounds anchor;

    // local knot vectors of this splines
    std::vector< std::vector<double> >  kv;

    // multiplicities of each distinct knot specified in kv.
    // This is later used for the evaluation of the TSpline
    // at given points
    std::vector< std::vector<int> >     mult;

    // A pointer to the father of the TSpline. This variable is not really used yet,
    // however, it is intended to also allow coarsening of the mesh, by deleting
    // respective tsplines. Since there has not been done any testing with it, and
    // coarsening of a TS_Triangulation is currently unavailable, this variable will
    // not be accessible from outside this class.
    std::shared_ptr< TSpline >          father1 = nullptr;

    // To sort the TSplines in a usefull order, we can use the barycenter of a TSpline.
    // Each coordinate is the sum of the local kv, and can thus be used to uniquely
    // sort every TSpline from, e.g. lower left to upper right, i.e. the lexicographical
    // order.
    Point<dim>    barycenter;

  // Some variables that will / will not change for each class.
  private:
    // The vector of polynomial degrees, with s_degree.size() == dim
    static std::vector<uint> s_degree;

    // static_id to define interal_id in each instantiation
    static int               s_id;

  // Public cal to set static variable
  public:
    // Set the degree for all TSpline instantiations
    static void setSplineData(std::vector<uint> degree){
      s_degree  =  degree;
    }

  // Constructors & Destructors
  public:
    TSplineFunction(
        const std::vector< std::vector< double > >& other_kv,
        const ControlPoint& other_cp,
        const std::shared_ptr<TSpline>& other_father1 = nullptr
        );

    TSplineFunction(
        const std::vector< std::vector< double > >& other_kv,
        const anchor_bounds& other_anchor,
        const ControlPoint& other_cp,
        const std::shared_ptr<TSpline>& other_father1
        );

    // delete default constructor
    TSplineFunction() = default;

    // Copy-Constructor
    TSplineFunction(const TSplineFunction<dim, spacedim>& TF) :
          Function<dim>(1){
      cp            = TF.cp;
      anchor        = TF.anchor;
      kv            = TF.kv;
      mult          = TF.mult;
      ind           = TF.ind;
      level_index   = TF.level_index;
      active        = TF.active;
      lev           = TF.lev;
    }

    ~TSplineFunction() = default;

    // Move-Constructor
    TSplineFunction<dim, spacedim>& operator=(TSplineFunction<dim, spacedim>&& other) noexcept;

    // Copy-Constructor
    TSplineFunction<dim, spacedim>& operator=(const TSplineFunction<dim, spacedim>& other) = default;

  // Public functions to access membaer variables
  public:
    // get control point of spline
    ControlPoint get_cp() const { return cp; }

    // get index over all levels of this splines
    int index() const {return ind;}

    // get current level of tspline
    int level() const {return lev;}

    // Check if this spline is active
    bool is_active() const {return active;}

    // de-activate or activate current spline
    void set_active(bool set) {this -> active = set;}

    // set level specific index of spline
    void set_level_index(unsigned int ind){this -> level_index = ind;}

    // get level specific index of spline
    unsigned int get_level_index(){return this -> level_index;}

    // return complete anchor
    anchor_bounds get_anchor() const { return anchor; }

    // return a pair that represents the anchor components in the specified direction
    std::pair<double, double> get_anchor(const int& d) const;

    // Return the multiplicity at a given direction
    std::vector< int > get_multiplicities(const int& d) const {
      return mult.at(d);
    }

    // clear father pointer (publicly accessible)
    void clear_father(){father1 = nullptr;}

    // return the local knot vector in dimension d /in [0, d) and store
    // it in a predefined vector
    void get_local_kv(std::vector<double>& kv_d,
               unsigned int d) const ;

    // return the local knot vector in dimension d
    std::vector<double> get_local_kv(int d) const {
      return kv.at(d); 
    }

    // Return complete set of knot vectors
    std::vector< std::vector <double> > get_kv() const {return kv; }

    // Return the complete barycenter of this spline
    Point<dim> get_barycenter() const {return barycenter;}

    // Return the d-th coordinate of the barycenter.
    double get_barycenter(int d) const {return barycenter(d);}

    // Return the size of the support of the current TSpline
    double get_support_volume();
  
  // Now, comes the part, where we simply evaluate the B-Splines at a given point;
  // by default value_list(...) uses a std::vector<double> to store the value at each point
  // in points. It does so by calling the value() method for each point. In order to also pass a
  // deal.ii Vector<double> we override the function accordingly. The standard method than can be called
  // in a similar way.
  //
  // Note, that for value and derivative evaluation of Splines at any specific point, we use
  // the algorithms specified in L. Piegl & W. Tiller, The NURBS Book, as the Tsplines are [locally]
  // just a tensor product of BSplines. The implementations are declared private, and come after
  // the current section
  public:
    virtual double value(
        const Point<dim>& p_d,
        const unsigned int component = 0
        ) const override ;

    virtual void value_list(
        const std::vector< Point<dim> >&  points,
        std::vector< double >&            values,
        const unsigned int                component = 0
        ) const override ;

      // Suppose, we are given the TSpline T(x) = B1(x1)*...*B{d}(x{d}),
      // where of course T:R^d -> R. The gradient is then given by
      //  grad[T(x)] = (B1'(x1)*B2(x2)* ... * B{d}(x{d}),
      //                B1(x1)*B2'(x2)* B3(x3)* ... * B{d}(x{d})
      //                          .
      //                          .
      //                          .
      //                B1(x1)* ... * B{d-1}(x{d-1}) * B{d}'(x{d}))
      virtual Tensor<1,dim> gradient(
          const Point<dim>& p,
          const unsigned int component = 0
          ) const override;

      // Compute the laplacian of the TSpline
      virtual double laplacian(
          const Point<dim>& p,
          const unsigned int component = 0
          ) const override;


      // Compute the hessian of the TSpline
      virtual SymmetricTensor<2, dim> hessian(
          const Point<dim>& p,
          const unsigned int component = 0
          ) const override;

  // Here follow the B-Spline recursion formula. We will need several functions for an efficient
  // implementation. We use the following source for implementation:
  //    L. Piegl, W. Tiller, "The NURBS Book", 2nd Ed., Algorithms A2.1 - A2.4
  // The first two algorithms might be unnecessary here, since the evaluation of a single Basis
  // function is done in the function BasisFunsSingle
  private:
    double BasisFunSingle1D(
        const double u,
        const unsigned int d) const;

    // Evaluate the k-th derivative of the BSpline specified by the d-th knot vector and evaluate at point u
    double BasisFunDerivSingle1D(
        const unsigned int k,
        const double u,
        const unsigned int d
        ) const;

  // Other functions to help us with our calculations
  public:
    // Returns wether or not this spline has support on a specified cell. Note, that this function
    // returns true, if the intersection of the TSpline with the cell lies on the boundary of the cell.
    bool has_support(const active_cell_iterator& cell) const ;
    bool has_support(const cell_iterator& cell) const ;

    // Same as above but with faces
    bool has_support(const active_face_iterator& face) const ;

    // Same as above but with points
    bool has_support(const Point<dim>& point) const ;

    // Checks if the specified point lies on the boundary of the support
    bool at_boundary(const Point<dim>& point) const;

    // Checks if a given knot is already given in the d-th local knot vector of this spline
    bool is_present(const double& knot, const int& d) const ;

    // Compute the extended kv of the d-th local knot vector, i.e.
    // insert the last and first knot until each have multiplicity s_degree[d] + 1.
    // The number of insertions to the front is stored in nf, and the number
    // of insertions to the end is stored in ne.
    std::vector<double> compute_extended_kv(
        const unsigned int& d,
        unsigned int& nf,
        unsigned int& ne
    ) const ;

    // Compute the bezier extraction operator row of the d-th local knot vector
    // for the cell specified by bezier_cell after inserting necessary interior_knots
    Vector<double> compute_be_operator(
        const std::vector< double >& interior_knots,
        const unsigned int& d,
        const active_cell_iterator& bezier_cell
    ) const ;

    void compute_be_operator(
        const std::vector< double >               &interior_knots,
        const unsigned int                        &d,
        const std::vector< active_cell_iterator > &cell_list,
              std::vector< Vector< double > >     &ops
    ) const ;

  // Some print routines
  public:
    void print_spline(const std::vector<Point<dim>>& Points,
                    const std::string& add = "") const;
  
  // The barycenter will be computed upon initialization
  private:
    void calc_barycenter();
  };

} // namespace dealt



#endif /* INC_TSPLINEFUNCTION_H_ */
