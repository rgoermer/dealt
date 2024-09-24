/*
 * IsoparametricFunction.h
 *
 *  Created on: Jul 1, 2020
 *      Author: goermer
 */

#ifndef ISOPARAMETRICFUNCTION_H_
#define ISOPARAMETRICFUNCTION_H_

// STL header
#include <iostream>
#include <iomanip>

#include <fstream>
#include <cmath>
#include <memory>
#include <algorithm>
#include <unordered_map>

// Deal.II Header
#include <deal.II/base/function.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/derivative_form.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>

#include <deal.II/grid/manifold.h>

#include <deal.II/fe/mapping_manifold.h>

//Own Header
#include <utilities.h>
#include <tspline_function.h>



namespace dealt {
  using namespace dealii;

  template<unsigned int dim, unsigned int spacedim>
  class IsoparametricFunction : public Function<dim>{
  private:
    using face_iterator           = TriaIterator<dealii::TriaAccessor<dim-1, dim, dim>>;
    using cell_iterator           = TriaIterator<dealii::CellAccessor<dim, dim>>;
    using active_cell_iterator    = TriaActiveIterator<dealii::CellAccessor<dim, dim>>;
    using vertex_iterator         = TriaIterator<dealii::TriaAccessor<0, dim, dim>>;
    using active_vertex_iterator  = TriaActiveIterator<dealii::TriaAccessor<0, dim, dim>>;

    using TSpline          = TSplineFunction<dim, spacedim>;
    using ControlPoint     = Point<spacedim + 1>;
    using ts_ptr           = std::shared_ptr< TSpline >;

  private:
    std::vector<TSpline> TSplineCollection;

  // Constructors:
  public:
    IsoparametricFunction(const std::vector< ts_ptr >& TSplineCollection);

    IsoparametricFunction() : Function<dim>(spacedim) {}

    IsoparametricFunction(const IsoparametricFunction<dim, spacedim>& other_IPF);

    ~IsoparametricFunction(){ TSplineCollection = {}; }

    IsoparametricFunction<dim, spacedim>& operator=(IsoparametricFunction<dim, spacedim>&& other);
    IsoparametricFunction<dim, spacedim>& operator=(const IsoparametricFunction<dim, spacedim>& other);
  //Evaluate function values
  public:
    virtual double value(
        const Point<dim>&     p_d,
        const unsigned int    component = 0
        ) const override;

    virtual void vector_value(
        const Point<dim>&   p_d,
        Vector< double >&   values
        ) const override;


    virtual void vector_value_list(
        const std::vector< Point<dim> >&    points,
        std::vector< Vector<double> >&      values
        ) const override;

    Point<spacedim> point_value(
        const Point<dim>& p_d
      ) const;
  // Evaluate Function derivatives
  public:
    virtual Tensor<1,dim> gradient(
        const Point<dim>&       p_d,
        const unsigned int      component = 0
        ) const override;

    virtual void vector_gradient(
        const Point<dim>&               p_d,
        std::vector< Tensor<1,dim> >&   gradients
        ) const override;

    virtual void gradient_list(
        const std::vector< Point<dim> >&    points,
        std::vector< Tensor<1,dim> >&      gradients,
        const unsigned int component = 0
        ) const override;

    virtual void vector_gradient_list(
        const std::vector< Point<dim> >&                points,
        std::vector< std::vector< Tensor<1,dim> > >&    gradients
        ) const override;

    void jacobian(
        const Point<dim>&     p_d,
        FullMatrix<double>&   J
        ) const;

    void jacobian_vector(
        const std::vector< Point<dim> >&     p_d,
        std::vector< FullMatrix<double> >&   J
        ) const;

  // Evaluate Function's second derivatives
    virtual SymmetricTensor<2, dim> hessian(
      const Point<dim>& p, 
      const unsigned int component = 0
    ) const override;


  // Additional Functions, such as print()
  public:
    void print(
        const std::vector< Point<dim> >& points,
        const std::string add_name = ""
        ) const ;

    void print_deriv(
        const std::vector< Point<dim> >& points,
        const std::string add_name = ""
        ) const ;
  
  };

  template<int dim, int spacedim = dim>
  class IsoparametricManifold : public ChartManifold<dim, spacedim> {
  private:
    using TSpline          = TSplineFunction<dim, spacedim>;
    using ts_ptr           = std::shared_ptr< TSpline >;

  private:
    IsoparametricFunction<dim, spacedim>    geometry;

  public: 
    IsoparametricManifold() = default; 

    IsoparametricManifold(const std::vector< ts_ptr >& TSplineCollection);
    IsoparametricManifold(const IsoparametricFunction<dim, spacedim>& other_geometry);

    virtual Point<spacedim> push_forward(
      const Point< dim >& p
    ) const override; 

    // Compute pull_back using newton's method solving 
    // phi(x) = \hat{x}
    // for x.
    virtual Point<dim> pull_back(
      const Point< spacedim >& p
    ) const override; 

    virtual DerivativeForm<1, dim, spacedim> push_forward_gradient(
      const Point< dim >& p
    ) const override;

    virtual std::unique_ptr<Manifold<dim, spacedim>> clone(
    ) const override;
  }; 


} // namespace dealt

#endif /* ISOPARAMETRICFUNCTION_H_ */
