/*
 * IsoParametricFunction.cc
 *
 *  Created on: Jul 1, 2020
 *      Author: goermer
 */


#include <isoparametric_function.h>

namespace dealt {
  template class IsoparametricFunction<1,1>;
  template class IsoparametricFunction<2,2>;
  template class IsoparametricFunction<3,3>;

  template class IsoparametricFunction<1,2>;
  template class IsoparametricFunction<1,3>;
  template class IsoparametricFunction<2,3>;

  template class IsoparametricFunction<2,4>;

  template class IsoparametricManifold<1,1>;
  template class IsoparametricManifold<2,2>;
  template class IsoparametricManifold<3,3>;

  template class IsoparametricManifold<1,2>;
  template class IsoparametricManifold<1,3>;
  template class IsoparametricManifold<2,3>;


  // Constructors:
  template<unsigned int dim, unsigned int spacedim>
  IsoparametricFunction<dim, spacedim>::IsoparametricFunction(
      const std::vector< ts_ptr >& TSplineCollection ) : Function<dim>(spacedim)
  {
#ifdef DEBUG
    unsigned int n_splines = TSplineCollection.size();
    Assert(n_splines > 0, ExcZero());
#endif
    for (const auto& ts : TSplineCollection)
      this -> TSplineCollection.push_back( *ts );
  }

  template<unsigned int dim, unsigned int spacedim>
  IsoparametricFunction<dim, spacedim>::IsoparametricFunction(
      const IsoparametricFunction<dim, spacedim>& other_IPF
      ) : Function<dim>(spacedim)
  {
    this -> TSplineCollection = other_IPF.TSplineCollection;
  }

  // move-assign:
  template<unsigned int dim, unsigned int spacedim>
  IsoparametricFunction<dim, spacedim>& IsoparametricFunction<dim, spacedim>::operator=(
                  IsoparametricFunction<dim, spacedim>&& other)
  {
    if (this == &other)
      return *this;

    TSplineCollection = other.TSplineCollection;

    other.TSplineCollection = {};

    return *this;
  }

  // copy-assign
  template<unsigned int dim, unsigned int spacedim>
  IsoparametricFunction<dim, spacedim>& IsoparametricFunction<dim, spacedim>::operator=(
                  const IsoparametricFunction<dim, spacedim>& other)
  {
    if (this == &other)
      return *this;

    TSplineCollection = other.TSplineCollection;
    return *this; 
  }

  // IPF-Values:
  template<unsigned int dim, unsigned int spacedim>
  double IsoparametricFunction<dim, spacedim>::value(
      const Point<dim>& p_d,
      const unsigned int component
      ) const
  {
    AssertIndexRange(component, spacedim);
    
    double out = 0;
    double D = 0;
    for (const auto& t : TSplineCollection){
      double currVal = t.value(p_d);
      const Point<spacedim+1>& currCP = t.get_cp();
      out += currVal*currCP[component];
      D += currVal*currCP[spacedim];

    }
    return out/D;
  }
  
  template<unsigned int dim, unsigned int spacedim>
  Point<spacedim> IsoparametricFunction<dim, spacedim>::point_value(
    const Point<dim>& p_d
  ) const {
    // Compute the value in the weighted space:
    Point<spacedim+1> val;
    for (const auto& t : TSplineCollection) {
      double t_val = t.value(p_d);
      const Point<spacedim + 1>& t_cp = t.get_cp();
      val += (t_val * t_cp);
    }

    // Use the perspective map to compute actual point in spacedim
    Point<spacedim> out;
    for (unsigned int d = 0; d < spacedim; d++)
      out(d) = val(d);

    out = out/val(spacedim);

    return out;
  }

  template<unsigned int dim, unsigned int spacedim>
  void IsoparametricFunction<dim, spacedim>::vector_value(
      const Point<dim>& p_d,
      Vector< double >& values
      ) const
  {
    values = Vector<double>(spacedim);
    double D = 0;

    for (const auto& t : TSplineCollection){
      double currVal = t.value(p_d);
      const Point<spacedim+1>&  currCP = t.get_cp();
      D += currVal*currCP[spacedim];

      for (unsigned int d = 0; d < spacedim; d++)
        values(d) += currVal*currCP[d];
    }

    Assert(D != 0, ExcDivideByZero());
    values /= D;
  }

  template<unsigned int dim, unsigned int spacedim>
  void IsoparametricFunction<dim, spacedim>::vector_value_list(
      const std::vector< Point<dim> >& points,
      std::vector< Vector<double> >& values
      ) const
  {
    AssertDimension(values.size(), points.size());
    for (unsigned int n = 0; n < points.size(); n++){
      vector_value(points[n], values[n]);
    }
  }


  // IPF-Derivatives:
  template<unsigned int dim, unsigned int spacedim>
  Tensor<1,dim> IsoparametricFunction<dim, spacedim>::gradient(
      const Point<dim>& p_d,
      const unsigned int component
      ) const
  {
    Tensor<1,dim> gradN;
    Tensor<1,dim> gradD;

    double D = 0;
    double N = 0;

    for (const auto& t : TSplineCollection){
      Tensor<1,dim> currGrad = t.gradient(p_d);
      double currVal = t.value(p_d);

      const Point<spacedim+1>& currCP = t.get_cp();

      gradN += currGrad*currCP[component];
      gradD += currGrad*currCP[spacedim];

      N += currVal*currCP[component];
      D += currVal*currCP[spacedim];
    }

    return (-N*gradD + D*gradN)/(D * D);
  }

  template<unsigned int dim, unsigned int spacedim>
  void IsoparametricFunction<dim, spacedim>::vector_gradient(
      const Point<dim>& p_d,
      std::vector< Tensor<1,dim> >& gradients
  ) const {
    std::vector< Tensor<1,dim> >  gradN(spacedim);
    Tensor<1, dim>                gradD;
    Vector<double>                N(spacedim);
    double                        D = 0;

    for (const auto& t : TSplineCollection){
      Tensor<1, dim> t_grad = t.gradient(p_d);
      double         t_val  = t.value(p_d);

      const Point<spacedim+1>& t_cp = t.get_cp();

      D     += t_val  * t_cp[spacedim];
      gradD += t_grad * t_cp[spacedim];
      for (unsigned int d = 0; d < spacedim; d++) {
        gradN[d]  += t_grad * t_cp[d];
        N(d)      += t_val  * t_cp[d];
      }
    }

    for(unsigned int d = 0; d < spacedim; d++)
      gradients[d] = ((-1.) * N(d) * gradD + D * gradN[d])/(D * D);
  }

  template<unsigned int dim, unsigned int spacedim>
  void IsoparametricFunction<dim, spacedim>::gradient_list(
      const std::vector< Point<dim> >&    points,
      std::vector< Tensor<1,dim> >&       gradients,
      const unsigned int component
      ) const
  {
    unsigned int pos = 0;
    for (const auto& p : points){
      gradients[pos++] = gradient(p, component);
    }
  }


  template<unsigned int dim, unsigned int spacedim>
  void IsoparametricFunction<dim, spacedim>::vector_gradient_list(
      const std::vector< Point<dim> >&                points,
      std::vector< std::vector< Tensor<1,dim> > >&    gradients
      ) const
  {
    unsigned int pos = 0;
    for (const auto& p : points){
      std::vector< Tensor<1,dim> > grad_p;
      vector_gradient(p, grad_p);
      gradients[pos++] = grad_p;
    }
  }


  template<unsigned int dim, unsigned int spacedim>
  void IsoparametricFunction<dim, spacedim>::jacobian(
      const Point<dim>&     p_d,
      FullMatrix<double>&   J
      ) const
  {
    J = FullMatrix<double>(spacedim, dim);
    std::vector< Tensor<1,dim> > gradients(spacedim);
    vector_gradient(p_d, gradients);

    for (unsigned int i = 0; i < spacedim; i++){
      for (unsigned int j = 0; j < dim; j++)
        J(i,j) = gradients[i][j];
    }
  }


  template<unsigned int dim, unsigned int spacedim>
  void IsoparametricFunction<dim, spacedim>::jacobian_vector(
      const std::vector< Point<dim> >&     p_d,
      std::vector< FullMatrix<double> >&   J
      ) const
  {
    unsigned int pos = 0;
    for (const auto& p : p_d){
      J[pos] = FullMatrix<double>(spacedim, dim);
      std::vector< Tensor<1, dim> > gradients(spacedim);

      vector_gradient(p, gradients);

      for (unsigned int i = 0; i < spacedim; i++){
        for (unsigned int j = 0; j < dim; j++)
          J[pos](i,j) = gradients[i][j];
      }
      pos++;
    }
  }

  template<unsigned int dim, unsigned int spacedim>
  SymmetricTensor<2, dim> IsoparametricFunction<dim, spacedim>::hessian(
    const Point<dim>&     p,
    const unsigned int    component
  ) const { 

    double D = 0;
    double N = 0;
    Tensor<1, dim>          gradN;
    Tensor<1, dim>          gradD;
    SymmetricTensor<2, dim> hessN;
    SymmetricTensor<2, dim> hessD;


    // build the necessary values:
    for (const auto& t : TSplineCollection){
      const double                  &currVal  = t.value(p);
      const Tensor<1, dim>          &currGrad = t.gradient(p);
      const SymmetricTensor<2, dim> &currHess = t.hessian(p);

      const Point<spacedim+1>& currCP = t.get_cp();

      N += currVal*currCP[component];
      D += currVal*currCP[spacedim];

      gradN += currGrad*currCP[component];
      gradD += currGrad*currCP[spacedim];

      hessN += currHess * currCP[component];
      hessD += currHess * currCP[spacedim];
    }
    const double D2 = D * D; 
    const double D4 = D2 * D2;

    SymmetricTensor<2, dim> out; 
    for (unsigned int i = 0; i < dim; i++)
      for (unsigned int j = i; j < dim; j++)
        out[i][j] = ( (hessN[j][i] * D 
                       + gradN[i] * gradD[j] 
                       - gradN[j] * gradD[i] 
                       - N * hessD[j][i] ) * D2 
                      - 2. * (gradN[i] * D - N * gradD[i]) *
                        gradD[j] * D ) / D4;
    return out;
  }

  template<unsigned int dim, unsigned int spacedim>
  void IsoparametricFunction<dim, spacedim>::print(
      const std::vector< Point<dim> >& points,
      const std::string add_name
      ) const 
  {
    unsigned int nPoints = points.size();
    std::vector< Vector<double> > values(nPoints, Vector<double>(dim));

    vector_value_list(points, values);

    FullMatrix<double> B(nPoints, spacedim);

    for (unsigned int i = 0; i < nPoints; i++){
      for (unsigned int j = 0; j < spacedim; j++){
        B(i, j) = values[i](j);
      }
    }

    std::filebuf f;
    std::string name = add_name + "_IPF";
    name += std::to_string(spacedim);// Dimension of image
    name += "d.dat";

    f.open(name.c_str(), std::ios::out);
    std::ostream out(&f);

    std::cout << "Values of IPFmap evaluation printed to " << name << "\n\n";

    B.print_formatted(out, 15, true, 1, "0");
    f.close();
  }


  template<unsigned int dim, unsigned int spacedim>
  void IsoparametricFunction<dim, spacedim>::print_deriv(
      const std::vector< Point<dim> >& points,
      const std::string add_name) const 
  {
    unsigned int nPoints = points.size();

    std::vector< FullMatrix<double> > J(nPoints);

    jacobian_vector(points, J);


    for (unsigned int i = 0; i < nPoints; i++){
      std::filebuf f;
      std::string name = "out/IPF";
      name += "_" + std::to_string(spacedim) + "d";
      name += "_deriv_";
      name += std::to_string(i);
      name += add_name + ".dat";


      f.open(name.c_str(), std::ios::out);
      std::ostream out(&f);

      J[i].print_formatted(out, 5, true, 1, "0");
      f.close();
    }
    std::ofstream gen_matlab;
    gen_matlab.open("out/loadIPF_deriv.m");
    gen_matlab << "spacedim = " << std::to_string(spacedim) << ";\n";
    gen_matlab << "dim      = " << std::to_string(dim) << ";\n";
    gen_matlab << "nPoints  = " << std::to_string(nPoints) << ";\n";
    gen_matlab << "J        = zeros(spacedim, dim, nPoints);\n\n";

    gen_matlab << "for i = 1:nPoints-1\n";
    gen_matlab << "  J(:,:,i) = load(['IPF_" << std::to_string(spacedim) << "d_deriv_' " << "num2str(i)" << " '.dat']);\n";
    gen_matlab << "  IPF(1,i) = J(1,1,i); \n";
    gen_matlab << "  IPF(2,i) = J(2,1,i);\n";
    gen_matlab << "end\n";
    gen_matlab << "plot(IPF(1,:), IPF(2,:))\n";

  }

  template<int dim, int spacedim>
  IsoparametricManifold<dim, spacedim>::IsoparametricManifold(
    const std::vector< ts_ptr >& TSplineCollection
  ) : ChartManifold<dim, spacedim>() {
    this -> geometry = IsoparametricFunction<dim, spacedim>(TSplineCollection);
  }

  template<int dim, int spacedim>
  IsoparametricManifold<dim, spacedim>::IsoparametricManifold(
    const IsoparametricFunction<dim, spacedim>& other_geometry
  ) : ChartManifold<dim, spacedim>() {
    this -> geometry = other_geometry;
  }

  template<int dim, int spacedim>
  Point<spacedim> IsoparametricManifold<dim, spacedim>::push_forward(
    const Point< dim >& p
  ) const {
    return geometry.point_value( p );
  }

  template<int dim, int spacedim>
  DerivativeForm<1, dim, spacedim> IsoparametricManifold<dim, spacedim>::
    push_forward_gradient(
    const Point< dim >& p
  ) const {
    DerivativeForm<1, dim, spacedim> out; 
    std::vector< Tensor<1, dim> > gradients(spacedim); 
    geometry.vector_gradient(p, gradients); 
    for (int d = 0; d < spacedim; d++)
      out[d] = gradients[d];

    return out;
  }

  template<int dim, int spacedim>
  Point<dim> IsoparametricManifold<dim, spacedim>::pull_back(
    const Point< spacedim >& p
  ) const {
    Point< dim > out;
    for (unsigned int d = 0; d < dim; d++)
      out(d) = p(d); // choose a starting point close to p

    long double norm = (geometry.point_value(out) + (-1.) * p).norm(); 
    unsigned int  k = 0;
    while (norm > 1e-10 && k < 100){
      const Tensor<1, spacedim>& F = this -> push_forward(out);
      const DerivativeForm<1, spacedim, dim>& DFI = 
              (this -> push_forward_gradient(out)).covariant_form().transpose(); 
      out = apply_transformation(DFI, F);
      norm = (geometry.point_value(out) + (-1.) * p).norm(); 
    }
    return out;
  }

  template<int dim, int spacedim>
  std::unique_ptr<Manifold<dim, spacedim>> 
    IsoparametricManifold<dim, spacedim>::clone(
  ) const {
    return std::make_unique<IsoparametricManifold<dim, spacedim>>(*this);
  }

} // namespace dealt
