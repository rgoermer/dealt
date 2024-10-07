#include <poisson_3d.h>
#include <residual_error_estimators.h>

namespace Poisson_Singularity_3D {
  Problem Singularity_RHS::prob = Problem::singularity; 
  Problem Singularity_SOL::prob = Problem::singularity; 

  // Declare the generator for the isogeometric mapping
  SingularityDataGenerator::SingularityDataGenerator(
  ) { 
    std::vector< Point< 3 > > cps = {
        Point<3>(-1., -1., -1.),
        Point<3>(+0., -1., -1.),
        Point<3>(+1., -1., -1.),
        Point<3>(-1., +0., -1.),
        Point<3>(+0., +0., -0.), // pinch a side
        Point<3>(+1., +0., -1.),
        Point<3>(-1., +1., -1.),
        Point<3>(+0., +1., -1.),
        Point<3>(+1., +1., -1.),
        Point<3>(-1., -1., +0.),
        Point<3>(+0., -1., +0.),
        Point<3>(+1., -1., +0.),
        Point<3>(-1., +0., +0.),
        Point<3>(+0., +0., +0.),
        Point<3>(+1., +0., +0.),
        Point<3>(-1., +1., +0.),
        Point<3>(+0., +1., +0.),
        Point<3>(+1., +1., +0.),
        Point<3>(-1., -1., +1.),
        Point<3>(+0., -1., +1.),
        Point<3>(+1., -1., +1.),
        Point<3>(-1., +0., +1.),
        Point<3>(+0., +0., +0.), // pinch a side
        Point<3>(+1., +0., +1.),
        Point<3>(-1., +1., +1.),
        Point<3>(+0., +1., +1.),
        Point<3>(+1., +1., +1.)
      };
    std::vector< double > w(cps.size(), +1.); 

    std::vector< std::vector< double > > kv = {
      {0, 0, 0, 1, 1, 1},
      {0, 0, 0, 1, 1, 1},
      {0, 0, 0, 1, 1, 1}
    };
    std::vector< unsigned int > deg = {2, 2, 2};
    data = IPF_Data<3, 3>(cps, w, kv, deg);
  } // SingularityDataGenerator

  double IPF_mapping::value(
    const Point<3>&     p,
    const unsigned int  component
  ) const { 
    double out;
    const double &x = p(0); 
    const double &y = p(1); 
    const double &z = p(2);
    const double xsq = p(0) * p(0);
    if (component == 0){
      out = -1. + 2. * x;
    } else if (component == 1) {
      out = -1. + 2. * y;
    } else {
      out = -1. + 2.*z 
                +   x * y * (+4. - 8. * z + y * (-4. +  8. * z))
                + xsq * y * (-4. + 8. * z + y * (+4. -  8. * z));
    }
    return out;
  }

  void IPF_mapping::vector_value(
    const Point<3>&       p, 
    Vector<double>&       values
  ) const {
    const double &x = p(0); 
    const double &y = p(1); 
    const double &z = p(2);
    const double xsq = p(0) * p(0);
    values(0) = -1. + 2. * x;
    values(1) = -1. + 2. * y;
    values(2) = -1. + 2.*z 
                +   x * y * (+4. - 8. * z + y * (-4. +  8. * z))
                + xsq * y * (-4. + 8. * z + y * (+4. -  8. * z));
  }

  Point<3> IPF_mapping::point_value(
    const Point<3>&       p
  ) const {
    Point<3> out; 
    const double &x = p(0); 
    const double &y = p(1); 
    const double &z = p(2);
    const double xsq = p(0) * p(0);
    out(0) = -1. + 2. * x;
    out(1) = -1. + 2. * y;
    out(2) = -1. + 2.*z 
                +   x * y * (+4. - 8. * z + y * (-4. +  8. * z))
                + xsq * y * (-4. + 8. * z + y * (+4. -  8. * z));
    return out;
  }

  Tensor<1, 3> IPF_mapping::tensor_value(
    const Tensor<1,3>&       p
  ) const {
    const double &x = p[0]; 
    const double &y = p[1]; 
    const double &z = p[2];
    const double xsq = p[0] * p[0];

    Tensor<1 ,3> values;
    values[0] = -1. + 2. * x;
    values[1] = -1. + 2. * y;
    values[2] = -1. + 2.*z 
                +   x * y * (+4. - 8. * z + y * (-4. +  8. * z))
                + xsq * y * (-4. + 8. * z + y * (+4. -  8. * z));
    return values;
  }

  Tensor<1, 3> IPF_mapping::gradient(
    const Point<3>&     p,
    const unsigned int  component
  ) const { 
    Tensor<1, 3> out; 
    const double &x = p(0); 
    const double &y = p(1); 
    const double &z = p(2);
    const double xsq = p(0) * p(0);
    const double ysq = p(1) * p(1); 

    if (component == 0) {
      out[0] = 2.;
      out[1] = 0.;
      out[2] = 0.;
    } else if (component == 1) {
      out[0] = 0.;
      out[1] = 2.;
      out[2] = 0.;
    } else {
      out[0] = 0. +   y * ( 4. - 8. * z + x * (-8. + 16. * z)) 
                  + ysq * (-4. + 8. * z + x * (+8. - 16. * z));
      out[1] = 0. +   x * ( 4. - 8. * z + y * (-8. + 16. * z)) 
                  + xsq * (-4. + 8. * z + y * (+8. - 16. * z));
      out[2] = 2. +   x * y * (-8. + 8. * y) 
                  + xsq * y * (+8. - 8. * y);
    }
    return out; 
  }

  void IPF_mapping::vector_gradient(
    const Point<3>&       p, 
    std::vector< Tensor<1, 3> >& gradients
  ) const {
    const double &x = p(0); 
    const double &y = p(1); 
    const double &z = p(2);
    const double xsq = p(0) * p(0);
    const double ysq = p(1) * p(1); 

    gradients[0][0] = 2.;
    gradients[0][1] = 0.;
    gradients[0][2] = 0.;
    
    gradients[1][0] = 0.;
    gradients[1][1] = 2.;
    gradients[1][2] = 0.;

    gradients[2][0] = 0. +   y * ( 4. - 8. * z + x * (-8. + 16. * z)) 
                         + ysq * (-4. + 8. * z + x * (+8. - 16. * z));
    gradients[2][1] = 0. +   x * ( 4. - 8. * z + y * (-8. + 16. * z)) 
                         + xsq * (-4. + 8. * z + y * (+8. - 16. * z));
    gradients[2][2] = 2. +   x * y * (-8. + 8. * y) 
                         + xsq * y * (+8. - 8. * y);

  }

  Tensor<2, 3> IPF_mapping::tensor_gradient(
    const Tensor<1, 3>&       p
  ) const {
    const double &x = p[0]; 
    const double &y = p[1]; 
    const double &z = p[2];
    const double xsq = p[0] * p[0];
    const double ysq = p[1] * p[1]; 

    Tensor<2, 3> gradients;

    gradients[0][0] = 2.;
    gradients[0][1] = 0.;
    gradients[0][2] = 0.;
    
    gradients[1][0] = 0.;
    gradients[1][1] = 2.;
    gradients[1][2] = 0.;

    gradients[2][0] = 0. +   y * ( 4. - 8. * z + x * (-8. + 16. * z)) 
                         + ysq * (-4. + 8. * z + x * (+8. - 16. * z));
    gradients[2][1] = 0. +   x * ( 4. - 8. * z + y * (-8. + 16. * z)) 
                         + xsq * (-4. + 8. * z + y * (+8. - 16. * z));
    gradients[2][2] = 2. +   x * y * (-8. + 8. * y) 
                         + xsq * y * (+8. - 8. * y);

    return gradients;
  }

  Tensor<1, 3> IPF_mapping::normal(
    const Point<3>&     p,
    const unsigned int  f
  ) const { 
    const Tensor<2, 3>& J = this -> tensor_gradient(p);

    Tensor<1, 3>  out;

    const int s = f % 2 == 0 ? -1  : +1; 
    unsigned int i, j; 
    switch (f / 2){
      case 0 : 
        i = 1; j = 2;
        break;
      case 1 : 
        i = 0; j = 2; 
        break; 
      case 2 :
        i = 0; j = 1; 
        break;
      default:
        i = -1; j = -1;
        break;
    } // switch


    out[0] = J[1][i] * J[2][j] - J[2][i] * J[1][j];
    out[1] = J[2][i] * J[0][j] - J[0][i] * J[2][j];
    out[2] = J[0][i] * J[1][j] - J[1][i] * J[0][j];

    if (out[f / 2] / std::fabs(out[f / 2]) != s ) {
      out *= -1.;

    }

    return out / out.norm();
  } 

  SymmetricTensor<2, 3> IPF_mapping::hessian(
    const Point<3>&     p,
    const unsigned int  component
  ) const { 
    SymmetricTensor<2, 3> out; 
    const double &x = p(0); 
    const double &y = p(1); 
    const double &z = p(2);
    if (component == 2) {
      out[0][0] = +0. + 2. * y * (-4. + y * (+4. - 8. * z) + 8. * z);
      out[0][1] = +4. - 8. * z + y * (-8. + 16. * z) + x * (-8. + y * (+16. -32. * z) + 16. * z);
      out[0][2] = +0. + y * (-8. + x * (+16. - 16. * y) +  8. * y);
      out[1][1] = +0. + 2. * x * (-4. + x * (4. - 8. * z) + 8. * z);
      out[1][2] = +0. + x * (-8. + x * (+ 8. - 16. * y) + 16. * y);
      out[2][2] = +0.;
    }
    return out; 
  }

  double IPF_mapping_simplified::value(
    const Point<3>&     p,
    const unsigned int  component
  ) const { 
    double out;
    const double &x = p(0); 
    const double &y = p(1); 
    const double &z = p(2);
    if (component == 0){
      out = -1. + 2. * x;
    } else if (component == 1) {
      out = -1. + 2. * y;
    } else {
      out = -1. + 2. * z;
    }
    return out;
  }

  void IPF_mapping_simplified::vector_value(
    const Point<3>&       p, 
    Vector<double>&       values
  ) const {
    const double &x = p(0); 
    const double &y = p(1); 
    const double &z = p(2);
    values(0) = -1. + 2. * x;
    values(1) = -1. + 2. * y;
    values(2) = -1. + 2. * z;
  }

  Point<3> IPF_mapping_simplified::point_value(
    const Point<3>&       p
  ) const {
    Point<3> out; 
    const double &x = p(0); 
    const double &y = p(1); 
    const double &z = p(2);
    out(0) = -1. + 2. * x;
    out(1) = -1. + 2. * y;
    out(2) = -1. + 2. * z;
    return out;
  }

  Tensor<1, 3> IPF_mapping_simplified::tensor_value(
    const Tensor<1,3>&       p
  ) const {
    const double &x = p[0]; 
    const double &y = p[1]; 
    const double &z = p[2];

    Tensor<1 ,3> values;
    values[0] = -1. + 2. * x;
    values[1] = -1. + 2. * y;
    values[2] = -1. + 2. * z;
    return values;
  }

  Tensor<1, 3> IPF_mapping_simplified::gradient(
    const Point<3>&     /* p */,
    const unsigned int  component
  ) const { 
    Tensor<1, 3> out; 
    if (component == 0) {
      out[0] = 2.;
      out[1] = 0.;
      out[2] = 0.;
    } else if (component == 1) {
      out[0] = 0.;
      out[1] = 2.;
      out[2] = 0.;
    } else {
      out[0] = 0.;
      out[1] = 0.;
      out[2] = 2.;
    }
    return out; 
  }

  void IPF_mapping_simplified::vector_gradient(
    const Point<3>&              /* p */, 
    std::vector< Tensor<1, 3> >& gradients
  ) const {

    gradients[0][0] = 2.;
    gradients[0][1] = 0.;
    gradients[0][2] = 0.;
    
    gradients[1][0] = 0.;
    gradients[1][1] = 2.;
    gradients[1][2] = 0.;

    gradients[2][0] = 0.; 
    gradients[2][1] = 0.; 
    gradients[2][2] = 2.; 
                         

  }

  Tensor<2, 3> IPF_mapping_simplified::tensor_gradient(
    const Tensor<1, 3>&       /* p */
  ) const {

    Tensor<2, 3> gradients;

    gradients[0][0] = 2.;
    gradients[0][1] = 0.;
    gradients[0][2] = 0.;
    
    gradients[1][0] = 0.;
    gradients[1][1] = 2.;
    gradients[1][2] = 0.;

    gradients[2][0] = 0.; 
    gradients[2][1] = 0.; 
    gradients[2][2] = 2.; 

    return gradients;
  }

  SymmetricTensor<2, 3> IPF_mapping_simplified::hessian(
    const Point<3>&     /* p */,
    const unsigned int  /* component */
  ) const { 
    SymmetricTensor<2, 3> out; 
    return out; 
  }

  Tensor<1, 3> IPF_mapping_simplified::normal(
    const Point<3>&     p,
    const unsigned int  f
  ) const { 
    const Tensor<2, 3>& J = this -> tensor_gradient(p);

    Tensor<1, 3> g0, g1, out;

    const int s = f % 2 == 0 ? -1  : +1; 
    unsigned int i, j; 
    switch (f / 2){
      case 0 : 
        i = 1; j = 2;
        break;
      case 1 : 
        i = 0; j = 2; 
        break; 
      case 2 :
        i = 0; j = 1; 
        break;
      default:
        i = -1; j = -1;
        break;
    } // switch


    g0[0] = J[0][i];
    g0[1] = J[1][i];
    g0[2] = J[2][i];

    g1[0] = J[0][j];
    g1[1] = J[1][j];
    g1[2] = J[2][j];

    out[0] = g0[1] * g1[2] - g0[2] * g1[1];
    out[1] = g0[2] * g1[0] - g0[0] * g1[2];
    out[2] = g0[0] * g1[1] - g0[1] * g1[0];

    if (out[f / 2] / std::fabs(out[f / 2]) != s )
      out *= s;

    return out / out.norm();
  }

  void IPF_inverse::vector_value(
    const Point<3>&     p, 
    Vector< double >&   values
  ) const {
    const IPF_mapping phi; 

    // Use Newton-Iteration to calculate the value: 
    Tensor<1, 3> xn, rhs;
    xn[0] = p(0);
    xn[1] = p(1);
    xn[2] = p(2);
    rhs = xn;

    DerivativeForm<1, 3, 3> Dphi(phi.tensor_gradient(xn));
    DerivativeForm<1, 3, 3> DphiInv = Dphi.covariant_form().transpose();
    
    Tensor<1, 3> xn1 = xn + apply_transformation(DphiInv, phi.tensor_value(xn) - rhs);

    const unsigned int kmax = 100; 
    unsigned int k = 0;
    while (k < kmax && (xn - xn1).norm() > 1e-15){
      xn = xn1; 
      Dphi             = phi.tensor_gradient(xn);
      DphiInv          = Dphi.covariant_form().transpose();
      xn1              = xn + apply_transformation(DphiInv, phi.tensor_value(xn) - rhs);
      k++;
    } // while 
    if (k == kmax){
      std::cout << "Maximum iteration steps reached for computation of inverse value. Terminated with an error of "
                << std::scientific << (xn - xn1).norm()
                << std::endl;
    } // if 

    values(0) = xn1[0];
    values(1) = xn1[1];
    values(2) = xn1[2];
  }

  void IPF_inverse::vector_hessian(
    const Point<3>&         p,
    std::vector< SymmetricTensor<2, 3> >& hessians
  ) const {
    // use these values to approximate second derivative
    SymmetricTensor<2, 3>& H0 = hessians[0];
    SymmetricTensor<2, 3>& H1 = hessians[1];
    SymmetricTensor<2, 3>& H2 = hessians[2];
    const double h = 1e-12; 
    const double hh = 4. * h * h;
    for (unsigned int i = 0; i < 3; i++){
      for (unsigned int j = i; j < 3; j++){
        const Point<3> p0 = p + h * Point<3>::unit_vector(i)
                              + h * Point<3>::unit_vector(j);
        const Point<3> p1 = p + h * Point<3>::unit_vector(i)
                              + (-1. * h) * Point<3>::unit_vector(j);
        const Point<3> p2 = p + (-1. * h) * Point<3>::unit_vector(i)
                              + h * Point<3>::unit_vector(j);
        const Point<3> p3 = p + (-1. * h) * Point<3>::unit_vector(i)
                              + (-1. * h) * Point<3>::unit_vector(j);

        Vector<double> v0(3), v1(3), v2(3), v3(3), V;
        vector_value(p0, v0);
        vector_value(p1, v1);
        vector_value(p2, v2);
        vector_value(p3, v3);
        
        V = v0; 
        V -= v1;
        V -= v2; 
        V += v3;
        H0[i][j] = V(0) / hh; 
        H1[i][j] = V(1) / hh;
        H2[i][j] = V(2) / hh;
      }
    }
   
  }
                  

  double Singularity_RHS::value(
    const Point<3>&     p,
    const unsigned int  /* component */
  ) const { 
    double out = nan("1"); 
    const double xsq = p(0) * p(0);
    const double ysq = p(1) * p(1);
    const double zsq = p(2) * p(2);

    switch (prob) {
      case Problem::simple: {
        out = 2. * p(2) * (2. - ysq - xsq);
        break;
      }
      case Problem::singularity: {
        const double sum = 25. * std::pow(xsq + ysq + zsq, 1.10 /* = 5/4 */);
        out = -0. + 50. * xsq * xsq + 50. * ysq * ysq 
                               -  4. * (1. + 25. * zsq)
                               + ysq * (-76.             + 50. * zsq)
                               + xsq * (-76. + 56. * ysq + 50. * zsq);
        out = -1. * out / sum;
        break;                           
      }
    }
    return out;
  }

  double Singularity_SOL::value(
    const Point<3>&     p,
    const unsigned int  /* component */
  ) const { 
    double out = nan("1");
    switch (prob){
      case Problem::simple : {
        out = (1. - p(0)) * (1. + p(0)) * (1. - p(1)) * (1. + p(1)) * p(2);
        break;
      }
      case Problem::singularity : {
        const double sum = std::pow(p(0) * p(0) + p(1) * p(1) + p(2) * p(2), 0.10);
        out =((1. - p(0)) * (1. + p(0)) * (1. - p(1)) * (1. + p(1))) / sum;
        break;
      }
    }

    return out;
  }

  Tensor<1, 3> Singularity_SOL::gradient(
    const Point<3>&     p,
    const unsigned int  /* component */
  ) const { 
    const double x   = p(0);
    const double y   = p(1);
    const double z   = p(2);
    const double xsq = p(0) * p(0);
    const double ysq = p(1) * p(1);
    const double zsq = p(2) * p(2);
    const double sum = 2. * std::pow(xsq + ysq + zsq, 1.10 /* = 5/4 */);

    Tensor<1, 3> out;

    switch (prob) {
      case Problem::simple : {
        out[0] = -2. * p(0) * (1. - ysq) * p(2);
        out[1] = -2. * p(1) * (1. - xsq) * p(2);
        out[2] = (1. - xsq) * (1. - ysq);
        break;
      }
      case Problem::singularity : {
        out[0] = x * (-1. + ysq) * (1. +  9. * xsq + 10. * ysq + 10. * zsq);
        out[1] = (-1. + xsq) * y * (1. + 10. * xsq +  9. * ysq + 10. * zsq);
        out[2] = -1. * (-1. + xsq) * (-1. + ysq) * z;
        out /= sum;
        break;
      }
      default : {
        out[0] = nan("1");
        out[1] = nan("1");
        out[2] = nan("1");
        break;
      }
    }
    return out;
  }

  double Singularity_NumSol::value(
    const Point<3>&     p, 
    const unsigned int  /* component */
  ) const {
    Assert(base.size() == solution.size(), 
            ExcDimensionMismatch(base.size(), solution.size()));

    const unsigned int N = base.size(); 
    double out = 0;
    for (unsigned int i = 0; i < N; i++)
      out += solution(i) * base[i].value(p);

    return out; 
  }

  void Singularity_NumSol::value_list(
    const std::vector< Point< 3 > >&  points,
    std::vector<double>&              values,
    const unsigned int                /* component = 0 */
  ) const {
    const unsigned int N = base.size(); 
    const unsigned int M = values.size();
    for (unsigned int i = 0; i < N; i++){
      std::vector<double> tmp(M); 
      base[i].value_list(points, tmp); 
      for (unsigned int j = 0; j < M; j++)
        values[j] += solution(i) * tmp[j];
    }
  } // value_list

  Tensor<1, 3> Singularity_NumSol::gradient(
    const Point<3>&     p, 
    const unsigned int  /* component */
  ) const {
    Assert(base.size() == solution.size(), 
            ExcDimensionMismatch(base.size(), solution.size()));

    const unsigned int N = base.size(); 
    Tensor<1, 3> out;
    for (unsigned int i = 0; i < N; i++)
      out += solution(i) * base[i].gradient(p);

    return out; 
  }

  double Singularity_NC::value(
    const Point<3>&     p,
    const unsigned int  /* component */
  ) const { 
    const Singularity_SOL u;
    const Tensor<1, 3> du = u.gradient(p);
    if (p(0) == -1.)
      return -du[0];
    else if (p(0) == +1.)
      return +du[0];
    else if (p(1) == -1.)
      return -du[1];
    else if (p(1) == +1.)
      return +du[1];


    Assert(false, ExcInternalError());
    throw ExcInternalError();
    return nan("1");
  }


  double Singularity_NC_Z0::value(
    const Point<3>&       p,
    const unsigned int    /* component */
  ) const {
    const Singularity_SOL u;
    const Tensor<1, 3> du = u.gradient(phi.point_value(p));
    const Tensor<1, 3> n({
                      p(1) * (-8. + p(0) * (16. - 16. * p(1)) +  8. * p(1)),
                      p(0) * (-8. + p(0) * ( 8. - 16. * p(1)) + 16. * p(1)),
                      4.
                    });
    return du * (-1. * n / n.norm()); 
  }

  Tensor<1, 3> Singularity_NC_Z0::normal(
    const Point<3>&     p
  ) const {
    const Tensor<1, 3> n({
                      p(1) * (-8. + p(0) * (16. - 16. * p(1)) +  8. * p(1)),
                      p(0) * (-8. + p(0) * ( 8. - 16. * p(1)) + 16. * p(1)),
                      4.
                    });
    return -1. * n / n.norm();
  }

  double Singularity_NC_Z1::value(
    const Point<3>&       p,
    const unsigned int    /* component */
  ) const {
    const Singularity_SOL u;
    const Tensor<1, 3> du = u.gradient(phi.point_value(p));
    const Tensor<1, 3> n({
                      p(1) * (8. -  8. * p(1) + p(0) * (-16. + 16. * p(1))),
                      p(0) * (8. - 16. * p(1) + p(0) * ( -8. + 16. * p(1))),
                      4.
                    });
    return du * (n / n.norm()); 
  }

  Poisson_Singularity::Poisson_Singularity(
      unsigned int ref, 
      unsigned int order
  ) : ref(ref),
    order(order), 
    data(SingularityDataGenerator().data),
    tria(data), 
    rhs_fcn(),
    sol_fcn(),
    neumann_bc(),
    neumann_bc_z0(),
    neumann_bc_z1()
  { 
    tria.degree_elevate_global(order); 

    out_init = "out/poisson_singularity_3d/";

    // set boundary indicators
    for (const auto& face : tria.active_face_iterators()){
      const Point<3>& c = face -> center();
      if (face -> at_boundary()){
        if (c(0) == +0)
          face -> set_boundary_id(Boundary::dirichleth_x0);
        else if (c(0) == +1)
          face -> set_boundary_id(Boundary::dirichleth_x1);
        else if (c(1) == +0)
          face -> set_boundary_id(Boundary::dirichleth_y0);
        else if (c(1) == +1)
          face -> set_boundary_id(Boundary::dirichleth_y1);
        else if (c(2) == +0)
          face -> set_boundary_id(Boundary::neumann_z0);
        else if (c(2) == +1)
          face -> set_boundary_id(Boundary::neumann_z1);
        else 
          face -> set_boundary_id(Boundary::none);
      } // if ( at_boundary )
    } // for ( face )

#ifdef DEBUG
    for (const auto& face : tria.active_face_iterators())
      Assert(face -> boundary_id() != Boundary::none, ExcInternalError());
#endif

    offset = 6; // ensures that [0, 0, 0] will not be evaluated
    tria.refine_global(offset);

    // Find boundary dofs
    tria.set_boundary_dofs();

    // switch mesh from parametric to bezier and compute
    // extraction operators
    tria.refine_bezier_elements();
    tria.compute_extraction_operators();


    // Initialize quantities of interest
    l2_error.reinit(ref+1);
    h1_error.reinit(ref+1);

    mesh_size.reinit(ref+1);
    dofs_per_level.reinit(ref+1);
  } // Constructor

  void Poisson_Singularity::print_grid(
      const bool bezier
  ) const { 
    std::string svg_name, wire_name;
    if (!bezier){
      svg_name = out_init + "00svg/o" 
              + std::to_string(data.max_degree() + order)
              + "/grid_l" 
              + std::to_string(tria.n_levels() - offset - 1);
      wire_name = out_init + "00svg/o"
              + std::to_string(data.max_degree() + order)
              + "/grid_l" 
              + std::to_string(tria.n_levels() - offset - 1);
    } else {
      svg_name = out_init + "00svg/o" 
              + std::to_string(data.max_degree() + order)
              + "/poisson_singularity_bezier_grid_l" 
              + std::to_string(tria.n_levels() - offset - 1);
      wire_name = out_init 
              + "poisson_singularity_bezier_grid_l" 
              + std::to_string(tria.n_levels() - offset - 1);
    }

    GridOutFlags::Eps<3> eps_flags;
    eps_flags.azimut_angle = 25;
    eps_flags.turn_angle = 45;
  
    std::ofstream out(svg_name + ".eps");
    GridOut       grid_out;
    grid_out.set_flags(eps_flags);
  
    grid_out.write_eps(tria, out);
    
    // tria.print_IPF_surface_wireframe(
    //         wire_name, 
    //         neumann_z1,
    //         15,
    //         5);

    // tria.print_IPF_surface_wireframe(
    //         wire_name, 
    //         dirichleth_x0,
    //         15,
    //         5);

    // tria.print_IPF_surface_wireframe(
    //         wire_name, 
    //         dirichleth_y0,
    //         15,
    //         5);
    
    std::map<unsigned int, Point<1>> solution_vertex_values; 
    const auto& splines          = tria.get_splines();
    const unsigned int n_splines = splines.size();

    for (auto it  = tria.begin_active_vertex(); 
              it != tria.end_vertex();
              it++){
      for (unsigned int i = 0; i < n_splines; i++){
        solution_vertex_values[it->index()](0)
                += splines[i] -> value(it -> vertex()) 
                    * solution(i);  
      }
    }

    print_surface_z0();

    tria.generate_mesh_file<0>(wire_name, true, 8);
    tria.generate_mesh_file<1>(wire_name, false, 8, solution_vertex_values);

    std::string name1 = "out/poisson_singularity_3d/poisson_o" + std::to_string(order + data.max_degree());
  
    std::string l2_err = name1 + "_l2e.dat";
    std::string h1_err = name1 + "_h1e.dat";
    std::string mesh   = name1 + "_h.dat";
    std::string dofs   = name1 + "_dofs.dat";
  
    std::filebuf l2e, h1e, h, d;
  
    l2e.open(l2_err.c_str(), std::ios::out);
    h1e.open(h1_err.c_str(), std::ios::out);
    h.open(mesh.c_str()    , std::ios::out);
    d.open(dofs.c_str()    , std::ios::out);
  
    std::ostream l2e_out(&l2e);
    std::ostream h1e_out(&h1e);
    std::ostream mesh_out(&h);
    std::ostream dofs_out(&d);
  
    l2_error.print(l2e_out, 16);
    h1_error.print(h1e_out, 16);
    mesh_size.print(mesh_out, 16);
    dofs_per_level.print(dofs_out, 16);
  
    l2e.close();
    h1e.close();
    h.close();
    d.close();
  } // print_grid

  void Poisson_Singularity::setup_system() {
    system_rhs.reinit(tria.n_active_splines());
    solution.reinit(tria.n_active_splines());
  
    const auto& IEN_array = tria.get_IEN_array();
    Assert(IEN_array.size() > 0, ExcInternalError());
  
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
    sparsity_pattern.reinit(
        tria.n_active_splines(),
        tria.n_active_splines(),
        tria.n_active_splines() );
  
    for (const auto& [_, arr] : IEN_array)
      for (unsigned int i : arr)
        for (unsigned int j : arr)
          sparsity_pattern.add(i, j);
  
    // Free superfluous space
    sparsity_pattern.compress();
    system_matrix.reinit(sparsity_pattern);
  
  } // setup_system

  void Poisson_Singularity::assemble_system(){
    std::cout << "Assembling system matrix ... " << std::endl;

    // Setup initial tables that store the bernstein values / grads / hessians.
    std::vector< unsigned int > degrees = tria.get_degree();
    degrees[0] = (degrees[0] * degrees[0])/2 + 1;
    degrees[1] = (degrees[1] * degrees[1])/2 + 1;
    degrees[2] = (degrees[2] * degrees[2])/2 + 1;

    TSValues<3> ts_values(
        &tria,
        degrees,
        update_values |
        update_gradients |
        update_quadrature_points |
        update_hessians |
        update_JxW_values);
    TSFaceValues<3> face_values(
        &tria,
        degrees,
        update_values |
        update_gradients |
        update_quadrature_points |
        update_normal_vectors |
        update_JxW_values);

    //print_reference_values(0.5);
  
  
    const unsigned int dofs_per_cell = ts_values.n_dofs_per_cell();
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);
  
    for (const auto& cell : tria.active_cell_iterators()){
      // Get the Bernstein values on the cell
      ts_values.reinit(cell);
  
      // Reset the cell matrix
      cell_matrix       = 0;
      cell_rhs          = 0;
  
      // Quadrature sum:
      for (const unsigned int q : ts_values.quadrature_point_indices()){
  
        // Build the cell matrix:
        double dx = ts_values.JxW(q);
        for (const unsigned int i : ts_values.dof_indices())
          for(const unsigned int j : ts_values.dof_indices())
            cell_matrix(i,j) += ( ts_values.shape_grad(i, q)
                                  * ts_values.shape_grad(j, q)
                                  * dx );
  
        // Map the quadrature point from real cell to parametric cell
        const Point<3>& mapped_q = ts_values.quadrature_point(q);
        const double rhs = rhs_fcn.value(mapped_q);
        AssertIsFinite(rhs);
        for (const unsigned int i : ts_values.dof_indices())
          cell_rhs(i) += ( ts_values.shape_value(i, q)
                            * rhs
                            * dx );
  
  
      } // for ( q )
  
      // Check for neumann conditions
      if (cell -> at_boundary()){
        for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; f++){
          if (cell -> face(f) -> at_boundary()
                && cell -> face(f) -> boundary_id() == Boundary::neumann_z0
              ){
            face_values.reinit(cell, f);
            for (const unsigned int q : face_values.quadrature_point_indices()){
              const long double g = sol_fcn.gradient(face_values.quadrature_point(q)) *
                                      face_values.normal_vector(q);
              for (const unsigned int i : face_values.dof_indices())
                cell_rhs(i) += g * face_values.shape_value(i, q)
                                 * face_values.JxW(q);
            }
          } else if (cell -> face(f) -> at_boundary()
                    && cell -> face(f) -> boundary_id() == Boundary::neumann_z1
              ){
            face_values.reinit(cell, f);
            for (const unsigned int q : face_values.quadrature_point_indices()){
              const long double g = sol_fcn.gradient(face_values.quadrature_point(q)) *
                                      face_values.normal_vector(q);
              for (const unsigned int i : face_values.dof_indices())
                cell_rhs(i) += g * face_values.shape_value(i, q)
                                 * face_values.JxW(q);
            }
          } // if ( boundary_id )
        } // for ( f )
      } // if
  
      // Add the values to the system
      std::vector< unsigned int > local_dof_indices =
              tria.get_IEN_array(cell);
  
      system_matrix.add(local_dof_indices, cell_matrix);
      system_rhs.add(local_dof_indices, cell_rhs);
  
    } // for ( cell )
  
  
    // Get a map that corresponds to boundary values
    const auto& boundary_dofs =
            tria.get_boundary_dofs();
  
    // Use the boundary dofs to eliminate rows and columns with 
    // dirichleth boundary conditions. In this example we impose
    // homogeneous boundary conditions.
    unsigned int n_global_dofs = tria.n_active_splines();
    std::vector< unsigned int > dirichleth = 
      {Boundary::dirichleth_x0, Boundary::dirichleth_x1,
       Boundary::dirichleth_y0, Boundary::dirichleth_y1};
    for (const unsigned int dc : dirichleth){
      for (const auto& dof : boundary_dofs.at(dc)){
        for (unsigned int j = 0; j < n_global_dofs; j++){
          system_matrix.set(dof, j, 0);
          system_matrix.set(j, dof, 0);
        }
        system_rhs(dof) = 0;
        system_matrix.set(dof, dof, 1);
      }
    }
  
    std::cout << " ... done!" << std::endl;
  } // assemble system


  void Poisson_Singularity::output_system()
  {
    std::cout << "Printing system matrix and rhs ... " << std::endl;
  
    const unsigned int level = tria.n_levels() - offset - 1;
    const unsigned int p = data.max_degree() + order;
    const std::string l = std::to_string(level);
    std::string name = out_init + "/singularity_o" + std::to_string(p) + "_l" + l; 
  
    std::string matrix = name + "_mat.dat" ;
    std::string vector = name + "_vec.dat" ;
    std::string soluti = name + "_sol.dat" ;
  
    std::filebuf mat, vec, sol;
    mat.open(matrix.c_str(), std::ios::out);
    vec.open(vector.c_str(), std::ios::out);
    sol.open(soluti.c_str(), std::ios::out);
  
    std::ostream mat_out(&mat);
    std::ostream vec_out(&vec);
    std::ostream sol_out(&sol);
  
    system_matrix.print_formatted(mat_out, 16, true, 1, "0");
    system_rhs.print(vec_out, 16);
    solution.print(sol_out, 16);
  
    mat.close();
    vec.close();
    sol.close();
  
    // tria.printIPF(name, 16, true, true);
  
    // Print the IPF wireframe:
    // tria.print_grid(name);
    tria.print_IPF_wireframe(name);
    tria.coarsen_bezier_elements();
    tria.print_IPF_wireframe(name);
    tria.refine_bezier_elements();
    
    std::cout << " ... output to files done!" << std::endl;
  } // output_system

  void Poisson_Singularity::solve(){
    std::cout << "Solving system ... " << std::endl;
  
    SolverControl            solver_control(750 * tria.n_active_splines(), 1e-12 * system_rhs.l2_norm());
    SolverCG<Vector<double>> solver(solver_control);
    solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());

    // print_pointwise_error();
  
    std::cout << " ... done!" << std::endl;
  }

  void Poisson_Singularity::compute_h1_error()
  {
    std::cout << "Computing H1-error at Level "
              << tria.n_levels() - offset - 1
              << " with degree "
              << data.max_degree() + order
              << " ... " << std::endl;
  
  
    const unsigned int r = tria.n_levels() - offset - 1;
    double h1 = 0;
    double l2 = 0;

    // Setup initial tables that store the bernstein values / grads / hessians.
    std::vector< unsigned int > degrees = tria.get_degree();
    degrees[0] = (degrees[0] * degrees[0]) + 1;
    degrees[1] = (degrees[1] * degrees[1]) + 1;
    degrees[2] = (degrees[2] * degrees[2]) + 1;

    TSValues<3> ts_values(
        &tria,
        degrees,
        update_values |
        update_gradients |
        update_quadrature_points |
        update_hessians |
        update_JxW_values);
  
    dofs_per_level(r) = tria.n_active_splines();
    for (const auto& cell : tria.active_cell_iterators()){
      // Get the Bernstein values on the cell
      ts_values.reinit(cell);
  
      std::vector<unsigned int> local_dof_indices = tria.get_IEN_array(cell);
  
      // Quadrature sum:
      for (const unsigned int q_index : ts_values.quadrature_point_indices()){
        // Map the quadrature point from real cell to parametric cell
        const Point<3> mapped_q = ts_values.quadrature_point(q_index);
  
        // Get the value of the approximation
        double u_diff = sol_fcn.value(mapped_q);
        for (const unsigned int i : ts_values.dof_indices())
          u_diff -= (solution(local_dof_indices[i]) *
                      ts_values.shape_value(i, q_index));
  
        // Build the gradient value at quadrature point:
        Tensor<1, 3> grad_u_diff = sol_fcn.gradient(mapped_q);
        for (const unsigned int i : ts_values.dof_indices())
          grad_u_diff -= ( solution(local_dof_indices[i]) *
                           ts_values.shape_grad(i, q_index));
  
        double dx = ts_values.JxW(q_index);
        h1 += (( grad_u_diff * grad_u_diff ) * dx) ;
        l2 += ((      u_diff * u_diff      ) * dx);
      } // for ( q_index )
    } // for ( cell )
    h1_error(r) = std::sqrt(l2 + h1);
    l2_error(r) = std::sqrt(l2);
  
    std::cout << " ... done!" << std::endl;
  } // compute_h1_error

  void Poisson_Singularity::estimate_and_refine(){
    // declare a container to store cells to be marked
    std::vector< TriaIterator<CellAccessor<3, 3>> > mark;
    const bool uniform = true; 
    std::vector< unsigned int > degrees = tria.get_degree();
    for (unsigned int& p : degrees)
      p = p * p + 1;
  
    if (uniform) {
      for (const auto& cell : tria.active_cell_iterators())
        mark.push_back(cell);

    } else {
      // special case, if there is only one cell
      if (tria.n_active_cells() == 1) {
        mark.push_back(tria.begin_active());
      } else {
        tria.get_IPF();
        Vector<double> local_residuals(tria.n_active_splines());
        std::map< types::boundary_id,
                  const Function<3>* >        
              neumann_data = 
                  {
                    {Boundary::neumann_z0, &neumann_bc_z0},
                    {Boundary::neumann_z1, &neumann_bc_z1}
                  };
        
        ResidualEstimators::Poisson<3>::estimate(
            &tria,
            degrees,
            solution,
            local_residuals,
            &rhs_fcn,
            neumann_data
        );

        tria.refine_fixed_number(local_residuals, 0.10);
        
      } // if ( special case )
    }
  
    // Finally prepare the next step, by coarsening the bezier elements
    tria.coarsen_bezier_elements();

    // And then setting the refinement flags
    tria.set_refine_flags(mark);

    // Perform refinement of grid:
    tria.execute_coarsening_and_refinement();

    // And set boundary dofs:
    tria.set_boundary_dofs();
  
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
  } // estimeate_and_refine

  void Poisson_Singularity::run(){
    unsigned int nlr = 0; 
    unsigned int old_level = 0;
    while (tria.n_levels() - offset - 1 < ref + 1){
      nlr++;
      this -> setup_system();
      this -> assemble_system();
      this -> solve();
      this -> print_solution();
      this -> print_grid();
      this -> output_system();
      this -> compute_h1_error();
      this -> estimate_and_refine();
      if (tria.n_levels() - offset - 1 != old_level) {
        nlr = 0;
      } else if (nlr == 5) {
        // Enforce global refinement
        tria.coarsen_bezier_elements();
        tria.refine_global();
        tria.set_boundary_dofs();
        tria.refine_bezier_elements();
        tria.compute_extraction_operators();
        nlr = 0; 
      }
    }
  }

  void Poisson_Singularity::print_solution(
  ){
    std::cout << "Printing numerical solution to files ... " << std::endl;
    const IPF_mapping                    Phi;
    const auto&                          base_ptr = tria.get_splines();
    std::vector< TSplineFunction<3, 3> > base(base_ptr.size());
    for (unsigned int i = 0; i < base_ptr.size(); i++)
      base[i] = *base_ptr[i];

    Singularity_NumSol num_fcn(solution, base);

    const unsigned int N = 100;
    const std::vector< double > zhat = {0.000, 0.125, 0.250,
                                        0.375, 0.500, 0.625,
                                        0.750, 0.875, 1.000};
    std::vector< double > xhat(N), yhat(N);
    FullMatrix<double> X(N, N), Y(N, N), Z(N, N), U(N, N);
    
    const unsigned int level = tria.n_levels() - offset - 1;
    const unsigned int p = data.max_degree() + order;
    const std::string l = std::to_string(level);

    for (unsigned int i = 0; i < N; i++) {
      xhat[i] = 0 + i / (N - 1.);
      yhat[i] = 0 + i / (N - 1.);
    }

    unsigned int ind = 0;
    for (const double zh : zhat){
      for (unsigned int i = 0; i < N; i++) {
        for (unsigned int j = 0; j < N; j++) {
          const Point<3> p(xhat[j], yhat[i], zh);
          X(i, j) = Phi.value(p, 0);
          Y(i, j) = Phi.value(p, 1);
          Z(i, j) = Phi.value(p, 2);
          U(i, j) = num_fcn.value(Point<3>(xhat[i], yhat[j], zh));
        }
      }

      const std::string name = out_init + "/singularity_o" + std::to_string(p) + "_l" + l; 
  
      const std::string name_x = name + "_X." + std::to_string(ind) + ".dat" ;
      const std::string name_y = name + "_Y." + std::to_string(ind) + ".dat" ;
      const std::string name_z = name + "_Z." + std::to_string(ind) + ".dat" ;
      const std::string name_u = name + "_U." + std::to_string(ind) + ".dat" ;
      ind++;


      std::filebuf fx, fy, fz, fu;
      fx.open(name_x.c_str(), std::ios::out);
      fy.open(name_y.c_str(), std::ios::out);
      fz.open(name_z.c_str(), std::ios::out);
      fu.open(name_u.c_str(), std::ios::out);

      std::ostream out_x(&fx), out_y(&fy), out_z(&fz), out_u(&fu);
      X.print_formatted(out_x, 16, true, 1, "0");
      Y.print_formatted(out_y, 16, true, 1, "0");
      Z.print_formatted(out_z, 16, true, 1, "0");
      U.print_formatted(out_u, 16, true, 1, "0");

      fx.close();
      fy.close();
      fz.close();
      fu.close();
    }

    std::cout << "... done!" << std::endl;

  } // print_solution


  void Poisson_Singularity::print_pointwise_error(
  ) const {
    const std::string log_name = "log/differences.txt"; 
    std::ofstream log(log_name);

    bool success = true;

    const unsigned int Nx = 50; 
    const unsigned int Ny = 25; 
    const unsigned int Nz = 15; 
    const unsigned int N  = Nx * Ny * Nz;

    const auto&                          mapping = tria.get_IPF();
    const auto&                          base_ptr = tria.get_splines();
    std::vector< TSplineFunction<3, 3> > base(base_ptr.size());
    for (unsigned int i = 0; i < base_ptr.size(); i++)
      base[i] = *base_ptr[i];

    Singularity_NumSol num_fcn(solution, base);

    std::vector< Point<3> > evals_parametric(Nx * Ny * Nz);
    std::vector< Point<3> > evals_physical(Nx * Ny * Nz);
    unsigned int ind = 0; 
    for (unsigned int z = 0; z < Nz; z++){
      for (unsigned int y = 0; y < Ny; y++){
        for (unsigned int x = 0; x < Nx; x++){
          const Point<3> p(+0. + x/(Nx - 1.), 
                           +0. + y/(Ny - 1.), 
                           +0. + z/(Nz - 1.));
          evals_parametric[ind] = p;
          evals_physical[ind]   = mapping.point_value(p);
          ind++;
        } // for ( x )
      } // for ( y )
    } // for ( z )


    std::vector< double > num_vals( N );
    std::vector< double > sol_vals( N );
    
    num_fcn.value_list(evals_parametric, num_vals);
    sol_fcn.value_list(evals_physical, sol_vals);

    log << std::scientific;
    for (unsigned int i = 0;i < N; i++){
      log << "p = (" << evals_parametric[i]
          << "); Phi(p) = (" << evals_physical[i] << ")"
          << std::endl;

      log << "u_h(Phi(p)) = " << num_vals[i] << std::endl
          << "u(Phi(p)) = " << sol_vals[i] << std::endl;
      
      log << "error = " << std::fabs(num_vals[i] - sol_vals[i]) << std::endl;

      if (std::fabs(num_vals[i] - sol_vals[i]) > 1e-14)
        success = false;
    } // for ( i )
#ifndef DEBUG
    if (!success)
      throw ExcInternalError();
#endif
    Assert(success, ExcInternalError());
  } // print_pointwise_error

  void Poisson_Singularity::print_surface_z0(
  ) const {

    const unsigned int N = 100; 
    const auto& map = tria.get_IPF(); 
    const std::string& l_str = std::to_string(tria.n_levels() - offset - 1); 
    const std::string& p_str = std::to_string(data.max_degree() + order);
    const std::string& out_name = out_init 
                                  + "surface_z0_l" 
                                  + l_str 
                                  + "_o"
                                  + p_str;
    const std::string& name_x = out_name + "_x.dat";
    const std::string& name_y = out_name + "_y.dat";
    const std::string& name_z = out_name + "_z.dat";

    FullMatrix<double> X(N, N), Y(N, N), Z(N, N);
    Vector<double> values(3);
    for (unsigned int j = 0; j < N; j++){
      for (unsigned int i = 0; i < N; i++){
        const Point<3> eval(0. + i/(N-1.), 0. + j/(N-1.), 0.);
        map.vector_value(eval, values); 

        X(i, j) = values(0); 
        Y(i, j) = values(1);
        Z(i, j) = values(2);
      }
    }

    std::filebuf fx, fy, fz;
    fx.open(name_x.c_str(), std::ios::out);
    fy.open(name_y.c_str(), std::ios::out);
    fz.open(name_z.c_str(), std::ios::out);

    std::ostream out_x(&fx), out_y(&fy), out_z(&fz);
    X.print_formatted(out_x, 16, true, 1, "0");
    Y.print_formatted(out_y, 16, true, 1, "0");
    Z.print_formatted(out_z, 16, true, 1, "0");

    fx.close();
    fy.close();
    fz.close();
  } // print_surface_z0

  void Poisson_Singularity::print_reference_values(
      const double z
  ) const {

    const std::string log_name = "log/active_splines.txt"; 
    std::ofstream log(log_name);
    log << "============================================================\n"
        << ">>>>>>>>>>>>>>> Active splines on level "
        << tria.n_levels() 
        << " <<<<<<<<<<<<<<<<<<<<\n"
        << "============================================================\n\n";
    const auto& T = tria.get_splines();
    for (const auto& t : T) {
      const Point< 3 + 1 >& cp = t -> get_cp();
      log << t -> get_level_index() << ": "; 
      for (unsigned int d = 0; d < 2; d++){
        std::vector< double > kvd = t -> get_local_kv(d);
        log << "[ "; 
        for (const double knot : kvd)
          log << knot << " ";
        log << "] x ";
      }

      {
        std::vector< double > kvd = t -> get_local_kv(2); 
        log << "[ ";
        for (const double knot : kvd)
          log << knot << " ";
        log << "]";
      }

      log << "; CP = (" << cp << ")";


      log << std::endl;
    } // for( t )
    log.close();

    const std::vector<unsigned int>& p = data.deg;
    const unsigned int N = 100;
    const unsigned int m = N*N;
    const unsigned int n_fcns = (p[0] + order + 1)*(p[1] + order + 1)*(p[2] + order + 1);
    std::vector< std::vector< Polynomials::Polynomial<double> >>
            base(3);
    std::vector< Point<3> > 
            reference_evals( m );
    FullMatrix<double> reference_values(m, n_fcns);
    bool success = true;

    // Get the values on the reference at z0
    unsigned int ind = 0;
    for (unsigned int j = 0; j < N; j++)
      for(unsigned int i = 0; i < N; i++)
        reference_evals[ind++] = Point<3>(0. + i / (N - 1.), 0. + j / (N - 1.), z);

    // Define the 3D basis: 
    base[0] = generate_complete_bernstein_basis<double>(p[0] + order);
    base[1] = generate_complete_bernstein_basis<double>(p[1] + order);
    base[2] = generate_complete_bernstein_basis<double>(p[2] + order);
    const AnisotropicPolynomials<3> abp(base);

    // Get values of Bersntein polynomials
    std::vector< double > vals(n_fcns); 
    std::vector< Tensor<1, 3> > grads;
    std::vector< Tensor<2, 3> > grad_grads;
    std::vector< Tensor<3, 3> > third_derivatives;
    std::vector< Tensor<4, 3> > fourth_derivatives;
    for (unsigned int j = 0; j < m; j++){
      abp.evaluate(reference_evals[j], vals, grads, grad_grads, third_derivatives, fourth_derivatives);

      for (unsigned int i = 0; i < n_fcns; i++)
        reference_values(j, i) = vals[i];
    }

    // Check if Bezier extraction is correct
    const auto& splines = tria.get_splines();
    for (const auto& cell : tria.active_cell_iterators()){
      std::cout << "On cell: " << cell -> level() << "." << cell -> index() << std::endl;
      // for each cell, get the list of splines that have support on this cell,
      // and the matrix of bezier extraction operators
      const std::vector<unsigned int>&  spline_ids = tria.get_IEN_array(cell); 
      const FullMatrix<double>&         C = tria.get_bezier_coefficients(cell);
      const unsigned int                n_splines = C.n();
      const Point<3>                    A = cell -> vertex(7) + (-1.) * cell -> vertex(0);
      const Point<3>                    b = cell -> vertex(0);
      
      // Transform each reference point on the cell
      std::vector< Point<3> >           evals(m);
      for (unsigned int i = 0; i < m; i++)
        evals[i] = Point<3>(b(0) + reference_evals[i](0)*A(0), 
                            b(1) + reference_evals[i](1)*A(1), 
                            b(2) + reference_evals[i](2)*A(2));
      
      for (unsigned int pnt = 0; pnt < m; pnt++){
        for (unsigned int j = 0; j < n_splines; j++){
          const unsigned int id = spline_ids[j];
          // for each spline with support on the cell
          const auto& T = *(splines[id]);

          // use the operator C to compare the value with the bernstein values
          const double value_exact = T.value(evals[pnt]);
          double value_approx      = 0;

          for (unsigned int k = 0; k < n_fcns; k++)
            value_approx += C(j, k) * reference_values(pnt, k);

          const double error = std::fabs(value_exact - value_approx);
          if (error > 1e-12) {
            std::cout << "id: " << id << std::endl;
            std::cout << "ref_eval = " << reference_evals[pnt] << std::endl;
            std::cout << "eval     = " << evals[pnt] << std::endl;
            std::cout << "exact = " << value_exact << std::endl;
            std::cout << "appro = " << value_approx << std::endl;
            std::cout << std::scientific << error << std::endl;
            success = false;
          }
        } // for ( j )
      } // for ( pnt )
    } // for (cell)

    IPF_mapping phi;
    std::vector< Tensor<1, 3> > phi_numeric(m);
    std::vector< Tensor<1, 3> > phi_analytic(m);
    for (unsigned int i = 0; i < m; i++){
      Point<4> lin_comb; 
      for (const auto& T : splines){
        const double val = T -> value(reference_evals[i]); 
        lin_comb += val * T -> get_cp();
      } // for ( T ) 
      phi_numeric[i][0] = lin_comb(0) / lin_comb(3);
      phi_numeric[i][1] = lin_comb(1) / lin_comb(3);
      phi_numeric[i][2] = lin_comb(2) / lin_comb(3);
      phi_analytic[i] = phi.tensor_value(reference_evals[i]);

      const Tensor<1, 3>& diff = (phi_analytic[i] - phi_numeric[i]);
      if (diff.norm() > 1e-10) {
        printf("pnt:      %.17g %.17g %.17g\n", reference_evals[i](0), reference_evals[i](1), reference_evals[i](2));
        printf("numeric:  %.17g %.17g %.17g\n", phi_numeric[i][0], phi_numeric[i][1], phi_numeric[i][2]);
        printf("analytic: %.17g %.17g %.17g\n", phi_analytic[i][0], phi_analytic[i][1], phi_analytic[i][2]);
        printf("diff:     %.17g %.17g %.17g\n", diff[0], diff[1], diff[2]);
        printf("error:    %.17g\n", diff.norm());
        success = false;
      }
      const double error = diff.norm();
#ifndef DEBUG
      if (error > 1e-15)
        throw ExcInternalError(); 
#endif
      Assert(error < 1e-15, ExcInternalError());
    } // for ( i )
#ifndef DEBUG
      if (!success)
        throw ExcInternalError(); 
#endif
    Assert(success, ExcInternalError());
  } // print_reference_values

} // namespace Poisson_Neumann
