/*
 * TSplineFunction.cc
 *
 *  Created on: Apr 2, 2020
 *      Author: goermer
 */


#include <tspline_function.h>
#include <fstream>



namespace dealt {

  template class TSplineFunction<1,1>;
  template class TSplineFunction<2,2>;
  template class TSplineFunction<3,3>;

  template class TSplineFunction<1,2>;
  template class TSplineFunction<1,3>;
  template class TSplineFunction<2,3>;

  template class TSplineFunction<2,4>;

  template<int dim, int spacedim>
  std::vector<unsigned int> TSplineFunction<dim, spacedim>::s_degree;

  template<int dim, int spacedim>
  int TSplineFunction<dim, spacedim>::s_id = 0;


  template<int dim, int spacedim>
  TSplineFunction<dim, spacedim>::TSplineFunction(
      const std::vector<std::vector<double>>& other_kv,
      const anchor_bounds& other_anchor,
      const ControlPoint& other_cp,
      const std::shared_ptr<TSpline>& other_father1
      ) : dealii::Function<dim>(1)
  {
    kv          = other_kv;
    anchor      = other_anchor;
    cp          = other_cp;
    ind = s_id++;
    father1 = other_father1;
    calc_barycenter();

    if (father1)
      lev = father1 -> level() + 1;
    else
      lev = 0; // ??

    // count multiplicities:
    for (int d = 0; d < dim; d++){
      std::vector<int> m({1});

      double k = kv[d][0];
      int j = 0;
      for (unsigned int i = 1; i < s_degree[d] + 2; i++){
        if (k == kv[d][i])
          m[j]++;
        else {
          m.push_back(1);
          j++;
        }
      }
    }
  }


  template<int dim, int spacedim>
  TSplineFunction<dim, spacedim>::TSplineFunction(
      const std::vector<std::vector<double>>& other_kv,
      const ControlPoint& other_cp,
      const std::shared_ptr<TSpline>& other_father1
      ) : dealii::Function<dim>(1)
  {

  //#ifdef DEBUG_TS
  //  std::cout << indent(2) <<  "Initializing TSpline " << s_id << " ... " << std::endl;
  //#endif

    kv          = other_kv;
    cp          = other_cp;
    ind = s_id++;
    father1 = other_father1;
    calc_barycenter();


    if (father1)
      lev = father1 -> level() + 1;
    else
      lev = 0; // ??

    // count multiplicities:
    std::vector< std::vector< int > > multiplicities;
    for (int d = 0; d < dim; d++){
      std::vector<int> m({1});

      double k = kv[d][0];
      int j = 0;
      for (unsigned int i = 1; i < s_degree[d] + 2; i++){
        if (k == kv[d][i]){
          m[j]++;
        } else {
          k = kv[d][i];
          m.push_back(1);
          j++;
        }
      }
      multiplicities.push_back(m);
    }
    mult = multiplicities;

  // #ifdef DEBUG_TS
  //   std::cout << indent(3) << "Local KV of this spline are: " << std::endl;
  //   for (int d = 0; d < dim; d++){
  //     std::cout << indent(3) << "d = " << d << ": ";
  //     for (const auto& k : kv[d])
  //       std::cout << k << " ";
  //
  //     std::cout << " | With multiplicities : ";
  //     for (const auto& m : mult[d])
  //       std::cout << m << " ";
  //
  //     std::cout << std::endl;
  //   }
  //
  //   std::cout << indent(3) << "ControlPoint is given as: ";
  //   printPoint<spacedim+1>(cp, ind);
  //
  //   if (father1)
  //     std::cout << indent(3) << "ID of father: " << father1 -> index() << std::endl;
  // #endif

    Point<dim> lower;
    Point<dim> upper;

    switch(dim){
      case 1:
        lower(0) = kv[0][std::floor((s_degree[0] + 1.)/2.)];
        upper(0) = kv[0][std::ceil ((s_degree[0] + 1.)/2.)];
        break;
      case 2:
        lower(0) = kv[0][std::floor((s_degree[0] + 1.)/2.)];
        lower(1) = kv[1][std::floor((s_degree[1] + 1.)/2.)];
        upper(0) = kv[0][std::ceil ((s_degree[0] + 1.)/2.)];
        upper(1) = kv[1][std::ceil ((s_degree[1] + 1.)/2.)];
        break;
      case 3:
        lower(0) = kv[0][std::floor((s_degree[0] + 1.)/2.)];
        lower(1) = kv[1][std::floor((s_degree[1] + 1.)/2.)];
        lower(2) = kv[2][std::floor((s_degree[2] + 1.)/2.)];
        upper(0) = kv[0][std::ceil ((s_degree[0] + 1.)/2.)];
        upper(1) = kv[1][std::ceil ((s_degree[1] + 1.)/2.)];
        upper(2) = kv[2][std::ceil ((s_degree[2] + 1.)/2.)];
        break;
    }
    anchor.first = lower;
    anchor.second = upper;
  }

  template<int dim, int spacedim>
  void TSplineFunction<dim, spacedim>::calc_barycenter()
  {
    for (int d = 0; d < dim; d++){
      barycenter(d) = 0;
      for (double k : kv[d])
        barycenter(d) += k;
    }
  }

  template<int dim, int spacedim>
  TSplineFunction<dim, spacedim>& TSplineFunction<dim, spacedim>::operator=(
        TSplineFunction<dim, spacedim>&& other) noexcept
  {
    if (this == &other)
      return *this;

    father1.reset();

    ind     = other.ind;
    mult            = other.mult;
    kv              = other.kv;
    anchor          = other.anchor;
    cp              = other.cp;
    std::swap(father1, other.father1);

    other.father1.reset();

    return *this;
  }


  template<int dim, int spacedim>
  double TSplineFunction<dim, spacedim>::value(
        const Point<dim>& p_d,
        const uint /* component */
  ) const {
    double out = 1;
    unsigned int d = 0;

    while (d < dim /* && std::fabs(out) > 1e-16 */ ) {
      out *= BasisFunSingle1D(p_d(d), d);
      d++;
    }

    // If the value is too small, set it to zero
    //out = std::fabs(out) < 1e-16 ? 0. : out;

    return out;
  }

  template< int dim, int spacedim>
  void TSplineFunction<dim, spacedim>::value_list(
    const std::vector< Point<dim> >&    points,
    std::vector< double  >&             values,
    const unsigned int /* component */
  ) const {
#ifdef DEBUG
    for (unsigned int i = 0; i < points.size(); i++)
      values.at(i) = value(points.at(i));
#else 
    for (unsigned int i = 0; i < points.size(); i++)
      values[i] = value(points[i]);
#endif
  } // vector_value_list


  template<int dim, int spacedim>
  Tensor<1,dim> TSplineFunction<dim,spacedim>::gradient(
      const Point<dim>& p,
      const uint
      ) const
  {
    Tensor<1,dim> out;
    for (uint d = 0; d < dim; d++){
      out[d] = BasisFunDerivSingle1D(1, p(d), d);
      for (uint d1 = 0; d1 < dim; d1++)
        out[d] *= (d1 == d ? 1 : BasisFunSingle1D(p(d1), d1));

      AssertIsFinite(out[d]);
    }
    return out;
  }

  template<int dim, int spacedim>
  double TSplineFunction<dim,spacedim>::laplacian(
      const Point<dim>& p,
      const unsigned int /* component = 0 */
  ) const
  {
    double out = 0;
    double helper = 1;
    Vector<double> BSplines(dim);

    for (unsigned int d = 0; d < dim; d++){
      BSplines(d) = BasisFunSingle1D(p(d), d);
      helper *= BSplines(d);
    }

    if (std::fabs(helper) < 1e-16)
      return 0;

    for (unsigned int d = 0; d < dim; d++){
      double tmp = BasisFunDerivSingle1D(2, p(d), d);
      out += tmp*helper / BSplines(d);
    }
    return out;
  }

  template<int dim, int spacedim>
  SymmetricTensor<2, dim> TSplineFunction<dim,spacedim>::hessian(
      const Point<dim>& p,
      const unsigned int /* i */
  ) const
  {
    SymmetricTensor<2,dim> out;

    double helper = 1;
    Vector<double> BSplines(dim);

    for (unsigned int d = 0; d < dim; d++){
      BSplines(d) = BasisFunSingle1D(p(d),d);
      helper *= BSplines(d); // Computes the function value
      out[d][d] = BasisFunDerivSingle1D(2, p(d), d);
    }
  
    if (helper < 1e-16)
      return out;

    for (unsigned int d1 = 0; d1 < dim; d1++){
      double tmp = BasisFunDerivSingle1D(1, p(d1), d1);
      for (unsigned int d2 = d1 + 1; d2 < dim; d2++){
        out[d1][d2] = tmp*BasisFunDerivSingle1D(1, p(d2), d2);
        out[d1][d2] *= helper/(BSplines(d1)*BSplines(d2));
      }
    }
    return out;
  }
  



  template<int dim, int spacedim>
  double TSplineFunction<dim,spacedim>::BasisFunSingle1D(
      const double u,
      const unsigned int d) const
  {

    std::vector<double> local_kv = kv[d];

    int p = s_degree[d];
    if (u == local_kv[0] && mult[d][0] == p + 1)
      return 1.;
    else if (u == local_kv[p+1] && mult[d][mult[d].size()-1] == p + 1) {
      return 1.;
    }

    if (u < local_kv[0] || u >= local_kv[p+1])
      return 0.;


    Vector<double> Spline_j(p+2);
    for (int i = 0; i < p + 1; i++){
      if (u >= local_kv[i] && u < local_kv[i+1])
        Spline_j(i) = 1.;
    }

    double saved;
    double temp;

    for (int k = 1; k < p + 1; k++){
      if (Spline_j(0) < 1e-16)
        saved = 0.;
      else
        saved = ((u - local_kv[0])*Spline_j[0])/(local_kv[k] - local_kv[0]);

      for (int i = 0; i < p - k + 1 ; i++) {
        if (Spline_j(i+1) < 1e-16){
          Spline_j(i) = saved;
          saved = 0.;
        } else {
          temp = Spline_j(i+1) / (local_kv[i+k+1] - local_kv[i+1]);
          Spline_j(i) = saved + (local_kv[i+k+1] - u)*temp;
          saved = (u - local_kv[i+1])*temp;
        }
      }
    }
    return Spline_j(0);

  } // BasisFunSingle1D


  template<int dim, int spacedim>
  double TSplineFunction<dim,spacedim>::BasisFunDerivSingle1D(
      const unsigned int n,
      const double u,
      const unsigned int d
  ) const {
    // Initialize some variables
    std::vector<double> local_kv = kv[d];
    unsigned int p = s_degree[d];
    std::vector<double> out(n);


    // the algorithm provided computes the derivative from the right.
    // Obviously, at the right end-point of the knot interval the
    // derivative from the right will be zero. Hence, for this special case
    // we compute the derivative from the left.
    if (u == local_kv[p+1]) {
      // initialize 0th degree basis functions
      FullMatrix<double> N(p + 1, p + 1);
      for (unsigned int i = 0; i < p + 1; i++)
        if ((u > local_kv[i] && u <= local_kv[i+1]))
          N(i, 0) = 1.;
        else
          N(i, 0) = 0.;


      // Compute the function values of the B-spline B_j as before.
      // However, we need the matrix used for this to compute the derivatives
      for (unsigned int k = 1; k < p + 1; k++){
        double saved;
        if (std::fabs(N(0, k-1)) < 1e-16)
          saved = 0.;
        else
          saved = ( (u - local_kv[0])*N(0, k-1) )/(local_kv[k] - local_kv[0]);

        AssertIsFinite(saved);

        for (unsigned int i = 0; i < p - k + 1; i++){
          if (std::fabs(N(i+1, k-1)) < 1e-16){
            N(i, k) = saved;
            saved = 0.;
          } else {
            double tmp2 = N(i+1, k-1)/(local_kv[i + k + 1] - local_kv[i + 1]);
            AssertIsFinite(tmp2);
            N(i, k) = saved + (local_kv[i + k + 1] - u)*tmp2;
            saved = (u - local_kv[i + 1])*tmp2;
          } // if (  )
        } // for ( i )
      } // for ( k )


      // Compute the derivatives;
      // loop has set ki (p+1)-times to the right, so we can just go one to the left
      out[0] = N(0, p);
      for (unsigned int k = 1; k < n + 1 ; k++){
        Vector<double> ND(k+1);

        // Load appropriate column
        for (unsigned int i = 0; i < k + 1; i++)
          ND(i) = N(i, p - k);


        for (unsigned int jj = 1; jj < k + 1; jj++){
          double saved;
          if (std::fabs(ND(0)) < 1e-16)
            saved = 0.;
          else
            saved = ND(0)/(local_kv[p - k + jj] - local_kv[0]);
  
          AssertIsFinite(saved);
  
          for (unsigned int i = 0; i < k - jj + 1; i++){
            if (std::fabs(ND(i+1)) < 1e-16) {
              ND(i) = (p - k + jj)*saved;
              saved = 0.;
            } else {
              double tmp2 = ND(i+1)/(local_kv[i + p + jj + 0] - local_kv[i + 0]);
              AssertIsFinite(tmp2);
              ND(i) = (p - k + jj)*(saved - tmp2);
              saved = tmp2;
            } // if ( )
          } // for ( i )
        } // for ( jj )
        out[k-1] = ND(0); /* k-th deivative just popped from heaven */
      } // for ( k )
    } else {
      // Local Property
      if ((u < local_kv[0] || u > local_kv[p+1] ) )
        return 0.;

      // initialize 0th degree basis functions
      FullMatrix<double> N(p + 1, p + 1);
      for (unsigned int i = 0; i < p + 1; i++)
        if ((u >= local_kv[i] && u < local_kv[i+1]))
          N(i, 0) = 1.;
        else
          N(i, 0) = 0.;
  

      // Compute the function values of the B-spline B_j as before.
      // However, we need the matrix used for this to compute the derivatives
      for (unsigned int k = 1; k < p + 1; k++){
        double saved;
        if (std::fabs(N(0, k-1)) < 1e-16)
          saved = 0.;
        else
          saved = ( (u - local_kv[0])*N(0, k-1) )/(local_kv[k] - local_kv[0]);

        AssertIsFinite(saved);

        for (unsigned int i = 0; i < p - k + 1; i++){
          if (std::fabs(N(i+1, k-1)) < 1e-16){
            N(i, k) = saved;
            saved = 0.;
          } else {
            double tmp2 = N(i+1, k-1)/(local_kv[i + k + 1] - local_kv[i + 1]);
            AssertIsFinite(tmp2);
            N(i, k) = saved + (local_kv[i + k + 1] - u)*tmp2;
            saved = (u - local_kv[i + 1])*tmp2;
          } // if (  )
        } // for ( i )
      } // for ( k )


      // Compute the derivatives;
      // loop has set ki (p+1)-times to the right, so we can just go one to the left
      out[0] = N(0, p);
      for (unsigned int k = 1; k < n + 1 ; k++){
        Vector<double> ND(k+1);

        // Load appropriate column
        for (unsigned int i = 0; i < k + 1; i++)
          ND(i) = N(i, p - k);


        for (unsigned int jj = 1; jj < k + 1; jj++){
          double saved;
          if (std::fabs(ND(0)) < 1e-16)
            saved = 0.;
          else
            saved = ND(0)/(local_kv[p - k + jj] - local_kv[0]);

          AssertIsFinite(saved);

          for (unsigned int i = 0; i < k - jj + 1; i++){
            if (std::fabs(ND(i+1)) < 1e-16) {
              ND(i) = (p - k + jj)*saved;
              saved = 0.;
            } else {
              double tmp2 = ND(i+1)/(local_kv[i + p + jj + 0] - local_kv[i + 0]);
              if (std::fabs(ND(i+1)) < 1e-16 && std::fabs(local_kv[i + p + jj + 0] - local_kv[i + 0]) < 1e-16)
                tmp2 = 0.;
              AssertIsFinite(tmp2);
              ND(i) = (p - k + jj)*(saved - tmp2);
              saved = tmp2;
            } // if ( )
          } // for ( i )
        } // for ( jj )
        out[k-1] = ND(0); /* k-th deivative just popped from heaven */
      } // for ( k )
    } // else
    return out[n-1];
  } // BasisFunDerivSingle1D



  template<int dim, int spacedim>
  void TSplineFunction<dim, spacedim>::get_local_kv(
      std::vector<double>& kv_d,
      unsigned int d) const
  {
    Assert(d < dim, ExcIndexRange(d, 0, dim));

    kv_d = kv[d];
  }


  template<int dim, int spacedim>
  bool TSplineFunction<dim, spacedim>::has_support(const cell_iterator& cell) const
  {
    bool support = false;
    for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell && !support; i++)
      support = support || has_support(cell -> vertex(i));

    return support;
  }

  template<int dim, int spacedim>
  bool TSplineFunction<dim, spacedim>::has_support(const active_cell_iterator& cell) const
  {
    bool support = false;
    for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell && !support; i++)
      support = support || has_support(cell -> vertex(i));

    return support;
  }

  template<int dim, int spacedim>
  bool TSplineFunction<dim, spacedim>::has_support(const active_face_iterator& face) const
  {
    bool support = false;
    for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_face && !support; i++){
      support = support || has_support(face -> vertex(i));
    }
  
    return support;
  }
  
  template<int dim, int spacedim>
  bool TSplineFunction<dim, spacedim>::has_support(const Point<dim>& point) const
  {
    const std::vector< unsigned int >& p = s_degree;
    for (int i = 0; i < dim; i++){
      if (point(i) < kv[i][0] || point(i) > kv[i][p[i] + 1])
        return false;
    }
    return true;
  }

  template<int dim, int spacedim>
  bool TSplineFunction<dim, spacedim>::at_boundary(const Point<dim>& point) const
  {
    bool is_at_boundary = true;
    for (int d = 0; d < dim && is_at_boundary; d++)
      is_at_boundary = (kv[d][0] == point(d) || kv[d][s_degree[d] + 1]) && is_at_boundary;

    return is_at_boundary;
  }

  template<int dim, int spacedim>
  auto TSplineFunction<dim, spacedim>::get_anchor(const int& d) const -> std::pair<double, double>
  {
    std::pair<double, double> out(anchor.first(d), anchor.second(d));
    return out;
  }

  template<int dim, int spacedim>
  bool TSplineFunction<dim, spacedim>::is_present(const double& knot, const int& d) const
  {
    for (const auto& k : kv[d])
      if (std::fabs(k - knot) < 1e-16){
        return true;
      }

    return false;
  }


  template<int dim, int spacedim>
  void TSplineFunction<dim, spacedim>::print_spline(
    const std::vector<Point<dim>>& Points,
    const std::string& add) const
  {
    int nPoints = Points.size();
    FullMatrix<double> B(nPoints, 1);

    for (int i = 0; i < nPoints; i++)
      B(i,0) = value(Points[i]);

    std::filebuf f;
    std::string name = "out/" + add + "ts_" + std::to_string(ind) + ".dat";

    f.open(name.c_str(), std::ios::out);
    std::ostream out(&f);

    B.print_formatted(out, 5, true, 1, "0");
    f.close();

    std::cout << indent(2) << " values of TSpline " << ind << " saved in " << name << std::endl;
  }


  template<int dim, int spacedim>
  std::vector<double >TSplineFunction<dim, spacedim>::compute_extended_kv(
      const unsigned int& d,
      unsigned int& nf,
      unsigned int& ne
  ) const {
    std::vector<double> out = kv[d];
    int n = out.size(), k;
    int p = s_degree[d];
    // First knot: count multiplicity
    for(k=0; kv[d][k] == kv[d][k+1]; k++);
    nf = p - k;

    // The first knot has multiplicity k+1.  Aim is multiplicity p+1.
    auto it = out.begin();
    if (nf > 0)
      out.insert( it, nf, kv[d][0] ); // insert missing entries

    // Last knot: count multiplicity
    int count = 1;
    for(k=n-1; kv[d][k] == kv[d][k-1]; k--)
      count++;

    ne = p - count + 1;

    // The last knot has multiplicity k+1.  Aim is multiplicity p+1.
    it = out.end();
    if (ne > 0)
      out.insert( it, ne, kv[d][p+1] ); // insert missing entries

    return out;
  }


  //template<int dim, int spacedim>
  //void TSplineFunction<dim, spacedim>::compute_be_operator()
  //{
  //  for (int d = 0; d < dim; d++)
  //    compute_be_operator(d);
  //}

  template<int dim, int spacedim>
  Vector<double> TSplineFunction<dim, spacedim>::compute_be_operator(
      const std::vector< double >& interior_knots,
      const unsigned int& d,
      const active_cell_iterator& bezier_cell
  ) const {
//    #ifdef DEBUG
//      std::cout << indent(1)
//                << "Computing bezier extraction operator for TSpline "
//                << this -> index()
//                << " in direction "
//                << d
//                << std::endl;
//      std::cout << indent(2)
//                << "KV is given by [";
//      for (const auto& knot : kv[d])
//        std::cout << knot << " ";
//      std::cout << " ];"
//                << std::endl;
//  
//      std::cout << indent(2)
//                << "Interior Knots are [";
//      for (const auto& knot : interior_knots)
//        std::cout << knot << " ";
//      std::cout << "];"
//                << std::endl;
//    #endif

    unsigned int nf, ne;
    auto extended_kv_d = compute_extended_kv(d, nf, ne);

    int n_elements = 0;
    int p = s_degree[d];

    int m = interior_knots.size();
    int mbar = ne + nf + m + 1;
    // count amount of elements / rows
    for (unsigned int k = 0; k < extended_kv_d.size() - 1; k++)
      if (extended_kv_d.at(k) != extended_kv_d.at(k+1))
        n_elements++;

    FullMatrix<double> local_be_operator(n_elements + m, p + 1);

    // initialize variables
    int a = p;
    int b = a+1;
    int nb = 0;
    int si = 0;
    local_be_operator(0, nf) = 1;
    while (b < mbar){
      // Some variables we need later
      int mult;
      int add = 0;

      // Count multiplicity of knot Ubar(b)
      if (si < m && !(extended_kv_d[b] < interior_knots[si] )){
        // knot from interior knot needs to be inserted before
        add = 1;
        mult = 0;  // This knot needs to be inserted p times, orginially it is not present
        extended_kv_d.insert(extended_kv_d.begin() + b, interior_knots[si]);
        Assert(extended_kv_d[b] == interior_knots[si], ExcInternalError());
        si++;
      } else {
        int i = b;
        for (; b < mbar && extended_kv_d[b+1] == extended_kv_d[b]; b++);
        mult = b - i + 1;
      }

      // Initialize the next extraction operator row with a 1 at the right place.
      // This is needed in some cases, but may be overwritten later
      local_be_operator(nb + 1, nf + 1 - (nb+1) + si - mult - add) = 1;

      // Perform knot insertion for bezier extraction
      if (mult < p){
        double numer = extended_kv_d.at(b) - extended_kv_d.at(a);
        std::vector< double > alphas(p - mult);

        for (int j = p; j > mult; j--){
          alphas.at((j-1) - mult) = numer/( extended_kv_d.at(a + j + add)
                                            - extended_kv_d.at(a) );
        }

        // repeat knot insertion procedure r times
        int r = p - mult;

        // Update matrix coefficients
        for (int j = 1; j < r+1; j++){
          int save = r - j;
          int s = mult + j; // This many new points

          for (int k = p; k > s - 1; k--){
            local_be_operator(nb, k) = alphas.at(k-s)*local_be_operator(nb, k) +
                                          (1. - alphas.at(k-s))*local_be_operator(nb, k-1);
          } // for( k )

          if (b < mbar /* && nb + 1 < n_elements + m */)
            local_be_operator(nb+1, save) = local_be_operator(nb, p);
        } // for( j )

        nb++; // finished with current row

        if (b < mbar){
          // Update indices for next iteration
          a = b;
          b++;
        }
      } else {
        // Special case:
        local_be_operator(nb + 1, 0) = 1.;
        a = b; 
        b++; nb++;
      }// if (mult)
    } // while ( b )

  // #ifdef DEBUG_TS
  //   std::cout << indent(3) << "Extraction Operator rows " << std::endl;
  //   int rows = local_be_operator.m();
  //   int cols = local_be_operator.n();
  //
  //   std::cout << indent(4) << "[";
  //   for (int i = 0; i < rows-1; i++){
  //     for (int j = 0; j < cols; j++)
  //       std::cout << local_be_operator(i, j) << " ";
  //
  //     std::cout << ";" << std::endl << indent(4);
  //   }
  //   for (int j = 0; j < cols; j++)
  //     std::cout << local_be_operator(rows-1,j) << " ";
  //   std::cout << "];" << std::endl;
  // #endif

    // Get the vector of singleton entries in Ubar = extended_kv_d
    std::vector<double> bezier_kv = extended_kv_d;
    auto it = bezier_kv.begin();
    for (; it != bezier_kv.end()-1 ;){
      if (*it == *(it+1))
        bezier_kv.erase(it);
      else
        ++it;
    }

    // Get the index where we hit lower
    const Point<dim>& lower = bezier_cell -> vertex(0);
    unsigned int i = 0;
    for (; bezier_kv[i] != lower(d); i++);

    Vector<double> out(p + 1);
    for (int j = 0; j < p+1; j++)
      out(j) = local_be_operator(i, j);

  // #ifdef DEBUG_TS
  //   std::cout << indent(3)
  //             << "Row used for bezier cell:"
  //             << bezier_cell -> vertex(0)
  //             << " x "
  //             << bezier_cell -> vertex(GeometryInfo<dim>::vertices_per_cell - 1)
  //             << std::endl;
  //   std::cout << indent(4) << "[";
  //   for (int i = 0; i < p + 1; i++)
  //     std::cout << out(i) << " ";
  //   std::cout << "];" << std::endl;
  // #endif

    return out;
  } // compute_be_operator

  template<int dim, int spacedim>
  void TSplineFunction<dim, spacedim>::compute_be_operator(
      const std::vector< double >               &interior_knots,
      const unsigned int                        &d,
      const std::vector< active_cell_iterator > &cell_list,
            std::vector< Vector< double > >     &ops
  ) const {
//    #ifdef DEBUG
//      std::cout << indent(1)
//                << "Computing bezier extraction operator for TSpline "
//                << this -> index()
//                << " in direction "
//                << d
//                << std::endl;
//      std::cout << indent(2)
//                << "KV is given by [";
//      for (const auto& knot : kv[d])
//        std::cout << knot << " ";
//      std::cout << " ];"
//                << std::endl;
//  
//      std::cout << indent(2)
//                << "Interior Knots are [";
//      for (const auto& knot : interior_knots)
//        std::cout << knot << " ";
//      std::cout << "];"
//                << std::endl;
//    #endif
    Assert(ops.size() == cell_list.size(),
             ExcDimensionMismatch(ops.size(), cell_list.size()));

    // Compute extraction operator as before:
    unsigned int nf, ne;
    auto extended_kv_d = compute_extended_kv(d, nf, ne);

    int n_elements = 0;
    int p = s_degree[d];

    int m = interior_knots.size();
    int mbar = ne + nf + m + 1;
    // count amount of elements / rows
    for (unsigned int k = 0; k < extended_kv_d.size() - 1; k++)
      if (extended_kv_d.at(k) != extended_kv_d.at(k+1))
        n_elements++;


    FullMatrix<double> local_be_operator(n_elements + m, p + 1);

    // initialize variables
    int a = p;
    int b = a+1;
    int nb = 0;
    int si = 0;
    local_be_operator(0, nf) = 1;
    while (b < mbar){
      // Some variables we need later
      int mult;
      int add = 0;

      // Count multiplicity of knot Ubar(b)
      if (si < m && !(extended_kv_d[b] < interior_knots[si] )){
        // knot from interior knot needs to be inserted before
        add = 1;
        mult = 0;  // This knot needs to be inserted p times, orginially it is not present
        extended_kv_d.insert(extended_kv_d.begin() + b, interior_knots[si]);
        Assert(extended_kv_d[b] == interior_knots[si], ExcInternalError());
        si++;
      } else {
        int i = b;
        for (; b < mbar && extended_kv_d[b+1] == extended_kv_d[b]; b++);
        mult = b - i + 1;
      }

      // Initialize the next extraction operator row with a 1 at the right place.
      // This is needed in some cases, but may be overwritten later
      local_be_operator(nb + 1, nf + 1 - (nb+1) + si - mult - add) = 1;

      // Perform knot insertion for bezier extraction
      if (mult < p){
        double numer = extended_kv_d.at(b) - extended_kv_d.at(a);
        std::vector< double > alphas(p - mult);

        for (int j = p; j > mult; j--){
          alphas.at((j-1) - mult) = numer/( extended_kv_d.at(a + j + add)
                                            - extended_kv_d.at(a) );
        }

        // repeat knot insertion procedure r times
        int r = p - mult;

        // Update matrix coefficients
        for (int j = 1; j < r+1; j++){
          int save = r - j;
          int s = mult + j; // This many new points

          for (int k = p; k > s - 1; k--){
            local_be_operator(nb, k) = alphas.at(k-s)*local_be_operator(nb, k) +
                                          (1. - alphas.at(k-s))*local_be_operator(nb, k-1);
          } // for( k )

          if (b < mbar /* && nb + 1 < n_elements + m */)
            local_be_operator(nb+1, save) = local_be_operator(nb, p);
        } // for( j )

        nb++; // finished with current row

        if (b < mbar){
          // Update indices for next iteration
          a = b;
          b++;
        }
      } else {
        // Special case:
        local_be_operator(nb + 1, 0) = 1.;
        a = b;
        b++; nb++;
      }// if (mult)
    } // while ( b )

    // Get the vector of singleton entries in Ubar = extended_kv_d
    std::vector<double> bezier_kv = extended_kv_d;
    auto it = bezier_kv.begin();
    for (; it != bezier_kv.end()-1 ;){
      if (*it == *(it+1))
        bezier_kv.erase(it);
      else
        ++it;
    }

    // after this is done, store the corresponding rows of the
    // operators for each cell in the Vector of ops.

    const unsigned int n_cells = cell_list.size();
    for (unsigned n = 0; n < n_cells; n++){
      const auto& cell    = cell_list[n];
      ops[n].reinit(p+1);

      // Get the index where we hit lower
      const double lower = (cell -> vertex(0)).operator()(d);
      unsigned int i = 0;
      for (; bezier_kv[i] != lower; i++);

      for (int j = 0; j < p+1; j++)
        ops[n](j) = local_be_operator(i, j);
    }
  } // compute_be_operator


  template<int dim, int spacedim>
  double TSplineFunction<dim, spacedim>::get_support_volume()
  {
    double out = 1;
    for (int d = 0; d < dim; d++)
      out *= (kv[d][s_degree[d] + 1] - kv[d][0]);

    return out;
  }


} // namespace dealt












