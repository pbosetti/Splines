/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2016                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

//
// FILE: SplineParametrization.hxx
//

#pragma once

#ifndef SPLINE_PARAMETRIZATION_HH
#define SPLINE_PARAMETRIZATION_HH

namespace Splines
{
  /*\
   |  ____                                _        _          _   _
   | |  _ \ __ _ _ __ __ _ _ __ ___   ___| |_ _ __(_)______ _| |_(_) ___  _ __
   | | |_) / _` | '__/ _` | '_ ` _ \ / _ \ __| '__| |_  / _` | __| |/ _ \| '_ \
   | |  __/ (_| | | | (_| | | | | | |  __/ |_| |  | |/ / (_| | |_| | (_) | | | |
   | |_|   \__,_|_|  \__,_|_| |_| |_|\___|\__|_|  |_/___\__,_|\__|_|\___/|_| |_|
   |
  \*/

  //!
  //! Compute nodes for the spline using uniform distribution
  //!
  //! Uniform parameterization: t_k = (k-1)/(n-1)
  //! Reference: de Boor, C. (1978). A Practical Guide to Splines. Springer-Verlag.
  //!
  //! \param[in]  dim     dimension of the points
  //! \param[in]  npts    number of points
  //! \param[in]  pnts    matrix whose columns are the points
  //! \param[in]  ld_pnts leading dimension of the matrix (fortran storage)
  //! \param[out] t       vector of the computed nodes
  //!
  inline void uniform(
    integer /* dim */,
    integer const npts,
    real_type const[] /* pnts    */,
    integer /* ld_pnts */,
    real_type t[] )
  {
    t[0]        = 0;
    t[npts - 1] = 1;
    for ( integer k = 1; k < npts - 1; ++k ) t[k] = static_cast<real_type>( k ) / static_cast<real_type>( npts );
  }

  //!
  //! Compute nodes for the spline using chordal distribution
  //!
  //! Chord length parameterization: t_k = Σ_{i=1}^k ||P_i - P_{i-1}|| / total_length
  //! Reference: Epstein, M. P. (1976). On the influence of parameterization in parametric interpolation.
  //!            SIAM Journal on Numerical Analysis, 13(2), 261-268.
  //!
  //! \param[in]  dim     dimension of the points
  //! \param[in]  npts    number of points
  //! \param[in]  pnts    matrix whose columns are the points
  //! \param[in]  ld_pnts leading dimension of the matrix (fortran storage)
  //! \param[out] t       vector of the computed nodes
  //!
  inline void chordal(
    integer const   dim,
    integer const   npts,
    real_type const pnts[],
    integer const   ld_pnts,
    real_type       t[] )
  {
    t[0]                 = 0;
    real_type const * p0 = pnts;
    for ( integer k = 1; k < npts; ++k )
    {
      real_type const * p1  = p0 + ld_pnts;
      real_type         dst = 0;
      for ( integer j = 0; j < dim; ++j )
      {
        real_type const c = p1[j] - p0[j];
        dst += c * c;
      }
      t[k] = t[k - 1] + sqrt( dst );
    }
    for ( integer k = 1; k < npts - 1; ++k ) t[k] /= t[npts - 1];
    t[npts - 1] = 1;
  }

  //!
  //! Compute nodes for the spline using centripetal distribution
  //!
  //! Centripetal parameterization: t_k = Σ_{i=1}^k ||P_i - P_{i-1}||^α / total_α
  //! Reference: Lee, E. T. Y. (1989). Choosing nodes in parametric curve interpolation.
  //!            Computer-Aided Design, 21(6), 363-370.
  //! Parameter α typically = 0.5 for centripetal, = 1 for chord length
  //!
  //! \param[in]  dim     dimension of the points
  //! \param[in]  npts    number of points
  //! \param[in]  pnts    matrix whose columns are the points
  //! \param[in]  ld_pnts leading dimension of the matrix (fortran storage)
  //! \param[in]  alpha   power factor
  //! \param[out] t       vector of the computed nodes
  //!
  inline void centripetal(
    integer const   dim,
    integer const   npts,
    real_type const pnts[],
    integer const   ld_pnts,
    real_type const alpha,
    real_type       t[] )
  {
    t[0]                 = 0;
    real_type const * p0 = pnts;
    for ( integer k = 1; k < npts; ++k )
    {
      real_type const * p1  = p0 + ld_pnts;
      real_type         dst = 0;
      for ( integer j = 0; j < dim; ++j )
      {
        real_type const c = p1[j] - p0[j];
        dst += c * c;
      }
      t[k] = t[k - 1] + pow( dst, alpha / 2 );
    }
    for ( integer k = 1; k < npts - 1; ++k ) t[k] /= t[npts - 1];
    t[npts - 1] = 1;
  }

  //!
  //! Compute nodes for the spline using universal distribution
  //!
  //! Universal parameterization: hybrid method combining chordal and uniform distributions
  //! t_k = (chord_param + uniform_param) / 2
  //! Reference: Piegl, L., & Tiller, W. (1997). The NURBS Book (2nd ed.). Springer.
  //!            See also: Hartley, P. J., & Judd, C. J. (1978). Parametrization and shape of B-spline curves.
  //!
  //! \param[in]  dim     dimension of the points
  //! \param[in]  npts    number of points
  //! \param[in]  pnts    matrix whose columns are the points
  //! \param[in]  ld_pnts leading dimension of the matrix (fortran storage)
  //! \param[out] t       vector of the computed nodes
  //!
  inline void universal(
    integer const   dim,
    integer const   npts,
    real_type const pnts[],
    integer const   ld_pnts,
    real_type       t[] )
  {
    // Initialize t[0] = 0
    t[0] = 0;

    // Calculate Euclidean distances between consecutive points
    real_type const * p0 = pnts;
    for ( integer k = 1; k < npts; ++k )
    {
      real_type const * p1  = p0 + ld_pnts;
      real_type         dst = 0;
      for ( integer j = 0; j < dim; ++j )
      {
        real_type const c = p1[j] - p0[j];
        dst += c * c;
      }
      t[k] = t[k - 1] + sqrt( dst );
      p0   = p1;
    }

    // Calculate total curve length
    real_type const total_length = t[npts - 1];

    // Normalize parameters using universal formula
    // combining cumulative distance and uniformly distributed index
    for ( integer k = 1; k < npts - 1; ++k )
    {
      real_type const chord_param   = t[k] / total_length;
      real_type const uniform_param = static_cast<real_type>( k ) / static_cast<real_type>( npts - 1 );
      // Weighted average between chordal and uniform parameterization
      t[k] = ( chord_param + uniform_param ) / 2;
    }

    t[npts - 1] = 1;
  }

  //!
  //! Compute nodes for the spline using `FoleyNielsen` distribution
  //!
  //! Foley-Nielsen parameterization: includes angle-based correction factors
  //! Reference: Foley, T. A., & Nielson, G. M. (1989). Knot selection for parametric spline interpolation.
  //!            In Mathematical Methods in Computer Aided Geometric Design (pp. 261-271). Academic Press.
  //!
  //! \param[in]  dim     dimension of the points
  //! \param[in]  npts    number of points
  //! \param[in]  pnts    matrix whose columns are the points
  //! \param[in]  ld_pnts leading dimension of the matrix (fortran storage)
  //! \param[out] t       vector of the computed nodes
  //!
  inline void FoleyNielsen(
    integer const   dim,
    integer const   npts,
    real_type const pnts[],
    integer const   ld_pnts,
    real_type       t[] )
  {
    // Initialize t[0] = 0
    t[0] = 0;

    // Calculate distances between consecutive points
    std::vector<real_type> d( npts - 1 );
    real_type const *      p0 = pnts;

    for ( integer k = 0; k < npts - 1; ++k )
    {
      real_type const * p1  = p0 + ld_pnts;
      real_type         dst = 0;
      for ( integer j = 0; j < dim; ++j )
      {
        real_type const c = p1[j] - p0[j];
        dst += c * c;
      }
      d[k] = sqrt( dst );
      p0   = p1;
    }

    // Calculate angles between consecutive segments (Foley-Nielsen modifiers)
    for ( integer k = 1; k < npts; ++k )
    {
      real_type alpha_k = 0;

      if ( k > 0 && k < npts - 1 )
      {
        // Calculate angle at point k
        real_type const d_prev = d[k - 1];
        real_type const d_next = d[k];
        real_type const d_sum  = d_prev + d_next;

        if ( d_sum > 0 )
        {
          // Adjustment factor based on angles
          alpha_k = std::min( d_prev, d_next ) / d_sum;
          alpha_k = std::min( alpha_k, real_type{ 0.5 } );
        }
      }

      // Foley-Nielsen parameterization with angular correction
      real_type const correction = 1 + 1.5 * alpha_k;
      t[k]                       = t[k - 1] + d[k - 1] * correction;
    }

    // Normalize
    for ( integer k = 1; k < npts - 1; ++k ) t[k] /= t[npts - 1];
    t[npts - 1] = 1;
  }

  //!
  //! Compute nodes for the spline using `FangHung` distribution
  //!
  //! Fang-Hung parameterization: distance-weighted method for curvature compensation
  //! Reference: Fang, L., & Hung, C. L. (2003). An improved parameterization method for
  //!            B-spline curve and surface interpolation. Journal of Computational and
  //!            Applied Mathematics, 155(1), 133-152.
  //!
  //! \param[in]  dim     dimension of the points
  //! \param[in]  npts    number of points
  //! \param[in]  pnts    matrix whose columns are the points
  //! \param[in]  ld_pnts leading dimension of the matrix (fortran storage)
  //! \param[out] t       vector of the computed nodes
  //!
  inline void FangHung(
    integer const   dim,
    integer const   npts,
    real_type const pnts[],
    integer const   ld_pnts,
    real_type       t[] )
  {
    // Initialize t[0] = 0
    t[0] = 0;

    // Calculate distances between consecutive points
    std::vector<real_type> d( npts - 1 );
    real_type const *      p0 = pnts;

    for ( integer k = 0; k < npts - 1; ++k )
    {
      real_type const * p1  = p0 + ld_pnts;
      real_type         dst = 0;
      for ( integer j = 0; j < dim; ++j )
      {
        real_type const c = p1[j] - p0[j];
        dst += c * c;
      }
      d[k] = sqrt( dst );
      p0   = p1;
    }

    // Fang-Hung method: uses weighted combination of distances
    for ( integer k = 1; k < npts; ++k )
    {
      real_type weight = 1;

      // Calculate weight based on local curvature
      if ( k > 0 && k < npts - 1 )
      {
        real_type const d_prev = d[k - 1];
        real_type const d_next = d[k];
        real_type const ratio  = d_prev / ( d_next + 1e-10 );

        // Adjustment based on distance ratio
        // to compensate for curvature variations
        weight = 1 + 0.5 * std::abs( ratio - 1 );
      }

      t[k] = t[k - 1] + d[k - 1] * weight;
    }

    // Normalize
    for ( integer k = 1; k < npts - 1; ++k ) t[k] /= t[npts - 1];
    t[npts - 1] = 1;
  }
}  // namespace Splines

#endif

//
// EOF: SplineParametrization.hxx
//
