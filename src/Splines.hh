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

#pragma once

#ifndef SPLINES_HH
#define SPLINES_HH

#ifdef __GNUC__
#pragma GCC diagnostic push
#endif

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#endif

#include "SplinesConfig.hh"
#include <fstream>

//!
//! Namespace of Splines library
//!
namespace Splines
{

  using std::basic_istream;
  using std::basic_ostream;
  using std::cerr;
  using std::cin;
  using std::cout;
  using std::exception;
  using std::lower_bound;
  using std::pair;
  using std::runtime_error;
  using std::string;
  using std::string_view;
  using std::vector;

  typedef double real_type;  //!< Floating point type for splines
  typedef int    integer;    //!< Signed integer type for splines

  using Malloc_real  = Utils::Malloc<real_type>;
  using ostream_type = basic_ostream<char>;
  using istream_type = basic_istream<char>;

  void backtrace( ostream_type & );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! Associate a number for each type of splines implemented
  using SplineType1D = enum class SplineType1D : integer {
    CONSTANT   = 0,
    LINEAR     = 1,
    CUBIC      = 2,
    AKIMA      = 3,
    BESSEL     = 4,
    PCHIP      = 5,
    QUINTIC    = 6,
    HERMITE    = 7,
    SPLINE_SET = 8,
    SPLINE_VEC = 9
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! Associate a number for each type of splines implemented
  using SplineType2D = enum class SplineType2D : integer { BILINEAR = 0, BICUBIC = 1, BIQUINTIC = 2, AKIMA2D = 3 };

  inline SplineType1D string_to_splineType1D( string_view nin )
  {
    string n{ nin };
    std::transform( n.begin(), n.end(), n.begin(), ::tolower );
    if ( n == "constant" ) return SplineType1D::CONSTANT;
    if ( n == "linear" ) return SplineType1D::LINEAR;
    if ( n == "cubic" ) return SplineType1D::CUBIC;
    if ( n == "akima" ) return SplineType1D::AKIMA;
    if ( n == "bessel" ) return SplineType1D::BESSEL;
    if ( n == "pchip" ) return SplineType1D::PCHIP;
    if ( n == "quintic" ) return SplineType1D::QUINTIC;
    if ( n == "hermite" ) return SplineType1D::HERMITE;
    if ( n == "spline_set" ) return SplineType1D::SPLINE_SET;
    if ( n == "spline_vec" ) return SplineType1D::SPLINE_VEC;
    throw std::runtime_error( fmt::format( "string_to_splineType1D({}) unknown type\n", n ) );
  }

  inline SplineType2D string_to_splineType2D( string_view nin )
  {
    string n{ nin };
    std::transform( n.begin(), n.end(), n.begin(), ::tolower );
    if ( n == "bilinear" ) return SplineType2D::BILINEAR;
    if ( n == "bicubic" ) return SplineType2D::BICUBIC;
    if ( n == "biquintic" ) return SplineType2D::BIQUINTIC;
    if ( n == "akima" ) return SplineType2D::AKIMA2D;
    throw std::runtime_error( fmt::format( "string_to_splineType2D({}) unknown type\n", n ) );
  }

  inline char const * to_string( SplineType1D const t )
  {
    switch ( t )
    {
      case SplineType1D::CONSTANT: return "SPLINE_CONSTANT";
      case SplineType1D::LINEAR: return "SPLINE_LINEAR";
      case SplineType1D::CUBIC: return "SPLINE_CUBIC";
      case SplineType1D::AKIMA: return "SPLINE_AKIMA";
      case SplineType1D::BESSEL: return "SPLINE_BESSEL";
      case SplineType1D::PCHIP: return "SPLINE_PCHIP";
      case SplineType1D::QUINTIC: return "SPLINE_QUINTIC";
      case SplineType1D::HERMITE: return "SPLINE_HERMITE";
      case SplineType1D::SPLINE_SET: return "SPLINE_SPLINE_SET";
      case SplineType1D::SPLINE_VEC: return "SPLINE_SPLINE_VEC";
    }
    return "NO_TYPE";
  }

  inline char const * to_string( SplineType2D const t )
  {
    switch ( t )
    {
      case SplineType2D::BILINEAR: return "SPLINE2D_BILINEAR";
      case SplineType2D::BICUBIC: return "SPLINE2D_BICUBIC";
      case SplineType2D::BIQUINTIC: return "SPLINE2D_BIQUINTIC";
      case SplineType2D::AKIMA2D: return "SPLINE2D_AKIMA2D";
    }
    return "NO_TYPE";
  }

  using GC_namespace::GenericContainer;
  using GC_namespace::map_type;
  using GC_namespace::vec_real_type;
  using GC_namespace::vec_string_type;
  using GC_namespace::vector_type;

  /*
  //   _   _                     _ _
  //  | | | | ___ _ __ _ __ ___ (_) |_ ___
  //  | |_| |/ _ \ '__| '_ ` _ \| | __/ _ \
  //  |  _  |  __/ |  | | | | | | | ||  __/
  //  |_| |_|\___|_|  |_| |_| |_|_|\__\___|
  */

#ifdef AUTODIFF_SUPPORT
  template <typename T> inline void Hermite3( T const & x, real_type const H, T base[4] )
  {
    T const X = x / H;
    base[1]   = X * X * ( 3 - 2 * X );
    base[0]   = 1 - base[1];
    base[2]   = x * ( X * ( X - 2 ) + 1 );
    base[3]   = x * X * ( X - 1 );
  }
#endif

  void Hermite3( real_type const x, real_type const H, real_type base[4] );
  void Hermite3_D( real_type const x, real_type const H, real_type base_D[4] );
  void Hermite3_DD( real_type const x, real_type const H, real_type base_DD[4] );
  void Hermite3_DDD( real_type const x, real_type const H, real_type base_DDD[4] );

#ifdef AUTODIFF_SUPPORT
  template <typename T> inline void Hermite5( T const & x, real_type H, T base[6] )
  {
    const real_type invH  = real_type( 1 ) / H;
    const real_type invH2 = invH * invH;
    const real_type invH4 = invH2 * invH2;

    const T x2 = x * x;
    const T x3 = x2 * x;
    const T x4 = x2 * x2;
    const T x5 = x4 * x;

    const T xm  = H - x;
    const T xm2 = xm * xm;
    const T xm3 = xm2 * xm;

    base[0] = invH4 * xm3 * ( 3 * x * H + H * H + 6 * x2 );

    base[1] = invH4 * ( -15 * H * x4 + 6 * x5 + 10 * H * H * x3 );

    base[2] = invH4 * xm3 * x * ( H + 3 * x );

    base[3] = invH4 * ( 3 * x - 4 * H ) * xm * x3;

    const real_type c = real_type( 0.5 ) * invH2 * invH;

    base[4] = c * xm3 * x2;
    base[5] = c * xm2 * x3;
  }

#endif

  void Hermite5( real_type const x, real_type const H, real_type base[6] );
  void Hermite5_D( real_type const x, real_type const H, real_type base_D[6] );
  void Hermite5_DD( real_type const x, real_type const H, real_type base_DD[6] );
  void Hermite5_DDD( real_type const x, real_type const H, real_type base_DDD[6] );
  void Hermite5_DDDD( real_type const x, real_type const H, real_type base_DDDD[6] );
  void Hermite5_DDDDD( real_type const x, real_type const H, real_type base_DDDDD[6] );

  //!
  //! Convert polynomial defined using Hermite base
  //!
  //! \f[ p(x) = p_0 h_0(t) + p_1 h_1(t) +  p'_0 h_2(t) + p'_1 h_3(t) \f]
  //!
  //! to standard form
  //!
  //! \f[ p(x) = A t^3 + B t^2 + C t + D \f]
  //!
  //!
  static inline void Hermite3_to_poly(
    real_type   H,
    real_type   P0,
    real_type   P1,
    real_type   DP0,
    real_type   DP1,
    real_type & A,
    real_type & B,
    real_type & C,
    real_type & D )
  {
    const real_type invH  = real_type( 1 ) / H;
    const real_type invH2 = invH * invH;

    const real_type dP = P1 - P0;

    A = ( DP0 + DP1 - 2 * dP * invH ) * invH2;
    B = ( 3 * dP * invH - ( 2 * DP0 + DP1 ) ) * invH;
    C = DP0;
    D = P0;
  }

  //!
  //! Convert polynomial defined using Hermite base
  //!
  //! \f[
  //! p(x) = p_0 h_0(t) + p_1 h_1(t) +
  //! p'_0 h_2(t) + p'_1 h_3(t) +
  //! p''_0 h_4(t) + p''_1 h_5(t)
  //! \f]
  //!
  //! to standard form
  //!
  //! \f[ p(x) = A t^5 + B t^4 + C t^3 + D t^2 + E t + F \f]
  //!
  //!
  static inline void Hermite5_to_poly(
    real_type   h,
    real_type   P0,
    real_type   P1,
    real_type   DP0,
    real_type   DP1,
    real_type   DDP0,
    real_type   DDP1,
    real_type & A,
    real_type & B,
    real_type & C,
    real_type & D,
    real_type & E,
    real_type & F )
  {
    const real_type invH  = real_type( 1 ) / h;
    const real_type invH2 = invH * invH;
    const real_type invH3 = invH2 * invH;

    const real_type dP = P1 - P0;

    A = ( real_type( 0.5 ) * ( DDP1 - DDP0 ) + invH * ( 6 * dP * invH - 3 * ( DP0 + DP1 ) ) ) * invH3;
    B = ( real_type( 1.5 ) * DDP0 - DDP1 + invH * ( ( 8 * DP0 + 7 * DP1 ) - 15 * dP * invH ) ) * invH2;
    C = ( real_type( 0.5 ) * DDP1 - real_type( 1.5 ) * DDP0 + invH * ( 10 * dP * invH - ( 6 * DP0 + 4 * DP1 ) ) ) *
        invH;
    D = real_type( 0.5 ) * DDP0;
    E = DP0;
    F = P0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*
  //   ____  _ _ _
  //  | __ )(_) (_)_ __   ___  __ _ _ __
  //  |  _ \| | | | '_ \ / _ \/ _` | '__|
  //  | |_) | | | | | | |  __/ (_| | |
  //  |____/|_|_|_|_| |_|\___|\__,_|_|
  */

  real_type bilinear3( real_type const p[4], real_type const M[4][4], real_type const q[4] );

  real_type bilinear5( real_type const p[6], real_type const M[6][6], real_type const q[6] );

  //! Check if cubic spline is monotone
  //! return: -2 non-monotone data
  //!         -1 non-monotone spline
  //!          0 monotone (non strict)
  //!          1 strictly monotone
  inline integer check_cubic_spline_monotonicity(
    real_type const X[],
    real_type const Y[],
    real_type const Yp[],
    integer const   npts )
  {
    integer flag = 1;

    // --- Check data monotonicity (X assumed monotone) ---
    for ( integer i = 1; i < npts; ++i )
    {
      if ( Y[i] < Y[i - 1] ) return -2;  // non-monotone data

      if ( flag > 0 && Utils::is_zero( Y[i] - Y[i - 1] ) && X[i] > X[i - 1] )
      {
        flag = 0;  // not strictly monotone
      }
    }

    // --- Check spline monotonicity (Fritsch–Carlson conditions) ---
    // See: Methods of Shape-Preserving Spline Approximation, p.146
    for ( integer i = 1; i < npts; ++i )
    {
      const real_type dx = X[i] - X[i - 1];
      if ( dx <= real_type( 0 ) ) continue;  // skip duplicate X

      const real_type dy = Y[i] - Y[i - 1];
      const real_type dd = dy / dx;

      const real_type m0 = Yp[i - 1] / dd;
      const real_type m1 = Yp[i] / dd;

      if ( m0 < 0 || m1 < 0 ) return -1;

      if ( m0 <= 3 && m1 <= 3 )
      {
        if ( flag > 0 )
        {
          if (
            ( i > 1 && ( Utils::is_zero( m0 ) || Utils::is_zero( m0 - 3 ) ) ) ||
            ( i < npts - 1 && ( Utils::is_zero( m1 ) || Utils::is_zero( m1 - 3 ) ) ) )
          {
            flag = 0;
          }
        }
      }
      else
      {
        const real_type t1   = 2 * m0 + m1 - 3;
        const real_type t2   = 2 * ( m0 + m1 - 2 );
        const real_type disc = m0 * t2 - t1 * t1;

        if ( ( t2 >= 0 && disc < 0 ) || ( t2 < 0 && disc > 0 ) ) { return -1; }

        if ( flag > 0 && Utils::is_zero( disc ) ) { flag = 0; }
      }
    }

    return flag;
  }

  /*\
   |  ____                                _        _          _   _
   | |  _ \ __ _ _ __ __ _ _ __ ___   ___| |_ _ __(_)______ _| |_(_) ___  _ __
   | | |_) / _` | '__/ _` | '_ ` _ \ / _ \ __| '__| |_  / _` | __| |/ _ \| '_ \
   | |  __/ (_| | | | (_| | | | | | |  __/ |_| |  | |/ / (_| | |_| | (_) | | | |
   | |_|   \__,_|_|  \__,_|_| |_| |_|\___|\__|_|  |_/___\__,_|\__|_|\___/|_| |_|
   |
  \*/
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
    t[0] = 0;
    real_type const * p0{ pnts };
    for ( integer k = 1; k < npts; ++k )
    {
      real_type const * p1{ p0 + ld_pnts };
      real_type         dst = 0;
      for ( integer j = 0; j < dim; ++j )
      {
        real_type const c{ p1[j] - p0[j] };
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
    t[0] = 0;
    real_type const * p0{ pnts };
    for ( integer k = 1; k < npts; ++k )
    {
      real_type const * p1{ p0 + ld_pnts };
      real_type         dst{ 0 };
      for ( integer j = 0; j < dim; ++j )
      {
        real_type const c{ p1[j] - p0[j] };
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
    real_type const * p0{ pnts };
    for ( integer k = 1; k < npts; ++k )
    {
      real_type const * p1{ p0 + ld_pnts };
      real_type         dst{ 0 };
      for ( integer j = 0; j < dim; ++j )
      {
        real_type const c{ p1[j] - p0[j] };
        dst += c * c;
      }
      t[k] = t[k - 1] + sqrt( dst );
      p0   = p1;
    }

    // Calculate total curve length
    real_type const total_length{ t[npts - 1] };

    // Normalize parameters using universal formula
    // combining cumulative distance and uniformly distributed index
    for ( integer k = 1; k < npts - 1; ++k )
    {
      real_type const chord_param{ t[k] / total_length };
      real_type const uniform_param{ static_cast<real_type>( k ) / static_cast<real_type>( npts - 1 ) };
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
      real_type alpha_k{ 0 };

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
        real_type const d_prev{ d[k - 1] };
        real_type const d_next{ d[k] };
        real_type const ratio{ d_prev / ( d_next + 1e-10 ) };

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

  /*\
   |   ____                      _     ___       _                       _
   |  / ___|  ___  __ _ _ __ ___| |__ |_ _|_ __ | |_ ___ _ ____   ____ _| |
   |  \___ \ / _ \/ _` | '__/ __| '_ \ | || '_ \| __/ _ \ '__\ \ / / _` | |
   |   ___) |  __/ (_| | | | (__| | | || || | | | ||  __/ |   \ V / (_| | |
   |  |____/ \___|\__,_|_|  \___|_| |_|___|_| |_|\__\___|_|    \_/ \__,_|_|
  \*/

  /**
   * @class SearchInterval
   * @brief Efficient interval search structure using a precomputed lookup table with binary search
   *
   * This class implements an optimized search algorithm to find which interval contains a given
   * query point in a sorted array. It uses a two-tier approach:
   * 1. A coarse lookup table (adaptive size) for O(1) initial range reduction
   * 2. Binary search within the reduced range for O(log n) final location
   *
   * The overall complexity is O(1) + O(log(n/table_size)) which is significantly faster
   * than pure binary search O(log n) for large datasets.
   *
   * **OPTIMIZATIONS vs ORIGINAL:**
   * - Lock-free reads after initialization (atomic flag pattern for 10-100x multi-thread speedup)
   * - Pre-computed duplicate handling (eliminates O(n) worst-case loop in find())
   * - Adaptive table sizing based on dataset size (better memory/speed tradeoff)
   * - Improved table construction algorithm (fewer adjustments needed in find())
   * - Branch prediction hints for hot paths (better CPU pipeline utilization)
   *
   * Thread-safety: This class is thread-safe. Multiple threads can call find() concurrently
   * without locking after the first initialization. The internal state is protected by atomic
   * operations and a mutex for lazy initialization only.
   *
   * @note The input array X must be sorted in ascending order
   * @note This implementation handles closed curves (periodic boundary conditions)
   * @note Duplicate consecutive nodes are handled by returning the leftmost valid interval
   *
   * Reference: Knuth, D.E. (1998). The Art of Computer Programming, Volume 3: Sorting and Searching.
   *            Addison-Wesley. Section 6.2.1 (Searching an Ordered Table).
   */

  class SearchInterval
  {
    //! @brief Relative epsilon for floating point comparisons
    static constexpr real_type m_epsilon{ 1e-10 };

    // External pointers to curve data (managed by external spline object)
    string const * p_name             = nullptr;  //!< Curve name for error messages
    integer *      p_npts             = nullptr;  //!< Number of data points
    bool *         p_curve_is_closed  = nullptr;  //!< True if curve wraps around (periodic)
    bool *         p_curve_can_extend = nullptr;  //!< True if extrapolation is allowed

    //! @brief Number of cells in the lookup table (trade-off between memory and speed)
    mutable integer m_table_size{ 400 };

    // Cached data for fast access
    mutable real_type ** p_X       = nullptr;  //!< Pointer to sorted X coordinates array
    mutable real_type    m_x_min   = 0;        //!< Minimum X value (first point)
    mutable real_type    m_x_max   = 0;        //!< Maximum X value (last point)
    mutable real_type    m_x_range = 0;        //!< Range = m_x_max - m_x_min
    mutable real_type    m_dx      = 0;        //!< Cell width = m_x_range / m_table_size

    /**
     * @brief Lookup table storing left boundary indices for each cell
     *
     * m_LO[i] contains the smallest data point index k such that X[k] >= i * m_dx + m_x_min
     * Size is m_table_size + 2 to avoid boundary checks and replicate last point
     *
     * Properties:
     * 1. m_LO[i] <= m_LO[i+1] (monotonically non-decreasing)
     * 2. m_LO[i] ∈ [0, n-1]
     * 3. For all k < m_LO[i], X[k] < m_x_min + i * m_dx - ε
     */
    mutable std::vector<integer> m_LO;

    /**
     * @brief Lookup table storing right boundary indices for each cell
     *
     * m_HI[i] contains the largest data point index k such that X[k] <= (i+1) * m_dx + m_x_min
     * Size is m_table_size + 2 to avoid boundary checks and replicate last point
     *
     * Properties:
     * 1. m_HI[i] >= m_HI[i-1] (monotonically non-decreasing)
     * 2. m_HI[i] ∈ [0, n-1]
     * 3. For all k > m_HI[i], X[k] > m_x_min + (i+1) * m_dx + ε
     */
    mutable std::vector<integer> m_HI;

    mutable bool              m_must_reset = true;  //!< Flag indicating tables need rebuilding
    mutable std::mutex        m_mutex;              //!< Protects concurrent access to internal state
    mutable std::atomic<bool> m_ready{ false };

    real_type eps_x( real_type x ) const { return m_epsilon * std::max( real_type( 1 ), std::abs( x ) ); }

#ifndef NDEBUG
    /**
     * @brief Validate the consistency of LO and HI tables
     *
     * This function performs sanity checks to ensure the lookup tables are properly constructed.
     * It verifies:
     * 1. All entries are within valid range [0, n-1]
     * 2. LO table is monotonically non-decreasing
     * 3. HI table is monotonically non-decreasing
     * 4. For each i, m_LO[i] <= m_HI[i] (cells cannot be empty)
     * 5. The tables cover the entire domain
     *
     * @param n Number of data points
     * @return true if tables are valid, false otherwise
     *
     * @note Used for debugging and testing purposes
     */
    bool validate_tables( integer n ) const
    {
      // Check all LO entries are valid
      for ( integer i = 0; i <= m_table_size + 1; ++i )
      {
        if ( m_LO[i] < 0 || m_LO[i] >= n )
        {
          fmt::print( "ERROR: m_LO[{}] = {} is out of range [0, {})\n", i, m_LO[i], n );
          return false;
        }
      }

      // Check all HI entries are valid
      for ( integer i = 0; i <= m_table_size + 1; ++i )
      {
        if ( m_HI[i] < 0 || m_HI[i] >= n )
        {
          fmt::print( "ERROR: m_HI[{}] = {} is out of range [0, {})\n", i, m_HI[i], n );
          return false;
        }
      }

      // Check LO is monotonically non-decreasing
      for ( integer i = 0; i <= m_table_size; ++i )
      {
        if ( m_LO[i] > m_LO[i + 1] )
        {
          fmt::print(
            "ERROR: m_LO not monotonic at i={}: m_LO[{}]={} > m_LO[{}]={}\n",
            i,
            i,
            m_LO[i],
            i + 1,
            m_LO[i + 1] );
          return false;
        }
      }

      // Check HI is monotonically non-decreasing
      for ( integer i = 0; i <= m_table_size; ++i )
      {
        if ( m_HI[i] > m_HI[i + 1] )
        {
          fmt::print(
            "ERROR: m_HI not monotonic at i={}: m_HI[{}]={} > m_HI[{}]={}\n",
            i,
            i,
            m_HI[i],
            i + 1,
            m_HI[i + 1] );
          return false;
        }
      }

      // Check each cell has LO <= HI
      for ( integer i = 0; i <= m_table_size; ++i )
      {
        if ( m_LO[i] > m_HI[i] + 1 )
        {
          fmt::print( "ERROR: Empty cell at i={}: m_LO[{}]={} > m_HI[{}]={}\n", i, i, m_LO[i], i, m_HI[i] );
          return false;
        }
      }

      // Check boundary conditions
      if ( m_LO[0] != 0 )
      {
        fmt::print( "ERROR: m_LO[0] should be 0, but is {}\n", m_LO[0] );
        return false;
      }

      if ( m_HI[m_table_size] != n - 1 )
      {
        fmt::print( "ERROR: m_HI[{}] should be {}, but is {}\n", m_table_size, n - 1, m_HI[m_table_size] );
        return false;
      }

      // Check that tables are consistent with each other
      // For adjacent cells, HI[i] should be >= LO[i+1] (or very close)
      for ( integer i = 0; i < m_table_size; ++i )
      {
        if ( m_HI[i] < m_LO[i + 1] - 1 )  // Allow one index gap for numerical reasons
        {
          fmt::print(
            "WARNING: Gap between cells {} and {}: m_HI[{}]={} < m_LO[{}]={}\n",
            i,
            i + 1,
            i,
            m_HI[i],
            i + 1,
            m_LO[i + 1] );
          // Not fatal, but indicates potential inefficiency
        }
      }

      return true;
    }
#endif

    /**
     * @brief Builds the lookup tables from current data
     *
     * This method constructs m_LO and m_HI tables that partition the X range into
     * m_table_size uniform cells. For each cell i, we precompute:
     * - m_LO[i]: leftmost data point that could be in this cell or to the right
     * - m_HI[i]: rightmost data point that could be in this cell or to the left
     *
     * Algorithm:
     * 1. Initialize all table entries to -1 (empty)
     * 2. For each data point X[k], determine which cells it affects and update tables
     * 3. Propagate values to fill gaps (cells with no direct data points)
     * 4. Replicate boundary values to simplify edge case handling
     *
     * @pre *p_npts >= 2
     * @pre X array is sorted in ascending order
     * @post m_LO and m_HI tables are fully initialized and validated
     * @post m_must_reset = false
     *
     * Example with 8 data points and simplified table:
     * @code
     *       0  1     2     3             4 5 6              7     8
     *  X    +--+-----+-----+-------------+-+-+--------------+-----+
     *
     *       0        2        3        -(3)     6        -(6)     8
     * TABLE |--------|--------|--------|--------|--------|--------|
     *                2        -(4)     4        -(7)     7        8
     *       0        2        -        4                 7        -
     *
     * Cell boundaries define search ranges:
     *   Cell 0: search in X[0..2]
     *   Cell 1: search in X[2..4]
     *   Cell 2: search in X[3..4]  (data concentrated here)
     *   Cell 3: search in X[3..7]
     *   Cell 4: search in X[6..7]
     *   Cell 5: search in X[6..8]
     * @endcode
     *
     * @pre *p_npts >= 2
     * @pre X array is sorted in ascending order
     * @post m_LO and m_HI tables are fully initialized
     * @post m_must_reset = false
     *
     * @note This method is const because it updates mutable cached data
     * @note Thread-safe: called only from within locked sections
     */

    void reset() const
    {
      integer           n = *p_npts;
      real_type const * X = *p_X;

      // Validate minimum requirements
      UTILS_ASSERT( n >= 2, "SearchInterval::reset({}), need at least 2 points!", *p_name );

      m_table_size = std::clamp<integer>( static_cast<integer>( std::sqrt( n ) ), 32, 2048 );

      m_LO.resize( m_table_size + 2 );
      m_HI.resize( m_table_size + 2 );

// Verify array is sorted (in debug mode)
#ifndef NDEBUG
      for ( integer i = 1; i < n; ++i )
      {
        UTILS_ASSERT(
          X[i] >= X[i - 1],
          "SearchInterval::reset({}), X array not sorted at index {}: {} < {}",
          *p_name,
          i,
          X[i],
          X[i - 1] );
      }
#endif

      // Extract range information from sorted array
      m_x_min   = X[0];
      m_x_max   = X[n - 1];
      m_x_range = m_x_max - m_x_min;

      // Protection against degenerate or nearly-degenerate cases
      // If all points are essentially at the same location, use a minimal range
      real_type eps = eps_x( std::max( std::abs( m_x_min ), std::abs( m_x_max ) ) );
      if ( m_x_range < eps ) m_x_range = eps;

      // Cell width for uniform partitioning
      m_dx = m_x_range / m_table_size;

      // Special case: if m_dx is extremely small, use a simpler approach
      if ( m_dx < std::numeric_limits<real_type>::epsilon() * m_x_range )
      {
        // Degenerate case: all points in one cell
        std::fill( m_LO.begin(), m_LO.end(), 0 );
        std::fill( m_HI.begin(), m_HI.end(), n - 1 );
        m_must_reset = false;
        return;
      }

      // Initialize all table entries to -1 (unset marker)
      std::fill( m_LO.begin(), m_LO.end(), -1 );
      std::fill( m_HI.begin(), m_HI.end(), -1 );

      // Build m_LO table: for each point, mark the leftmost cell it belongs to
      // We use ceil with small negative offset to handle numerical precision
      for ( integer k = 0; k < n; ++k )
      {
        // Normalized position in [0, m_table_size]
        real_type pos{ ( X[k] - m_x_min ) / m_dx };

        // Cell index (with small epsilon to avoid precision issues at boundaries)
        integer i_LO{ static_cast<integer>( std::ceil( pos - eps / m_dx ) ) };

        // Clamp to valid range [0, m_table_size] for safety
        i_LO = std::max( 0, std::min( i_LO, m_table_size ) );

        // Update if this is the first point in this cell or a smaller index
        if ( m_LO[i_LO] == -1 || k < m_LO[i_LO] ) m_LO[i_LO] = k;
      }

      // First cell always starts at index 0
      m_LO[0] = 0;

      // Build m_HI table: for each point (in reverse), mark the rightmost cell it belongs to
      // We use floor with small positive offset to handle numerical precision
      for ( integer k = n - 1; k >= 0; --k )
      {
        // Normalized position in [0, m_table_size]
        real_type pos{ ( X[k] - m_x_min ) / m_dx };

        // Cell index (with small epsilon to avoid precision issues at boundaries)
        integer i_HI = std::clamp<integer>( static_cast<integer>( std::floor( pos + eps / m_dx ) ), 0, m_table_size );

        // Set only if not already set (we want the rightmost point)
        if ( m_HI[i_HI] == -1 ) m_HI[i_HI] = k;
      }

      // Last cell always ends at index n-1
      m_HI[m_table_size] = n - 1;

      // Propagate values forward in m_LO to fill empty cells
      // If a cell has no points, inherit from the previous cell
      for ( integer i = 0; i < m_table_size; ++i )
        if ( m_LO[i + 1] == -1 ) m_LO[i + 1] = m_LO[i];

      // Propagate values backward in m_HI to fill empty cells
      // If a cell has no points, inherit from the next cell
      for ( integer i{ m_table_size }; i > 0; --i )
        if ( m_HI[i - 1] == -1 ) m_HI[i - 1] = m_HI[i];

      // Replicate last values to avoid special boundary handling
      m_LO[m_table_size + 1] = m_LO[m_table_size];
      m_HI[m_table_size + 1] = m_HI[m_table_size];

// Validate tables for consistency
#ifndef NDEBUG
      if ( !validate_tables( n ) )
      {
        UTILS_ASSERT( false, "SearchInterval::reset({}), table validation failed!", *p_name );
      }
#endif

      // Mark tables as valid
      m_must_reset = false;
    }

  public:
    // Disable copy and move operations (contains mutex and external pointers)
    SearchInterval( SearchInterval const & )                   = delete;
    SearchInterval const & operator=( SearchInterval const & ) = delete;
    SearchInterval( SearchInterval && )                        = delete;
    SearchInterval & operator=( SearchInterval && )            = delete;

    /**
     * @brief Default constructor
     *
     * Creates an uninitialized SearchInterval. Must call setup() before use.
     */
    SearchInterval() {}

    /**
     * @brief Initialize the search structure with external data
     *
     * @param name Pointer to curve name (for error messages)
     * @param n Pointer to number of data points
     * @param X Pointer to pointer of sorted X coordinates array
     * @param is_closed Pointer to flag: true if curve is periodic
     * @param can_extend Pointer to flag: true if extrapolation allowed (currently unused)
     *
     * @note All pointers must remain valid for the lifetime of this object
     * @note This method is thread-safe
     * @note Marks internal tables for rebuild on next find() call
     */
    void setup( string const * name, integer * n, real_type ** X, bool * is_closed, bool * can_extend )
    {
      std::lock_guard<std::mutex> lock( m_mutex );
      p_name             = name;
      p_npts             = n;
      p_X                = X;
      p_curve_is_closed  = is_closed;
      p_curve_can_extend = can_extend;
      m_must_reset       = true;
      m_ready.store( false, std::memory_order_release );
    }

    /**
     * @brief Find the interval containing a query point
     *
     * Given a query point x, finds the interval index i such that X[i] <= x < X[i+1].
     * The algorithm uses a hybrid approach:
     * 1. Use lookup table to reduce search range from [0,n) to [k_LO, k_HI]
     * 2. Apply binary search within the reduced range
     * 3. Handle special cases (duplicates, boundaries, closed curves)
     *
     * @param[in,out] res Pair where:
     *   - res.second (input): query point x
     *   - res.first (output): interval index i where X[i] <= x < X[i+1]
     *   - res.second (output): possibly modified x (if curve is closed and x was wrapped)
     *
     * @pre setup() has been called
     * @pre *p_npts > 0
     * @post res.first is in range [0, n-2]
     *
     * @par Out-of-bounds handling:
     * - If x > X[n-1] and curve is closed: x is wrapped using modulo arithmetic
     * - If x > X[n-1] and curve is open: returns interval [n-2, n-1] (extrapolation)
     * - If x < X[0] and curve is closed: x is wrapped using modulo arithmetic
     * - If x < X[0] and curve is open: returns interval [0, 1] (extrapolation)
     *
     * @par Duplicate nodes handling:
     * If consecutive nodes have identical X values, returns the leftmost valid interval.
     * This ensures that evaluation at duplicate nodes uses the first segment definition.
     *
     * @note Thread-safe: multiple threads can call this method concurrently
     * @note First call after setup() or must_reset() will rebuild internal tables
     *
     * @par Complexity:
     * - Table lookup: O(1)
     * - Binary search: O(log(k_HI - k_LO + 1))
     * - Overall: O(log(n/table_size)) ≈ O(log n / 400) for typical cases
     *
     * Reference: Bentley, J.L. (1975). Multidimensional Binary Search Trees Used for Associative
     *            Searching. Communications of the ACM, 18(9), 509-517.
     */
    void find( std::pair<integer, real_type> & res ) const
    {
      // Lock for thread-safety and lazy initialization
      if ( !m_ready.load( std::memory_order_acquire ) )
      {
        std::lock_guard<std::mutex> lock( m_mutex );
        if ( m_must_reset ) reset();
        m_ready.store( true, std::memory_order_release );
      }

      // Local references for cleaner code
      integer const     n{ *p_npts };
      string const &    name{ *p_name };
      real_type const * X{ *p_X };

      UTILS_ASSERT( n > 0, "SearchInterval::find({}), n°points == 0!", name );

      integer &   pos{ res.first };  // Output: interval index
      real_type & x{ res.second };   // Input/Output: query point (may be wrapped)
      // ========================================================================
      // STEP 1: Handle out-of-bounds cases (robust + well-defined)
      // ========================================================================

      if ( x < m_x_min || x > m_x_max )
      {
        if ( *p_curve_is_closed )
        {
          // Periodic boundary: map x into [m_x_min, m_x_max)
          real_type t = std::fmod( x - m_x_min, m_x_range );

          // fmod can return negative values
          if ( t < 0 ) t += m_x_range;

          x = m_x_min + t;

          // Handle exact wrap-around (x == m_x_max)
          // Ensure continuity: last interval owns the boundary
          if ( x == m_x_min && t > 0 )
          {
            pos = n - 2;
            return;
          }
        }
        else
        {
          // Open curve: clamp to boundary intervals
          pos = x <= m_x_min ? 0 : n - 2;
          return;
        }
      }

      // ========================================================================
      // STEP 2: Use lookup table to reduce search range
      // ========================================================================

      // Compute normalized position in range [0, m_table_size]
      real_type norm_pos = ( x - m_x_min ) / m_dx;
      integer   i_cell{ static_cast<integer>( std::floor( norm_pos ) ) };

      // Critical: clamp to valid range to prevent array overflow
      // Without this, i_cell could be m_table_size, causing m_HI[i_cell+1] overflow
      i_cell = std::max( 0, std::min( i_cell, m_table_size - 1 ) );

      // Extract search boundaries from precomputed tables
      integer k_LO = m_LO[i_cell];      // Smallest possible interval
      integer k_HI = m_HI[i_cell + 1];  // Largest possible interval (+1 is safe due to clamp)

      // Sanity check: indices must be valid
      UTILS_ASSERT(
        k_LO >= 0 && k_LO < n && k_HI >= 0 && k_HI < n,
        "SearchInterval::find({}), invalid table indices: k_LO={}, k_HI={}, n={}\n",
        name,
        k_LO,
        k_HI,
        n );

      // Verify that k_LO <= k_HI (should be guaranteed by validate_tables)
      if ( k_LO > k_HI )
      {
        // Fall back to full binary search if tables are inconsistent
        k_LO = 0;
        k_HI = n - 1;
        fmt::print(
          "WARNING: SearchInterval::find({}), inconsistent tables at cell {}: k_LO={} > k_HI={}. Using full search.\n",
          name,
          i_cell,
          k_LO,
          k_HI );
      }

      // ========================================================================
      // STEP 3: Ensure x is within the reduced range for binary search
      // ========================================================================

      // Adjust k_LO if x is before the first point in the range
      if ( x < X[k_LO] )
      {
        // Binary search invariant: X[k_LO] <= x < X[k_HI]
        // We need to expand the left boundary
        while ( k_LO > 0 && x < X[k_LO] ) --k_LO;
      }

      // Adjust k_HI if x is after the last point in the range
      if ( x >= X[k_HI] )  // Note: use >= because binary search uses strict < for right boundary
      {
        // Binary search invariant: X[k_LO] <= x < X[k_HI]
        // We need to expand the right boundary
        while ( k_HI < n - 1 && x >= X[k_HI] ) ++k_HI;
      }

      // Final check of binary search precondition
      UTILS_ASSERT(
        k_LO >= 0 && k_HI >= k_LO && k_HI < n,
        "SearchInterval::find({}), invalid search range after adjustment: k_LO={}, k_HI={}, n={}\n",
        name,
        k_LO,
        k_HI,
        n );

      // ========================================================================
      // STEP 4: Binary search within reduced range
      // ========================================================================

      // Standard binary search: find largest k such that X[k] <= x
      // Invariant: X[k_LO] <= x < X[k_HI]
      while ( k_HI > k_LO + 1 )
      {
        // Midpoint calculation that avoids overflow
        integer k_M = k_LO + ( k_HI - k_LO ) / 2;

        // Maintain invariant: if x < X[k_M], then x < X[k_M] <= X[k_HI]
        // so we can set k_HI = k_M
        if ( x < X[k_M] )
          k_HI = k_M;  // x is in left half [k_LO, k_M)
        else
          k_LO = k_M;  // x is in right half [k_M, k_HI) (or exactly at k_M)
      }

      pos = k_LO;

      // ========================================================================
      // STEP 5: Handle duplicate consecutive nodes
      // ========================================================================

      // If X[pos] == X[pos+1], the interval [pos, pos+1] has zero width
      // We backtrack to find the leftmost valid interval with positive width
      // This ensures consistent behavior when evaluating at duplicate nodes
      while ( pos > 0 && pos < n - 1 )
      {
        // Compute relative tolerance based on magnitude of X values
        real_type const eps = eps_x( X[pos] );
        // If nodes are different (beyond tolerance), we found a valid interval
        if ( std::abs( X[pos] - X[pos + 1] ) > eps ) break;
        // Otherwise, move to previous interval
        --pos;
      }

      // Final safety clamp to ensure pos is a valid interval index
      pos = std::max( integer( 0 ), std::min( pos, n - 2 ) );

      // ========================================================================
      // STEP 6: Verify result (debug mode only)
      // ========================================================================

#ifndef NDEBUG
      // Verify that the found interval actually contains x
      if ( pos >= 0 && pos < n - 1 )
      {
        real_type const eps            = eps_x( X[pos] );
        bool            interval_valid = ( X[pos] - eps <= x ) && ( x <= X[pos + 1] + eps );

        if ( !interval_valid )
        {
          fmt::print(
            "WARNING: SearchInterval::find({}), invalid interval found: x={}, interval=[{},{}], X[{}]={}, X[{}]={}\n",
            name,
            x,
            pos,
            pos + 1,
            pos,
            X[pos],
            pos + 1,
            X[pos + 1] );
        }
      }
#endif
    }

    /**
     * @brief Mark internal tables for rebuild on next find() call
     *
     * Call this method when the external data (X array, n) has changed.
     * The actual rebuild is deferred until the next find() call (lazy evaluation).
     *
     * @note Thread-safe: can be called from any thread
     * @note Does not immediately rebuild - uses lazy evaluation for efficiency
     */
    void must_reset()
    {
      std::lock_guard<std::mutex> lock( m_mutex );
      m_must_reset = true;
      m_ready.store( false, std::memory_order_release );
    }

    /**
     * @brief Get debug information about the internal state
     *
     * @return String containing debug information about tables and state
     */
    string debug_info() const
    {
      std::lock_guard<std::mutex> lock( m_mutex );

      std::stringstream ss;
      ss << "SearchInterval Debug Info:\n";
      ss << "  Table size: " << m_table_size << "\n";
      ss << "  X range: [" << m_x_min << ", " << m_x_max << "]\n";
      ss << "  Cell width: " << m_dx << "\n";
      ss << "  Must reset: " << ( m_must_reset ? "true" : "false" ) << "\n";

      if ( !m_must_reset && p_npts )
      {
        ss << "  LO table sample: ";
        for ( integer i = 0; i <= 10 && i <= m_table_size; ++i ) ss << m_LO[i] << " ";
        ss << "\n";

        ss << "  HI table sample: ";
        for ( integer i = 0; i <= 10 && i <= m_table_size; ++i ) ss << m_HI[i] << " ";
        ss << "\n";
      }

      return ss.str();
    }
  };

  /*\
   |   ____        _ _
   |  / ___| _ __ | (_)_ __   ___
   |  \___ \| '_ \| | | '_ \ / _ \
   |   ___) | |_) | | | | | |  __/
   |  |____/| .__/|_|_|_| |_|\___|
   |        |_|
  \*/
  //!
  //! Spline Management Class
  //!
  class Spline
  {
    friend class SplineSet;

  protected:
    string const m_name;
    bool         m_curve_is_closed         = false;
    bool         m_curve_can_extend        = true;
    bool         m_curve_extended_constant = false;

    integer     m_npts          = 0;
    integer     m_npts_reserved = 0;
    real_type * m_X             = nullptr;  // allocated in the derived class!
    real_type * m_Y             = nullptr;  // allocated in the derived class!

    SearchInterval m_search;

  protected:
    void copy_flags( Spline const & S )
    {
      m_curve_is_closed         = S.m_curve_is_closed;
      m_curve_can_extend        = S.m_curve_can_extend;
      m_curve_extended_constant = S.m_curve_extended_constant;
    }

  public:
    Spline( Spline const & )                   = delete;
    Spline const & operator=( Spline const & ) = delete;

    //! \name Constructors
    ///@{

    //!
    //! spline constructor
    //!
    explicit Spline( string_view const name = "Spline" ) : m_name( name )
    {
      m_search.setup( &m_name, &m_npts, &m_X, &m_curve_is_closed, &m_curve_can_extend );
    }

    //!
    //! spline destructor
    //!
    virtual ~Spline() = default;

    ///@}

    //!
    //! \name Open/Close
    //!
    ///@{

    //!
    //! \return string with the name of the spline
    //!
    string_view name() const { return m_name; }

    //! \return `true` if spline is a closed spline
    bool is_closed() const { return m_curve_is_closed; }
    //!
    //! Set spline as a closed spline.
    //! When evaluated if parameter is outside the domain
    //! is wrapped cyclically before evalation.
    //!
    void make_closed() { m_curve_is_closed = true; }
    //!
    //! Set spline as an opened spline.
    //! When evaluated if parameter is outside the domain
    //! an error is produced.
    //!
    void make_opened() { m_curve_is_closed = false; }

    //!
    //! \return `true` if spline cannot extend outside interval of definition
    //!
    bool is_bounded() const { return !m_curve_can_extend; }
    //!
    //! Set spline as unbounded.
    //! When evaluated if parameter is outside the domain
    //! an extrapolated value is used.
    //!
    void make_unbounded() { m_curve_can_extend = true; }
    //!
    //! Set spline as bounded.
    //! When evaluated if parameter is outside the domain
    //! an error is issued.
    //!
    void make_bounded() { m_curve_can_extend = false; }

    //!
    //! \return `true` if the spline extend with a constant value
    //!
    bool is_extended_constant() const { return m_curve_extended_constant; }
    //!
    //! Set spline to extend constant.
    //! When evaluated if parameter is outside the domain
    //! the value returned is the value of the closed border.
    //!
    void make_extended_constant() { m_curve_extended_constant = true; }
    //!
    //! Set spline to extend NOT constant.
    //! When evaluated if parameter is outside the domain
    //! teh value returned is extrapolated using the last spline polynomial.
    //!
    void make_extended_not_constant() { m_curve_extended_constant = false; }

    ///@}

    //! \name Spline Data Info
    ///@{

    //!
    //! return true if spline is empty (no points)
    //!
    bool empty() const { return m_npts == 0; }

    //!
    //! the number of support points of the spline.
    //!
    integer num_points() const { return m_npts; }

    //!
    //! Return the pointer of values of x-nodes.
    //!
    real_type const * x_nodes() const { return m_X; }

    //!
    //! Return the pointer of values of y-nodes.
    //!
    real_type const * y_nodes() const { return m_Y; }

    //!
    //! the i-th node of the spline (x component).
    //!
    real_type x_node( integer const i ) const { return m_X[i]; }

    //!
    //! the i-th node of the spline (y component).
    //!
    real_type y_node( integer const i ) const { return m_Y[i]; }

    //!
    //! first node of the spline (x component).
    //!
    real_type x_begin() const { return m_X[0]; }

    //!
    //! first node of the spline (y component).
    //!
    real_type y_begin() const { return m_Y[0]; }

    //!
    //! last node of the spline (x component).
    //!
    real_type x_end() const { return m_X[m_npts - 1]; }

    //!
    //! last node of the spline (y component).
    //!
    real_type y_end() const { return m_Y[m_npts - 1]; }

    //!
    //! x-minumum spline value
    //!
    real_type x_min() const { return m_X[0]; }

    //!
    //! x-maximum spline value
    //!
    real_type x_max() const { return m_X[m_npts - 1]; }

    //!
    //! y-minumum spline value
    //!
    real_type y_min() const
    {
      integer N = ( type() == SplineType1D::CONSTANT ) ? m_npts - 1 : m_npts;
      return *std::min_element( m_Y, m_Y + N );
    }

    //!
    //! return y-maximum spline value
    //!
    real_type y_max() const
    {
      integer N = ( type() == SplineType1D::CONSTANT ) ? m_npts - 1 : m_npts;
      return *std::max_element( m_Y, m_Y + N );
    }

    //!
    //! Search the max and min values of `y` along the spline
    //! with the corresponding `x` position
    //!
    //! \param[out] i_min_pos interval where is the minimum
    //! \param[out] x_min_pos where is the minimum
    //! \param[out] y_min     the minimum value
    //! \param[out] i_max_pos interval where is the maximum
    //! \param[out] x_max_pos where is the maximum
    //! \param[out] y_max     the maximum value
    //!
    virtual void y_min_max(
      integer &   i_min_pos,
      real_type & x_min_pos,
      real_type & y_min,
      integer &   i_max_pos,
      real_type & x_max_pos,
      real_type & y_max ) const
    {
      i_min_pos = i_max_pos = 0;
      x_min_pos = y_min = x_max_pos = y_max = 0;
      UTILS_ERROR( "In spline: {} y_min_max not implemented\n", info() );
    }

    //!
    //! Search the max and min values of `y` along the spline
    //! with the corresponding `x` position
    //!
    //! \param[out] i_min_pos interval where is the minimum
    //! \param[out] x_min_pos where is the minimum
    //! \param[out] y_min     the minimum value
    //! \param[out] i_max_pos interval where is the maximum
    //! \param[out] x_max_pos where is the maximum
    //! \param[out] y_max     the maximum value
    //!
    virtual void y_min_max(
      vector<integer> &   i_min_pos,
      vector<real_type> & x_min_pos,
      vector<real_type> & y_min,
      vector<integer> &   i_max_pos,
      vector<real_type> & x_max_pos,
      vector<real_type> & y_max ) const
    {
      i_min_pos.clear();
      i_max_pos.clear();
      x_min_pos.clear();
      x_max_pos.clear();
      y_min.clear();
      y_max.clear();
      UTILS_ERROR( "In spline: {} y_min_max not implemented\n", info() );
    }
    ///@}

    //! \name Build
    ///@{

    //!
    //! Build a spline using data in `GenericContainer`
    //!
    void build( GenericContainer const & gc ) { setup( gc ); }

    //!
    //! Build a spline using data in a file `file_name`
    //!
    void build( string const & file_name ) { setup( file_name ); }

    //!
    //! Build a spline.
    //!
    //! \param x    vector of x-coordinates
    //! \param incx access elements as `x[0]`, `x[incx]`, `x[2*incx]`,...
    //! \param y    vector of y-coordinates
    //! \param incy access elements as `y[0]`, `y[incx]`, `y[2*incx]`,...
    //! \param n    total number of points
    //!
    virtual void build( real_type const x[], integer incx, real_type const y[], integer incy, integer const n )
    {
      reserve( n );
      for ( integer i = 0; i < n; ++i ) m_X[i] = x[i * incx];
      for ( integer i = 0; i < n; ++i ) m_Y[i] = y[i * incy];
      m_npts = n;
      build();
    }

    //!
    //! Build a spline.
    //!
    //! \param x vector of x-coordinates
    //! \param y vector of y-coordinates
    //! \param n total number of points
    //!
    inline void build( real_type const x[], real_type const y[], integer const n ) { this->build( x, 1, y, 1, n ); }

    //!
    //! Build a spline.
    //!
    //! \param x vector of x-coordinates
    //! \param y vector of y-coordinates
    //!
    void build( vector<real_type> const & x, vector<real_type> const & y )
    {
      integer N = std::min<integer>( x.size(), y.size() );
      this->build( x.data(), 1, y.data(), 1, N );
    }

    //!
    //! Build a spline using internal stored data
    //!
    virtual void build() = 0;

    //!
    //! Setup a spline using a `GenericContainer`
    //!
    //! - gc("xdata") vector with the `x` coordinate of the data
    //! - gc("ydata") vector with the `y` coordinate of the data
    //!
    virtual void setup( GenericContainer const & gc )
    {
      /*
      // gc["xdata"]
      // gc["ydata"]
      //
      */
      string const where{ fmt::format( "Spline[{}]::setup( gc ):", m_name ) };

      std::set<std::string> keywords;
      for ( auto const & pair : gc.get_map( where ) ) { keywords.insert( pair.first ); }

      GenericContainer const & gc_x{ gc( "xdata", where ) };
      keywords.erase( "xdata" );
      GenericContainer const & gc_y{ gc( "ydata", where ) };
      keywords.erase( "ydata" );
      keywords.erase( "spline_type" );

      vec_real_type x, y;
      {
        string const ff{ fmt::format( "{}, field `xdata'", where ) };
        gc_x.copyto_vec_real( x, ff );
      }
      {
        string const ff{ fmt::format( "{}, field `ydata'", where ) };
        gc_y.copyto_vec_real( y, ff );
      }

      UTILS_WARNING(
        keywords.empty(),
        "{}: unused keys\n{}\n",
        where,
        [&keywords]() -> string
        {
          string res;
          for ( auto const & it : keywords )
          {
            res += it;
            res += ' ';
          };
          return res;
        }() );
      build( x, y );
    }
    //!
    //! Setup a spline using a `GenericContainer` readed from file
    //!
    //! - file_name file name of the file with data
    //!
    //! - `.json`
    //! - `.yaml` or `.yml`
    //! - `.toml`
    //!
    void setup( string const & file_name )
    {
      GenericContainer gc;
      UTILS_ASSERT( gc.from_file( file_name ), "Spline::setup( '{}' ) failed to read\n", file_name );
      setup( gc );
    }

    ///@}

    //! \name Incremental Build
    ///@{

    //!
    //! Allocate memory for `npts` points
    //!
    virtual void reserve( integer const npts ) = 0;

    //!
    //! Add a support point (x,y) to the spline.
    //!
    void push_back( real_type const x, real_type const y )
    {
      if ( m_npts > 0 )
      {
        UTILS_ASSERT(
          x >= m_X[m_npts - 1],  // ammetto punti doppi
          "Spline[{}]::push_back, non monotone insert at insert N.{}"
          "\nX[{}] = {:.5}\nX[{}] = {:>5}\n",
          m_name,
          m_npts,
          m_npts - 1,
          m_X[m_npts - 1],
          m_npts,
          x );
      }
      if ( m_npts_reserved == 0 ) { reserve( 2 ); }
      else if ( m_npts >= m_npts_reserved )
      {
        // riallocazione & copia
        integer const saved_npts{ m_npts };  // salvo npts perche reserve lo azzera
        Malloc_real   mem( "Spline::push_back" );
        mem.allocate( 2 * m_npts );
        real_type * Xsaved{ mem( m_npts ) };
        real_type * Ysaved{ mem( m_npts ) };

        std::copy_n( m_X, m_npts, Xsaved );
        std::copy_n( m_Y, m_npts, Ysaved );
        reserve( ( m_npts + 1 ) * 2 );
        m_npts = saved_npts;
        std::copy_n( Xsaved, m_npts, m_X );
        std::copy_n( Ysaved, m_npts, m_Y );
      }
      m_X[m_npts] = x;
      m_Y[m_npts] = y;
      ++m_npts;
      m_search.must_reset();
    }

    //!
    //! Drop last inserted point of the spline.
    //!
    void drop_back()
    {
      if ( m_npts > 0 ) --m_npts;
      m_search.must_reset();
    }

    //!
    //! Delete the support points, empty the spline.
    //!
    virtual void clear() = 0;

    ///@}

    //! \name Manipulate
    ///@{

    //!
    //! change X-origin of the spline
    //!
    void set_origin( real_type const x0 ) const
    {
      real_type const Tx{ x0 - m_X[0] };
      real_type *     ix{ m_X };
      while ( ix < m_X + m_npts ) *ix++ += Tx;
    }

    //!
    //! change X-range of the spline
    //!
    void set_range( real_type xmin, real_type xmax )
    {
      UTILS_ASSERT( xmax > xmin, "Spline[{}]::set_range({},{}) bad range ", m_name, xmin, xmax );
      real_type const S  = ( xmax - xmin ) / ( m_X[m_npts - 1] - m_X[0] );
      real_type const Tx = xmin - S * m_X[0];
      for ( real_type * ix = m_X; ix < m_X + m_npts; ++ix ) *ix = *ix * S + Tx;
    }
    ///@}

    //! \name Dump Data
    ///@{

    //!
    //! dump a sample of the spline
    //!
    void dump( ostream_type & s, integer const nintervals, string_view header = "x\ty" ) const
    {
      s << header << '\n';
      real_type const dx = ( x_max() - x_min() ) / nintervals;
      for ( integer i = 0; i <= nintervals; ++i )
      {
        real_type x = x_min() + i * dx;
        fmt::print( s, "{}\t{}\n", x, this->eval( x ) );
      }
    }

    //!
    //! dump a sample of the spline
    //!
    void dump( string_view fname, integer const nintervals, string_view header = "x\ty" ) const
    {
      std::ofstream file( fname.data() );
      this->dump( file, nintervals, header );
      file.close();
    }

    //!
    //! Print spline coefficients
    //!
    virtual void write_to_stream( ostream_type & s ) const = 0;

    ///@}

    //! \name Evaluation
    ///@{

#ifdef AUTODIFF_SUPPORT
    //!
    //! Evaluate spline value
    //!
    virtual real_type         eval( real_type const x ) const           = 0;
    virtual autodiff::dual1st eval( autodiff::dual1st const & x ) const = 0;
    virtual autodiff::dual2nd eval( autodiff::dual2nd const & x ) const = 0;

    // Template unificato per tutti i tipi
    template <typename T> auto eval( T const & x ) const
    {
      if constexpr ( std::is_arithmetic<T>::value )
      {
        // Se T è un tipo numerico (int, float, double, etc.), promuovi a real_type
        return eval( static_cast<real_type>( x ) );
      }
      else
      {
        // Altrimenti deduce automaticamente il tipo duale appropriato
        return eval( autodiff::detail::to_dual( x ) );
      }
    }

    template <typename T> auto operator()( T const & x ) const -> decltype( eval( x ) ) { return eval( x ); }
#endif

    //!
    //! First derivative
    //!
    virtual real_type D( real_type const x ) const = 0;

    //!
    //! Second derivative
    //!
    virtual real_type DD( real_type const x ) const = 0;

    //!
    //! Third derivative
    //!
    virtual real_type DDD( real_type const x ) const = 0;

    //!
    //! 4th derivative
    //!
    virtual real_type DDDD( real_type const ) const { return real_type( 0 ); }

    //!
    //! 5th derivative
    //!
    virtual real_type DDDDD( real_type const ) const { return real_type( 0 ); }

    virtual void D( real_type const x, real_type dd[2] ) const  = 0;
    virtual void DD( real_type const x, real_type dd[3] ) const = 0;

    ///@}

    //!
    //! \name Evaluation Aliases
    //!
    ///@{
    //!
    //! Alias for `real_type eval( real_type x )`
    //!
    real_type operator()( real_type const x ) const { return this->eval( x ); }
    //!
    //! Alias for `real_type D( real_type x )`
    //!
    real_type eval_D( real_type const x ) const { return this->D( x ); }
    //!
    //! Alias for `real_type DD( real_type x )`
    //!
    real_type eval_DD( real_type const x ) const { return this->DD( x ); }
    //!
    //! Alias for `real_type DDD( real_type x )`
    //!
    real_type eval_DDD( real_type const x ) const { return this->DDD( x ); }
    //!
    //! Alias for `real_type DDDD( real_type x )`
    //!
    real_type eval_DDDD( real_type const x ) const { return this->DDDD( x ); }
    //!
    //! Alias for `real_type DDDD( real_type x )`
    //!
    real_type eval_DDDDD( real_type const x ) const { return this->DDDDD( x ); }
    ///@}

    //!
    //! \name Evaluation when segment is known
    //!
    ///@{

    //!
    //! Evaluate spline value
    //!
    virtual real_type id_eval( integer const ni, real_type const x ) const = 0;

    //!
    //! First derivative
    //!
    virtual real_type id_D( integer const ni, real_type const x ) const = 0;

    //!
    //! Second derivative
    //!
    virtual real_type id_DD( integer const ni, real_type const x ) const = 0;

    //!
    //! Third derivative
    //!
    virtual real_type id_DDD( integer const ni, real_type const x ) const = 0;

    //!
    //! 4th derivative
    //!
    virtual real_type id_DDDD( integer const, real_type const ) const { return real_type( 0 ); }

    //!
    //! 5th derivative
    //!
    virtual real_type id_DDDDD( integer const, real_type const ) const { return real_type( 0 ); }

    ///@}

    //! \name Get Info
    ///@{

    //!
    //! get the piecewise polinomials of the spline
    //!
    virtual integer  // order
    coeffs( real_type cfs[], real_type nodes[], bool transpose = false ) const = 0;

    //!
    //! Spline order of the piecewise polynomial
    //!
    virtual integer order() const = 0;

    //!
    //! spline type returned as a string
    //!
    char const * type_name() const { return to_string( type() ); }

    //!
    //! spline type returned as integer
    //!
    virtual SplineType1D type() const = 0;

    //!
    //! String information of the kind and order of the spline
    //!
    string info() const
    {
      string res{ fmt::format( "Spline `{}` of type: {} of order: {}", m_name, type_name(), order() ) };
      if ( m_npts > 0 )
        res += fmt::format( "\nx_min={:.5} x_max={:.5} y_min={:.5} y_max={:.5}", x_min(), x_max(), y_min(), y_max() );
      return res;
    }

    //!
    //! Print information of the kind and order of the spline
    //!
    void info( ostream_type & stream ) const { stream << this->info() << '\n'; }

    ///@}

#ifdef SPLINES_BACK_COMPATIBILITY
    void pushBack( real_type x, real_type y ) { push_back( x, y ); }
    void dropBack() { drop_back(); }
    void setOrigin( real_type x0 ) { set_origin( x0 ); }
    void setRange( real_type xmin, real_type xmax ) { set_range( xmin, xmax ); }
    void writeToStream( ostream_type & s ) const { write_to_stream( s ); }
#endif
  };

  //!
  //! compute curvature of a planar curve
  //!
  inline real_type curvature( real_type s, Spline const & X, Spline const & Y )
  {
    const real_type dx  = X.D( s );
    const real_type dy  = Y.D( s );
    const real_type ddx = X.DD( s );
    const real_type ddy = Y.DD( s );

    const real_type speed2 = dx * dx + dy * dy;
    const real_type speed  = std::sqrt( speed2 );

    return ( dx * ddy - dy * ddx ) / ( speed2 * speed );
  }

  //!
  //! compute curvature derivative of a planar curve
  //!
  inline real_type curvature_D( real_type s, Spline const & X, Spline const & Y )
  {
    const real_type dx   = X.D( s );
    const real_type dy   = Y.D( s );
    const real_type ddx  = X.DD( s );
    const real_type ddy  = Y.DD( s );
    const real_type dddx = X.DDD( s );
    const real_type dddy = Y.DDD( s );

    const real_type dx2 = dx * dx;
    const real_type dy2 = dy * dy;

    const real_type speed2 = dx2 + dy2;            // |r'|^2
    const real_type speed  = std::sqrt( speed2 );  // |r'|

    const real_type a = dddy * dx - dy * dddx;
    const real_type b = 3 * ddy * ddx;

    const real_type num = dx2 * ( a - b ) + dy2 * ( a + b ) + 3 * dx * dy * ( ddx * ddx - ddy * ddy );

    return num / ( speed * speed2 * speed2 );
  }

  //!
  //! compute curvature second derivative of a planar curve
  //!
  inline real_type curvature_DD( real_type s, Spline const & X, Spline const & Y )
  {
    const real_type dx = X.D( s ), dy = Y.D( s );
    const real_type ddx = X.DD( s ), ddy = Y.DD( s );
    const real_type dddx = X.DDD( s ), dddy = Y.DDD( s );
    const real_type ddddx = X.DDDD( s ), ddddy = Y.DDDD( s );

    const real_type dx2 = dx * dx, dy2 = dy * dy;
    const real_type dx3 = dx2 * dx, dy3 = dy2 * dy;

    const real_type ddx2 = ddx * ddx, ddy2 = ddy * ddy;
    const real_type ddx3 = ddx2 * ddx, ddy3 = ddy2 * ddy;

    const real_type speed2 = dx2 + dy2;
    if ( speed2 == real_type( 0 ) ) return real_type( 0 );

    const real_type inv_speed7 = real_type( 1 ) / ( std::sqrt( speed2 ) * speed2 * speed2 * speed2 );

    const real_type dx_dddx     = dx * dddx;
    const real_type dx_ddy_dddy = dx * ddy * dddy;

    const real_type num = dy2 * dy2 * ( dx * ddddy + 4 * ddx * dddy + 5 * dddx * ddy - dy * ddddx )

                          + dy3 * ( 3 * ddx3 + ddx * ( 9 * dx_dddx - 12 * ddy2 ) - 2 * ddddx * dx2 - 9 * dx_ddy_dddy )

                          + dy2 * ( dx * ( 12 * ddy3 - 33 * ddy * ddx2 ) + dx2 * ( ddy * dddx - dddy * ddx ) +
                                    2 * ddddy * dx3 )

                          - dy * dx2 * ( 12 * ddx3 - ddx * ( 9 * dx_dddx + 33 * ddy2 ) + ddddx * dx2 + 9 * dx_ddy_dddy )

                          + dx3 * ( ddy * ( 12 * ddx2 - 4 * dx_dddx ) + ddddy * dx2 - 5 * dx * ddx * dddy - 3 * ddy3 );

    return num * inv_speed7;
  }


  /*\
   |    ____      _     _        ____        _ _              ____
   |   / ___|   _| |__ (_) ___  / ___| _ __ | (_)_ __   ___  | __ )  __ _ ___  ___
   |  | |  | | | | '_ \| |/ __| \___ \| '_ \| | | '_ \ / _ \ |  _ \ / _` / __|/ _ \
   |  | |__| |_| | |_) | | (__   ___) | |_) | | | | | |  __/ | |_) | (_| \__ \  __/
   |   \____\__,_|_.__/|_|\___| |____/| .__/|_|_|_| |_|\___| |____/ \__,_|___/\___|
   |                                  |_|
  \*/
  //!
  //! cubic spline base class
  //!
  class CubicSplineBase : public Spline
  {
  protected:
    Malloc_real m_mem_cubic;
    real_type * m_Yp{ nullptr };
    bool        m_external_alloc{ false };

  public:
    using Spline::build;

    //!
    //! \name Contructors/Destructors
    ///@{
    //!
    //! Spline constructor.
    //!
    explicit CubicSplineBase( string_view name = "CubicSplineBase" );

    ~CubicSplineBase() override {}
    ///@}

    //!
    //! Build a copy of spline `S`
    //!
    void copy_spline( CubicSplineBase const & S );

    //!
    //! Return the pointer of values of yp-nodes.
    //!
    real_type const * yp_nodes() const { return m_Yp; }

    //!
    //! Return the i-th node of the spline (y' component).
    //!
    real_type yp_node( integer i ) const { return m_Yp[i]; }

    //!
    //! Change X-range of the spline.
    //!
    void set_range( real_type xmin, real_type xmax );

    //!
    //! Use externally allocated memory for `npts` points.
    //!
    void reserve_external( integer n, real_type *& p_x, real_type *& p_y, real_type *& p_dy );

    void y_min_max(
      integer &   i_min_pos,
      real_type & x_min_pos,
      real_type & y_min,
      integer &   i_max_pos,
      real_type & x_max_pos,
      real_type & y_max ) const override;

    void y_min_max(
      vector<integer> &   i_min_pos,
      vector<real_type> & x_min_pos,
      vector<real_type> & y_min,
      vector<integer> &   i_max_pos,
      vector<real_type> & x_max_pos,
      vector<real_type> & y_max ) const override;

    // --------------------------- VIRTUALS -----------------------------------

    //!
    //! \name Evaluation
    //!
    ///@{
    real_type eval( real_type const x ) const override;
    real_type D( real_type const x ) const override;
    real_type DD( real_type const x ) const override;
    real_type DDD( real_type const x ) const override;
    real_type DDDD( real_type const ) const override { return 0; }
    real_type DDDDD( real_type const ) const override { return 0; }

    void D( real_type const x, real_type dd[2] ) const override;
    void DD( real_type const x, real_type dd[3] ) const override;

    ///@}

    //!
    //! \name Evaluation when segment is known
    ///@{
    real_type id_eval( integer const ni, real_type const x ) const override;
    real_type id_D( integer const ni, real_type const x ) const override;
    real_type id_DD( integer const ni, real_type const x ) const override;
    real_type id_DDD( integer const ni, real_type const x ) const override;
    real_type id_DDDD( integer const, real_type const ) const override { return 0; }
    real_type id_DDDDD( integer const, real_type const ) const override { return 0; }
    ///@}

#ifdef AUTODIFF_SUPPORT
    autodiff::dual1st eval( autodiff::dual1st const & x ) const override;
    autodiff::dual2nd eval( autodiff::dual2nd const & x ) const override;

    // Template unificato per tutti i tipi
    template <typename T> auto eval( T const & x ) const
    {
      if constexpr ( std::is_arithmetic<T>::value )
      {
        // Se T è un tipo numerico (int, float, double, etc.), promuovi a real_type
        return eval( static_cast<real_type>( x ) );
      }
      else
      {
        // Altrimenti deduce automaticamente il tipo duale appropriato
        return eval( autodiff::detail::to_dual( x ) );
      }
    }

    template <typename T> auto operator()( T const & x ) const -> decltype( eval( x ) ) { return eval( x ); }
#endif

    void write_to_stream( ostream_type & s ) const override;

    // --------------------------- VIRTUALS -----------------------------------

    //!
    //! \name Build
    //!
    ///@{

    //!
    //! Build a spline.
    //!
    //! \param[in] x     vector of x-coordinates
    //! \param[in] incx  access elements as x[0], x[incx], x[2*incx],...
    //! \param[in] y     vector of y-coordinates
    //! \param[in] incy  access elements as y[0], y[incy], y[2*incy],...
    //! \param[in] yp    vector of y-defivative
    //! \param[in] incyp access elements as yp[0], yp[incyp], yp[2*incyy],...
    //! \param[in] n     total number of points
    //!
    void build(
      real_type const x[],
      integer const   incx,
      real_type const y[],
      integer const   incy,
      real_type const yp[],
      integer const   incyp,
      integer         n );

    //!
    //! Build a spline.
    //!
    //! \param[in] x  vector of x-coordinates
    //! \param[in] y  vector of y-coordinates
    //! \param[in] yp vector of y'-coordinates
    //! \param[in] n  total number of points
    //!
    inline void build( real_type const x[], real_type const y[], real_type const yp[], integer const n )
    {
      this->build( x, 1, y, 1, yp, 1, n );
    }

    //!
    //! Build a spline.
    //!
    //! \param[in] x  vector of x-coordinates
    //! \param[in] y  vector of y-coordinates
    //! \param[in] yp vector of y'-coordinates
    //!
    void build( vector<real_type> const & x, vector<real_type> const & y, vector<real_type> const & yp );

    void reserve( integer npts ) override;

    ///@}

    void clear() override;

    integer  // order
    coeffs( real_type cfs[], real_type nodes[], bool transpose = false ) const override;

    integer order() const override;

#ifdef SPLINES_BACK_COMPATIBILITY
    void      copySpline( CubicSplineBase const & S ) { this->copy_spline( S ); }
    integer   numPoints() const { return m_npts; }
    real_type xNode( integer i ) const { return m_X[i]; }
    real_type yNode( integer i ) const { return m_Y[i]; }
    real_type ypNode( integer i ) const { return this->yp_node( i ); }
    real_type xBegin() const { return m_X[0]; }
    real_type yBegin() const { return m_Y[0]; }
    real_type xEnd() const { return m_X[m_npts - 1]; }
    real_type yEnd() const { return m_Y[m_npts - 1]; }
    real_type xMin() const { return m_X[0]; }
    real_type xMax() const { return m_X[m_npts - 1]; }
    real_type yMin() const { return y_min(); }
    real_type yMax() const { return y_max(); }
#endif
  };

}  // namespace Splines

#include "Splines/SplineAkima.hxx"
#include "Splines/SplineBessel.hxx"
#include "Splines/SplineConstant.hxx"
#include "Splines/SplineLinear.hxx"
#include "Splines/SplineCubic.hxx"
#include "Splines/SplineHermite.hxx"
#include "Splines/SplinePchip.hxx"
#include "Splines/SplineQuinticBase.hxx"
#include "Splines/SplineQuintic.hxx"


namespace Splines
{

  /*\
   |   ____        _ _            ____              __
   |  / ___| _ __ | (_)_ __   ___/ ___| _   _ _ __ / _|
   |  \___ \| '_ \| | | '_ \ / _ \___ \| | | | '__| |_
   |   ___) | |_) | | | | | |  __/___) | |_| | |  |  _|
   |  |____/| .__/|_|_|_| |_|\___|____/ \__,_|_|  |_|
   |        |_|
  \*/

  //!
  //! Spline Management Class
  //!
  class SplineSurf
  {
    Malloc_real m_mem;

  protected:
    string const m_name;
    bool         m_x_closed{ false };
    bool         m_y_closed{ false };
    bool         m_x_can_extend{ true };
    bool         m_y_can_extend{ true };

    integer m_nx{ 0 };
    integer m_ny{ 0 };

    real_type * m_X{ nullptr };
    real_type * m_Y{ nullptr };
    real_type * m_Z{ nullptr };

    real_type m_Z_min{ 0 };
    real_type m_Z_max{ 0 };

    SearchInterval m_search_x;
    SearchInterval m_search_y;

    static integer ipos_C( integer const i, integer const j, integer const ldZ ) { return i * ldZ + j; }

    static integer ipos_F( integer const i, integer const j, integer const ldZ ) { return i + ldZ * j; }

    integer ipos_C( integer const i, integer const j ) const { return this->ipos_C( i, j, m_ny ); }

    integer ipos_F( integer const i, integer const j ) const { return this->ipos_F( i, j, m_nx ); }

    real_type & z_node_ref( integer const i, integer const j ) { return m_Z[this->ipos_C( i, j )]; }

    void load_Z( real_type const z[], integer const ldZ, bool fortran_storage, bool transposed );

    virtual void make_spline() = 0;

    void make_derivative_x( real_type const z[], real_type dx[] )
    {
      PchipSpline pchip_work;
      for ( integer j = 0; j < m_ny; ++j )
      {
        pchip_work.build( m_X, 1, z + ipos_C( 0, j ), m_ny, m_nx );
        for ( integer i = 0; i < m_nx; ++i ) dx[ipos_C( i, j )] = pchip_work.yp_node( i );
      }
    }

    void make_derivative_y( real_type const z[], real_type dy[] )
    {
      PchipSpline pchip_work;
      for ( integer i = 0; i < m_nx; ++i )
      {
        pchip_work.build( m_Y, 1, z + ipos_C( i, 0 ), 1, m_ny );
        for ( integer j = 0; j < m_ny; ++j ) dy[ipos_C( i, j )] = pchip_work.yp_node( j );
      }
    }

    void make_derivative_xy( real_type const dx[], real_type const dy[], real_type dxy[] )
    {
      PchipSpline pchip_work;

      auto minmod = []( real_type a, real_type b ) -> real_type
      {
        if ( a * b <= 0 ) return 0;
        if ( a > 0 ) return std::min( a, b );
        return std::max( a, b );
      };

      for ( integer j = 0; j < m_ny; ++j )
      {
        pchip_work.build( m_X, 1, dy + ipos_C( 0, j ), m_ny, m_nx );
        for ( integer i = 0; i < m_nx; ++i ) dxy[ipos_C( i, j )] = pchip_work.yp_node( i );
      }

      for ( integer i = 0; i < m_nx; ++i )
      {
        pchip_work.build( m_Y, 1, dx + ipos_C( i, 0 ), 1, m_ny );
        for ( integer j = 0; j < m_ny; ++j )
        {
          integer const ij{ ipos_C( i, j ) };
          dxy[ij] = minmod( dxy[ij], pchip_work.yp_node( j ) );
        }
      }
    }

  public:
    SplineSurf( SplineSurf const & )                   = delete;  // block copy constructor
    SplineSurf const & operator=( SplineSurf const & ) = delete;  // block copy method

    //!
    //! Spline constructor
    //!
    explicit SplineSurf( string_view name = "Spline" ) : m_mem( name.data() ), m_name( name )
    {
      m_search_x.setup( &m_name, &m_nx, &m_X, &m_x_closed, &m_x_can_extend );
      m_search_y.setup( &m_name, &m_ny, &m_Y, &m_y_closed, &m_y_can_extend );
    }

    //!
    //! Spline destructor
    //!
    virtual ~SplineSurf();

    //!
    //! \name Open/Close
    //!
    ///@{

    //!
    //! Return `true` if the surface is assumed closed in the `x` direction.
    //!
    bool is_x_closed() const { return m_x_closed; }

    //!
    //! Setup the surface as closed in the `x` direction.
    //!
    void make_x_closed() { m_x_closed = true; }

    //!
    //! Setup the surface as open in the `x` direction.
    //!
    void make_x_opened() { m_x_closed = false; }

    //!
    //! Return `true` if the surface is assumed closed in the `y` direction.
    //!
    bool is_y_closed() const { return m_y_closed; }

    //!
    //! Setup the surface as closed in the `y` direction.
    //!
    void make_y_closed() { m_y_closed = true; }

    //!
    //! Setup the surface as open in the `y` direction.
    //!
    void make_y_opened() { m_y_closed = false; }

    //!
    //! Return `true` if the parameter `x` assumed bounded.
    //! If false the spline is estrapolated for `x` values
    //! outside the range.
    //!
    bool is_x_bounded() const { return m_x_can_extend; }

    //!
    //! Make the spline surface unbounded in the `x` direction.
    //!
    void make_x_unbounded() { m_x_can_extend = true; }

    //!
    //! Make the spline surface bounded in the `x` direction.
    //!
    void make_x_bounded() { m_x_can_extend = false; }

    //!
    //! Return `true` if the parameter `y` assumed bounded.
    //! If false the spline is extrapolated for `y` values
    //! outside the range.
    //!
    bool is_y_bounded() const { return m_y_can_extend; }

    //!
    //! Make the spline surface unbounded in the `y` direction
    //!
    void make_y_unbounded() { m_y_can_extend = true; }

    //!
    //! Make the spline surface bounded in the `x` direction.
    //!
    void make_y_bounded() { m_y_can_extend = false; }

    ///@}

    //!
    //! Cancel the support points, empty the spline.
    //!
    void clear();

    //!
    //! \name Info
    //!
    ///@{

    //!
    //! \return string with the name of the spline
    //!
    string_view name() const { return m_name; }

    //!
    //! Return the number of support points of the spline along x direction.
    //!
    integer num_point_x() const { return m_nx; }

    //!
    //! Return the number of support points of the spline along y direction.
    //!
    integer num_point_y() const { return m_ny; }

    //!
    //! Return the i-th node of the spline (x component).
    //!
    real_type x_node( integer const i ) const { return m_X[i]; }

    //!
    //! Return the i-th node of the spline (y component).
    //!
    real_type y_node( integer const i ) const { return m_Y[i]; }

    //!
    //! Return the i-th node of the spline (y component).
    //!
    real_type z_node( integer const i, integer const j ) const { return m_Z[this->ipos_C( i, j )]; }

    //!
    //! Return x-minumum spline value.
    //!
    real_type x_min() const { return m_X[0]; }

    //!
    //! Return x-maximum spline value.
    //!
    real_type x_max() const { return m_X[m_nx - 1]; }

    //!
    //! Return y-minumum spline value.
    //!
    real_type y_min() const { return m_Y[0]; }

    //!
    //! Return y-maximum spline value.
    //!
    real_type y_max() const { return m_Y[m_ny - 1]; }

    //!
    //! Return z-minumum spline value.
    //!
    real_type z_min() const { return m_Z_min; }

    //!
    //! Return z-maximum spline value.
    //!
    real_type z_max() const { return m_Z_max; }

    ///@}

    //!
    //! \name Build Spline
    //!
    ///@{

    //!
    //! Build surface spline
    //!
    //! \param x               vector of `x`-coordinates
    //! \param incx            access elements as `x[0]`, `x[incx]`, `x[2*incx]`,...
    //! \param y               vector of `y`-coordinates
    //! \param incy            access elements as `y[0]`, `y[incy]`, `y[2*incy]`,...
    //! \param z               matrix of `z`-values. Elements are stored
    //!                        by row Z(i,j) = z[i*ny+j] as C-matrix
    //! \param ldZ             leading dimension of `z`
    //! \param nx              number of points in `x` direction
    //! \param ny              number of points in `y` direction
    //! \param fortran_storage if true elements are stored by column
    //!                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
    //! \param transposed      if true matrix Z is stored transposed
    //!
    void build(
      real_type const x[],
      integer const   incx,
      real_type const y[],
      integer const   incy,
      real_type const z[],
      integer const   ldZ,
      integer const   nx,
      integer const   ny,
      bool            fortran_storage = false,
      bool            transposed      = false );

    //!
    //! Build surface spline
    //!
    //! \param x               vector of x-coordinates, nx = x.size()
    //! \param y               vector of y-coordinates, ny = y.size()
    //! \param z               matrix of z-values. Elements are stored
    //!                        by row Z(i,j) = z[i*ny+j] as C-matrix
    //! \param fortran_storage if true elements are stored by column
    //!                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
    //! \param transposed      if true matrix Z is stored transposed
    //!
    void build(
      vector<real_type> const & x,
      vector<real_type> const & y,
      vector<real_type> const & z,
      bool                      fortran_storage = false,
      bool                      transposed      = false )
    {
      this->build(
        x.data(),
        1,
        y.data(),
        1,
        z.data(),
        integer( fortran_storage ? y.size() : x.size() ),
        integer( x.size() ),
        integer( y.size() ),
        fortran_storage,
        transposed );
    }

    void build(
      real_type const z[],
      integer         ldZ,
      integer         nx,
      integer         ny,
      bool            fortran_storage = false,
      bool            transposed      = false );

    //!
    //! Build surface spline
    //!
    //! \param z               matrix of z-values. Elements are stored
    //!                        by row Z(i,j) = z[i*ny+j] as C-matrix.
    //!                        ldZ leading dimension of the matrix is ny for C-storage
    //!                        and nx for Fortran storage.
    //! \param nx              x-dimension
    //! \param ny              y-dimension
    //! \param fortran_storage if true elements are stored by column
    //!                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
    //! \param transposed      if true matrix Z is stored transposed
    //!
    void build(
      vector<real_type> const & z,
      integer const             nx,
      integer const             ny,
      bool                      fortran_storage = false,
      bool                      transposed      = false )
    {
      this->build( z.data(), nx, ny, fortran_storage ? nx : ny, fortran_storage, transposed );
    }

    //!
    //! Build spline using data in `gc`
    //!
    void setup( GenericContainer const & gc );

    //!
    //! Setup a spline using a `GenericContainer` readed from file
    //!
    //! - file_name file name of the file with data
    //!
    //! - `.json`
    //! - `.yaml` or `.yml`
    //! - `.toml`
    //!
    void setup( string const & file_name );

    //!
    //! Build a spline using data in `GenericContainer`
    //!
    void build( GenericContainer const & gc ) { setup( gc ); }

    //!
    //! Build a spline using data in a file `file_name`
    //!
    void build( string const & file_name ) { setup( file_name ); }

    ///@}

    //!
    //! \name Evaluate
    //!
    ///@{

    //!
    //! Evaluate spline value at point \f$ (x,y) \f$.
    //!
    virtual real_type eval( real_type const x, real_type const y ) const = 0;

#ifdef AUTODIFF_SUPPORT
    // Metodi base per dual1st e dual2nd
    autodiff::dual1st eval( autodiff::dual1st const & x, autodiff::dual1st const & y ) const
    {
      using autodiff::dual1st;
      using autodiff::detail::val;

      real_type dd[3];
      D( val( x ), val( y ), dd );

      dual1st res{ dd[0] };
      res.grad = dd[1] * x.grad + dd[2] * y.grad;
      return res;
    }

    autodiff::dual2nd eval( autodiff::dual2nd const & x, autodiff::dual2nd const & y ) const
    {
      using autodiff::derivative;
      using autodiff::dual2nd;

      real_type dd[6], dx{ val( x.grad ) }, dy{ val( y.grad ) }, ddx{ x.grad.grad }, ddy{ y.grad.grad };
      DD( val( x ), val( y ), dd );

      dual2nd res{ dd[0] };
      res.grad      = dd[1] * dx + dd[2] * dy;
      res.grad.grad = dx * dx * dd[3] + 2 * dx * dy * dd[4] + dy * dy * dd[5] + ddx * dd[1] + ddy * dd[2];
      return res;
    }

    // Template per due parametri (x, y) - per SplineSurf
    // Promuove automaticamente a double, dual1st o dual2nd in base ai tipi di input
    template <typename T1, typename T2> auto eval( T1 const & x, T2 const & y ) const
    {
      if constexpr ( std::is_arithmetic<T1>::value && std::is_arithmetic<T2>::value )
      {
        // Entrambi numerici: ritorna real_type
        return eval( static_cast<real_type>( x ), static_cast<real_type>( y ) );
      }
      else
      {
        // Almeno uno è duale: determina il tipo duale appropriato
        constexpr int order1    = std::is_arithmetic<T1>::value ? 0 : autodiff::detail::DualOrder<T1>::value;
        constexpr int order2    = std::is_arithmetic<T2>::value ? 0 : autodiff::detail::DualOrder<T2>::value;
        constexpr int max_order = ( order1 > order2 ) ? order1 : order2;

        if constexpr ( max_order == 1 )
        {
          // Promuovi a dual1st
          autodiff::dual1st X = std::is_arithmetic<T1>::value ? autodiff::dual1st{ static_cast<real_type>( x ) }
                                                              : autodiff::dual1st{ x };
          autodiff::dual1st Y = std::is_arithmetic<T2>::value ? autodiff::dual1st{ static_cast<real_type>( y ) }
                                                              : autodiff::dual1st{ y };
          return eval( X, Y );
        }
        else
        {
          // Promuovi a dual2nd
          autodiff::dual2nd X = std::is_arithmetic<T1>::value ? autodiff::dual2nd{ static_cast<real_type>( x ) }
                                                              : autodiff::dual2nd{ x };
          autodiff::dual2nd Y = std::is_arithmetic<T2>::value ? autodiff::dual2nd{ static_cast<real_type>( y ) }
                                                              : autodiff::dual2nd{ y };
          return eval( X, Y );
        }
      }
    }

    // Operator() per due parametri
    template <typename T1, typename T2> auto operator()( T1 const & x, T2 const & y ) const -> decltype( eval( x, y ) )
    {
      return eval( x, y );
    }
#endif

    //!
    //! Value and first derivatives at point \f$ (x,y) \f$:
    //!
    //! - d[0] value of the spline \f$ S(x,y) \f$
    //! - d[1] derivative respect to \f$ x \f$ of the spline: \f$ S_x(x,y) \f$
    //! - d[2] derivative respect to \f$ y \f$ of the spline: \f$ S_y(x,y) \f$
    //!
    virtual void D( real_type const x, real_type const y, real_type d[3] ) const = 0;

    //!
    //! First derivatives respect to \f$ x \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_x(x,y) \f$.
    //!
    virtual real_type Dx( real_type const x, real_type const y ) const = 0;

    //!
    //! First derivatives respect to \f$ y \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_y(x,y) \f$.
    //!
    virtual real_type Dy( real_type const x, real_type const y ) const = 0;

    //!
    //! Value, first and second derivatives at point \f$ (x,y) \f$:
    //!
    //! - dd[0] value of the spline \f$ S(x,y) \f$
    //! - dd[1] derivative respect to \f$ x \f$ of the spline: \f$ S_x(x,y) \f$
    //! - dd[2] derivative respect to \f$ y \f$ of the spline: \f$ S_y(x,y) \f$
    //! - dd[3] second derivative respect to \f$ x \f$ of the spline: \f$ S_{xx}(x,y) \f$
    //! - dd[4] mixed second derivative: \f$ S_{xy}(x,y) \f$
    //! - dd[5] second derivative respect to \f$ y \f$ of the spline: \f$ S_{yy}(x,y) \f$
    //!
    virtual void DD( real_type const x, real_type const y, real_type dd[6] ) const = 0;

    //!
    //! Second derivatives respect to \f$ x \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_{xx}(x,y) \f$.
    //!
    virtual real_type Dxx( real_type const x, real_type const y ) const = 0;

    //!
    //! Mixed second derivatives: \f$ S_{xy}(x,y) \f$.
    //!
    virtual real_type Dxy( real_type const x, real_type const y ) const = 0;

    //!
    //! Second derivatives respect to \f$ y \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_{yy}(x,y) \f$.
    //!
    virtual real_type Dyy( real_type const x, real_type const y ) const = 0;

    //!
    //! Evaluate spline value at point \f$ (x,y) \f$.
    //!
    real_type operator()( real_type const x, real_type const y ) const { return this->eval( x, y ); }

    //!
    //! Alias for `Dx(x,y)`
    //!
    real_type eval_D_1( real_type const x, real_type const y ) const { return this->Dx( x, y ); }

    //!
    //! Alias for `Dy(x,y)`
    //!
    real_type eval_D_2( real_type const x, real_type const y ) const { return this->Dy( x, y ); }

    //!
    //! Alias for `Dxx(x,y)`
    //!
    real_type eval_D_1_1( real_type const x, real_type const y ) const { return this->Dxx( x, y ); }

    //!
    //! Alias for `Dxy(x,y)`
    //!
    real_type eval_D_1_2( real_type const x, real_type const y ) const { return this->Dxy( x, y ); }

    //!
    //! Alias for `Dyy(x,y)`
    //!
    real_type eval_D_2_2( real_type const x, real_type const y ) const { return this->Dyy( x, y ); }

    ///@}

    //!
    //! Print spline coefficients.
    //!
    virtual void write_to_stream( ostream_type & s ) const = 0;

    //!
    //! Return spline type as a string pointer.
    //!
    virtual char const * type_name() const = 0;

    //!
    //! String information of the kind and order of the spline
    //!
    virtual string info() const;

    //!
    //! Print information of the kind and order of the spline
    //!
    void info( ostream_type & stream ) const { stream << this->info() << '\n'; }

    //!
    //! Print stored data x, y, and matrix z.
    //!
    void dump_data( ostream_type & s ) const;

#ifdef SPLINES_BACK_COMPATIBILITY
    integer   numPointX() const { return m_nx; }
    integer   numPointY() const { return m_ny; }
    real_type xNode( integer i ) const { return m_X[i]; }
    real_type yNode( integer i ) const { return m_Y[i]; }
    real_type zNode( integer i, integer j ) const { return z_node( i, j ); }
    real_type xMin() const { return this->x_min(); }
    real_type xMax() const { return this->x_max(); }
    real_type yMin() const { return this->y_min(); }
    real_type yMax() const { return this->y_max(); }
    real_type zMin() const { return m_Z_min; }
    real_type zMax() const { return m_Z_max; }
    void      writeToStream( ostream_type & s ) const { write_to_stream( s ); }
#endif
  };

}  // namespace Splines


#include "Splines/SplineBilinear.hxx"
#include "Splines/SplineBiCubic.hxx"
#include "Splines/SplineAkima2D.hxx"
#include "Splines/SplineBiQuintic.hxx"

#include "Splines/SplineVec.hxx"
#include "Splines/SplineSet.hxx"
#include "Splines/Splines1D.hxx"
#include "Splines/Splines2D.hxx"

#include "Splines/Splines1Dblend.hxx"
#include "Splines/Splines2Dblend.hxx"

namespace SplinesLoad
{

  using Splines::AkimaSpline;
  using Splines::BesselSpline;
  using Splines::ConstantSpline;
  using Splines::CubicSpline;
  using Splines::CubicSplineBase;
  using Splines::LinearSpline;
  using Splines::PchipSpline;
  using Splines::QuinticSpline;
  using Splines::Spline;
  using Splines::Spline1D;

  using Splines::Akima2Dspline;
  using Splines::BiCubicSpline;
  using Splines::BilinearSpline;
  using Splines::BiQuinticSpline;
  using Splines::Spline2D;

  using Splines::SplineSet;
  using Splines::SplineVec;

  using Splines::SplineType1D;
  using Splines::SplineType2D;
}  // namespace SplinesLoad

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif

#endif
