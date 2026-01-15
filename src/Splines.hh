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
    T const X{ x / H };
    base[1] = X * X * ( 3 - 2 * X );
    base[0] = 1 - base[1];
    base[2] = x * ( X * ( X - 2 ) + 1 );
    base[3] = x * X * ( X - 1 );
  }
#endif

  void Hermite3( real_type const x, real_type const H, real_type base[4] );
  void Hermite3_D( real_type const x, real_type const H, real_type base_D[4] );
  void Hermite3_DD( real_type const x, real_type const H, real_type base_DD[4] );
  void Hermite3_DDD( real_type const x, real_type const H, real_type base_DDD[4] );

#ifdef AUTODIFF_SUPPORT
  template <typename T> inline void Hermite5( T const & x, real_type const H, T base[6] )
  {
    auto const t1{ H * H };
    auto const t4{ x * x };
    auto const t7{ H - x };
    auto const t8{ t7 * t7 };
    auto const t9{ t8 * t7 };
    auto const t11{ t1 * t1 };
    auto const t2{ 1 / t11 };
    auto const t3{ 1 / H };
    auto const t13{ t3 * t2 };
    auto const t14{ t4 * x };
    auto const t17{ t4 * t4 };
    base[0]             = t13 * t9 * ( 3.0 * x * H + t1 + 6.0 * t4 );
    base[1]             = t13 * ( -15.0 * H * t17 + 6.0 * t17 * x + 10.0 * t1 * t14 );
    base[2]             = t2 * t9 * x * ( H + 3 * x );
    base[3]             = t2 * ( 3 * x - 4 * H ) * t7 * t14;
    real_type const t36 = t3 / t1 / 2;
    base[4]             = t36 * t9 * t4;
    base[5]             = t36 * t8 * t14;
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
    real_type const H,
    real_type const P0,
    real_type const P1,
    real_type const DP0,
    real_type const DP1,
    real_type &     A,
    real_type &     B,
    real_type &     C,
    real_type &     D )
  {
    real_type H2{ H * H };
    real_type P10{ P1 - P0 };
    A = ( DP0 + DP1 - 2 * P10 / H ) / H2;
    B = ( 3 * P10 / H - ( 2 * DP0 + DP1 ) ) / H;
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
    real_type const h,
    real_type const P0,
    real_type const P1,
    real_type const DP0,
    real_type const DP1,
    real_type const DDP0,
    real_type const DDP1,
    real_type &     A,
    real_type &     B,
    real_type &     C,
    real_type &     D,
    real_type &     E,
    real_type &     F )
  {
    real_type h2{ h * h };
    real_type h3{ h * h2 };
    real_type P10{ P1 - P0 };
    A = ( ( DDP1 - DDP0 ) / 2 + ( 6 * P10 / h - 3 * ( DP0 + DP1 ) ) / h ) / h3;
    B = ( ( 1.5 * DDP0 - DDP1 ) + ( ( 8 * DP0 + 7 * DP1 ) - 15 * P10 / h ) / h ) / h2;
    C = ( 0.5 * DDP1 - 1.5 * DDP0 + ( 10 * P10 / h - ( 6 * DP0 + 4 * DP1 ) ) / h ) / h;
    D = DDP0 / 2;
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

  //! Check if cubic spline with this data is monotone, return -1 no, 0 yes, 1 strictly monotone
  inline integer check_cubic_spline_monotonicity(
    real_type const X[],
    real_type const Y[],
    real_type const Yp[],
    integer const   npts )
  {
    // check monotonicity of data: (assuming X monotone)
    integer flag{ 1 };
    for ( integer i{ 1 }; i < npts; ++i )
    {
      if ( Y[i - 1] > Y[i] ) return -2;                                      // non monotone data
      if ( Utils::is_zero( Y[i - 1] - Y[i] ) && X[i - 1] < X[i] ) flag = 0;  // non strict monotone
    }
    // pag 146 Methods of Shape-Preserving Spline Approximation, K
    for ( integer i{ 1 }; i < npts; ++i )
    {
      if ( X[i] <= X[i - 1] ) continue;  // skip duplicate points
      real_type const dd = ( Y[i] - Y[i - 1] ) / ( X[i] - X[i - 1] );
      real_type const m0 = Yp[i - 1] / dd;
      real_type const m1 = Yp[i] / dd;
      if ( m0 < 0 || m1 < 0 ) return -1;  // non monotone
      if ( m0 <= 3 && m1 <= 3 )
      {
        if ( flag > 0 && i > 1 && ( Utils::is_zero( m0 ) || Utils::is_zero( m0 - 3 ) ) ) flag = 0;
        if ( flag > 0 && i < npts - 1 && ( Utils::is_zero( m1 ) || Utils::is_zero( m1 - 3 ) ) ) flag = 0;
      }
      else
      {
        real_type const tmp1 = 2 * m0 + m1 - 3;
        real_type const tmp2 = 2 * ( m0 + m1 - 2 );
        real_type const tmp3 = m0 * tmp2 - ( tmp1 * tmp1 );
        if ( tmp2 >= 0 )
        {
          if ( tmp3 < 0 ) return -1;  // non monotone spline
        }
        else
        {
          if ( tmp3 > 0 ) return -1;
        }
        if ( Utils::is_zero( tmp3 ) ) flag = 0;
      }
    }
    return flag;  // passed all check
  }
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
    for ( integer k{ 1 }; k < npts - 1; ++k ) t[k] = static_cast<real_type>( k ) / static_cast<real_type>( npts );
  }

  //!
  //! Compute nodes for the spline using chordal distribution
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
    for ( integer k{ 1 }; k < npts; ++k )
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
    for ( integer k{ 1 }; k < npts - 1; ++k ) t[k] /= t[npts - 1];
    t[npts - 1] = 1;
  }

  //!
  //! Compute nodes for the spline using centripetal distribution
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
    for ( integer k{ 1 }; k < npts; ++k )
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
    for ( integer k{ 1 }; k < npts - 1; ++k ) t[k] /= t[npts - 1];
    t[npts - 1] = 1;
  }

#if 0

  void
  universal(
    integer         dim,
    integer         npts,
    real_type const pnts[],
    integer         ld_pnts,
    real_type       t[]
  ); // to be done

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  FoleyNielsen(
    integer         dim,
    integer         npts,
    real_type const pnts[],
    integer         ld_pnts,
    real_type       t[]
  ); // to be done

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  FangHung(
    integer         dim,
    integer         npts,
    real_type const pnts[],
    integer         ld_pnts,
    real_type       t[]
  ); // to be done

#endif

  //!
  //! Compute nodes for the spline using universal distribution
  //!
  //! \param[in]  dim     dimension of the points
  //! \param[in]  npts    number of points
  //! \param[in]  pnts    matrix whose columns are the points
  //! \param[in]  ld_pnts leading dimension of the matrix (fortran storage)
  //! \param[out] t       vector of the computed nodes
  //!
  void universal( integer const dim, integer const npts, real_type const pnts[], integer const ld_pnts, real_type t[] );

  //!
  //! Compute nodes for the spline using `FoleyNielsen` distribution
  //!
  //! \param[in]  dim     dimension of the points
  //! \param[in]  npts    number of points
  //! \param[in]  pnts    matrix whose columns are the points
  //! \param[in]  ld_pnts leading dimension of the matrix (fortran storage)
  //! \param[out] t       vector of the computed nodes
  //!
  void FoleyNielsen(
    integer const   dim,
    integer const   npts,
    real_type const pnts[],
    integer const   ld_pnts,
    real_type       t[] );

  //!
  //! Compute nodes for the spline using `FangHung` distribution
  //!
  //! \param[in]  dim     dimension of the points
  //! \param[in]  npts    number of points
  //! \param[in]  pnts    matrix whose columns are the points
  //! \param[in]  ld_pnts leading dimension of the matrix (fortran storage)
  //! \param[out] t       vector of the computed nodes
  //!
  void FangHung( integer const dim, integer const npts, real_type const pnts[], integer const ld_pnts, real_type t[] );

  /*\
   |   ____                      _     ___       _                       _
   |  / ___|  ___  __ _ _ __ ___| |__ |_ _|_ __ | |_ ___ _ ____   ____ _| |
   |  \___ \ / _ \/ _` | '__/ __| '_ \ | || '_ \| __/ _ \ '__\ \ / / _` | |
   |   ___) |  __/ (_| | | | (__| | | || || | | | ||  __/ |   \ V / (_| | |
   |  |____/ \___|\__,_|_|  \___|_| |_|___|_| |_|\__\___|_|    \_/ \__,_|_|
  \*/

  //!
  //! Manage Search intervals
  //!
  class SearchInterval
  {
    static integer const m_table_size{ 400 };

    string const * p_name{ nullptr };
    integer *      p_npts{ nullptr };
    bool *         p_curve_is_closed{ nullptr };
    bool *         p_curve_can_extend{ nullptr };

    mutable real_type ** p_X{ nullptr };
    mutable real_type    m_x_min{ 0 };
    mutable real_type    m_x_max{ 0 };
    mutable real_type    m_x_range{ 0 };
    mutable real_type    m_dx{ 0 };
    mutable integer      m_LO[m_table_size + 2];  // to avoid overflow and replicate last point
    mutable integer      m_HI[m_table_size + 2];  // to avoid overflow and replicate last point
    mutable bool         m_must_reset{ true };
    mutable std::mutex   m_mutex;

    void reset() const
    {
      integer           n{ *p_npts };
      real_type const * X{ *p_X };

      m_x_min   = X[0];
      m_x_max   = X[n - 1];
      m_x_range = m_x_max - m_x_min;
      m_dx      = m_x_range / m_table_size;

      //
      //
      //       0  1     2     3             4 5 6              7     8
      //  X    +--+-----+-----+-------------+-+-+--------------+-----+
      //
      //       0        2        3        -(3)     6        -(6)     8
      // TABLE |--------|--------|--------|--------|--------|--------|
      //                2        -(4)     4        -(7)     7        8
      //       0        2        -        4                 7        -
      //
      //       +--------+                                              [0..2]
      //                +-------------------+                          [2..4]
      //                       +------------+                          [3..4]
      //                       +---------------------------------+     [3..7]
      //                                        +----------------+     [6..7]
      //                                        +--------------------+ [6..8]
      //
      std::fill_n( m_LO, m_table_size + 1, -1 );
      std::fill_n( m_HI, m_table_size + 1, -1 );
      for ( integer k = 0; k < n; ++k )
      {
        real_type pos{ ( X[k] - m_x_min ) / m_dx };
        integer   i_LO{ static_cast<integer>( std::ceil( pos + 1e-6 ) ) };
        m_LO[i_LO] = std::min( k, n - 1 );
      }
      m_LO[0] = 0;
      for ( integer k{ n - 1 }; k >= 0; --k )
      {
        real_type pos{ ( X[k] - m_x_min ) / m_dx };
        integer   i_HI{ static_cast<integer>( std::floor( pos - 1e-6 ) ) };
        if ( i_HI < 0 ) i_HI = 0;
        if ( m_HI[i_HI] == -1 ) m_HI[i_HI] = k;
      }
      m_HI[m_table_size] = n - 1;

      for ( integer i = 0; i < m_table_size; ++i )
        if ( m_LO[i + 1] == -1 ) m_LO[i + 1] = m_LO[i];
      for ( integer i{ m_table_size }; i > 0; --i )
        if ( m_HI[i - 1] == -1 ) m_HI[i - 1] = m_HI[i];
      m_LO[m_table_size + 1] = m_LO[m_table_size];  // replica ultimo nodo
      m_HI[m_table_size + 1] = m_HI[m_table_size];  // replica ultimo nodo

      m_must_reset = false;
    }

  public:
    SearchInterval( SearchInterval const & )                   = delete;
    SearchInterval const & operator=( SearchInterval const & ) = delete;
    SearchInterval( SearchInterval && )                        = delete;
    SearchInterval & operator=( SearchInterval && )            = delete;

    SearchInterval() {}

    void setup( string const * name, integer * n, real_type ** X, bool * is_closed, bool * can_extend )
    {
      p_name             = name;
      p_npts             = n;
      p_X                = X;
      p_curve_is_closed  = is_closed;
      p_curve_can_extend = can_extend;
      m_must_reset       = true;
    }

    //!
    //! Find interval containing `res.second` using binary search.
    //! Return result in `res.first`
    //!
    void find( std::pair<integer, real_type> & res ) const
    {
      {
        std::lock_guard<std::mutex> lock( m_mutex );
        if ( m_must_reset ) this->reset();
      }

      integer const &   n{ *p_npts };
      string const &    name{ *p_name };
      real_type const * X{ *p_X };
      UTILS_ASSERT( n > 0, "in SearchInterval::find({}), n⁰points == 0!", name );

      integer &   pos{ res.first };
      real_type & x{ res.second };

#if 1
      // casi out of bound
      if ( x > m_x_max )
      {
        if ( *p_curve_is_closed ) { x -= m_x_range * std::floor( ( x - m_x_min ) / m_x_range ); }
        else
        {
          pos = n - 2;
          return;
        }
      }
      else if ( x < m_x_min )
      {
        if ( *p_curve_is_closed ) { x -= m_x_range * std::floor( ( x - m_x_min ) / m_x_range ); }
        else
        {
          pos = 0;
          return;
        }
      }

      // uso table
      integer i_cell{ static_cast<integer>( std::floor( ( x - m_x_min ) / m_dx ) ) };
      integer k_LO{ m_LO[i_cell] };
      integer k_HI{ m_HI[i_cell + 1] };

      UTILS_ASSERT(
        x >= X[k_LO] && x <= X[k_HI],
        "Spline::SearchInterval, x={}, ipos={}, dx={}, X[{}]={}, X[{}]={}, range=[{},{}]\n",
        x,
        i_cell,
        m_dx,
        k_LO,
        X[k_LO],
        k_HI,
        X[k_HI],
        m_x_min,
        m_x_max );
#else
      integer k_LO = 0;
      integer k_HI = n;
#endif

      // binary search
      while ( k_HI > k_LO + 1 )
      {
        integer k_M{ k_LO + ( k_HI - k_LO ) / 2 };
        if ( x < X[k_M] )
          k_HI = k_M;
        else
          k_LO = k_M;
      }

      pos = k_LO;
      if ( Utils::is_zero( X[pos] - X[pos + 1] ) ) --pos;  // caso nodi ripetuti
    }

    void must_reset() { m_must_reset = true; }
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
    bool         m_curve_is_closed{ false };
    bool         m_curve_can_extend{ true };
    bool         m_curve_extended_constant{ false };

    integer     m_npts{ 0 };
    integer     m_npts_reserved{ 0 };
    real_type * m_X{ nullptr };  // allocated in the derived class!
    real_type * m_Y{ nullptr };  // allocated in the derived class!

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
      integer N{ m_npts };
      if ( type() == SplineType1D::CONSTANT ) --N;
      return *std::min_element( m_Y, m_Y + N );
    }

    //!
    //! return y-maximum spline value
    //!
    real_type y_max() const
    {
      integer N{ m_npts };
      if ( type() == SplineType1D::CONSTANT ) --N;
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
      integer N{ integer( x.size() ) };
      if ( N > integer( y.size() ) ) N = integer( y.size() );
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
      real_type const dx{ ( x_max() - x_min() ) / nintervals };
      for ( integer i = 0; i <= nintervals; ++i )
      {
        real_type x{ x_min() + i * dx };
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
