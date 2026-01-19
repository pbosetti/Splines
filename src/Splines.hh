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
#include "Utils_search_intervals2.hh"

//!
//! Namespace of Splines library
//!
namespace Splines
{

  using std::basic_istream;
  using std::basic_ostream;
  using std::cerr;
  using std::cin;
  using std::copy_n;
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

  using GC_namespace::GC_type;
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
    const real_type dP    = P1 - P0;

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
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#include "Splines/SplineParametrization.hxx"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

namespace Splines {

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

    Utils::SearchInterval<real_type, integer> m_search;

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
      string const where = fmt::format( "Spline[{}]::setup( gc ):", m_name );

      std::set<std::string> keywords;
      for ( auto const & pair : gc.get_map( where ) ) { keywords.insert( pair.first ); }

      GenericContainer const & gc_x = gc( "xdata", where );
      keywords.erase( "xdata" );
      GenericContainer const & gc_y = gc( "ydata", where );
      keywords.erase( "ydata" );
      keywords.erase( "spline_type" );

      vec_real_type x, y;
      {
        string const ff = fmt::format( "{}, field `xdata'", where );
        gc_x.copyto_vec_real( x, ff );
      }
      {
        string const ff = fmt::format( "{}, field `ydata'", where );
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
        integer const saved_npts = m_npts;  // salvo npts perche reserve lo azzera
        Malloc_real   mem( "Spline::push_back" );
        mem.allocate( 2 * m_npts );
        real_type * Xsaved = mem( m_npts );
        real_type * Ysaved = mem( m_npts );

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
      real_type const Tx = x0 - m_X[0];
      real_type *     ix = m_X;
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
      string res = fmt::format( "Spline `{}` of type: {} of order: {}", m_name, type_name(), order() );
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


  inline void Hermite3( real_type const x, real_type const H, real_type base[4] )
  {
    real_type const X = x / H;
    base[1]           = X * X * ( 3 - 2 * X );
    base[0]           = 1 - base[1];
    base[2]           = x * ( X * ( X - 2 ) + 1 );
    base[3]           = x * X * ( X - 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline void Hermite3_D( real_type const x, real_type const H, real_type base_D[4] )
  {
    real_type const X = x / H;
    base_D[0]         = 6.0 * X * ( X - 1.0 ) / H;
    base_D[1]         = -base_D[0];
    base_D[2]         = ( ( 3 * X - 4 ) * X + 1 );
    base_D[3]         = X * ( 3 * X - 2 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline void Hermite3_DD( real_type const x, real_type const H, real_type base_DD[4] )
  {
    real_type const X = x / H;
    base_DD[0]        = ( 12 * X - 6 ) / ( H * H );
    base_DD[1]        = -base_DD[0];
    base_DD[2]        = ( 6 * X - 4 ) / H;
    base_DD[3]        = ( 6 * X - 2 ) / H;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline void Hermite3_DDD( real_type, real_type const H, real_type base_DDD[4] )
  {
    base_DDD[0] = 12 / ( H * H * H );
    base_DDD[1] = -base_DDD[0];
    base_DDD[2] = 6 / ( H * H );
    base_DDD[3] = base_DDD[2];
  }

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

}  // namespace Splines

#include "Splines/SplineCubicBase.hxx"
#include "Splines/SplineHermite.hxx"
#include "Splines/SplineAkima.hxx"
#include "Splines/SplineBessel.hxx"
#include "Splines/SplineConstant.hxx"
#include "Splines/SplineLinear.hxx"
#include "Splines/SplineCubic.hxx"
#include "Splines/SplinePchip.hxx"
#include "Splines/SplineQuinticBase.hxx"
#include "Splines/SplineQuintic.hxx"

#include "Splines/SplineSurf.hxx"


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
