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

#ifndef SPLINE_BESSEL_HXX
#define SPLINE_BESSEL_HXX

/*\
 |   ____                     _ ____        _ _
 |  | __ )  ___  ___ ___  ___| / ___| _ __ | (_)_ __   ___
 |  |  _ \ / _ \/ __/ __|/ _ \ \___ \| '_ \| | | '_ \ / _ \
 |  | |_) |  __/\__ \__ \  __/ |___) | |_) | | | | | |  __/
 |  |____/ \___||___/___/\___|_|____/| .__/|_|_|_| |_|\___|
 |                                   |_|
\*/

namespace Splines
{

  /**
   * @brief Bessel cubic derivative reconstruction (true order 3).
   *
   * First derivatives computed as derivatives of the local cubic
   * interpolating polynomial (4 points).
   *
   * Requirements:
   *  - X strictly monotonic
   *  - npts >= 4
   */
  inline void Bessel_build( real_type const X[], real_type const Y[], real_type Yp[], integer const npts )
  {
    UTILS_ASSERT( npts >= 2, "Bessel_build: npts={} >= 2 required\n", npts );

    if ( npts == 2 )
    {  // solo 2 punti, niente da fare
      Yp[0] = Yp[1] = ( Y[1] - Y[0] ) / ( X[1] - X[0] );
      return;
    }

    // ---- left boundary ----
    if ( npts <= 3 )
    {
      Yp[0] = Utils::first_derivative_3p( X[0], Y[0], X[1], Y[1], X[2], Y[2] );
      Yp[1] = Utils::first_derivative_3p( X[1], Y[1], X[2], Y[2], X[0], Y[0] );
    }
    else
    {
      Yp[0] = Utils::first_derivative_4p( X[0], Y[0], X[1], Y[1], X[2], Y[2], X[3], Y[3] );
      Yp[1] = Utils::first_derivative_4p( X[1], Y[1], X[2], Y[2], X[3], Y[3], X[0], Y[0] );
    }

    // ---- interior points ----
    for ( integer i = 2; i < npts - 2; ++i )
      Yp[i] = Utils::first_derivative_5p(
        X[i - 0],
        Y[i - 0],
        X[i + 1],
        Y[i + 1],
        X[i - 1],
        Y[i - 1],
        X[i + 2],
        Y[i + 2],
        X[i - 2],
        Y[i - 2] );

    // ---- right boundary ----
    if ( npts <= 3 )
    {
      Yp[npts - 1] =
        Utils::first_derivative_3p( X[npts - 1], Y[npts - 1], X[npts - 2], Y[npts - 2], X[npts - 3], Y[npts - 3] );
      Yp[npts - 2] =
        Utils::first_derivative_3p( X[npts - 2], Y[npts - 2], X[npts - 3], Y[npts - 3], X[npts - 1], Y[npts - 1] );
    }
    else
    {
      Yp[npts - 1] = Utils::first_derivative_4p(
        X[npts - 1],
        Y[npts - 1],
        X[npts - 2],
        Y[npts - 2],
        X[npts - 3],
        Y[npts - 3],
        X[npts - 4],
        Y[npts - 4] );
      Yp[npts - 2] = Utils::first_derivative_4p(
        X[npts - 2],
        Y[npts - 2],
        X[npts - 3],
        Y[npts - 3],
        X[npts - 4],
        Y[npts - 4],
        X[npts - 1],
        Y[npts - 1] );
    }
  }

  //!
  //! Bessel spline class
  //!
  class BesselSpline : public CubicSplineBase
  {
  public:
    using CubicSplineBase::build;
    using CubicSplineBase::reserve;

    //!
    //! Build an empty spline of `BesselSpline` type
    //!
    //! \param name the name of the spline
    //!
    explicit BesselSpline( string_view name = "BesselSpline" ) : CubicSplineBase( name ) {}

    //!
    //! spline destructor
    //!
    ~BesselSpline() override {}

    //!
    //! Return spline type (as number)
    //!
    SplineType1D type() const override { return SplineType1D::BESSEL; }

    // --------------------------- VIRTUALS -----------------------------------

    void build() override
    {
      string msg{ fmt::format( "BesselSpline[{}]::build():", m_name ) };
      UTILS_ASSERT( m_npts > 1, "{} npts={} not enought points\n", msg, m_npts );
      Utils::check_NaN( m_X, msg + " X", m_npts, __LINE__, __FILE__ );
      Utils::check_NaN( m_Y, msg + " Y", m_npts, __LINE__, __FILE__ );
      integer ibegin{ 0 };
      integer iend{ 0 };
      do
      {
        // cerca intervallo monotono strettamente crescente
        for ( ++iend; iend < m_npts && m_X[iend - 1] < m_X[iend]; ++iend ) {}
        Bessel_build( m_X + ibegin, m_Y + ibegin, m_Yp + ibegin, iend - ibegin );
        ibegin = iend;
      } while ( iend < m_npts );

      Utils::check_NaN( m_Yp, msg + " Yp", m_npts, __LINE__, __FILE__ );
      m_search.must_reset();
    }

    void setup( GenericContainer const & gc ) override
    {
      /*
      // gc["xdata"]
      // gc["ydata"]
      //
      */
      string const where{ fmt::format( "BesselSpline[{}]::setup( gc ):", m_name ) };

      std::set<std::string> keywords;
      for ( auto const & pair : gc.get_map( where ) ) { keywords.insert( pair.first ); }
      keywords.erase( "spline_type" );

      GenericContainer const & gc_x{ gc( "xdata", where ) };
      keywords.erase( "xdata" );
      GenericContainer const & gc_y{ gc( "ydata", where ) };
      keywords.erase( "ydata" );

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

      this->build( x, y );
    }
  };

}  // namespace Splines

#endif
