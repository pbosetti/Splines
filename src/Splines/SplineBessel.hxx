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
