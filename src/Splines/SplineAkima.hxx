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

#ifndef SPLINE_AKIMA_HXX
#define SPLINE_AKIMA_HXX

/*\
 |      _    _    _                   ____        _ _
 |     / \  | | _(_)_ __ ___   __ _  / ___| _ __ | (_)_ __   ___
 |    / _ \ | |/ / | '_ ` _ \ / _` | \___ \| '_ \| | | '_ \ / _ \
 |   / ___ \|   <| | | | | | | (_| |  ___) | |_) | | | | | |  __/
 |  /_/   \_\_|\_\_|_| |_| |_|\__,_| |____/| .__/|_|_|_| |_|\___|
 |                                         |_|
\*/

namespace Splines
{

  //!
  //! Smooth Curve Fitting Based on Local Procedures
  //!
  //! *Reference*
  //!
  //! - *Hiroshi Akima*, Journal of the ACM, Vol.17, No. 4, 589-602, 1970.
  //!
  class AkimaSpline : public CubicSplineBase
  {
  public:
    using CubicSplineBase::build;
    using CubicSplineBase::reserve;

    //!
    //! Build an empty spline of `AkimaSpline` type
    //!
    //! \param name the name of the spline
    //!
    explicit AkimaSpline( string_view name = "AkimaSpline" ) : CubicSplineBase( name ) {}

    //!
    //! Spline destructor.
    //!
    ~AkimaSpline() override {}

    //!
    //! Return spline type (as number).
    //!
    SplineType1D type() const override { return SplineType1D::AKIMA; }

    // --------------------------- VIRTUALS -----------------------------------

    void build() override
    {
      string const msg = fmt::format( "AkimaSpline[{}]::build():", m_name );
      UTILS_ASSERT( m_npts > 1, "{} npts = {} not enought points\n", msg, m_npts );
      Utils::check_NaN( m_X, msg + " X ", m_npts, __LINE__, __FILE__ );
      Utils::check_NaN( m_Y, msg + " Y ", m_npts, __LINE__, __FILE__ );
      integer ibegin = 0;
      integer iend   = 0;

      Malloc_real m_mem( "AkimaSpline::work memory" );
      real_type * m_work{ m_mem.malloc( m_npts ) };

      do
      {
        // cerca intervallo monotono strettamente crescente
        for ( ++iend; iend < m_npts && m_X[iend - 1] < m_X[iend]; ++iend ) {}
        Akima_build( m_X + ibegin, m_Y + ibegin, m_Yp + ibegin, m_work, iend - ibegin );
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
      string const where = fmt::format( "AkimaSpline[{}]::setup():", m_name );

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

//
// EOF: SplineAkima.hxx
//
