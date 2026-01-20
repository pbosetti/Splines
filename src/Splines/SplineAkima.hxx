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


  inline void Akima_build(
    real_type const X[],  // puntatori ai dati
    real_type const Y[],
    real_type       Yp[],
    real_type       m[],
    int             N )
  {
    assert( N >= 2 && "Akima_build requires at least 2 points" );

    if ( N == 2 )
    {
      Yp[0] = Yp[1] = ( Y[1] - Y[0] ) / ( X[1] - X[0] );
      return;
    }

    const int n = N - 1;
    
    for ( integer i = 0; i < n; ++i )
      m[i] = (Y[i+1]-Y[i])/(X[i+1]-X[i]);

    auto slope = [&]( int i ) -> real_type
    {
      if ( i == -2 ) return 2 * m[0] - m[1];
      if ( i == -1 ) return 2 * ( 2 * m[0] - m[1] ) - m[0];
      if ( i >= 0 && i <= n - 1 ) return m[i];
      if ( i == n ) return 2 * m[n - 1] - m[n - 2];
      if ( i == n + 1 ) return 2 * ( 2 * m[n - 1] - m[n - 2] ) - m[n - 1];
      return 0.0;
    };

    real_type m_i   = slope( 0 );
    real_type m_im1 = slope( 0 - 1 );
    real_type m_im2 = slope( 0 - 2 );
    real_type m_ip1 = slope( 0 + 1 );
    for ( integer i = 0; i <= n; ++i )
    {
      real_type w1    = std::abs( m_ip1 - m_i ) + 0.5 * std::abs( m_ip1 + m_i );
      real_type w2    = std::abs( m_im1 - m_im2 ) + 0.5 * std::abs( m_im1 + m_im2 );
      real_type sum_w = w1 + w2;

      if ( sum_w > 1e-12 )
        Yp[i] = ( w1 * m_im1 + w2 * m_i ) / sum_w;
      else
        Yp[i] = 0.5 * ( m_im1 + m_i );  // fallback stabile
      m_im2 = m_im1;
      m_im1 = m_i;
      m_i   = m_ip1;
      m_ip1 = slope( i + 2 );
    }
  }

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
