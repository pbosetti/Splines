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

#ifndef SPLINE_HERMITE_HXX
#define SPLINE_HERMITE_HXX

namespace Splines
{

  /*\
   |    _   _                     _ _       ____        _ _
   |   | | | | ___ _ __ _ __ ___ (_) |_ ___/ ___| _ __ | (_)_ __   ___
   |   | |_| |/ _ \ '__| '_ ` _ \| | __/ _ \___ \| '_ \| | | '_ \ / _ \
   |   |  _  |  __/ |  | | | | | | | ||  __/___) | |_) | | | | | |  __/
   |   |_| |_|\___|_|  |_| |_| |_|_|\__\___|____/| .__/|_|_|_| |_|\___|
   |                                             |_|
  \*/

  //! Hermite Spline Management Class
  class HermiteSpline : public CubicSplineBase
  {
  public:
    using CubicSplineBase::build;
    using CubicSplineBase::reserve;

    //!
    //! Build an empty spline of `HermiteSpline` type
    //!
    //! \param name the name of the spline
    //!
    HermiteSpline( string_view name = "HermiteSpline" ) : CubicSplineBase( name ) {}

    //!
    //! Spline destructor.
    //!
    ~HermiteSpline() override {}

    //! Return spline type (as number)
    SplineType1D type() const override { return SplineType1D::HERMITE; }

    // --------------------------- VIRTUALS -----------------------------------

    void build() override { m_search.must_reset(); }

    // block method!
    void build( real_type const[], integer, real_type const[], integer, integer ) override
    {
      UTILS_ERROR( "HermiteSpline[{}]::build(x,incx,y,incy,n) cannot be used\n", m_name );
    }

    //!
    //!
    //! Setup a spline using a `GenericContainer`
    //!
    //! - gc("xdata")  vector with the `x` coordinate of the data
    //! - gc("ydata")  vector with the `y` coordinate of the data
    //! - gc("ypdata") vector with the `y` derivative of the data
    //!
    void setup( GenericContainer const & gc ) override
    {
      /*
      // gc["xdata"]
      // gc["ydata"]
      //
      */
      string const where{ fmt::format( "HermiteSpline[{}]::setup( gc ):", m_name ) };

      std::set<std::string> keywords;
      for ( auto const & pair : gc.get_map( where ) ) { keywords.insert( pair.first ); }
      keywords.erase( "spline_type" );

      GenericContainer const & gc_x{ gc( "xdata", where ) };
      keywords.erase( "xdata" );
      GenericContainer const & gc_y{ gc( "ydata", where ) };
      keywords.erase( "ydata" );
      GenericContainer const & gc_yp{ gc( "ypdata", where ) };
      keywords.erase( "ypdata" );

      vec_real_type x, y, yp;
      {
        string const ff{ fmt::format( "{}, field `xdata'", where ) };
        gc_x.copyto_vec_real( x, ff );
      }
      {
        string const ff{ fmt::format( "{}, field `ydata'", where ) };
        gc_y.copyto_vec_real( y, ff );
      }
      {
        string const ff{ fmt::format( "{}, field `ypdata'", where ) };
        gc_yp.copyto_vec_real( yp, ff );
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

      this->build( x, y, yp );
    }
  };

}  // namespace Splines

#endif  // SPLINE_HERMITE_HXX
