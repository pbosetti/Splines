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
   |   _   _                     _ _
   |  | | | | ___ _ __ _ __ ___ (_) |_ ___
   |  | |_| |/ _ \ '__| '_ ` _ \| | __/ _ \
   |  |  _  |  __/ |  | | | | | | | ||  __/
   |  |_| |_|\___|_|  |_| |_| |_|_|\__\___|
  \*/

  // --------------------------------------------------------------------------

inline void Hermite5(real_type x, real_type H, real_type base[6])
{
  // Potenze di base
  const real_type x2 = x * x;
  const real_type x3 = x2 * x;
  const real_type x4 = x2 * x2;

  const real_type H2 = H * H;
  const real_type H4 = H2 * H2;

  // (H - x)
  const real_type dx  = H - x;
  const real_type dx2 = dx * dx;
  const real_type dx3 = dx2 * dx;

  // Inversi
  const real_type invH  = real_type(1) / H;
  const real_type invH4 = real_type(1) / H4;
  const real_type invH5 = invH * invH4;

  // Base functions
  base[0] = invH5 * dx3 * (3 * x * H + H2 + 6 * x2);

  base[1] = invH5 * (
              -15 * H * x4 +
                6 * x4 * x +
               10 * H2 * x3
            );

  base[2] = invH4 * dx3 * x * (H + 3 * x);

  base[3] = invH4 * (3 * x - 4 * H) * dx * x3;

  const real_type c = invH / (2 * H2);

  base[4] = c * dx3 * x2;
  base[5] = c * dx2 * x3;
}


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline void Hermite5_D(real_type x, real_type H, real_type base_D[6])
{
  // Differenza
  const real_type dx  = H - x;
  const real_type dx2 = dx * dx;

  // Potenze di x
  const real_type x2 = x * x;
  const real_type x4 = x2 * x2;

  // Potenze di H
  const real_type H2 = H * H;
  const real_type H4 = H2 * H2;

  // Inversi
  const real_type invH  = real_type(1) / H;
  const real_type invH4 = real_type(1) / H4;

  // Termine comune
  const real_type common = 30 * x2 * dx2 * invH * invH4;

  // Derivate delle basi
  base_D[0] = -common;
  base_D[1] =  common;

  base_D[2] = invH4 * (H - 3 * x) * (H + 5 * x) * dx2;

  base_D[3] = invH4 * (
                x2 * (28 * x * H - 12 * H2)
                - 15 * x4
              );

  const real_type c = invH / (2 * H2);

  base_D[4] = c * (2 * H - 5 * x) * dx2 * x;
  base_D[5] = c * (3 * H - 5 * x) * x2 * dx;
}

inline void Hermite5_DD(real_type x, real_type H, real_type base_DD[6])
{
  // Differenza
  const real_type dx = H - x;

  // Potenze di x
  const real_type x2 = x * x;

  // Potenze di H
  const real_type H2 = H * H;
  const real_type H4 = H2 * H2;

  // Inversi
  const real_type invH  = real_type(1) / H;
  const real_type invH4 = real_type(1) / H4;

  // Termine comune
  const real_type common = 60 * (H - 2 * x) * x * dx * invH * invH4;

  base_DD[0] = -common;
  base_DD[1] =  common;

  base_DD[2] = 12 * invH4 * dx * (5 * x - 3 * H) * x;

  base_DD[3] = 12 * invH4 * x * dx * (5 * x - 2 * H);

  const real_type c = invH / H2;

  base_DD[4] = c * (10 * x2 + H2 - 8 * H * x) * dx;

  base_DD[5] = c * (x2 * (10 * x - 12 * H) + 3 * x * H2);
}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
inline void Hermite5_DDD(real_type x, real_type H, real_type base_DDD[6])
{
  // Potenze
  const real_type x2 = x * x;
  const real_type H2 = H * H;
  const real_type H4 = H2 * H2;

  // Prodotti utili
  const real_type Hx = H * x;

  // Inversi
  const real_type invH  = real_type(1) / H;
  const real_type invH4 = real_type(1) / H4;

  // Termine comune
  const real_type common = invH * invH4 * (360 * Hx - 60 * H2 - 360 * x2);

  base_DDD[0] =  common;
  base_DDD[1] = -common;

  base_DDD[2] = invH4 * (192 * Hx - 36 * H2 - 180 * x2);
  base_DDD[3] = invH4 * (168 * Hx - 24 * H2 - 180 * x2);

  const real_type c = invH / H2;

  base_DDD[4] = c * (36 * Hx - 9 * H2 - 30 * x2);
  base_DDD[5] = c * (3 * H2 - 24 * Hx + 30 * x2);
}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
inline void Hermite5_DDDD(real_type x, real_type H, real_type base_DDDD[6])
{
  // Potenze di H
  const real_type H2 = H * H;
  const real_type H4 = H2 * H2;

  // Inversi
  const real_type invH  = real_type(1) / H;
  const real_type invH4 = real_type(1) / H4;

  // Termine comune
  const real_type common = invH * invH4 * (360 * H - 720 * x);

  base_DDDD[0] =  common;
  base_DDDD[1] = -common;

  base_DDDD[2] = invH4 * (192 * H - 360 * x);
  base_DDDD[3] = invH4 * (168 * H - 360 * x);

  const real_type c = invH / H2;

  base_DDDD[4] = c * (36 * H - 60 * x);
  base_DDDD[5] = c * (60 * x - 24 * H);
}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
inline void Hermite5_DDDDD(real_type /*x*/, real_type H, real_type base_DDDDD[6])
{
  // Potenze di H
  const real_type H2 = H * H;
  const real_type H4 = H2 * H2;

  // Inversi
  const real_type invH  = real_type(1) / H;
  const real_type invH4 = real_type(1) / H4;

  // Costanti
  const real_type c0 = 720 * invH * invH4;
  const real_type c1 = 360 * invH4;
  const real_type c2 = 60  * invH / H2;

  base_DDDDD[0] = -c0;
  base_DDDDD[1] =  c0;

  base_DDDDD[2] = -c1;
  base_DDDDD[3] = -c1;

  base_DDDDD[4] = -c2;
  base_DDDDD[5] =  c2;
}

  // --------------------------------------------------------------------------

  inline real_type bilinear3( real_type const p[4], real_type const M[4][4], real_type const q[4] )
  {
    return p[0] * ( M[0][0] * q[0] + M[0][1] * q[1] + M[0][2] * q[2] + M[0][3] * q[3] ) +
           p[1] * ( M[1][0] * q[0] + M[1][1] * q[1] + M[1][2] * q[2] + M[1][3] * q[3] ) +
           p[2] * ( M[2][0] * q[0] + M[2][1] * q[1] + M[2][2] * q[2] + M[2][3] * q[3] ) +
           p[3] * ( M[3][0] * q[0] + M[3][1] * q[1] + M[3][2] * q[2] + M[3][3] * q[3] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline real_type bilinear5( real_type const p[6], real_type const M[6][6], real_type const q[6] )
  {
    return p[0] *
             ( M[0][0] * q[0] + M[0][1] * q[1] + M[0][2] * q[2] + M[0][3] * q[3] + M[0][4] * q[4] + M[0][5] * q[5] ) +
           p[1] *
             ( M[1][0] * q[0] + M[1][1] * q[1] + M[1][2] * q[2] + M[1][3] * q[3] + M[1][4] * q[4] + M[1][5] * q[5] ) +
           p[2] *
             ( M[2][0] * q[0] + M[2][1] * q[1] + M[2][2] * q[2] + M[2][3] * q[3] + M[2][4] * q[4] + M[2][5] * q[5] ) +
           p[3] *
             ( M[3][0] * q[0] + M[3][1] * q[1] + M[3][2] * q[2] + M[3][3] * q[3] + M[3][4] * q[4] + M[3][5] * q[5] ) +
           p[4] *
             ( M[4][0] * q[0] + M[4][1] * q[1] + M[4][2] * q[2] + M[4][3] * q[3] + M[4][4] * q[4] + M[4][5] * q[5] ) +
           p[5] *
             ( M[5][0] * q[0] + M[5][1] * q[1] + M[5][2] * q[2] + M[5][3] * q[3] + M[5][4] * q[4] + M[5][5] * q[5] );
  }

  // --------------------------------------------------------------------------

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
