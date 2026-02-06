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

namespace Splines
{

  /*\
   |   ____  _ _ _                       ____        _ _
   |  | __ )(_) (_)_ __   ___  __ _ _ __/ ___| _ __ | (_)_ __   ___
   |  |  _ \| | | | '_ \ / _ \/ _` | '__\___ \| '_ \| | | '_ \ / _ \
   |  | |_) | | | | | | |  __/ (_| | |   ___) | |_) | | | | | |  __/
   |  |____/|_|_|_|_| |_|\___|\__,_|_|  |____/| .__/|_|_|_| |_|\___|
   |                                          |_|
  \*/

  //! bilinear spline base class
  class BilinearSpline : public SplineSurf
  {
    using SplineSurf::m_nx;
    using SplineSurf::m_ny;
    using SplineSurf::mX;
    using SplineSurf::mY;
    using SplineSurf::mZ;

    void make_spline() override
    {
      m_search_x.must_reset();
      m_search_y.must_reset();
    }

  public:
    using SplineSurf::eval;

    //!
    //! Build an empty spline of `BilinearSpline` type
    //!
    //! \param name the name of the spline
    //!
    explicit BilinearSpline( string_view name = "BilinearSpline" ) : SplineSurf( name ) {}

    //!
    //! Spline destructor.
    //!
    ~BilinearSpline() override = default;

    //! Evaluate the spline at (x,y)
    real_type eval( real_type const x, real_type const y ) const override
    {
      auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

      real_type const u  = dx / DX;
      real_type const v  = dy / DY;
      real_type const u1 = 1 - u;
      real_type const v1 = 1 - v;

      Mat2x2 Z   = mZ.block<2, 2>( i, j ).matrix();
      auto   Z00 = Z.coeff( 0, 0 );
      auto   Z01 = Z.coeff( 0, 1 );
      auto   Z10 = Z.coeff( 1, 0 );
      auto   Z11 = Z.coeff( 1, 1 );
      return u1 * ( Z00 * v1 + Z01 * v ) + u * ( Z10 * v1 + Z11 * v );
    }

    //! Compute x-derivative at (x,y)
    real_type Dx( real_type const x, real_type const y ) const override
    {
      auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

      real_type const v = dy / DY;

      Mat2x2 Z   = mZ.block<2, 2>( i, j ).matrix();
      auto   Z00 = Z.coeff( 0, 0 );
      auto   Z01 = Z.coeff( 0, 1 );
      auto   Z10 = Z.coeff( 1, 0 );
      auto   Z11 = Z.coeff( 1, 1 );
      return ( ( Z10 - Z00 ) * ( 1 - v ) + ( Z11 - Z01 ) * v ) / DX;
    }

    //! Compute y-derivative at (x,y)
    real_type Dy( real_type const x, real_type const y ) const override
    {
      auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

      real_type const u = dx / DX;

      Mat2x2 Z   = mZ.block<2, 2>( i, j ).matrix();
      auto   Z00 = Z.coeff( 0, 0 );
      auto   Z01 = Z.coeff( 0, 1 );
      auto   Z10 = Z.coeff( 1, 0 );
      auto   Z11 = Z.coeff( 1, 1 );
      return ( ( Z01 - Z00 ) * ( 1 - u ) + ( Z11 - Z10 ) * u ) / DY;
    }

    //! Compute value and first derivatives at (x,y)
    void D( real_type const x, real_type const y, real_type d[3] ) const override
    {
      auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

      // Ottimizzazione: Calcolo inversi
      real_type const invDX = 1.0 / DX;
      real_type const invDY = 1.0 / DY;

      real_type const u  = dx * invDX;
      real_type const v  = dy * invDY;
      real_type const u1 = 1 - u;
      real_type const v1 = 1 - v;

      Mat2x2 Z   = mZ.block<2, 2>( i, j ).matrix();
      auto   Z00 = Z.coeff( 0, 0 );
      auto   Z01 = Z.coeff( 0, 1 );
      auto   Z10 = Z.coeff( 1, 0 );
      auto   Z11 = Z.coeff( 1, 1 );

      // Value
      // Riordino per minimizzare operazioni: interpolazione bilineare standard
      real_type const c0 = Z00 * v1 + Z01 * v;
      real_type const c1 = Z10 * v1 + Z11 * v;
      d[0]               = u1 * c0 + u * c1;

      // Dx
      d[1] = ( c1 - c0 ) * invDX;
      // Nota: (c1 - c0) equivale a v1*(Z10-Z00) + v*(Z11-Z01)

      // Dy
      // Derivata rispetto a Y interpolata lungo X
      real_type const dZ_dy_0 = Z01 - Z00;
      real_type const dZ_dy_1 = Z11 - Z10;
      d[2]                    = ( u1 * dZ_dy_0 + u * dZ_dy_1 ) * invDY;
    }

    //! Compute value and all derivatives up to second order at (x,y)
    void DD( real_type const x, real_type const y, real_type dd[6] ) const override
    {
      this->D( x, y, dd );

      auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

      Mat2x2 Z   = mZ.block<2, 2>( i, j ).matrix();
      auto   Z00 = Z.coeff( 0, 0 );
      auto   Z01 = Z.coeff( 0, 1 );
      auto   Z10 = Z.coeff( 1, 0 );
      auto   Z11 = Z.coeff( 1, 1 );

      dd[3] = 0;                                        // Dxx
      dd[4] = ( Z11 - Z10 - Z01 + Z00 ) / ( DX * DY );  // Dxy
      dd[5] = 0;                                        // Dyy
    }

    //! Second derivatives
    real_type Dxx( real_type const, real_type const ) const override { return 0; }

    // FIX MATEMATICO
    real_type Dxy( real_type const x, real_type const y ) const override
    {
      auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

      Mat2x2 Z   = mZ.block<2, 2>( i, j ).matrix();
      auto   Z00 = Z.coeff( 0, 0 );
      auto   Z01 = Z.coeff( 0, 1 );
      auto   Z10 = Z.coeff( 1, 0 );
      auto   Z11 = Z.coeff( 1, 1 );

      return ( Z11 - Z10 - Z01 + Z00 ) / ( DX * DY );
    }

    real_type Dyy( real_type const, real_type const ) const override { return 0; }

#ifdef AUTODIFF_SUPPORT
    //!
    //! \name Autodiff
    //!
    ///@{
    autodiff::dual1st eval( autodiff::dual1st const & x, autodiff::dual1st const & y ) const override
    {
      auto [i, j, dx, dy, DX, DY] = find_patch( x.val, y.val );

      autodiff::dual1st dx1st = x;
      dx1st.val               = dx;
      autodiff::dual1st dy1st = y;
      dy1st.val               = dy;

      autodiff::dual1st u  = dx1st / DX;
      autodiff::dual1st v  = dy1st / DY;
      autodiff::dual1st u1 = 1 - u;
      autodiff::dual1st v1 = 1 - v;

      Mat2x2 Z   = mZ.block<2, 2>( i, j ).matrix();
      auto   Z00 = Z.coeff( 0, 0 );
      auto   Z01 = Z.coeff( 0, 1 );
      auto   Z10 = Z.coeff( 1, 0 );
      auto   Z11 = Z.coeff( 1, 1 );
      return u1 * ( Z00 * v1 + Z01 * v ) + u * ( Z10 * v1 + Z11 * v );
    }

    autodiff::dual2nd eval( autodiff::dual2nd const & x, autodiff::dual2nd const & y ) const override
    {
      auto [i, j, dx, dy, DX, DY] = find_patch( x.val.val, y.val.val );

      autodiff::dual2nd dx2nd = x;
      dx2nd.val.val           = dx;
      autodiff::dual2nd dy2nd = y;
      dy2nd.val.val           = dy;

      autodiff::dual2nd u  = dx2nd / DX;
      autodiff::dual2nd v  = dy2nd / DY;
      autodiff::dual2nd u1 = 1 - u;
      autodiff::dual2nd v1 = 1 - v;

      Mat2x2 Z   = mZ.block<2, 2>( i, j ).matrix();
      auto   Z00 = Z.coeff( 0, 0 );
      auto   Z01 = Z.coeff( 0, 1 );
      auto   Z10 = Z.coeff( 1, 0 );
      auto   Z11 = Z.coeff( 1, 1 );
      return u1 * ( Z00 * v1 + Z01 * v ) + u * ( Z10 * v1 + Z11 * v );
    }

    template <typename T1, typename T2>
    autodiff::HigherOrderDual<autodiff::detail::DualOrder<T1, T2>::value, real_type> eval( T1 const & x, T2 const & y )
      const
    {
      autodiff::HigherOrderDual<autodiff::detail::DualOrder<T1, T2>::value, real_type> X{ x }, Y{ y };
      return eval( X, Y );
    }
///@}
#endif

    //! Write spline information to output stream
    void write_to_stream( ostream_type & s ) const override
    {
      fmt::print( s, "Nx={} Ny={}\n", m_nx, m_ny );
      for ( integer i = 1; i < m_nx; ++i )
      {
        for ( integer j = 1; j < m_ny; ++j )
        {
          fmt::print(
            s,
            "patch ({},{})\n"
            "  DX    = {:<12.4}  DY    = {:<12.4}\n"
            "  Z00   = {:<12.4}  Z10   = {:<12.4}\n"
            "  Z01   = {:<12.4}  Z11   = {:<12.4}\n",
            i,
            j,
            mX.coeff( i ) - mX.coeff( i - 1 ),
            mY.coeff( j ) - mY.coeff( j - 1 ),
            mZ.coeff( i - 1, j - 1 ),
            mZ.coeff( i, j - 1 ),
            mZ.coeff( i - 1, j ),
            mZ.coeff( i, j ) );
        }
      }
    }

    //! Return spline type name
    char const * type_name() const override { return "bilinear"; }
  };

}  // namespace Splines
