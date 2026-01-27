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
    using SplineSurf::ipos_C;  // Accesso esplicito per chiarezza
    using SplineSurf::m_nx;
    using SplineSurf::m_ny;
    using SplineSurf::m_X;
    using SplineSurf::m_Y;
    using SplineSurf::m_Z;

    void make_spline() override
    {
      m_search_x.must_reset();
      m_search_y.must_reset();
    }

    // Helper per ottenere i 4 valori Z della cella corrente
    // Ottimizzabile ulteriormente se si conosce lo stride di memoria di m_Z
    inline void get_Z_patch( integer i, integer j, real_type & Z00, real_type & Z10, real_type & Z01, real_type & Z11 )
      const
    {
      Z00 = m_Z[ipos_C( i, j )];
      Z10 = m_Z[ipos_C( i + 1, j )];
      Z01 = m_Z[ipos_C( i, j + 1 )];
      Z11 = m_Z[ipos_C( i + 1, j + 1 )];
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
    ~BilinearSpline() override {}

    //! Evaluate the spline at (x,y)
    real_type eval( real_type const x, real_type const y ) const override
    {
      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i = X.first;
      integer const j = Y.first;

      // Ottimizzazione: Moltiplicazione inversa invece di divisione
      real_type const DX = m_X[i + 1] - m_X[i];
      real_type const DY = m_Y[j + 1] - m_Y[j];
      real_type const u  = ( X.second - m_X[i] ) / DX;
      real_type const v  = ( Y.second - m_Y[j] ) / DY;

      real_type const u1 = 1 - u;
      real_type const v1 = 1 - v;

      real_type Z00, Z10, Z01, Z11;
      get_Z_patch( i, j, Z00, Z10, Z01, Z11 );

      return u1 * ( Z00 * v1 + Z01 * v ) + u * ( Z10 * v1 + Z11 * v );
    }

    //! Compute x-derivative at (x,y)
    real_type Dx( real_type const x, real_type const y ) const override
    {
      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i = X.first;
      integer const j = Y.first;

      real_type const DX = m_X[i + 1] - m_X[i];
      real_type const DY = m_Y[j + 1] - m_Y[j];
      real_type const v  = ( Y.second - m_Y[j] ) / DY;

      real_type Z00, Z10, Z01, Z11;
      get_Z_patch( i, j, Z00, Z10, Z01, Z11 );

      return ( ( Z10 - Z00 ) * ( 1 - v ) + ( Z11 - Z01 ) * v ) / DX;
    }

    //! Compute y-derivative at (x,y)
    real_type Dy( real_type const x, real_type const y ) const override
    {
      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i = X.first;
      integer const j = Y.first;

      real_type const DX = m_X[i + 1] - m_X[i];
      real_type const DY = m_Y[j + 1] - m_Y[j];
      real_type const u  = ( X.second - m_X[i] ) / DX;

      real_type Z00, Z10, Z01, Z11;
      get_Z_patch( i, j, Z00, Z10, Z01, Z11 );

      return ( ( Z01 - Z00 ) * ( 1 - u ) + ( Z11 - Z10 ) * u ) / DY;
    }

    //! Compute value and first derivatives at (x,y)
    void D( real_type const x, real_type const y, real_type d[3] ) const override
    {
      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i = X.first;
      integer const j = Y.first;

      real_type const DX = m_X[i + 1] - m_X[i];
      real_type const DY = m_Y[j + 1] - m_Y[j];

      // Ottimizzazione: Calcolo inversi
      real_type const invDX = 1.0 / DX;
      real_type const invDY = 1.0 / DY;

      real_type const u  = ( X.second - m_X[i] ) * invDX;
      real_type const v  = ( Y.second - m_Y[j] ) * invDY;
      real_type const u1 = 1 - u;
      real_type const v1 = 1 - v;

      real_type Z00, Z10, Z01, Z11;
      get_Z_patch( i, j, Z00, Z10, Z01, Z11 );

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

      // Recupero indici per il calcolo del termine misto
      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );
      integer const i = X.first;
      integer const j = Y.first;

      real_type const DX = m_X[i + 1] - m_X[i];
      real_type const DY = m_Y[j + 1] - m_Y[j];

      real_type Z00, Z10, Z01, Z11;
      get_Z_patch( i, j, Z00, Z10, Z01, Z11 );

      dd[3] = 0;                                        // Dxx
      dd[4] = ( Z11 - Z10 - Z01 + Z00 ) / ( DX * DY );  // Dxy
      dd[5] = 0;                                        // Dyy
    }

    //! Second derivatives
    real_type Dxx( real_type const, real_type const ) const override { return 0; }

    // FIX MATEMATICO
    real_type Dxy( real_type const x, real_type const y ) const override
    {
      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );
      integer const i = X.first;
      integer const j = Y.first;

      real_type const DX = m_X[i + 1] - m_X[i];
      real_type const DY = m_Y[j + 1] - m_Y[j];

      real_type Z00, Z10, Z01, Z11;
      get_Z_patch( i, j, Z00, Z10, Z01, Z11 );

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
      using autodiff::dual1st;
      using autodiff::detail::val;

      real_type dd[3];
      D( val( x ), val( y ), dd );

      dual1st res{ dd[0] };
      res.grad = dd[1] * x.grad + dd[2] * y.grad;

      return res;
    }

    autodiff::dual2nd eval( autodiff::dual2nd const & x, autodiff::dual2nd const & y ) const override
    {
      using autodiff::derivative;
      using autodiff::dual2nd;

      real_type dd[6];
      DD( val( x ), val( y ), dd );

      dual2nd res{ dd[0] };
      res.grad      = dd[1] * x.grad + dd[2] * y.grad;
      res.grad.grad = 2 * dd[4] * val( x.grad ) * val( y.grad ) + dd[1] * x.grad.grad + dd[2] * y.grad.grad;

      return res;
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
          integer const i00{ ipos_C( i - 1, j - 1 ) };
          integer const i10{ ipos_C( i, j - 1 ) };
          integer const i01{ ipos_C( i - 1, j ) };
          integer const i11{ ipos_C( i, j ) };
          fmt::print(
            s,
            "patch ({},{})\n"
            "  DX    = {:<12.4}  DY    = {:<12.4}\n"
            "  Z00   = {:<12.4}  Z10   = {:<12.4}\n"
            "  Z01   = {:<12.4}  Z11   = {:<12.4}\n",
            i,
            j,
            m_X[i] - m_X[i - 1],
            m_Y[j] - m_Y[j - 1],
            m_Z[i00],
            m_Z[i10],
            m_Z[i01],
            m_Z[i11] );
        }
      }
    }

    //! Return spline type name
    char const * type_name() const override { return "bilinear"; }
  };

}  // namespace Splines
