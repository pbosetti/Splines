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

/*\
 |   ____  _  ___        _       _   _      ____        _ _            ____
 |  | __ )(_)/ _ \ _   _(_)_ __ | |_(_) ___/ ___| _ __ | (_)_ __   ___| __ )  __ _ ___  ___
 |  |  _ \| | | | | | | | | '_ \| __| |/ __\___ \| '_ \| | | '_ \ / _ \  _ \ / _` / __|/ _ \
 |  | |_) | | |_| | |_| | | | | | |_| | (__ ___) | |_) | | | | | |  __/ |_) | (_| \__ \  __/
 |  |____/|_|\__\_\\__,_|_|_| |_|\__|_|\___|____/| .__/|_|_|_| |_|\___|____/ \__,_|___/\___|
 |                                               |_|
\*/

#pragma once

#ifndef SPLINE_BIQUINTIC_SPLINE_BASE_HXX
#define SPLINE_BIQUINTIC_SPLINE_BASE_HXX

namespace Splines
{

  //! Bi-quintic spline base class
  class BiQuinticSplineBase : public SplineSurf
  {
  protected:
    Malloc_real m_mem_biquintic{ "BiQuinticSplineBase" };

    real_type * m_DX = nullptr;
    real_type * m_DY = nullptr;

    real_type * m_DXX = nullptr;
    real_type * m_DYY = nullptr;
    real_type * m_DXY = nullptr;

    real_type * m_DXYY  = nullptr;
    real_type * m_DXXY  = nullptr;
    real_type * m_DXXYY = nullptr;

    Spline_sub_type m_sub_type;

    using SplineSurf::m_nx;
    using SplineSurf::m_ny;

    using SplineSurf::m_X;
    using SplineSurf::m_Y;
    using SplineSurf::m_Z;

    real_type & Dx_node_ref( integer const i, integer const j ) { return m_DX[ipos_C( i, j )]; }
    real_type & Dy_node_ref( integer const i, integer const j ) { return m_DY[ipos_C( i, j )]; }

    real_type & Dxx_node_ref( integer const i, integer const j ) { return m_DXX[ipos_C( i, j )]; }
    real_type & Dyy_node_ref( integer const i, integer const j ) { return m_DYY[ipos_C( i, j )]; }
    real_type & Dxy_node_ref( integer const i, integer const j ) { return m_DXY[ipos_C( i, j )]; }

    real_type & Dxyy_node_ref( integer const i, integer const j ) { return m_DXYY[ipos_C( i, j )]; }
    real_type & Dxxy_node_ref( integer const i, integer const j ) { return m_DXXY[ipos_C( i, j )]; }
    real_type & Dxxyy_node_ref( integer const i, integer const j ) { return m_DXXYY[ipos_C( i, j )]; }

    void load( integer const i, integer const j, real_type bili5[6][6] ) const
    {
      integer const i00{ ipos_C( i, j ) };
      integer const i01{ ipos_C( i, j + 1 ) };
      integer const i10{ ipos_C( i + 1, j ) };
      integer const i11{ ipos_C( i + 1, j + 1 ) };

      //
      //  1    3
      //
      //  0    2
      //
      // H0, H1, dH0, dH1, ddH0, ddH1

      // + + + + + +
      // + + + + + +
      // + + + + . .
      // + + + + . .
      // + + . . . .
      // + + . . . .

      // P
      bili5[0][0] = m_Z[i00];
      bili5[0][1] = m_Z[i01];
      bili5[1][0] = m_Z[i10];
      bili5[1][1] = m_Z[i11];

      // DX
      bili5[2][0] = m_DX[i00];
      bili5[2][1] = m_DX[i01];
      bili5[3][0] = m_DX[i10];
      bili5[3][1] = m_DX[i11];

      // DXX
      bili5[4][0] = m_DXX[i00];
      bili5[4][1] = m_DXX[i01];
      bili5[5][0] = m_DXX[i10];
      bili5[5][1] = m_DXX[i11];

      // DY
      bili5[0][2] = m_DY[i00];
      bili5[0][3] = m_DY[i01];
      bili5[1][2] = m_DY[i10];
      bili5[1][3] = m_DY[i11];

      // DYY
      bili5[0][4] = m_DYY[i00];
      bili5[0][5] = m_DYY[i01];
      bili5[1][4] = m_DYY[i10];
      bili5[1][5] = m_DYY[i11];

      // DXY
      bili5[2][2] = m_DXY[i00];
      bili5[2][3] = m_DXY[i01];
      bili5[3][2] = m_DXY[i10];
      bili5[3][3] = m_DXY[i11];

      // DXXY
      bili5[4][2] = m_DXXY[i00];
      bili5[4][3] = m_DXXY[i01];
      bili5[5][2] = m_DXXY[i10];
      bili5[5][3] = m_DXXY[i11];

      // DXYY
      bili5[2][4] = m_DXYY[i00];
      bili5[2][5] = m_DXYY[i01];
      bili5[3][4] = m_DXYY[i10];
      bili5[3][5] = m_DXYY[i11];

      // DXXYY
      bili5[4][4] = m_DXXYY[i00];
      bili5[4][5] = m_DXXYY[i01];
      bili5[5][4] = m_DXXYY[i10];
      bili5[5][5] = m_DXXYY[i11];
    }

  public:
    using SplineSurf::eval;

    //! spline constructor
    explicit BiQuinticSplineBase( Spline_sub_type sub_type, string_view name = "Spline" )
      : SplineSurf( name ), m_sub_type( sub_type )
    {
    }

    ~BiQuinticSplineBase() override { m_mem_biquintic.free(); }

    //!
    //! \name Estimated derivatives at interpolation nodes
    //!
    ///@{

    //!
    //! Estimated `x` derivatives at node `(i,j)`
    //!
    real_type Dx_node( integer const i, integer const j ) const { return m_DX[ipos_C( i, j )]; }

    //!
    //! Estimated `y` derivatives at node `(i,j)`
    //!
    real_type Dy_node( integer const i, integer const j ) const { return m_DY[ipos_C( i, j )]; }

    //!
    //! Estimated second `x` second derivatives at node `(i,j)`
    //!
    real_type Dxx_node( integer const i, integer const j ) const { return m_DXX[ipos_C( i, j )]; }

    //!
    //! Estimated second`y` derivatives at node `(i,j)`
    //!
    real_type Dyy_node( integer const i, integer const j ) const { return m_DYY[ipos_C( i, j )]; }

    //!
    //! Estimated mixed `xy` derivatives at node `(i,j)`
    //!
    real_type Dxy_node( integer const i, integer const j ) const { return m_DXY[ipos_C( i, j )]; }

    //!
    //! Estimated `xxy` second derivatives at node `(i,j)`
    //!
    real_type Dxxy_node( integer const i, integer const j ) const { return m_DXXY[ipos_C( i, j )]; }

    //!
    //! Estimated `xyy` derivatives at node `(i,j)`
    //!
    real_type Dxyy_node( integer const i, integer const j ) const { return m_DXYY[ipos_C( i, j )]; }

    //!
    //! Estimated mixed `xxyy` derivatives at node `(i,j)`
    //!
    real_type Dxxyy_node( integer const i, integer const j ) const { return m_DXXYY[ipos_C( i, j )]; }

    ///@}

    //!
    //! \name Evaluate
    //!
    ///@{
    //!
    //! Evaluate spline at point \f$ (x,y) \f$
    //!
    real_type eval( real_type x, real_type y ) const override
    {
      real_type bili5[6][6], u[6], v[6];

      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i{ X.first };
      integer const j{ Y.first };

      real_type const dx{ X.second - m_X[i] };
      real_type const dy{ Y.second - m_Y[j] };
      real_type const DX{ m_X[i + 1] - m_X[i] };
      real_type const DY{ m_Y[j + 1] - m_Y[j] };

      Hermite5( dx, DX, u );
      Hermite5( dy, DY, v );

      load( i, j, bili5 );

      return bilinear5( u, bili5, v );
    }


    //!
    //! Evaluate spline with derivative at point \f$ (x,y) \f$
    //!
    //! - `d[0]` the value of the spline
    //! - `d[1]` the value of the spline `x` derivative
    //! - `d[2]` the value of the spline `y` derivative
    //!
    void D( real_type const x, real_type const y, real_type d[3] ) const override
    {
      real_type bili5[6][6], u[6], u_D[6], v[6], v_D[6];

      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i{ X.first };
      integer const j{ Y.first };

      real_type const dx{ X.second - m_X[i] };
      real_type const dy{ Y.second - m_Y[j] };
      real_type const DX{ m_X[i + 1] - m_X[i] };
      real_type const DY{ m_Y[j + 1] - m_Y[j] };

      Hermite5( dx, DX, u );
      Hermite5_D( dx, DX, u_D );
      Hermite5( dy, DY, v );
      Hermite5_D( dy, DY, v_D );

      load( i, j, bili5 );

      d[0] = bilinear5( u, bili5, v );
      d[1] = bilinear5( u_D, bili5, v );
      d[2] = bilinear5( u, bili5, v_D );
    }

    //!
    //! Evaluate spline `x`  derivative at point \f$ (x,y) \f$
    //!
    real_type Dx( real_type const x, real_type const y ) const override
    {
      real_type bili5[6][6], u_D[6], v[6];

      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i{ X.first };
      integer const j{ Y.first };

      real_type const dx{ X.second - m_X[i] };
      real_type const dy{ Y.second - m_Y[j] };
      real_type const DX{ m_X[i + 1] - m_X[i] };
      real_type const DY{ m_Y[j + 1] - m_Y[j] };

      Hermite5_D( dx, DX, u_D );
      Hermite5( dy, DY, v );

      load( i, j, bili5 );

      return bilinear5( u_D, bili5, v );
    }

    //!
    //! Evaluate spline `y`  derivative at point \f$ (x,y) \f$
    //!
    real_type Dy( real_type const x, real_type const y ) const override
    {
      real_type bili5[6][6], u[6], v_D[6];

      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i{ X.first };
      integer const j{ Y.first };

      real_type const dx{ X.second - m_X[i] };
      real_type const dy{ Y.second - m_Y[j] };
      real_type const DX{ m_X[i + 1] - m_X[i] };
      real_type const DY{ m_Y[j + 1] - m_Y[j] };

      Hermite5( dx, DX, u );
      Hermite5_D( dy, DY, v_D );

      load( i, j, bili5 );

      return bilinear5( u, bili5, v_D );
    }

    //!
    //! Evaluate spline with derivative at point \f$ (x,y) \f$
    //!
    //! - `d[0]` the value of the spline
    //! - `d[1]` the value of the spline `x` derivative
    //! - `d[2]` the value of the spline `y` derivative
    //! - `d[3]` the value of the spline `x` second derivative
    //! - `d[4]` the value of the spline `y` second derivative
    //! - `d[5]` the value of the spline `xy` mixed derivative
    //!
    void DD( real_type const x, real_type const y, real_type dd[6] ) const override
    {
      real_type bili5[6][6], u[6], u_D[6], u_DD[6], v[6], v_D[6], v_DD[6];

      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i{ X.first };
      integer const j{ Y.first };

      real_type const dx{ X.second - m_X[i] };
      real_type const dy{ Y.second - m_Y[j] };
      real_type const DX{ m_X[i + 1] - m_X[i] };
      real_type const DY{ m_Y[j + 1] - m_Y[j] };

      Hermite5( dx, DX, u );
      Hermite5_D( dx, DX, u_D );
      Hermite5_DD( dx, DX, u_DD );
      Hermite5( dy, DY, v );
      Hermite5_D( dy, DY, v_D );
      Hermite5_DD( dy, DY, v_DD );

      load( i, j, bili5 );

      dd[0] = bilinear5( u, bili5, v );
      dd[1] = bilinear5( u_D, bili5, v );
      dd[2] = bilinear5( u, bili5, v_D );
      dd[3] = bilinear5( u_DD, bili5, v );
      dd[4] = bilinear5( u_D, bili5, v_D );
      dd[5] = bilinear5( u, bili5, v_DD );
    }

    //!
    //! Evaluate spline `x` second derivative at point \f$ (x,y) \f$
    //!
    real_type Dxx( real_type const x, real_type const y ) const override
    {
      real_type bili5[6][6], u_DD[6], v[6];

      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i{ X.first };
      integer const j{ Y.first };

      real_type const dx{ X.second - m_X[i] };
      real_type const dy{ Y.second - m_Y[j] };
      real_type const DX{ m_X[i + 1] - m_X[i] };
      real_type const DY{ m_Y[j + 1] - m_Y[j] };

      Hermite5_DD( dx, DX, u_DD );
      Hermite5( dy, DY, v );

      load( i, j, bili5 );

      return bilinear5( u_DD, bili5, v );
    }

    //!
    //! Evaluate spline `xy` mixed derivative at point \f$ (x,y) \f$
    //!
    real_type Dxy( real_type const x, real_type const y ) const override
    {
      real_type bili5[6][6], u_D[6], v_D[6];

      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i{ X.first };
      integer const j{ Y.first };

      real_type const dx{ X.second - m_X[i] };
      real_type const dy{ Y.second - m_Y[j] };
      real_type const DX{ m_X[i + 1] - m_X[i] };
      real_type const DY{ m_Y[j + 1] - m_Y[j] };

      Hermite5_D( dx, DX, u_D );
      Hermite5_D( dy, DY, v_D );

      load( i, j, bili5 );

      return bilinear5( u_D, bili5, v_D );
    }

    //!
    //! Evaluate spline `y` second derivative at point \f$ (x,y) \f$
    //!
    real_type Dyy( real_type const x, real_type const y ) const override
    {
      real_type bili5[6][6], u[6], v_DD[6];

      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i{ X.first };
      integer const j{ Y.first };

      real_type const dx{ X.second - m_X[i] };
      real_type const dy{ Y.second - m_Y[j] };
      real_type const DX{ m_X[i + 1] - m_X[i] };
      real_type const DY{ m_Y[j + 1] - m_Y[j] };

      Hermite5( dx, DX, u );
      Hermite5_DD( dy, DY, v_DD );

      load( i, j, bili5 );

      return bilinear5( u, bili5, v_DD );
    }

    ///@}

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
      real_type dx  = val( x.grad );
      real_type dy  = val( y.grad );
      real_type ddx = x.grad.grad;
      real_type ddy = y.grad.grad;

      DD( val( x ), val( y ), dd );

      dual2nd res{ dd[0] };
      res.grad      = dd[1] * dx + dd[2] * dy;
      res.grad.grad = dd[3] * dx * dx + 2 * dx * dy * dd[4] + dy * dy * dd[5] + ddx * dd[1] + ddy * dd[2];
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
  };

}  // namespace Splines

#endif

//
// EOF: SplineBiQuintic.hxx
//
