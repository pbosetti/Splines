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
    MatC mDX;
    MatC mDY;

    MatC mDXX;
    MatC mDYY;
    MatC mDXY;

    MatC mDXYY;
    MatC mDXXY;
    MatC mDXXYY;

    Spline_sub_type m_sub_type;

    using SplineSurf::m_nx;
    using SplineSurf::m_ny;

    using SplineSurf::mX;
    using SplineSurf::mY;
    using SplineSurf::mZ;

    real_type & Dx_node_ref( integer const i, integer const j ) { return mDX.coeffRef( i, j ); }
    real_type & Dy_node_ref( integer const i, integer const j ) { return mDY.coeffRef( i, j ); }

    real_type & Dxx_node_ref( integer const i, integer const j ) { return mDXX.coeffRef( i, j ); }
    real_type & Dyy_node_ref( integer const i, integer const j ) { return mDYY.coeffRef( i, j ); }
    real_type & Dxy_node_ref( integer const i, integer const j ) { return mDXY.coeffRef( i, j ); }

    real_type & Dxyy_node_ref( integer const i, integer const j ) { return mDXYY.coeffRef( i, j ); }
    real_type & Dxxy_node_ref( integer const i, integer const j ) { return mDXXY.coeffRef( i, j ); }

    real_type & Dxxyy_node_ref( integer const i, integer const j ) { return mDXXYY.coeffRef( i, j ); }

    void load( integer const i, integer const j, Mat6x6 & bili5 ) const
    {
      bili5.block<2, 2>( 0, 0 ) = mZ.block( i, j, 2, 2 );
      bili5.block<2, 2>( 2, 0 ) = mDX.block( i, j, 2, 2 );
      bili5.block<2, 2>( 4, 0 ) = mDXX.block( i, j, 2, 2 );

      bili5.block<2, 2>( 0, 2 ) = mDY.block( i, j, 2, 2 );
      bili5.block<2, 2>( 2, 2 ) = mDXY.block( i, j, 2, 2 );
      bili5.block<2, 2>( 4, 2 ) = mDXXY.block( i, j, 2, 2 );

      bili5.block<2, 2>( 0, 4 ) = mDYY.block( i, j, 2, 2 );
      bili5.block<2, 2>( 2, 4 ) = mDXYY.block( i, j, 2, 2 );
      bili5.block<2, 2>( 4, 4 ) = mDXXYY.block( i, j, 2, 2 );
    }

  public:
    using SplineSurf::eval;

    //! spline constructor
    explicit BiQuinticSplineBase( Spline_sub_type sub_type, string_view name = "Spline" )
      : SplineSurf( name ), m_sub_type( sub_type )
    {
    }

    ~BiQuinticSplineBase() override = default;

    //!
    //! \name Estimated derivatives at interpolation nodes
    //!
    ///@{

    //!
    //! Estimated `x` derivatives at node `(i,j)`
    //!
    real_type Dx_node( integer const i, integer const j ) const { return mDX.coeff( i, j ); }

    //!
    //! Estimated `y` derivatives at node `(i,j)`
    //!
    real_type Dy_node( integer const i, integer const j ) const { return mDY.coeff( i, j ); }

    //!
    //! Estimated second `x` second derivatives at node `(i,j)`
    //!
    real_type Dxx_node( integer const i, integer const j ) const { return mDXX.coeff( i, j ); }

    //!
    //! Estimated second`y` derivatives at node `(i,j)`
    //!
    real_type Dyy_node( integer const i, integer const j ) const { return mDYY.coeff( i, j ); }

    //!
    //! Estimated mixed `xy` derivatives at node `(i,j)`
    //!
    real_type Dxy_node( integer const i, integer const j ) const { return mDXY.coeff( i, j ); }

    //!
    //! Estimated `xxy` second derivatives at node `(i,j)`
    //!
    real_type Dxxy_node( integer const i, integer const j ) const { return mDXXY.coeff( i, j ); }

    //!
    //! Estimated `xyy` derivatives at node `(i,j)`
    //!
    real_type Dxyy_node( integer const i, integer const j ) const { return mDXYY.coeff( i, j ); }

    //!
    //! Estimated mixed `xxyy` derivatives at node `(i,j)`
    //!
    real_type Dxxyy_node( integer const i, integer const j ) const { return mDXXYY.coeff( i, j ); }

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
      Mat6x6 M;
      Vec6   u, v;

      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i = X.first;
      integer const j = Y.first;

      real_type const dx = X.second - mX.coeff( i );
      real_type const dy = Y.second - mY.coeff( j );
      real_type const DX = mX.coeff( i + 1 ) - mX.coeff( i );
      real_type const DY = mY.coeff( j + 1 ) - mY.coeff( j );

      Hermite5( dx, DX, u.data() );
      Hermite5( dy, DY, v.data() );

      load( i, j, M );

      return u.dot( M * v );
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
      Mat6x6 M;
      Vec6   u, u_D, v, v_D;

      auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

      Hermite5( dx, DX, u.data() );
      Hermite5_D( dx, DX, u_D.data() );
      Hermite5( dy, DY, v.data() );
      Hermite5_D( dy, DY, v_D.data() );

      load( i, j, M );

      auto Mv = M * v;

      d[0] = u.dot( Mv );
      d[1] = u_D.dot( Mv );
      d[2] = u.dot( M * v_D );
    }

    //!
    //! Evaluate spline `x`  derivative at point \f$ (x,y) \f$
    //!
    real_type Dx( real_type const x, real_type const y ) const override
    {
      Mat6x6 M;
      Vec6   u_D, v;

      auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

      Hermite5_D( dx, DX, u_D.data() );
      Hermite5( dy, DY, v.data() );

      load( i, j, M );

      return u_D.dot( M * v );
    }

    //!
    //! Evaluate spline `y`  derivative at point \f$ (x,y) \f$
    //!
    real_type Dy( real_type const x, real_type const y ) const override
    {
      Mat6x6 M;
      Vec6   u, v_D;

      auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

      Hermite5( dx, DX, u.data() );
      Hermite5_D( dy, DY, v_D.data() );

      load( i, j, M );

      return u.dot( M * v_D );
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
      Mat6x6 M;
      Vec6   u, u_D, u_DD, v, v_D, v_DD;

      auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

      Hermite5( dx, DX, u.data() );
      Hermite5_D( dx, DX, u_D.data() );
      Hermite5_DD( dx, DX, u_DD.data() );
      Hermite5( dy, DY, v.data() );
      Hermite5_D( dy, DY, v_D.data() );
      Hermite5_DD( dy, DY, v_DD.data() );

      load( i, j, M );

      auto Mv   = M * v;
      auto Mv_D = M * v_D;

      dd[0] = u.dot( Mv );
      dd[1] = u_D.dot( Mv );
      dd[2] = u.dot( Mv_D );
      dd[3] = u_DD.dot( Mv );
      dd[4] = u_D.dot( Mv_D );
      dd[5] = u.dot( M * v_DD );
    }

    //!
    //! Evaluate spline `x` second derivative at point \f$ (x,y) \f$
    //!
    real_type Dxx( real_type const x, real_type const y ) const override
    {
      Mat6x6 M;
      Vec6   u_DD, v;

      auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

      Hermite5_DD( dx, DX, u_DD.data() );
      Hermite5( dy, DY, v.data() );

      load( i, j, M );

      return u_DD.dot( M * v );
    }

    //!
    //! Evaluate spline `xy` mixed derivative at point \f$ (x,y) \f$
    //!
    real_type Dxy( real_type const x, real_type const y ) const override
    {
      Mat6x6 M;
      Vec6   u_D, v_D;

      auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

      Hermite5_D( dx, DX, u_D.data() );
      Hermite5_D( dy, DY, v_D.data() );

      load( i, j, M );

      return u_D.dot( M * v_D );
    }

    //!
    //! Evaluate spline `y` second derivative at point \f$ (x,y) \f$
    //!
    real_type Dyy( real_type const x, real_type const y ) const override
    {
      Mat6x6 M;
      Vec6   u, v_DD;

      auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

      Hermite5( dx, DX, u.data() );
      Hermite5_DD( dy, DY, v_DD.data() );

      load( i, j, M );

      return u.dot( M * v_DD );
    }

    ///@}

#ifdef AUTODIFF_SUPPORT
    //!
    //! \name Autodiff
    //!
    ///@{
    autodiff::dual1st eval( autodiff::dual1st const & x, autodiff::dual1st const & y ) const override
    {
      Mat6x6            M;
      autodiff::dual1st u[6], v[6], Mv[6];

      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i = X.first;
      integer const j = Y.first;

      autodiff::dual1st XX, YY;
      XX.val  = X.second;
      XX.grad = x.grad;
      YY.val  = Y.second;
      YY.grad = y.grad;

      autodiff::dual1st const dx = XX - mX.coeff( i );
      autodiff::dual1st const dy = YY - mY.coeff( j );
      real_type const         DX = mX.coeff( i + 1 ) - mX.coeff( i );
      real_type const         DY = mY.coeff( j + 1 ) - mY.coeff( j );

      Hermite5<autodiff::dual1st>( dx, DX, u );
      Hermite5<autodiff::dual1st>( dy, DY, v );

      load( i, j, M );

      Mv[0] = M.coeff( 0, 0 ) * v[0] + M.coeff( 0, 1 ) * v[1] + M.coeff( 0, 2 ) * v[2] + M.coeff( 0, 3 ) * v[3] +
              M.coeff( 0, 4 ) * v[4] + M.coeff( 0, 5 ) * v[5];

      Mv[1] = M.coeff( 1, 0 ) * v[0] + M.coeff( 1, 1 ) * v[1] + M.coeff( 1, 2 ) * v[2] + M.coeff( 1, 3 ) * v[3] +
              M.coeff( 1, 4 ) * v[4] + M.coeff( 1, 5 ) * v[5];

      Mv[2] = M.coeff( 2, 0 ) * v[0] + M.coeff( 2, 1 ) * v[1] + M.coeff( 2, 2 ) * v[2] + M.coeff( 2, 3 ) * v[3] +
              M.coeff( 2, 4 ) * v[4] + M.coeff( 2, 5 ) * v[5];

      Mv[3] = M.coeff( 3, 0 ) * v[0] + M.coeff( 3, 1 ) * v[1] + M.coeff( 3, 2 ) * v[2] + M.coeff( 3, 3 ) * v[3] +
              M.coeff( 3, 4 ) * v[4] + M.coeff( 3, 5 ) * v[5];

      Mv[4] = M.coeff( 4, 0 ) * v[0] + M.coeff( 4, 1 ) * v[1] + M.coeff( 4, 2 ) * v[2] + M.coeff( 4, 3 ) * v[3] +
              M.coeff( 4, 4 ) * v[4] + M.coeff( 4, 5 ) * v[5];

      Mv[5] = M.coeff( 5, 0 ) * v[0] + M.coeff( 5, 1 ) * v[1] + M.coeff( 5, 2 ) * v[2] + M.coeff( 5, 3 ) * v[3] +
              M.coeff( 5, 4 ) * v[4] + M.coeff( 5, 5 ) * v[5];

      return u[0] * Mv[0] + u[1] * Mv[1] + u[2] * Mv[2] + u[3] * Mv[3] + u[4] * Mv[4] + u[5] * Mv[5];
    }

    autodiff::dual2nd eval( autodiff::dual2nd const & x, autodiff::dual2nd const & y ) const override
    {
      Mat6x6            M;
      autodiff::dual2nd u[6], v[6], Mv[6];

      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i = X.first;
      integer const j = Y.first;

      autodiff::dual2nd XX, YY;

      XX.val.val   = X.second;
      XX.val.grad  = x.val.grad;
      XX.grad.val  = x.grad.val;
      XX.grad.grad = x.grad.grad;

      YY.val.val   = Y.second;
      YY.val.grad  = y.val.grad;
      YY.grad.val  = y.grad.val;
      YY.grad.grad = y.grad.grad;

      autodiff::dual2nd const dx = XX - mX.coeff( i );
      autodiff::dual2nd const dy = YY - mY.coeff( j );
      real_type const         DX = mX.coeff( i + 1 ) - mX.coeff( i );
      real_type const         DY = mY.coeff( j + 1 ) - mY.coeff( j );

      Hermite5<autodiff::dual2nd>( dx, DX, u );
      Hermite5<autodiff::dual2nd>( dy, DY, v );

      load( i, j, M );

      Mv[0] = M.coeff( 0, 0 ) * v[0] + M.coeff( 0, 1 ) * v[1] + M.coeff( 0, 2 ) * v[2] + M.coeff( 0, 3 ) * v[3] +
              M.coeff( 0, 4 ) * v[4] + M.coeff( 0, 5 ) * v[5];

      Mv[1] = M.coeff( 1, 0 ) * v[0] + M.coeff( 1, 1 ) * v[1] + M.coeff( 1, 2 ) * v[2] + M.coeff( 1, 3 ) * v[3] +
              M.coeff( 1, 4 ) * v[4] + M.coeff( 1, 5 ) * v[5];

      Mv[2] = M.coeff( 2, 0 ) * v[0] + M.coeff( 2, 1 ) * v[1] + M.coeff( 2, 2 ) * v[2] + M.coeff( 2, 3 ) * v[3] +
              M.coeff( 2, 4 ) * v[4] + M.coeff( 2, 5 ) * v[5];

      Mv[3] = M.coeff( 3, 0 ) * v[0] + M.coeff( 3, 1 ) * v[1] + M.coeff( 3, 2 ) * v[2] + M.coeff( 3, 3 ) * v[3] +
              M.coeff( 3, 4 ) * v[4] + M.coeff( 3, 5 ) * v[5];

      Mv[4] = M.coeff( 4, 0 ) * v[0] + M.coeff( 4, 1 ) * v[1] + M.coeff( 4, 2 ) * v[2] + M.coeff( 4, 3 ) * v[3] +
              M.coeff( 4, 4 ) * v[4] + M.coeff( 4, 5 ) * v[5];

      Mv[5] = M.coeff( 5, 0 ) * v[0] + M.coeff( 5, 1 ) * v[1] + M.coeff( 5, 2 ) * v[2] + M.coeff( 5, 3 ) * v[3] +
              M.coeff( 5, 4 ) * v[4] + M.coeff( 5, 5 ) * v[5];

      return u[0] * Mv[0] + u[1] * Mv[1] + u[2] * Mv[2] + u[3] * Mv[3] + u[4] * Mv[4] + u[5] * Mv[5];
    }

///@}
#endif
  };

}  // namespace Splines

#endif

//
// EOF: SplineBiQuintic.hxx
//
