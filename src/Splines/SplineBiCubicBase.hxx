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

//
// FILE: SplineBiCubicBase.hxx
//

#pragma once

#ifndef SPLINE_BICUBIC_BASE_HXX
#define SPLINE_BICUBIC_BASE_HXX

namespace Splines
{

  /*\
   |   ____  _  ____      _     _      ____        _ _            ____
   |  | __ )(_)/ ___|   _| |__ (_) ___/ ___| _ __ | (_)_ __   ___| __ )  __ _ ___  ___
   |  |  _ \| | |  | | | | '_ \| |/ __\___ \| '_ \| | | '_ \ / _ \  _ \ / _` / __|/ _ \
   |  | |_) | | |__| |_| | |_) | | (__ ___) | |_) | | | | | |  __/ |_) | (_| \__ \  __/
   |  |____/|_|\____\__,_|_.__/|_|\___|____/| .__/|_|_|_| |_|\___|____/ \__,_|___/\___|
   |                                        |_|
  \*/

  //!
  //! Bi-cubic spline base class
  //!
  class BiCubicSplineBase : public SplineSurf
  {
  protected:
    MatC mDX;
    MatC mDY;
    MatC mDXY;

    Spline_sub_type m_sub_type;

    using SplineSurf::m_nx;
    using SplineSurf::m_ny;

    using SplineSurf::mX;
    using SplineSurf::mY;
    using SplineSurf::mZ;

    void load( integer const i, integer const j, Mat4x4 & bili3 ) const
    {
      // Copia a blocchi
      // Invece di copiare scalarmente, copiamo 4 blocchi 2x2.
      // Eigen ottimizzerà queste operazioni usando istruzioni SIMD.

      // Quadrante in alto a sinistra: Z
      bili3.topLeftCorner<2, 2>() = mZ.block<2, 2>( i, j );

      // Quadrante in alto a destra: DY
      bili3.topRightCorner<2, 2>() = mDY.block<2, 2>( i, j );

      // Quadrante in basso a sinistra: DX
      bili3.bottomLeftCorner<2, 2>() = mDX.block<2, 2>( i, j );

      // Quadrante in basso a destra: DXY
      bili3.bottomRightCorner<2, 2>() = mDXY.block<2, 2>( i, j );
    }

    real_type & Dx_node_ref( integer const i, integer const j ) { return mDX.coeffRef( i, j ); }
    real_type & Dy_node_ref( integer const i, integer const j ) { return mDY.coeffRef( i, j ); }
    real_type & Dxy_node_ref( integer const i, integer const j ) { return mDXY.coeffRef( i, j ); }

  public:
    using SplineSurf::eval;

    //! spline constructor
    explicit BiCubicSplineBase(
      Spline_sub_type sub_type = Spline_sub_type::PCHIP,
      string_view     name     = "BiCubicSplineBase" )
      : SplineSurf( name ), m_sub_type( sub_type )
    {
    }

    ~BiCubicSplineBase() override {}

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
    //! Estimated mixed `xy` derivatives at node `(i,j)`
    //!
    real_type Dxy_node( integer const i, integer const j ) const { return mDXY.coeff( i, j ); }

    ///@}

    //!
    //! \name Evaluate
    //!
    ///@{
    //!
    //! Evaluate spline at point \f$ (x,y) \f$
    //!
    real_type eval( real_type const x, real_type const y ) const override
    {
      Mat4x4 bili3;
      Vec4   u, v;

      auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

      Hermite3( dx, DX, u.data() );
      Hermite3( dy, DY, v.data() );

      load( i, j, bili3 );

      return u.dot( bili3 * v );
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
      Mat4x4 bili3;
      Vec4   u, u_D, v, v_D;

      auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

      Hermite3( dx, DX, u.data() );
      Hermite3_D( dx, DX, u_D.data() );

      Hermite3( dy, DY, v.data() );
      Hermite3_D( dy, DY, v_D.data() );

      load( i, j, bili3 );
      
      auto Bv = bili3 * v;

      d[0] = u.dot( Bv );
      d[1] = u_D.dot( Bv );
      d[2] = u.dot( bili3 * v_D );
    }

    //!
    //! Evaluate spline `x`  derivative at point \f$ (x,y) \f$
    //!
    real_type Dx( real_type const x, real_type const y ) const override
    {
      Mat4x4 bili3;
      Vec4   u_D, v;

      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

      Hermite3_D( dx, DX, u_D.data() );
      Hermite3( dy, DY, v.data() );

      load( i, j, bili3 );

      return u_D.dot( bili3 * v );
    }

    //!
    //! Evaluate spline `y`  derivative at point \f$ (x,y) \f$
    //!
    real_type Dy( real_type const x, real_type const y ) const override
    {
      Mat4x4 bili3;
      Vec4   u, v_D;

      auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

      Hermite3( dx, DX, u.data() );
      Hermite3_D( dy, DY, v_D.data() );

      load( i, j, bili3 );

      return u.dot( bili3 * v_D );
    }

    //!
    //! Evaluate spline with derivative at point \f$ (x,y) \f$
    //!
    //! - `d[0]` the value of the spline
    //! - `d[1]` the value of the spline `x` derivative
    //! - `d[2]` the value of the spline `y` derivative
    //! - `d[3]` the value of the spline `x` second derivative
    //! - `d[4]` the value of the spline `xy` mixed derivative
    //! - `d[5]` the value of the spline `y` second derivative
    //!
    void DD( real_type const x, real_type const y, real_type dd[6] ) const override
    {
      Mat4x4 bili3;
      Vec4   u, u_D, u_DD, v, v_D, v_DD;

      auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

      Hermite3( dx, DX, u.data() );
      Hermite3_D( dx, DX, u_D.data() );
      Hermite3_DD( dx, DX, u_DD.data() );

      Hermite3( dy, DY, v.data() );
      Hermite3_D( dy, DY, v_D.data() );
      Hermite3_DD( dy, DY, v_DD.data() );

      load( i, j, bili3 );
      
      auto Bv   = bili3 * v;
      auto Bv_D = bili3 * v_D;

      dd[0] = u.dot( Bv );
      dd[1] = u_D.dot( Bv );
      dd[2] = u.dot( Bv_D );
      dd[3] = u_DD.dot( Bv );
      dd[4] = u_D.dot( Bv_D );
      dd[5] = u.dot( bili3 * v_DD );
    }

    //!
    //! Evaluate spline `x` second derivative at point \f$ (x,y) \f$
    //!
    real_type Dxx( real_type const x, real_type const y ) const override
    {
      Mat4x4 bili3;
      Vec4   u_DD, v;

      auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

      Hermite3_DD( dx, DX, u_DD.data() );
      Hermite3( dy, DY, v.data() );

      load( i, j, bili3 );

      return u_DD.dot( bili3 * v );
    }

    //!
    //! Evaluate spline `xy` mixed derivative at point \f$ (x,y) \f$
    //!
    real_type Dxy( real_type const x, real_type const y ) const override
    {
      Mat4x4 bili3;
      Vec4   u_D, v_D;

      auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

      Hermite3_D( dx, DX, u_D.data() );
      Hermite3_D( dy, DY, v_D.data() );

      load( i, j, bili3 );

      return u_D.dot( bili3 * v_D );
    }

    //!
    //! Evaluate spline `y` second derivative at point \f$ (x,y) \f$
    //!
    real_type Dyy( real_type const x, real_type const y ) const override
    {
      Mat4x4 bili3;
      Vec4   u, v_DD;

      auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

      Hermite3( dx, DX, u.data() );
      Hermite3_DD( dy, DY, v_DD.data() );

      load( i, j, bili3 );

      return u.dot( bili3 * v_DD );
    }
    ///@}
  };

}  // namespace Splines

#endif

//
// EOF: SplineBiCubicBase.hxx
//
