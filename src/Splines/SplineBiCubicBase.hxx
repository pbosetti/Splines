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
    Malloc_real m_mem_bicubic;

    real_type * m_DX  = nullptr;
    real_type * m_DY  = nullptr;
    real_type * m_DXY = nullptr;

    using SplineSurf::m_nx;
    using SplineSurf::m_ny;

    using SplineSurf::m_X;
    using SplineSurf::m_Y;
    using SplineSurf::m_Z;

    void load( integer const i, integer const j, real_type bili3[4][4] ) const
    {
      //
      //  1    3
      //
      //  0    2
      //
      integer const i0 = ipos_C( i, j );
      integer const i1 = ipos_C( i, j + 1 );
      integer const i2 = ipos_C( i + 1, j );
      integer const i3 = ipos_C( i + 1, j + 1 );

      bili3[0][0] = m_Z[i0];
      bili3[0][1] = m_Z[i1];
      bili3[0][2] = m_DY[i0];
      bili3[0][3] = m_DY[i1];

      bili3[1][0] = m_Z[i2];
      bili3[1][1] = m_Z[i3];
      bili3[1][2] = m_DY[i2];
      bili3[1][3] = m_DY[i3];

      bili3[2][0] = m_DX[i0];
      bili3[2][1] = m_DX[i1];
      bili3[2][2] = m_DXY[i0];
      bili3[2][3] = m_DXY[i1];

      bili3[3][0] = m_DX[i2];
      bili3[3][1] = m_DX[i3];
      bili3[3][2] = m_DXY[i2];
      bili3[3][3] = m_DXY[i3];
    }


    real_type & Dx_node_ref( integer const i, integer const j ) { return m_DX[this->ipos_C( i, j )]; }
    real_type & Dy_node_ref( integer const i, integer const j ) { return m_DY[this->ipos_C( i, j )]; }
    real_type & Dxy_node_ref( integer const i, integer const j ) { return m_DXY[this->ipos_C( i, j )]; }

  public:
    using SplineSurf::eval;

    //! spline constructor
    explicit BiCubicSplineBase( string_view name = "BiCubicSplineBase" )
      : SplineSurf( name ), m_mem_bicubic( fmt::format( "BiCubicSplineBase[{}]", name ) )
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
    real_type Dx_node( integer const i, integer const j ) const { return m_DX[this->ipos_C( i, j )]; }

    //!
    //! Estimated `y` derivatives at node `(i,j)`
    //!
    real_type Dy_node( integer const i, integer const j ) const { return m_DY[this->ipos_C( i, j )]; }

    //!
    //! Estimated mixed `xy` derivatives at node `(i,j)`
    //!
    real_type Dxy_node( integer const i, integer const j ) const { return m_DXY[this->ipos_C( i, j )]; }

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
      real_type bili3[4][4], u[4], v[4];

      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i = X.first;
      integer const j = Y.first;

      real_type const dx = X.second - m_X[i];
      real_type const dy = Y.second - m_Y[j];
      real_type const DX = m_X[i + 1] - m_X[i];
      real_type const DY = m_Y[j + 1] - m_Y[j];

      Hermite3( dx, DX, u );
      Hermite3( dy, DY, v );

      load( i, j, bili3 );

      return bilinear3( u, bili3, v );
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
      real_type bili3[4][4], u[4], u_D[4], v[4], v_D[4];

      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i = X.first;
      integer const j = Y.first;

      real_type const dx = X.second - m_X[i];
      real_type const dy = Y.second - m_Y[j];
      real_type const DX = m_X[i + 1] - m_X[i];
      real_type const DY = m_Y[j + 1] - m_Y[j];

      Hermite3( dx, DX, u );
      Hermite3_D( dx, DX, u_D );

      Hermite3( dy, DY, v );
      Hermite3_D( dy, DY, v_D );

      load( i, j, bili3 );

      d[0] = bilinear3( u, bili3, v );
      d[1] = bilinear3( u_D, bili3, v );
      d[2] = bilinear3( u, bili3, v_D );
    }

    //!
    //! Evaluate spline `x`  derivative at point \f$ (x,y) \f$
    //!
    real_type Dx( real_type const x, real_type const y ) const override
    {
      real_type bili3[4][4], u_D[4], v[4];

      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i = X.first;
      integer const j = Y.first;

      real_type const dx = X.second - m_X[i];
      real_type const dy = Y.second - m_Y[j];
      real_type const DX = m_X[i + 1] - m_X[i];
      real_type const DY = m_Y[j + 1] - m_Y[j];

      Hermite3_D( dx, DX, u_D );
      Hermite3( dy, DY, v );

      load( i, j, bili3 );

      return bilinear3( u_D, bili3, v );
    }

    //!
    //! Evaluate spline `y`  derivative at point \f$ (x,y) \f$
    //!
    real_type Dy( real_type const x, real_type const y ) const override
    {
      real_type bili3[4][4], u[4], v_D[4];

      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i = X.first;
      integer const j = Y.first;

      real_type const dx = X.second - m_X[i];
      real_type const dy = Y.second - m_Y[j];
      real_type const DX = m_X[i + 1] - m_X[i];
      real_type const DY = m_Y[j + 1] - m_Y[j];

      Hermite3( dx, DX, u );
      Hermite3_D( dy, DY, v_D );

      load( i, j, bili3 );

      return bilinear3( u, bili3, v_D );
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
      real_type bili3[4][4], u[4], u_D[4], u_DD[4], v[4], v_D[4], v_DD[4];

      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i = X.first;
      integer const j = Y.first;

      real_type const dx = X.second - m_X[i];
      real_type const dy = Y.second - m_Y[j];
      real_type const DX = m_X[i + 1] - m_X[i];
      real_type const DY = m_Y[j + 1] - m_Y[j];

      Hermite3( dx, DX, u );
      Hermite3_D( dx, DX, u_D );
      Hermite3_DD( dx, DX, u_DD );

      Hermite3( dy, DY, v );
      Hermite3_D( dy, DY, v_D );
      Hermite3_DD( dy, DY, v_DD );

      load( i, j, bili3 );

      dd[0] = bilinear3( u, bili3, v );
      dd[1] = bilinear3( u_D, bili3, v );
      dd[2] = bilinear3( u, bili3, v_D );
      dd[3] = bilinear3( u_DD, bili3, v );
      dd[4] = bilinear3( u_D, bili3, v_D );
      dd[5] = bilinear3( u, bili3, v_DD );
    }

    //!
    //! Evaluate spline `x` second derivative at point \f$ (x,y) \f$
    //!
    real_type Dxx( real_type const x, real_type const y ) const override
    {
      real_type bili3[4][4], u_DD[4], v[4];

      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i = X.first;
      integer const j = Y.first;

      real_type const dx = X.second - m_X[i];
      real_type const dy = Y.second - m_Y[j];
      real_type const DX = m_X[i + 1] - m_X[i];
      real_type const DY = m_Y[j + 1] - m_Y[j];

      Hermite3_DD( dx, DX, u_DD );
      Hermite3( dy, DY, v );

      load( i, j, bili3 );

      return bilinear3( u_DD, bili3, v );
    }

    //!
    //! Evaluate spline `xy` mixed derivative at point \f$ (x,y) \f$
    //!
    real_type Dxy( real_type const x, real_type const y ) const override
    {
      real_type bili3[4][4], u_D[4], v_D[4];

      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i = X.first;
      integer const j = Y.first;

      real_type const dx = X.second - m_X[i];
      real_type const dy = Y.second - m_Y[j];
      real_type const DX = m_X[i + 1] - m_X[i];
      real_type const DY = m_Y[j + 1] - m_Y[j];

      Hermite3_D( dx, DX, u_D );
      Hermite3_D( dy, DY, v_D );

      load( i, j, bili3 );

      return bilinear3( u_D, bili3, v_D );
    }

    //!
    //! Evaluate spline `y` second derivative at point \f$ (x,y) \f$
    //!
    real_type Dyy( real_type const x, real_type const y ) const override
    {
      real_type bili3[4][4], u[4], v_DD[4];

      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i = X.first;
      integer const j = Y.first;

      real_type const dx = X.second - m_X[i];
      real_type const dy = Y.second - m_Y[j];
      real_type const DX = m_X[i + 1] - m_X[i];
      real_type const DY = m_Y[j + 1] - m_Y[j];

      Hermite3( dx, DX, u );
      Hermite3_DD( dy, DY, v_DD );

      load( i, j, bili3 );

      return bilinear3( u, bili3, v_DD );
    }
    ///@}
  };

}  // namespace Splines

#endif

//
// EOF: SplineBiCubicBase.hxx
//
