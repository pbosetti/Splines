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

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

#include "Splines.hh"
#include "Utils_fmt.hh"

#include <cmath>
#include <iomanip>
#include <set>

namespace Splines
{

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void BiCubicSplineBase::load( integer const i, integer const j, real_type bili3[4][4] ) const
  {
    //
    //  1    3
    //
    //  0    2
    //
    integer const i0{ ipos_C( i, j ) };
    integer const i1{ ipos_C( i, j + 1 ) };
    integer const i2{ ipos_C( i + 1, j ) };
    integer const i3{ ipos_C( i + 1, j + 1 ) };

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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type BiCubicSplineBase::eval( real_type const x, real_type const y ) const
  {
    real_type bili3[4][4], u[4], v[4];

    std::pair<integer, real_type> X( 0, x ), Y( 0, y );
    m_search_x.find( X );
    m_search_y.find( Y );

    integer const i{ X.first };
    integer const j{ Y.first };

    real_type const dx{ X.second - m_X[i] };
    real_type const dy{ Y.second - m_Y[j] };
    real_type const DX{ m_X[i + 1] - m_X[i] };
    real_type const DY{ m_Y[j + 1] - m_Y[j] };

    Hermite3( dx, DX, u );
    Hermite3( dy, DY, v );

    load( i, j, bili3 );

    return bilinear3( u, bili3, v );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type BiCubicSplineBase::Dx( real_type const x, real_type const y ) const
  {
    real_type bili3[4][4], u_D[4], v[4];

    std::pair<integer, real_type> X( 0, x ), Y( 0, y );
    m_search_x.find( X );
    m_search_y.find( Y );

    integer const i{ X.first };
    integer const j{ Y.first };

    real_type const dx{ X.second - m_X[i] };
    real_type const dy{ Y.second - m_Y[j] };
    real_type const DX{ m_X[i + 1] - m_X[i] };
    real_type const DY{ m_Y[j + 1] - m_Y[j] };

    Hermite3_D( dx, DX, u_D );
    Hermite3( dy, DY, v );

    load( i, j, bili3 );

    return bilinear3( u_D, bili3, v );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type BiCubicSplineBase::Dy( real_type const x, real_type const y ) const
  {
    real_type bili3[4][4], u[4], v_D[4];

    std::pair<integer, real_type> X( 0, x ), Y( 0, y );
    m_search_x.find( X );
    m_search_y.find( Y );

    integer const i{ X.first };
    integer const j{ Y.first };

    real_type const dx{ X.second - m_X[i] };
    real_type const dy{ Y.second - m_Y[j] };
    real_type const DX{ m_X[i + 1] - m_X[i] };
    real_type const DY{ m_Y[j + 1] - m_Y[j] };

    Hermite3( dx, DX, u );
    Hermite3_D( dy, DY, v_D );

    load( i, j, bili3 );

    return bilinear3( u, bili3, v_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type BiCubicSplineBase::Dxy( real_type const x, real_type const y ) const
  {
    real_type bili3[4][4], u_D[4], v_D[4];

    std::pair<integer, real_type> X( 0, x ), Y( 0, y );
    m_search_x.find( X );
    m_search_y.find( Y );

    integer const i{ X.first };
    integer const j{ Y.first };

    real_type const dx{ X.second - m_X[i] };
    real_type const dy{ Y.second - m_Y[j] };
    real_type const DX{ m_X[i + 1] - m_X[i] };
    real_type const DY{ m_Y[j + 1] - m_Y[j] };

    Hermite3_D( dx, DX, u_D );
    Hermite3_D( dy, DY, v_D );

    load( i, j, bili3 );

    return bilinear3( u_D, bili3, v_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type BiCubicSplineBase::Dxx( real_type const x, real_type const y ) const
  {
    real_type bili3[4][4], u_DD[4], v[4];

    std::pair<integer, real_type> X( 0, x ), Y( 0, y );
    m_search_x.find( X );
    m_search_y.find( Y );

    integer const i{ X.first };
    integer const j{ Y.first };

    real_type const dx{ X.second - m_X[i] };
    real_type const dy{ Y.second - m_Y[j] };
    real_type const DX{ m_X[i + 1] - m_X[i] };
    real_type const DY{ m_Y[j + 1] - m_Y[j] };

    Hermite3_DD( dx, DX, u_DD );
    Hermite3( dy, DY, v );

    load( i, j, bili3 );

    return bilinear3( u_DD, bili3, v );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type BiCubicSplineBase::Dyy( real_type const x, real_type const y ) const
  {
    real_type bili3[4][4], u[4], v_DD[4];

    std::pair<integer, real_type> X( 0, x ), Y( 0, y );
    m_search_x.find( X );
    m_search_y.find( Y );

    integer const i{ X.first };
    integer const j{ Y.first };

    real_type const dx{ X.second - m_X[i] };
    real_type const dy{ Y.second - m_Y[j] };
    real_type const DX{ m_X[i + 1] - m_X[i] };
    real_type const DY{ m_Y[j + 1] - m_Y[j] };

    Hermite3( dx, DX, u );
    Hermite3_DD( dy, DY, v_DD );

    load( i, j, bili3 );

    return bilinear3( u, bili3, v_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void BiCubicSplineBase::D( real_type const x, real_type const y, real_type d[3] ) const
  {
    real_type bili3[4][4], u[4], u_D[4], v[4], v_D[4];

    std::pair<integer, real_type> X( 0, x ), Y( 0, y );
    m_search_x.find( X );
    m_search_y.find( Y );

    integer const i{ X.first };
    integer const j{ Y.first };

    real_type const dx{ X.second - m_X[i] };
    real_type const dy{ Y.second - m_Y[j] };
    real_type const DX{ m_X[i + 1] - m_X[i] };
    real_type const DY{ m_Y[j + 1] - m_Y[j] };

    Hermite3( dx, DX, u );
    Hermite3_D( dx, DX, u_D );

    Hermite3( dy, DY, v );
    Hermite3_D( dy, DY, v_D );

    load( i, j, bili3 );

    d[0] = bilinear3( u, bili3, v );
    d[1] = bilinear3( u_D, bili3, v );
    d[2] = bilinear3( u, bili3, v_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void BiCubicSplineBase::DD( real_type const x, real_type const y, real_type d[6] ) const
  {
    real_type bili3[4][4], u[4], u_D[4], u_DD[4], v[4], v_D[4], v_DD[4];

    std::pair<integer, real_type> X( 0, x ), Y( 0, y );
    m_search_x.find( X );
    m_search_y.find( Y );

    integer const i{ X.first };
    integer const j{ Y.first };

    real_type const dx{ X.second - m_X[i] };
    real_type const dy{ Y.second - m_Y[j] };
    real_type const DX{ m_X[i + 1] - m_X[i] };
    real_type const DY{ m_Y[j + 1] - m_Y[j] };

    Hermite3( dx, DX, u );
    Hermite3_D( dx, DX, u_D );
    Hermite3_DD( dx, DX, u_DD );
    Hermite3( dy, DY, v );
    Hermite3_D( dy, DY, v_D );
    Hermite3_DD( dy, DY, v_DD );

    load( i, j, bili3 );

    d[0] = bilinear3( u, bili3, v );
    d[1] = bilinear3( u_D, bili3, v );
    d[2] = bilinear3( u, bili3, v_D );
    d[3] = bilinear3( u_DD, bili3, v );
    d[4] = bilinear3( u_D, bili3, v_D );
    d[5] = bilinear3( u, bili3, v_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void BiQuinticSplineBase::load( integer const i, integer const j, real_type bili5[6][6] ) const
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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type BiQuinticSplineBase::eval( real_type const x, real_type const y ) const
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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type BiQuinticSplineBase::Dx( real_type const x, real_type const y ) const
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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type BiQuinticSplineBase::Dy( real_type const x, real_type const y ) const
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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type BiQuinticSplineBase::Dxy( real_type const x, real_type const y ) const
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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type BiQuinticSplineBase::Dxx( real_type const x, real_type const y ) const
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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type BiQuinticSplineBase::Dyy( real_type const x, real_type const y ) const
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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void BiQuinticSplineBase::D( real_type const x, real_type const y, real_type d[3] ) const
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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void BiQuinticSplineBase::DD( real_type const x, real_type const y, real_type d[6] ) const
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

    d[0] = bilinear5( u, bili5, v );
    d[1] = bilinear5( u_D, bili5, v );
    d[2] = bilinear5( u, bili5, v_D );
    d[3] = bilinear5( u_DD, bili5, v );
    d[4] = bilinear5( u_D, bili5, v_D );
    d[5] = bilinear5( u, bili5, v_DD );
  }

}  // namespace Splines
