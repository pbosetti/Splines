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
#endif

#include "Splines.hh"
#include "SplinesUtils.hh"
#include "Utils_fmt.hh"

#include <cmath>
#include <set>

using namespace std;  // load standard namspace

namespace Splines
{

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*\
     Sistema lineare da risolvere

     D U
     L D U
       L D U
         L D U
           .....
              L D U
                L D U
                  L D

  \*/

  static void QuinticSpline_Yppp_continuous(
    real_type const X[],
    real_type const Y[],
    real_type const Yp[],
    real_type       Ypp[],
    integer const   npts,
    bool const      setbc )
  {
    UTILS_ASSERT( npts >= 2, "QuinticSpline_Yppp_continuous, npts={} must be >= 2\n", npts );

    integer const n{ npts - 1 };

    Malloc_real mem( "QuinticSpline_Yppp_continuous" );
    mem.allocate( 3 * ( n + 1 ) );
    real_type * L{ mem( n + 1 ) };
    real_type * D{ mem( n + 1 ) };
    real_type * U{ mem( n + 1 ) };
    real_type * Z{ Ypp };

    for ( integer i = 1; i < n; ++i )
    {
      real_type const hL{ X[i] - X[i - 1] };
      real_type const hL2{ hL * hL };
      real_type const hL3{ hL * hL2 };
      real_type const hR{ X[i + 1] - X[i] };
      real_type const hR2{ hR * hR };
      real_type const hR3{ hR * hR2 };
      real_type const DL{ 60 * ( Y[i] - Y[i - 1] ) / hL3 };
      real_type const DR{ 60 * ( Y[i + 1] - Y[i] ) / hR3 };
      real_type const DDL{ ( 36 * Yp[i] + 24 * Yp[i - 1] ) / hL2 };
      real_type const DDR{ ( 36 * Yp[i] + 24 * Yp[i + 1] ) / hR2 };
      L[i] = -3 / hL;
      D[i] = 9 / hL + 9 / hR;
      U[i] = -3 / hR;
      Z[i] = DR - DL + DDL - DDR;
    }
    L[0] = U[0] = 0;
    D[0]        = 1;
    L[n] = U[n] = 0;
    D[n]        = 1;
    if ( setbc )
    {
      {
        real_type const hL{ X[1] - X[0] };
        real_type const hR{ X[2] - X[1] };
        real_type const SL{ ( Y[1] - Y[0] ) / hL };
        real_type const SR{ ( Y[2] - Y[1] ) / hR };
        real_type const dp0{ Yp[1] };
        real_type const dpL{ Yp[0] };
        real_type const dpR{ Yp[2] };
        Z[0] = second_deriv3p_L( SL, hL, SR, hR, dpL, dp0, dpR );
      }
      {
        real_type const hL{ X[n - 1] - X[n - 2] };
        real_type const hR{ X[n] - X[n - 1] };
        real_type const SL{ ( Y[n - 1] - Y[n - 2] ) / hL };
        real_type const SR{ ( Y[n] - Y[n - 1] ) / hR };
        real_type const dp0{ Yp[n - 1] };
        real_type const dpL{ Yp[n - 2] };
        real_type const dpR{ Yp[n] };
        Z[n] = second_deriv3p_R( SL, hL, SR, hR, dpL, dp0, dpR );
      }
    }

    integer i = 0;
    do
    {
      Z[i] /= D[i];
      U[i] /= D[i];
      D[i + 1] -= L[i + 1] * U[i];
      Z[i + 1] -= L[i + 1] * Z[i];
    } while ( ++i < n );

    Z[i] /= D[i];

    do
    {
      --i;
      Z[i] -= U[i] * Z[i + 1];
    } while ( i > 0 );
  }

  static void QuinticSpline_Ypp_build(
    real_type const X[],
    real_type const Y[],
    real_type const Yp[],
    real_type       Ypp[],
    integer const   npts )
  {
    UTILS_ASSERT( npts >= 2, "QuinticSpline_Ypp_build, npts={} must be >= 2\n", npts );

    integer const n{ npts - 1 };

    if ( n == 1 )
    {
      Ypp[0] = Ypp[1] = 0;
      return;
    }

    {
      real_type const hL{ X[1] - X[0] };
      real_type const hR{ X[2] - X[1] };
      real_type const SL{ ( Y[1] - Y[0] ) / hL };
      real_type const SR{ ( Y[2] - Y[1] ) / hR };
      // Ypp[0] = (2*SL-SR)*al+SL*be;
      // Ypp[0] = second_deriv3p_L( SL, hL, SR, hR, Yp[0] );
      Ypp[0] = second_deriv3p_L( SL, hL, SR, hR, Yp[0], Yp[1], Yp[2] );
    }
    {
      real_type const hL{ X[n - 1] - X[n - 2] };
      real_type const hR{ X[n] - X[n - 1] };
      real_type const SL{ ( Y[n - 1] - Y[n - 2] ) / hL };
      real_type const SR{ ( Y[n] - Y[n - 1] ) / hR };
      // Ypp[n] = (2*SR-SL)*be+SR*al;
      // Ypp[n] = second_deriv3p_R( SL, hL, SR, hR, Yp[n] );
      Ypp[n] = second_deriv3p_R( SL, hL, SR, hR, Yp[n - 2], Yp[n - 1], Yp[n] );
    }

    for ( integer i = 1; i < n; ++i )
    {
      real_type const hL{ X[i] - X[i - 1] };
      real_type const hR{ X[i + 1] - X[i] };
      real_type const SL{ ( Y[i] - Y[i - 1] ) / hL };
      real_type const SR{ ( Y[i + 1] - Y[i] ) / hR };
      // Ypp[i] = second_deriv3p_C( SL, hL, SR, hR, Yp[i] );
      real_type ddC = second_deriv3p_C( SL, hL, SR, hR, Yp[i - 1], Yp[i], Yp[i + 1] );
      if ( i > 1 )
      {
        real_type const hLL{ X[i - 1] - X[i - 2] };
        real_type const SLL{ ( Y[i - 1] - Y[i - 2] ) / hLL };
        // real_type dd  = second_deriv3p_R( SLL, hLL, SL, hR, Yp[i] );
        real_type ddL{ second_deriv3p_R( SLL, hLL, SL, hR, Yp[i - 2], Yp[i - 1], Yp[i] ) };
        if ( ddL * ddC < 0 )
          ddC = 0;
        else if ( std::abs( ddL ) < std::abs( ddC ) )
          ddC = ddL;
      }
      if ( i < n - 1 )
      {
        real_type const hRR{ X[i + 2] - X[i + 1] };
        real_type const SRR{ ( Y[i + 2] - Y[i + 1] ) / hRR };
        // real_type dd = second_deriv3p_L( SR, hR, SRR, hRR, Yp[i] );
        real_type ddR{ second_deriv3p_L( SR, hR, SRR, hRR, Yp[i], Yp[i + 1], Yp[i + 2] ) };
        if ( ddR * ddC < 0 )
          ddC = 0;
        else if ( std::abs( ddR ) < std::abs( ddC ) )
          ddC = ddR;
      }
      Ypp[i] = ddC;
    }
  }

  /*\
   |    ___        _       _   _      ____        _ _
   |   / _ \ _   _(_)_ __ | |_(_) ___/ ___| _ __ | (_)_ __   ___
   |  | | | | | | | | '_ \| __| |/ __\___ \| '_ \| | | '_ \ / _ \
   |  | |_| | |_| | | | | | |_| | (__ ___) | |_) | | | | | |  __/
   |   \__\_\\__,_|_|_| |_|\__|_|\___|____/| .__/|_|_|_| |_|\___|
   |                                       |_|
   |
  \*/

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void Quintic_build(
    QuinticSpline_sub_type const q_sub_type,
    real_type const              X[],
    real_type const              Y[],
    real_type                    Yp[],
    real_type                    Ypp[],
    integer const                npts )
  {
    UTILS_ASSERT( npts >= 2, "Quintic_build, npts={} must be >= 2\n", npts );

    switch ( q_sub_type )
    {
      case QuinticSpline_sub_type::CUBIC:
        CubicSpline_build( X, Y, Yp, npts, CubicSpline_BC::EXTRAPOLATE, CubicSpline_BC::EXTRAPOLATE );
        break;
      case QuinticSpline_sub_type::PCHIP: Pchip_build( X, Y, Yp, npts ); break;
      case QuinticSpline_sub_type::AKIMA:
      {
        Malloc_real mem( "Quintic_build::work memory" );
        Akima_build( X, Y, Yp, mem.malloc( npts ), npts );
        mem.free();
      }
      break;
      case QuinticSpline_sub_type::BESSEL: BesselSpline::build( X, Y, Yp, npts ); break;
    }
    QuinticSpline_Ypp_build( X, Y, Yp, Ypp, npts );
  }


}  // namespace Splines
