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
// FILE: SplinePatch.hxx
//

#pragma once

#ifndef SPLINE_PATCH_HH
#define SPLINE_PATCH_HH

namespace Splines
{

  /*
  //   ____  _ _ _
  //  | __ )(_) (_)_ __   ___  __ _ _ __
  //  |  _ \| | | | '_ \ / _ \/ _` | '__|
  //  | |_) | | | | | | |  __/ (_| | |
  //  |____/|_|_|_|_| |_|\___|\__,_|_|
  */

  real_type bilinear3( real_type const p[4], real_type const M[4][4], real_type const q[4] );

  real_type bilinear5( real_type const p[6], real_type const M[6][6], real_type const q[6] );

}  // namespace Splines

#endif

//
// EOF: SplinePatch.hxx
//
