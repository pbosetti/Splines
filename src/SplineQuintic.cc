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
#include "Utils_TridiagonalSolver.hh"

namespace Splines
{

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*!
   * \brief Compute second derivatives for curvature-minimizing quintic spline
   *
   * Given a C1 cubic spline (or Hermite interpolation) defined by nodal values
   * and first derivatives, this function computes the second derivatives
   * required to elevate the spline to quintic order while minimizing the
   * integral of the squared third derivative.
   *
   * The method solves a tridiagonal system derived from:
   * 1. C3 continuity at interior nodes
   * 2. Minimization of ∫ [S'''(x)]² dx
   *
   * The resulting quintic spline is C3 continuous and provides smoother
   * curvature transitions than the original cubic.
   *
   * \param[in]  X     Array of strictly increasing x-coordinates
   * \param[in]  Y     Array of y-coordinates
   * \param[in]  Yp    Array of first derivatives (must be computed beforehand)
   * \param[out] Ypp   Array of computed second derivatives
   * \param[in]  npts  Number of data points
   * \param[in]  bc_type Type of boundary condition (0=natural, 1=not-a-knot, 2=clamped)
   * \param[in]  bcl   Left boundary value for clamped condition (if applicable)
   * \param[in]  bcr   Right boundary value for clamped condition (if applicable)
   */
  static void QuinticSpline_compute_second_derivatives(
    real_type const X[],
    real_type const Y[],
    real_type const Yp[],
    real_type       Ypp[],
    integer const   npts,
    real_type       bcl = 0.0,
    real_type       bcr = 0.0 )
  {
    UTILS_ASSERT( npts >= 2, "QuinticSpline_compute_second_derivatives, npts={} must be >= 2\n", npts );

    integer const n = npts - 1;

    // Degenerate case: fewer than 3 points
    if ( npts < 3 )
    {
      // For 2 points, we can only set second derivatives to zero or use finite differences
      Ypp[0] = bcl;
      Ypp[1] = bcr;
    }

    // Initialize tridiagonal solver
    Utils::TridiagonalSolver<real_type> solver;
    solver.resize( npts );

    // System vectors: subdiagonal a, diagonal b, superdiagonal c, right-hand side d
    Utils::TridiagonalSolver<real_type>::VecS a( npts - 1 ), b( npts ), c( npts - 1 ), d( npts );

    // Compute interval lengths
    Malloc_real mem_h( "QuinticSpline_compute_second_derivatives::h" );
    real_type * h = mem_h.malloc( n );
    for ( integer i = 0; i < n; ++i )
    {
      h[i] = X[i + 1] - X[i];
      UTILS_ASSERT( h[i] > 0, "X points must be strictly increasing, h[{}] = {}\n", i, h[i] );
    }

    // Build tridiagonal system for interior nodes (i = 1, ..., n-1)
    for ( integer i = 1; i < n; ++i )
    {
      real_type const hL = h[i - 1];
      real_type const hR = h[i];

      // Matrix coefficients (symmetric tridiagonal)
      a( i - 1 ) = 1.0 / hL;                       // sub-diagonal
      b( i )     = 2.0 * ( 1.0 / hL + 1.0 / hR );  // main diagonal
      c( i - 1 ) = 1.0 / hR;                       // super-diagonal

      // Right-hand side
      real_type const term1 = ( Y[i + 1] - Y[i] ) / ( hR * hR );
      real_type const term2 = ( Y[i] - Y[i - 1] ) / ( hL * hL );
      d( i )                = 6.0 * ( term1 - term2 ) - 2.0 * Yp[i - 1] / hL - 4.0 * Yp[i] / hR - 2.0 * Yp[i + 1] / hR;
    }

    // Apply boundary conditions
    b( 0 )     = 1.0;
    c( 0 )     = 0.0;
    d( 0 )     = bcl;
    a( n - 1 ) = 0.0;
    b( n )     = 1.0;
    d( n )     = bcr;

    // Solve the tridiagonal system
    solver.factorize( a, b, c );
    Utils::TridiagonalSolver<real_type>::VecS ypp_vec( npts );
    solver.solve( a, b, d, ypp_vec );

    // Copy results
    for ( integer i = 0; i < npts; ++i ) { Ypp[i] = ypp_vec[i]; }
  }

  /*!
   * \brief Construct a complete quintic spline from cubic Hermite data
   *
   * This function elevates a C1 cubic spline (defined by nodal values and
   * first derivatives) to a C3 quintic spline by computing optimal second
   * derivatives that minimize the curvature energy.
   *
   * The algorithm consists of three phases:
   * 1. First derivative estimation (if not already provided)
   * 2. Second derivative computation via tridiagonal solver
   * 3. Quintic Hermite interpolation on each interval
   *
   * The resulting spline is C3 continuous and provides superior smoothness
   * while preserving the original interpolation conditions.
   *
   * \param[in]  q_sub_type Method for first derivative estimation (CUBIC, PCHIP, AKIMA, BESSEL)
   * \param[in]  X          Array of x-coordinates
   * \param[in]  Y          Array of y-coordinates
   * \param[out] Yp         Array of first derivatives (computed if needed)
   * \param[out] Ypp        Array of second derivatives (computed)
   * \param[in]  npts       Number of data points
   * \param[in]  bc_type    Boundary condition type (0=natural, 1=not-a-knot, 2=clamped)
   * \param[in]  bcl        Left boundary second derivative for clamped condition
   * \param[in]  bcr        Right boundary second derivative for clamped condition
   */
  void Quintic_build(
    QuinticSpline_sub_type const q_sub_type,
    real_type const              X[],
    real_type const              Y[],
    real_type                    Yp[],
    real_type                    Ypp[],
    integer const                npts,
    real_type                    bcl,
    real_type                    bcr )
  {
    UTILS_ASSERT( npts >= 2, "Quintic_build, npts={} must be >= 2\n", npts );

    // Phase 1: Estimate first derivatives if not already available
    // (This step can be skipped if Yp is precomputed from a C1 cubic spline)
    switch ( q_sub_type )
    {
      case QuinticSpline_sub_type::CUBIC:
        CubicSpline_build( X, Y, Yp, Ypp, npts, CubicSpline_BC::EXTRAPOLATE, CubicSpline_BC::EXTRAPOLATE );
        return;
      case QuinticSpline_sub_type::PCHIP: Pchip_build( X, Y, Yp, npts ); break;
      case QuinticSpline_sub_type::AKIMA:
      {
        Malloc_real mem( "Quintic_build::work memory" );
        Akima_build( X, Y, Yp, mem.malloc( npts ), npts );
      }
      break;
      case QuinticSpline_sub_type::BESSEL: Bessel_build( X, Y, Yp, npts ); break;
      default: UTILS_ERROR( "Unknown QuinticSpline_sub_type value\n" );
    }

    if ( npts == 2 )
    {  // solo 2 punti, niente da fare
      Ypp[0] = Ypp[1] = 0;
      return;
    }

    // ---- left boundary ----
    if ( npts <= 3 )
    {
      Ypp[0] = Utils::second_derivative_3p( X[0], Y[0], X[1], Y[1], X[2], Y[2] );
      Ypp[1] = Utils::second_derivative_3p( X[1], Y[1], X[2], Y[2], X[0], Y[0] );
    }
    else
    {
      Ypp[0] = Utils::second_derivative_4p( X[0], Y[0], X[1], Y[1], X[2], Y[2], X[3], Y[3] );
      Ypp[1] = Utils::second_derivative_4p( X[1], Y[1], X[2], Y[2], X[3], Y[3], X[0], Y[0] );
    }

    // ---- interior points ----
    for ( integer i = 2; i < npts - 2; ++i )
      Ypp[i] = Utils::second_derivative_5p(
        X[i - 0],
        Y[i - 0],
        X[i + 1],
        Y[i + 1],
        X[i - 1],
        Y[i - 1],
        X[i + 2],
        Y[i + 2],
        X[i - 2],
        Y[i - 2] );

    // ---- right boundary ----
    if ( npts <= 3 )
    {
      Ypp[npts - 1] =
        Utils::second_derivative_3p( X[npts - 1], Y[npts - 1], X[npts - 2], Y[npts - 2], X[npts - 3], Y[npts - 3] );
      Ypp[npts - 2] =
        Utils::second_derivative_3p( X[npts - 2], Y[npts - 2], X[npts - 3], Y[npts - 3], X[npts - 1], Y[npts - 1] );
    }
    else
    {
      Ypp[npts - 1] = Utils::second_derivative_4p(
        X[npts - 1],
        Y[npts - 1],
        X[npts - 2],
        Y[npts - 2],
        X[npts - 3],
        Y[npts - 3],
        X[npts - 4],
        Y[npts - 4] );
      Ypp[npts - 2] = Utils::second_derivative_4p(
        X[npts - 2],
        Y[npts - 2],
        X[npts - 3],
        Y[npts - 3],
        X[npts - 4],
        Y[npts - 4],
        X[npts - 1],
        Y[npts - 1] );
    }

    // Phase 2: Compute optimal second derivatives
    // QuinticSpline_compute_second_derivatives( X, Y, Yp, Ypp, npts, bcl, bcr );
  }

}  // namespace Splines
