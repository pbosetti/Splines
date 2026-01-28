/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2026                                                      |
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

#include "Splines.hh"
#include "Utils_fmt.hh"
#include "Utils_eigen.hh"

#include <random>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

namespace SplinesTest
{

  using namespace Splines;
  using namespace std;
  using Eigen::Map;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  // ===========================================================================
  // TEST FUNCTIONS - Well-behaved smooth functions for spline testing
  // ===========================================================================

  // Function 1: Simple 2D polynomial (for basic testing)
  real_type poly_function( real_type x, real_type y )
  {
    return 1.0 + 0.5 * x - 0.3 * y + 0.2 * x * x - 0.1 * x * y + 0.05 * y * y;
  }

  real_type poly_dx( real_type x, real_type y )
  {
    return 0.5 + 0.4 * x - 0.1 * y;
  }

  real_type poly_dy( real_type x, real_type y )
  {
    return -0.3 - 0.1 * x + 0.1 * y;
  }

  real_type poly_dxx( real_type x, real_type y )
  {
    return 0.4;
  }

  real_type poly_dyy( real_type x, real_type y )
  {
    return 0.1;
  }

  real_type poly_dxy( real_type x, real_type y )
  {
    return -0.1;
  }

  // Function 2: Smooth exponential function (for more challenging tests)
  real_type exp_function( real_type x, real_type y )
  {
    return exp( -0.1 * ( x * x + y * y ) ) * cos( x ) * sin( y );
  }

  real_type exp_dx( real_type x, real_type y )
  {
    real_type g  = exp( -0.1 * ( x * x + y * y ) );
    real_type dg = -0.2 * x * g;
    return dg * cos( x ) * sin( y ) - g * sin( x ) * sin( y );
  }

  real_type exp_dy( real_type x, real_type y )
  {
    real_type g  = exp( -0.1 * ( x * x + y * y ) );
    real_type dg = -0.2 * y * g;
    return dg * cos( x ) * sin( y ) + g * cos( x ) * cos( y );
  }

  real_type exp_dxx( real_type x, real_type y )
  {
    real_type g   = exp( -0.1 * ( x * x + y * y ) );
    real_type dg  = -0.2 * x * g;
    real_type d2g = ( 0.04 * x * x - 0.2 ) * g;
    return d2g * cos( x ) * sin( y ) - 2.0 * dg * sin( x ) * sin( y ) - g * cos( x ) * sin( y );
  }

  real_type exp_dyy( real_type x, real_type y )
  {
    real_type g   = exp( -0.1 * ( x * x + y * y ) );
    real_type dg  = -0.2 * y * g;
    real_type d2g = ( 0.04 * y * y - 0.2 ) * g;
    return d2g * cos( x ) * sin( y ) + 2.0 * dg * cos( x ) * cos( y ) - g * cos( x ) * sin( y );
  }

  real_type exp_dxy( real_type x, real_type y )
  {
    real_type g    = exp( -0.1 * ( x * x + y * y ) );
    real_type dgx  = -0.2 * x * g;
    real_type dgy  = -0.2 * y * g;
    real_type dgxy = 0.04 * x * y * g;

    return dgxy * cos( x ) * sin( y ) + dgx * cos( x ) * cos( y ) - dgy * sin( x ) * sin( y ) - g * sin( x ) * cos( y );
  }

  // Function 3: Gentle sine-cosine product (for accuracy tests)
  real_type trig_function( real_type x, real_type y )
  {
    return sin( 0.5 * x ) * cos( 0.5 * y );
  }

  real_type trig_dx( real_type x, real_type y )
  {
    return 0.5 * cos( 0.5 * x ) * cos( 0.5 * y );
  }

  real_type trig_dy( real_type x, real_type y )
  {
    return -0.5 * sin( 0.5 * x ) * sin( 0.5 * y );
  }

  real_type trig_dxx( real_type x, real_type y )
  {
    return -0.25 * sin( 0.5 * x ) * cos( 0.5 * y );
  }

  real_type trig_dyy( real_type x, real_type y )
  {
    return -0.25 * sin( 0.5 * x ) * cos( 0.5 * y );
  }

  real_type trig_dxy( real_type x, real_type y )
  {
    return -0.25 * cos( 0.5 * x ) * sin( 0.5 * y );
  }

  // ===========================================================================
  // UTILITY FUNCTIONS
  // ===========================================================================

  // Finite differences using Eigen-friendly interface
  template <typename Func> real_type finite_diff_dx( Func f, real_type x, real_type y, real_type h = 1e-6 )
  {
    return ( f( x + h, y ) - f( x - h, y ) ) / ( 2 * h );
  }

  template <typename Func> real_type finite_diff_dy( Func f, real_type x, real_type y, real_type h = 1e-6 )
  {
    return ( f( x, y + h ) - f( x, y - h ) ) / ( 2 * h );
  }

  template <typename Func> real_type finite_diff_dxx( Func f, real_type x, real_type y, real_type h = 1e-6 )
  {
    return ( f( x + h, y ) - 2 * f( x, y ) + f( x - h, y ) ) / ( h * h );
  }

  template <typename Func> real_type finite_diff_dyy( Func f, real_type x, real_type y, real_type h = 1e-6 )
  {
    return ( f( x, y + h ) - 2 * f( x, y ) + f( x, y - h ) ) / ( h * h );
  }

  template <typename Func> real_type finite_diff_dxy( Func f, real_type x, real_type y, real_type h = 1e-6 )
  {
    real_type h2 = 0.5 * h;
    return ( f( x + h2, y + h2 ) - f( x + h2, y - h2 ) - f( x - h2, y + h2 ) + f( x - h2, y - h2 ) ) / ( h * h );
  }

  // Statistics structure con punto di errore massimo
  struct DerivativeStats
  {
    string    name;
    real_type max_abs_error = 0;
    real_type rms_error     = 0;
    real_type max_rel_error = 0;
    integer   count         = 0;
    bool      is_derivative = false;

    // Punti dove si verificano gli errori massimi
    real_type max_error_x   = 0;
    real_type max_error_y   = 0;
    real_type exact_at_max  = 0;
    real_type approx_at_max = 0;

    void update( real_type x, real_type y, real_type exact, real_type approx )
    {
      real_type abs_err = abs( exact - approx );

      // Aggiorna punto di errore massimo
      if ( abs_err > max_abs_error )
      {
        max_abs_error = abs_err;
        max_error_x   = x;
        max_error_y   = y;
        exact_at_max  = exact;
        approx_at_max = approx;
      }

      real_type rel_err = 0.0;
      if ( abs( exact ) > 1e-12 ) { rel_err = abs_err / abs( exact ); }
      max_rel_error = max( max_rel_error, rel_err );

      rms_error += abs_err * abs_err;
      ++count;
    }

    void finalize()
    {
      if ( count > 0 ) { rms_error = sqrt( rms_error / static_cast<real_type>( count ) ); }
    }

    void print_row() const
    {
      fmt::color color;
      if ( max_abs_error < 1e-6 ) { color = fmt::color::green; }
      else if ( max_abs_error < 1e-3 ) { color = fmt::color::yellow; }
      else
      {
        color = fmt::color::red;
      }

      if ( is_derivative )
      {
        fmt::print(
          fg( color ),
          "│ {:<14} │ {:^12.2e} │ {:^12.2e} │ {:^12.2e} │ {:18} │\n",
          name,
          max_abs_error,
          rms_error,
          max_rel_error,
          fmt::format( "({:.3f},{:.3f})", max_error_x, max_error_y ) );
      }
      else
      {
        fmt::print(
          fg( color ),
          "│ {:<14} │ {:^12.2e} │ {:^12.2e} │ {:^12.2e} │ {:18} │\n",
          name,
          max_abs_error,
          rms_error,
          max_rel_error,
          fmt::format( "({:.3f},{:.3f})", max_error_x, max_error_y ) );
      }
    }

    void print_max_error_details() const
    {
      if ( count > 0 )
      {
        fmt::print(
          fg( fmt::color::light_gray ),
          "    Max error at ({:.4f}, {:.4f}): exact={:.6e}, approx={:.6e}, diff={:.2e}\n",
          max_error_x,
          max_error_y,
          exact_at_max,
          approx_at_max,
          max_abs_error );
      }
    }
  };

  struct SplineTestResults
  {
    string          spline_name;
    DerivativeStats value{ "f(x,y)", 0, 0, 0, 0, false };
    DerivativeStats dx{ "∂f/∂x", 0, 0, 0, 0, true };
    DerivativeStats dy{ "∂f/∂y", 0, 0, 0, 0, true };
    DerivativeStats dxx{ "∂²f/∂x²", 0, 0, 0, 0, true };
    DerivativeStats dyy{ "∂²f/∂y²", 0, 0, 0, 0, true };
    DerivativeStats dxy{ "∂²f/∂x∂y", 0, 0, 0, 0, true };
    DerivativeStats dx_fd{ "∂f/∂x (FD)", 0, 0, 0, 0, true };  // Consistency with FD
    DerivativeStats dy_fd{ "∂f/∂y (FD)", 0, 0, 0, 0, true };
    DerivativeStats dxx_fd{ "∂²f/∂x² (FD)", 0, 0, 0, 0, true };
    DerivativeStats dyy_fd{ "∂²f/∂y² (FD)", 0, 0, 0, 0, true };
    DerivativeStats dxy_fd{ "∂²f/∂x∂y (FD)", 0, 0, 0, 0, true };

    // Statistiche per interpolazione sui punti di griglia
    DerivativeStats grid_interp{ "Grid Interp", 0, 0, 0, 0, false };
    real_type       grid_max_interp_error{ 0 };
    real_type       grid_max_interp_x{ 0 };
    real_type       grid_max_interp_y{ 0 };

    SplineTestResults( const string & name ) : spline_name( name ) {}

    void finalize()
    {
      value.finalize();
      dx.finalize();
      dy.finalize();
      dxx.finalize();
      dyy.finalize();
      dxy.finalize();
      dx_fd.finalize();
      dy_fd.finalize();
      dxx_fd.finalize();
      dyy_fd.finalize();
      dxy_fd.finalize();
      grid_interp.finalize();
    }

    void print_table() const
    {
      fmt::print(
        fg( fmt::color::cyan ),
        "\n┌{0:─^{2}}┐\n"
        "│{1: ^{2}}│\n"
        "└{0:─^{2}}┘\n",
        "",
        spline_name,
        78 );

      fmt::print(
        "┌────────────────┬──────────────┬──────────────┬──────────────┬────────────────────┐\n"
        "│ Derivative     │ Max Abs Err  │ RMS Error    │ Max Rel Err  │    Max Err At      │\n"
        "├────────────────┼──────────────┼──────────────┼──────────────┼────────────────────┤\n" );

      value.print_row();
      dx.print_row();
      dy.print_row();

      if ( dxx.count > 0 )
      {
        dxx.print_row();
        dyy.print_row();
        dxy.print_row();
      }

      fmt::print(
        "├────────────────┼──────────────┼──────────────┼──────────────┼────────────────────┤\n"
        "│ FD Consistency │              │              │              │                    │\n" );

      dx_fd.print_row();
      dy_fd.print_row();

      if ( dxx_fd.count > 0 )
      {
        dxx_fd.print_row();
        dyy_fd.print_row();
        dxy_fd.print_row();
      }

      // Aggiungi riga per interpolazione sui punti di griglia
      fmt::print( "├────────────────┼──────────────┼──────────────┼──────────────┼────────────────────┤\n" );
      grid_interp.print_row();

      fmt::print( "└────────────────┴──────────────┴──────────────┴──────────────┴────────────────────┘\n" );

      // Print summary con dettagli errori massimi
      fmt::print( "\n📋 Summary for {}:\n", spline_name );

      // Dettagli errori massimi per ogni tipo
      fmt::print( "  Max error details:\n" );
      value.print_max_error_details();
      dx.print_max_error_details();
      dy.print_max_error_details();

      if ( dxx.count > 0 )
      {
        dxx.print_max_error_details();
        dyy.print_max_error_details();
        dxy.print_max_error_details();
      }

      // Dettagli per interpolazione griglia
      if ( grid_interp.count > 0 )
      {
        fmt::print( "  Grid interpolation:\n" );
        grid_interp.print_max_error_details();
      }

      bool all_good = true;
      if ( value.max_abs_error > 1e-4 )
      {
        fmt::print( fg( fmt::color::yellow ), "  ⚠️  Function error: {:.2e}\n", value.max_abs_error );
        all_good = false;
      }

      if ( grid_interp.max_abs_error > 1e-10 )
      {
        fmt::print( fg( fmt::color::yellow ), "  ⚠️  Grid interpolation error: {:.2e}\n", grid_interp.max_abs_error );
        all_good = false;
      }

      if ( dx_fd.max_abs_error > 1e-6 )
      {
        fmt::print( fg( fmt::color::red ), "  ❌ ∂f/∂x FD consistency: {:.2e}\n", dx_fd.max_abs_error );
        all_good = false;
      }

      if ( dy_fd.max_abs_error > 1e-6 )
      {
        fmt::print( fg( fmt::color::red ), "  ❌ ∂f/∂y FD consistency: {:.2e}\n", dy_fd.max_abs_error );
        all_good = false;
      }

      if ( dxx.count > 0 )
      {
        if ( dxx_fd.max_abs_error > 1e-4 )
        {
          fmt::print( fg( fmt::color::yellow ), "  ⚠️  ∂²f/∂x² FD consistency: {:.2e}\n", dxx_fd.max_abs_error );
          all_good = false;
        }

        if ( dyy_fd.max_abs_error > 1e-4 )
        {
          fmt::print( fg( fmt::color::yellow ), "  ⚠️  ∂²f/∂y² FD consistency: {:.2e}\n", dyy_fd.max_abs_error );
          all_good = false;
        }

        if ( dxy_fd.max_abs_error > 1e-4 )
        {
          fmt::print( fg( fmt::color::yellow ), "  ⚠️  ∂²f/∂x∂y FD consistency: {:.2e}\n", dxy_fd.max_abs_error );
          all_good = false;
        }
      }

      if ( all_good ) { fmt::print( fg( fmt::color::green ), "  ✅ All tests passed!\n" ); }
    }
  };

  // ===========================================================================
  // TEST FUNCTION
  // ===========================================================================

  void test_spline_type(
    SplineType2D     spline_type,
    const VectorXd & x_grid,
    const VectorXd & y_grid,
    const MatrixXd & z_data,
    integer          n_test_points = 500 )
  {
    string type_name;
    switch ( spline_type )
    {
      case SplineType2D::BILINEAR: type_name = "Bilinear"; break;
      case SplineType2D::BICUBIC_CUBIC: type_name = "Bicubic"; break;
      case SplineType2D::BICUBIC_AKIMA: type_name = "Bicubic[Akima]"; break;
      case SplineType2D::BICUBIC_BESSEL: type_name = "Bicubic[Bessel]"; break;
      case SplineType2D::BICUBIC_PCHIP: type_name = "Bicubic[Pchip]"; break;
      case SplineType2D::BIQUINTIC_CUBIC: type_name = "Biquintic"; break;
      case SplineType2D::BIQUINTIC_AKIMA: type_name = "Biquintic[Akima]"; break;
      case SplineType2D::BIQUINTIC_BESSEL: type_name = "Biquintic[Bessel]"; break;
      case SplineType2D::BIQUINTIC_PCHIP: type_name = "Biquintic[Pchip]"; break;
      default: type_name = "Unknown"; break;
    }

    // Create and build spline
    Spline2D spline( type_name );
    integer  nx = static_cast<integer>( x_grid.size() );
    integer  ny = static_cast<integer>( y_grid.size() );
    spline.build( spline_type, x_grid.data(), 1, y_grid.data(), 1, z_data.data(), nx, nx, ny, true, false );

    // Get bounds
    real_type x_min = spline.x_min();
    real_type x_max = spline.x_max();
    real_type y_min = spline.y_min();
    real_type y_max = spline.y_max();

    // Test interpolation at grid points

    SplineTestResults results( type_name );

    // Prima testiamo l'interpolazione sui punti di griglia
    fmt::print( fg( fmt::color::light_blue ), "\n  Testing interpolation at grid points for {}...\n", type_name );

    for ( integer i = 0; i < nx; ++i )
    {
      for ( integer j = 0; j < ny; ++j )
      {
        real_type x          = x_grid[i];
        real_type y          = y_grid[j];
        real_type exact_val  = z_data( i, j );
        real_type spline_val = spline.eval( x, y );
        results.grid_interp.update( x, y, exact_val, spline_val );
      }
    }

    fmt::print( fg( fmt::color::light_green ), "    Grid points tested: {} × {} = {}\n", nx, ny, nx * ny );

    // Generate test points using Eigen (avoid boundaries)
    std::mt19937                              rng( 1000 + static_cast<unsigned>( spline_type ) );
    std::uniform_real_distribution<real_type> dist_x( x_min * 0.9, x_max * 0.9 );
    std::uniform_real_distribution<real_type> dist_y( y_min * 0.9, y_max * 0.9 );

    // Use trigonometric function for testing (smooth and well-behaved)
    auto & func     = trig_function;
    auto & func_dx  = trig_dx;
    auto & func_dy  = trig_dy;
    auto & func_dxx = trig_dxx;
    auto & func_dyy = trig_dyy;
    auto & func_dxy = trig_dxy;

    // Test points - use Eigen for generating test points
    VectorXd x_test = VectorXd::NullaryExpr( n_test_points, [&]() { return dist_x( rng ); } );
    VectorXd y_test = VectorXd::NullaryExpr( n_test_points, [&]() { return dist_y( rng ); } );

    for ( integer i_test = 0; i_test < n_test_points; ++i_test )
    {
      real_type x = x_test( i_test );
      real_type y = y_test( i_test );

      // Avoid points too close to boundaries
      real_type margin = 0.05 * ( x_max - x_min );
      if ( x < x_min + margin || x > x_max - margin || y < y_min + margin || y > y_max - margin ) continue;

      // Exact values
      real_type exact_val = func( x, y );
      real_type exact_dx  = func_dx( x, y );
      real_type exact_dy  = func_dy( x, y );
      real_type exact_dxx = func_dxx( x, y );
      real_type exact_dyy = func_dyy( x, y );
      real_type exact_dxy = func_dxy( x, y );

      // Spline values
      real_type spline_val = spline.eval( x, y );
      real_type spline_dx  = spline.Dx( x, y );
      real_type spline_dy  = spline.Dy( x, y );
      real_type spline_dxx = spline.Dxx( x, y );
      real_type spline_dyy = spline.Dyy( x, y );
      real_type spline_dxy = spline.Dxy( x, y );

      // Update accuracy stats
      results.value.update( x, y, exact_val, spline_val );
      results.dx.update( x, y, exact_dx, spline_dx );
      results.dy.update( x, y, exact_dy, spline_dy );

      if ( spline_type != SplineType2D::BILINEAR )
      {
        results.dxx.update( x, y, exact_dxx, spline_dxx );
        results.dyy.update( x, y, exact_dyy, spline_dyy );
        results.dxy.update( x, y, exact_dxy, spline_dxy );
      }

      // Finite difference consistency check
      auto spline_eval = [&]( real_type xx, real_type yy ) { return spline.eval( xx, yy ); };

      real_type fd_dx  = finite_diff_dx( spline_eval, x, y );
      real_type fd_dy  = finite_diff_dy( spline_eval, x, y );
      real_type fd_dxx = finite_diff_dxx( spline_eval, x, y );
      real_type fd_dyy = finite_diff_dyy( spline_eval, x, y );
      real_type fd_dxy = finite_diff_dxy( spline_eval, x, y );

      results.dx_fd.update( x, y, spline_dx, fd_dx );
      results.dy_fd.update( x, y, spline_dy, fd_dy );

      if ( spline_type != SplineType2D::BILINEAR )
      {
        results.dxx_fd.update( x, y, spline_dxx, fd_dxx );
        results.dyy_fd.update( x, y, spline_dyy, fd_dyy );
        results.dxy_fd.update( x, y, spline_dxy, fd_dxy );
      }
    }

    results.finalize();
    results.print_table();
  }

  // Test interpolation at grid points
  void test_interpolation_accuracy( const VectorXd & x_grid, const VectorXd & y_grid, const MatrixXd & z_data )
  {
    fmt::print(
      fg( fmt::color::light_blue ) | fmt::emphasis::bold,
      "\n"
      "╔══════════════════════════════════════════════════════════╗\n"
      "║          Interpolation Accuracy at Grid Points           ║\n"
      "╚══════════════════════════════════════════════════════════╝\n" );

    // Test with Bicubic spline (should interpolate exactly for all types)
    Spline2D spline( "TestSpline" );
    integer  nx = static_cast<integer>( x_grid.size() );
    integer  ny = static_cast<integer>( y_grid.size() );
    spline.build(
      SplineType2D::BIQUINTIC_CUBIC,
      x_grid.data(),
      1,
      y_grid.data(),
      1,
      z_data.data(),
      nx,
      nx,
      ny,
      true,
      false );

    real_type max_error   = 0.0;
    real_type max_error_x = 0.0;
    real_type max_error_y = 0.0;
    real_type avg_error   = 0.0;
    integer   count       = 0;

    for ( integer i = 0; i < nx; ++i )
    {
      for ( integer j = 0; j < ny; ++j )
      {
        real_type exact  = z_data( i, j );
        real_type interp = spline.eval( x_grid[i], y_grid[j] );
        real_type error  = abs( interp - exact );

        if ( error > max_error )
        {
          max_error   = error;
          max_error_x = x_grid[i];
          max_error_y = y_grid[j];
        }

        avg_error += error;
        ++count;
      }
    }

    avg_error /= static_cast<real_type>( count );

    fmt::print( "\nGrid interpolation test ({}×{} points):\n", nx, ny );

    fmt::print( "  Maximum error: {:12.2e} at ({:.4f}, {:.4f})\n", max_error, max_error_x, max_error_y );
    fmt::print( "  Average error: {:12.2e}\n", avg_error );

    if ( max_error < 1e-10 ) { fmt::print( fg( fmt::color::green ), "  ✅ Interpolation is exact at grid points\n" ); }
    else if ( max_error < 1e-6 )
    {
      fmt::print( fg( fmt::color::yellow ), "  ⚠️  Small interpolation errors (numerical precision)\n" );
    }
    else
    {
      fmt::print( fg( fmt::color::red ), "  ❌ Large interpolation errors - check implementation!\n" );
    }
  }

  // Test derivative method consistency
  void test_derivative_consistency( Spline2D & spline, const VectorXd & x_grid, const VectorXd & y_grid )
  {
    fmt::print(
      fg( fmt::color::light_blue ) | fmt::emphasis::bold,
      "\n"
      "╔══════════════════════════════════════════════════════════╗\n"
      "║          Derivative Method Consistency Test              ║\n"
      "╚══════════════════════════════════════════════════════════╝\n" );

    std::mt19937                              rng( 1234 );
    std::uniform_real_distribution<real_type> dist_x( x_grid( 0 ) * 0.8, x_grid.tail( 1 )( 0 ) * 0.8 );
    std::uniform_real_distribution<real_type> dist_y( y_grid( 0 ) * 0.8, y_grid.tail( 1 )( 0 ) * 0.8 );

    real_type max_D_error   = 0.0;
    real_type max_DD_error  = 0.0;
    real_type max_D_error_x = 0.0, max_D_error_y = 0.0;
    real_type max_DD_error_x = 0.0, max_DD_error_y = 0.0;
    integer   test_count = 0;

    // Generate test points using Eigen
    VectorXd x_test = VectorXd::NullaryExpr( 100, [&]() { return dist_x( rng ); } );
    VectorXd y_test = VectorXd::NullaryExpr( 100, [&]() { return dist_y( rng ); } );

    for ( integer i = 0; i < 100; ++i )
    {
      real_type x = x_test( i );
      real_type y = y_test( i );

      // Skip boundary regions
      real_type margin = 0.1 * ( x_grid.tail( 1 )( 0 ) - x_grid( 0 ) );
      if (
        x < x_grid( 0 ) + margin || x > x_grid.tail( 1 )( 0 ) - margin || y < y_grid( 0 ) + margin ||
        y > y_grid.tail( 1 )( 0 ) - margin )
        continue;

      // Test D() method
      real_type d[3];
      spline.D( x, y, d );

      real_type dx_indiv = spline.Dx( x, y );
      real_type dy_indiv = spline.Dy( x, y );

      real_type d_error = max( abs( d[1] - dx_indiv ), abs( d[2] - dy_indiv ) );
      if ( d_error > max_D_error )
      {
        max_D_error   = d_error;
        max_D_error_x = x;
        max_D_error_y = y;
      }

      // Test DD() method
      real_type dd[6];
      spline.DD( x, y, dd );

      real_type dxx_indiv = spline.Dxx( x, y );
      real_type dyy_indiv = spline.Dyy( x, y );
      real_type dxy_indiv = spline.Dxy( x, y );

      real_type dd_error = max( max( abs( dd[3] - dxx_indiv ), abs( dd[5] - dyy_indiv ) ), abs( dd[4] - dxy_indiv ) );
      if ( dd_error > max_DD_error )
      {
        max_DD_error   = dd_error;
        max_DD_error_x = x;
        max_DD_error_y = y;
      }

      ++test_count;
    }

    fmt::print( "\nMethod consistency test ({} points):\n", test_count );
    fmt::print(
      "  D() method max inconsistency:  {:12.2e} at ({:.4f}, {:.4f})\n",
      max_D_error,
      max_D_error_x,
      max_D_error_y );
    fmt::print(
      "  DD() method max inconsistency: {:12.2e} at ({:.4f}, {:.4f})\n",
      max_DD_error,
      max_DD_error_x,
      max_DD_error_y );

    if ( max_D_error < 1e-12 && max_DD_error < 1e-12 )
    {
      fmt::print( fg( fmt::color::green ), "  ✅ All derivative methods are perfectly consistent\n" );
    }
    else if ( max_D_error < 1e-8 && max_DD_error < 1e-8 )
    {
      fmt::print( fg( fmt::color::yellow ), "  ⚠️  Minor inconsistencies (numerical precision)\n" );
    }
    else
    {
      fmt::print( fg( fmt::color::red ), "  ❌ Significant inconsistencies found!\n" );
    }
  }

}  // namespace SplinesTest

int main()
{
  using namespace SplinesTest;

  fmt::print(
    fg( fmt::color::light_blue ) | fmt::emphasis::bold,
    "\n"
    "╔══════════════════════════════════════════════════════════╗\n"
    "║          2D Spline Library - Validation Test             ║\n"
    "║                  (Well-behaved functions)                ║\n"
    "╚══════════════════════════════════════════════════════════╝\n\n" );

  // Create a moderately sized uniform grid using Eigen
  integer nx = 51;  // Reduced for faster testing
  integer ny = 63;

  // Domain for trigonometric function (gentle oscillations)
  real_type XMIN = -M_PI;
  real_type XMAX = M_PI;
  real_type YMIN = -M_PI;
  real_type YMAX = M_PI;

  // Create uniform grid using Eigen
  VectorXd x_grid = VectorXd::LinSpaced( nx, XMIN, XMAX );
  VectorXd y_grid = VectorXd::LinSpaced( ny, YMIN, YMAX );

  // Generate data using trigonometric function (well-behaved)
  MatrixXd z_data = MatrixXd::NullaryExpr(
    nx,
    ny,
    [&]( Eigen::Index i, Eigen::Index j ) { return trig_function( x_grid( i ), y_grid( j ) ); } );

  fmt::print(
    "Test configuration:\n"
    "  Grid: {} × {} points\n"
    "  Domain: x ∈ [{:.2f}, {:.2f}], y ∈ [{:.2f}, {:.2f}]\n"
    "  Test function: sin(0.5x) * cos(0.5y) (smooth, well-behaved)\n"
    "  Test points: 500 per spline type\n\n",
    nx,
    ny,
    x_grid.minCoeff(),
    x_grid.maxCoeff(),
    y_grid.minCoeff(),
    y_grid.maxCoeff() );

  // Test all spline types
  vector<SplineType2D> spline_types = { SplineType2D::BILINEAR,        SplineType2D::BICUBIC_CUBIC,
                                        SplineType2D::BICUBIC_AKIMA,   SplineType2D::BICUBIC_BESSEL,
                                        SplineType2D::BICUBIC_PCHIP,   SplineType2D::BIQUINTIC_CUBIC,
                                        SplineType2D::BIQUINTIC_AKIMA, SplineType2D::BIQUINTIC_BESSEL,
                                        SplineType2D::BIQUINTIC_PCHIP };

  for ( auto type : spline_types ) { test_spline_type( type, x_grid, y_grid, z_data, 500 ); }

  // Test interpolation accuracy
  test_interpolation_accuracy( x_grid, y_grid, z_data );

  // Test derivative consistency
  Spline2D test_spline( "ConsistencyTest" );

  test_spline
    .build( SplineType2D::BIQUINTIC_CUBIC, x_grid.data(), 1, y_grid.data(), 1, z_data.data(), nx, nx, ny, true, false );
  test_derivative_consistency( test_spline, x_grid, y_grid );

  // Final summary
  fmt::print(
    fg( fmt::color::light_green ) | fmt::emphasis::bold,
    "\n\n"
    "════════════════════════════════════════════════════════════\n"
    "                    TEST COMPLETED                          \n"
    "════════════════════════════════════════════════════════════\n" );

  fmt::print(
    "\nInterpretation guide:\n"
    "  ✅ Excellent (error < 1e-6)\n"
    "  ⚠️  Acceptable (error < 1e-3)\n"
    "  ❌ Problematic (error ≥ 1e-3)\n"
    "\n"
    "Note: Finite Difference (FD) consistency errors measure how well\n"
    "the spline's analytical derivatives match numerical derivatives\n"
    "of the spline itself. This is a key consistency check.\n"
    "\n"
    "Grid interpolation test: Checks if the spline exactly matches\n"
    "the input data at the grid points (as expected for interpolating splines).\n" );

  return 0;
}
