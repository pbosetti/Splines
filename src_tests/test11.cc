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

#ifdef __clang__
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

namespace SplinesTest
{

  using namespace Splines;
  using namespace std;

  // Test function: f(x,y) = sin(x)*cos(y) + exp(-(x^2+y^2)/10)
  // This function is smooth but non-trivial
  real_type test_function( real_type x, real_type y )
  {
    return sin( x ) * cos( y ) + exp( -( x * x + y * y ) / 10.0 );
  }

  // Analytical derivatives of test function
  real_type test_function_dx( real_type x, real_type y )
  {
    return cos( x ) * cos( y ) - ( x / 5.0 ) * exp( -( x * x + y * y ) / 10.0 );
  }

  real_type test_function_dy( real_type x, real_type y )
  {
    return -sin( x ) * sin( y ) - ( y / 5.0 ) * exp( -( x * x + y * y ) / 10.0 );
  }

  real_type test_function_dxx( real_type x, real_type y )
  {
    real_type exp_term = exp( -( x * x + y * y ) / 10.0 );
    return -sin( x ) * cos( y ) + ( ( x * x / 25.0 ) - 1.0 / 5.0 ) * exp_term;
  }

  real_type test_function_dyy( real_type x, real_type y )
  {
    real_type exp_term = exp( -( x * x + y * y ) / 10.0 );
    return -sin( x ) * cos( y ) + ( ( y * y / 25.0 ) - 1.0 / 5.0 ) * exp_term;
  }

  real_type test_function_dxy( real_type x, real_type y )
  {
    real_type exp_term = exp( -( x * x + y * y ) / 10.0 );
    return -cos( x ) * sin( y ) + ( x * y / 25.0 ) * exp_term;
  }

  // Find interval index for a sorted array
  integer find_interval( const vector<real_type> & X, real_type x )
  {
    const integer n = static_cast<integer>( X.size() );
    if ( n <= 1 ) return 0;

    if ( x <= X.front() ) return 0;
    if ( x >= X.back() ) return n - 2;

    auto it = std::upper_bound( X.begin(), X.end(), x );
    return static_cast<integer>( std::distance( X.begin(), it ) ) - 1;
  }

  inline real_type safe_step( real_type x, const vector<real_type> & grid, integer i, real_type base_h )
  {
    const real_type left  = x - grid[i];
    const real_type right = grid[i + 1] - x;
    return std::min( base_h, real_type( 0.49 ) * std::min( left, right ) );
  }

  // Adaptive finite difference that respects node boundaries
  template <typename Func> real_type finite_diff_dx_adaptive(
    Func                      f,
    real_type                 x,
    real_type                 y,
    const vector<real_type> & X_grid,
    const vector<real_type> & Y_grid,
    real_type                 base_h = 1e-6 )
  {
    const integer   i = find_interval( X_grid, x );
    const real_type h = safe_step( x, X_grid, i, base_h );

    return ( f( x + h, y ) - f( x - h, y ) ) / ( 2 * h );
  }

  template <typename Func> real_type finite_diff_dy_adaptive(
    Func                      f,
    real_type                 x,
    real_type                 y,
    const vector<real_type> & X_grid,
    const vector<real_type> & Y_grid,
    real_type                 base_h = 1e-6 )
  {
    const integer   j = find_interval( Y_grid, y );
    const real_type h = safe_step( y, Y_grid, j, base_h );

    return ( f( x, y + h ) - f( x, y - h ) ) / ( 2 * h );
  }

  template <typename Func> real_type finite_diff_dxx_adaptive(
    Func                      f,
    real_type                 x,
    real_type                 y,
    const vector<real_type> & X_grid,
    const vector<real_type> & Y_grid,
    real_type                 base_h = 1e-6 )
  {
    const integer   i = find_interval( X_grid, x );
    const real_type h = safe_step( x, X_grid, i, base_h );

    return ( f( x + h, y ) - 2 * f( x, y ) + f( x - h, y ) ) / ( h * h );
  }

  template <typename Func> real_type finite_diff_dyy_adaptive(
    Func                      f,
    real_type                 x,
    real_type                 y,
    const vector<real_type> & X_grid,
    const vector<real_type> & Y_grid,
    real_type                 base_h = 1e-6 )
  {
    const integer   j = find_interval( Y_grid, y );
    const real_type h = safe_step( y, Y_grid, j, base_h );

    return ( f( x, y + h ) - 2 * f( x, y ) + f( x, y - h ) ) / ( h * h );
  }

  template <typename Func> real_type finite_diff_dxy_adaptive(
    Func                      f,
    real_type                 x,
    real_type                 y,
    const vector<real_type> & X_grid,
    const vector<real_type> & Y_grid,
    real_type                 base_h = 1e-6 )
  {
    const integer i = find_interval( X_grid, x );
    const integer j = find_interval( Y_grid, y );

    const real_type hx = safe_step( x, X_grid, i, base_h );
    const real_type hy = safe_step( y, Y_grid, j, base_h );
    const real_type h  = std::min( hx, hy );

    return ( f( x + h, y + h ) - f( x + h, y - h ) - f( x - h, y + h ) + f( x - h, y - h ) ) / ( 4 * h * h );
  }


  // Statistics structure
  struct Stats
  {
    real_type max_error{ 0 };
    real_type avg_error{ 0 };
    real_type rms_error{ 0 };
    integer   count{ 0 };

    void update( real_type error )
    {
      max_error = max( max_error, abs( error ) );
      avg_error += error;
      rms_error += error * error;
      ++count;
    }

    void finalize()
    {
      if ( count > 0 )
      {
        avg_error /= static_cast<real_type>( count );
        rms_error = sqrt( rms_error / static_cast<real_type>( count ) );
      }
    }
  };

  // Color formatting helpers
  string format_error( real_type error )
  {
    auto color = error < 1e-8 ? fmt::color::green : ( error < 1e-4 ? fmt::color::yellow : fmt::color::red );
    return fmt::format( fmt::fg( color ), "{:12.2e}", error );
  }

  string format_ratio( real_type ratio )
  {
    auto color = ratio < 1.1 ? fmt::color::green : ratio < 2.0 ? fmt::color::yellow : fmt::color::red;
    return fmt::format( fmt::fg( color ), "{:30.4g}", ratio );
  }

  // Test a specific spline type
  void test_spline_type(
    SplineType2D              spline_type,
    const vector<real_type> & x_grid,
    const vector<real_type> & y_grid,
    const vector<real_type> & z_data,
    integer                   n_test_points = 1000 )
  {
    string type_name;
    switch ( spline_type )
    {
      case SplineType2D::BILINEAR: type_name = "Bilinear"; break;
      case SplineType2D::BICUBIC: type_name = "Bicubic"; break;
      case SplineType2D::BIQUINTIC: type_name = "Biquintic"; break;
      case SplineType2D::AKIMA2D: type_name = "Akima2D"; break;
      default: type_name = "Unknown"; break;
    }

    fmt::print(
      fmt::fg( fmt::color::cyan ),
      "\n┌{0:─^{2}}┐\n"
      "│{1: ^{2}}│\n"
      "└{0:─^{2}}┘\n",
      "",
      fmt::format( "Testing {} Spline", type_name ),
      50 );

    // Create and build spline
    Spline2D spline( type_name );
    spline.build( spline_type, x_grid, y_grid, z_data );

    // Get bounds
    real_type x_min = spline.x_min();
    real_type x_max = spline.x_max();
    real_type y_min = spline.y_min();
    real_type y_max = spline.y_max();

    // Generate random test points
    mt19937                              rng( 1000 );
    uniform_real_distribution<real_type> dist_x( x_min, x_max );
    uniform_real_distribution<real_type> dist_y( y_min, y_max );

    // Statistics
    Stats stats_value, stats_dx, stats_dy, stats_dxx, stats_dyy, stats_dxy;
    Stats stats_dx_fd, stats_dy_fd, stats_dxx_fd, stats_dyy_fd, stats_dxy_fd;

    // Test points
    integer skipped_points = 0;
    for ( integer i_test = 0; i_test < n_test_points; ++i_test )
    {
      real_type x = dist_x( rng );
      real_type y = dist_y( rng );

      // Skip points too close to boundaries
      real_type h_min = 1e-6;
      if ( x <= x_min + h_min || x >= x_max - h_min || y <= y_min + h_min || y >= y_max - h_min )
      {
        ++skipped_points;
        continue;
      }

      // Find intervals to check if we're too close to nodes
      integer i_interval = find_interval( x_grid, x );
      integer j_interval = find_interval( y_grid, y );

      size_t i_size     = static_cast<size_t>( i_interval );
      size_t i_plus_one = i_size + 1;
      size_t j_size     = static_cast<size_t>( j_interval );
      size_t j_plus_one = j_size + 1;

      real_type dist_to_left   = x - x_grid[i_size];
      real_type dist_to_right  = x_grid[i_plus_one] - x;
      real_type dist_to_bottom = y - y_grid[j_size];
      real_type dist_to_top    = y_grid[j_plus_one] - y;

      // Skip if too close to any node for accurate finite differences
      real_type min_dist_threshold = 5e-4;
      if (
        dist_to_left < min_dist_threshold || dist_to_right < min_dist_threshold ||
        dist_to_bottom < min_dist_threshold || dist_to_top < min_dist_threshold )
      {
        ++skipped_points;
        continue;
      }

      // Exact values
      real_type exact_val = test_function( x, y );
      real_type exact_dx  = test_function_dx( x, y );
      real_type exact_dy  = test_function_dy( x, y );
      real_type exact_dxx = test_function_dxx( x, y );
      real_type exact_dyy = test_function_dyy( x, y );
      real_type exact_dxy = test_function_dxy( x, y );

      // Spline values
      real_type spline_val = spline.eval( x, y );
      real_type spline_dx  = spline.Dx( x, y );
      real_type spline_dy  = spline.Dy( x, y );
      real_type spline_dxx = spline.Dxx( x, y );
      real_type spline_dyy = spline.Dyy( x, y );
      real_type spline_dxy = spline.Dxy( x, y );

      // Adaptive finite difference approximations of spline
      auto spline_eval = [&]( real_type xx, real_type yy ) { return spline.eval( xx, yy ); };

      real_type fd_dx  = finite_diff_dx_adaptive( spline_eval, x, y, x_grid, y_grid );
      real_type fd_dy  = finite_diff_dy_adaptive( spline_eval, x, y, x_grid, y_grid );
      real_type fd_dxx = finite_diff_dxx_adaptive( spline_eval, x, y, x_grid, y_grid );
      real_type fd_dyy = finite_diff_dyy_adaptive( spline_eval, x, y, x_grid, y_grid );
      real_type fd_dxy = finite_diff_dxy_adaptive( spline_eval, x, y, x_grid, y_grid );

      // Update statistics
      stats_value.update( std::abs( spline_val - exact_val ) );
      stats_dx.update( std::abs( spline_dx - exact_dx ) );
      stats_dy.update( std::abs( spline_dy - exact_dy ) );
      stats_dxx.update( std::abs( spline_dxx - exact_dxx ) );
      stats_dyy.update( std::abs( spline_dyy - exact_dyy ) );
      stats_dxy.update( std::abs( spline_dxy - exact_dxy ) );

      stats_dx_fd.update( std::abs( spline_dx - fd_dx ) );
      stats_dy_fd.update( std::abs( spline_dy - fd_dy ) );
      stats_dxx_fd.update( std::abs( spline_dxx - fd_dxx ) );
      stats_dyy_fd.update( std::abs( spline_dyy - fd_dyy ) );
      stats_dxy_fd.update( std::abs( spline_dxy - fd_dxy ) );
    }

    // Finalize statistics
    stats_value.finalize();
    stats_dx.finalize();
    stats_dy.finalize();
    stats_dxx.finalize();
    stats_dyy.finalize();
    stats_dxy.finalize();
    stats_dx_fd.finalize();
    stats_dy_fd.finalize();
    stats_dxx_fd.finalize();
    stats_dyy_fd.finalize();
    stats_dxy_fd.finalize();

    fmt::print( "  Test points: {}, Skipped (near nodes/boundaries): {}\n", n_test_points, skipped_points );

    // Print results table
    fmt::print(
      "\n"
      "┌─────────────────────────────────────────────────────────────────────┐\n"
      "│ {:^67} │\n"
      "├────────────────────────┬──────────────┬──────────────┬──────────────┤\n"
      "│      Derivative        │  Max Error   │  Avg Error   │  RMS Error   │\n"
      "├────────────────────────┼──────────────┼──────────────┼──────────────┤\n",
      fmt::format( "{} Spline Error Analysis", type_name ) );

    // Function value
    fmt::print(
      "│ {:<22} │ {} │ {} │ {} │\n",
      "f(x,y)",
      format_error( stats_value.max_error ),
      format_error( stats_value.avg_error ),
      format_error( stats_value.rms_error ) );

    // First derivatives
    fmt::print(
      "│ {:<22} │ {} │ {} │ {} │\n",
      "∂f/∂x (vs exact)",
      format_error( stats_dx.max_error ),
      format_error( stats_dx.avg_error ),
      format_error( stats_dx.rms_error ) );

    fmt::print(
      "│ {:<22} │ {} │ {} │ {} │\n",
      "∂f/∂x (vs FD)",
      format_error( stats_dx_fd.max_error ),
      format_error( stats_dx_fd.avg_error ),
      format_error( stats_dx_fd.rms_error ) );

    fmt::print(
      "│ {:<22} │ {} │ {} │ {} │\n",
      "∂f/∂y (vs exact)",
      format_error( stats_dy.max_error ),
      format_error( stats_dy.avg_error ),
      format_error( stats_dy.rms_error ) );

    fmt::print(
      "│ {:<22} │ {} │ {} │ {} │\n",
      "∂f/∂y (vs FD)",
      format_error( stats_dy_fd.max_error ),
      format_error( stats_dy_fd.avg_error ),
      format_error( stats_dy_fd.rms_error ) );

    // Second derivatives (if applicable)
    if ( spline_type != SplineType2D::BILINEAR )
    {
      fmt::print(
        "│ {:<22} │ {} │ {} │ {} │\n",
        "∂²f/∂x² (vs exact)",
        format_error( stats_dxx.max_error ),
        format_error( stats_dxx.avg_error ),
        format_error( stats_dxx.rms_error ) );

      fmt::print(
        "│ {:<22} │ {} │ {} │ {} │\n",
        "∂²f/∂x² (vs FD)",
        format_error( stats_dxx_fd.max_error ),
        format_error( stats_dxx_fd.avg_error ),
        format_error( stats_dxx_fd.rms_error ) );

      fmt::print(
        "│ {:<22} │ {} │ {} │ {} │\n",
        "∂²f/∂y² (vs exact)",
        format_error( stats_dyy.max_error ),
        format_error( stats_dyy.avg_error ),
        format_error( stats_dyy.rms_error ) );

      fmt::print(
        "│ {:<22} │ {} │ {} │ {} │\n",
        "∂²f/∂y² (vs FD)",
        format_error( stats_dyy_fd.max_error ),
        format_error( stats_dyy_fd.avg_error ),
        format_error( stats_dyy_fd.rms_error ) );

      fmt::print(
        "│ {:<22} │ {} │ {} │ {} │\n",
        "∂²f/∂x∂y (vs exact)",
        format_error( stats_dxy.max_error ),
        format_error( stats_dxy.avg_error ),
        format_error( stats_dxy.rms_error ) );

      fmt::print(
        "│ {:<22} │ {} │ {} │ {} │\n",
        "∂²f/∂x∂y (vs FD)",
        format_error( stats_dxy_fd.max_error ),
        format_error( stats_dxy_fd.avg_error ),
        format_error( stats_dxy_fd.rms_error ) );
    }
    else
    {
      // For bilinear spline, derivatives are piecewise constant
      fmt::print( "│ {:<22} │ {:>12} │ {:>12} │ {:>12} │\n", "∂²f/∂x² (vs exact)", "N/A", "N/A", "N/A" );
      fmt::print( "│ {:<22} │ {:>12} │ {:>12} │ {:>12} │\n", "∂²f/∂y² (vs exact)", "N/A", "N/A", "N/A" );
      fmt::print( "│ {:<22} │ {:>12} │ {:>12} │ {:>12} │\n", "∂²f/∂x∂y (vs exact)", "N/A", "N/A", "N/A" );
    }

    fmt::print( "└────────────────────────┴──────────────┴──────────────┴──────────────┘\n" );

    // Print consistency ratios (Spline vs Finite Difference)
    fmt::print(
      "\n"
      "┌─────────────────────────────────────────────────────────┐\n"
      "│     Consistency Check (Spline vs Finite Difference)     │\n"
      "├────────────────────────┬────────────────────────────────┤\n"
      "│      Derivative        │   Error Ratio (Spline/FD)      │\n"
      "├────────────────────────┼────────────────────────────────┤\n" );

    real_type ratio_dx = ( stats_dx_fd.rms_error > 0 ) ? stats_dx.rms_error / stats_dx_fd.rms_error : 0;
    real_type ratio_dy = ( stats_dy_fd.rms_error > 0 ) ? stats_dy.rms_error / stats_dy_fd.rms_error : 0;

    fmt::print( "│ {:<22} │ {} │\n", "∂f/∂x", format_ratio( ratio_dx ) );

    fmt::print( "│ {:<22} │ {} │\n", "∂f/∂y", format_ratio( ratio_dy ) );

    if ( spline_type != SplineType2D::BILINEAR )
    {
      real_type ratio_dxx = ( stats_dxx_fd.rms_error > 0 ) ? stats_dxx.rms_error / stats_dxx_fd.rms_error : 0;
      real_type ratio_dyy = ( stats_dyy_fd.rms_error > 0 ) ? stats_dyy.rms_error / stats_dyy_fd.rms_error : 0;
      real_type ratio_dxy = ( stats_dxy_fd.rms_error > 0 ) ? stats_dxy.rms_error / stats_dxy_fd.rms_error : 0;

      fmt::print( "│ {:<22} │ {:>30} │\n", "∂²f/∂x²", format_ratio( ratio_dxx ) );

      fmt::print( "│ {:<22} │ {:>30} │\n", "∂²f/∂y²", format_ratio( ratio_dyy ) );

      fmt::print( "│ {:<22} │ {:>30} │\n", "∂²f/∂x∂y", format_ratio( ratio_dxy ) );
    }

    fmt::print( "└────────────────────────┴────────────────────────────────┘\n" );
  }

}  // namespace SplinesTest

int main()
{
  using namespace SplinesTest;

  fmt::print(
    fmt::fg( fmt::color::light_blue ) | fmt::emphasis::bold,
    "\n"
    "╔══════════════════════════════════════════════════════════╗\n"
    "║         2D Spline Library - Comprehensive Test           ║\n"
    "╚══════════════════════════════════════════════════════════╝\n\n" );

  // Create a non-uniform grid for testing
  integer nx = 113;
  integer ny = 121;

  // Use size_t for vector sizes to avoid warnings
  size_t nx_size = static_cast<size_t>( nx );
  size_t ny_size = static_cast<size_t>( ny );

  vector<real_type> x_grid( nx_size );
  vector<real_type> y_grid( ny_size );
  vector<real_type> z_data( nx_size * ny_size );

  real_type XMIN = -3;
  real_type XMAX = 3;
  real_type YMIN = -3;
  real_type YMAX = 3;

  // Create a slightly non-uniform grid
  for ( integer i = 0; i < nx; ++i )
  {
    size_t    i_idx = static_cast<size_t>( i );
    real_type t     = real_type( i ) / ( nx - 1 );
    x_grid[i_idx]   = XMIN + ( XMAX - XMIN ) * ( t + 0.1 * sin( 2.0 * M_PI * t ) );
  }

  for ( integer j = 0; j < ny; ++j )
  {
    size_t    j_idx = static_cast<size_t>( j );
    real_type t     = real_type( j ) / ( ny - 1 );
    y_grid[j_idx]   = YMIN + ( YMAX - YMIN ) * ( t + 0.1 * cos( 2.0 * M_PI * t ) );
  }

  // Sort grids to ensure monotonicity
  sort( x_grid.begin(), x_grid.end() );
  sort( y_grid.begin(), y_grid.end() );

  // Generate test data
  for ( integer i = 0; i < nx; ++i )
  {
    size_t i_idx = static_cast<size_t>( i );
    for ( integer j = 0; j < ny; ++j )
    {
      size_t j_idx = static_cast<size_t>( j );
      size_t idx   = i_idx * ny_size + j_idx;
      z_data[idx]  = test_function( x_grid[i_idx], y_grid[j_idx] );
    }
  }

  fmt::print(
    "Test grid: {} x {} points, domain: x∈[{:.2f}, {:.2f}], y∈[{:.2f}, {:.2f}]\n\n",
    nx,
    ny,
    x_grid.front(),
    x_grid.back(),
    y_grid.front(),
    y_grid.back() );

  // Test all spline types
  vector<SplineType2D> spline_types = { SplineType2D::BILINEAR,
                                        SplineType2D::BICUBIC,
                                        SplineType2D::BIQUINTIC,
                                        SplineType2D::AKIMA2D };

  for ( auto type : spline_types ) { test_spline_type( type, x_grid, y_grid, z_data, 2000 ); }

  // Additional test: Evaluate at grid points (should be exact for interpolating splines)
  fmt::print(
    fmt::fg( fmt::color::light_blue ) | fmt::emphasis::bold,
    "\n\n"
    "╔══════════════════════════════════════════════════════════╗\n"
    "║          Interpolation Accuracy at Grid Points           ║\n"
    "╚══════════════════════════════════════════════════════════╝\n" );

  Spline2D test_spline( "TestSpline" );
  test_spline.build( SplineType2D::BIQUINTIC, x_grid, y_grid, z_data );

  real_type max_interp_error = 0.0;
  real_type avg_interp_error = 0.0;
  integer   count            = 0;

  for ( integer i = 0; i < nx; ++i )
  {
    size_t i_idx = static_cast<size_t>( i );
    for ( integer j = 0; j < ny; ++j )
    {
      size_t    j_idx  = static_cast<size_t>( j );
      size_t    idx    = i_idx * ny_size + j_idx;
      real_type exact  = z_data[idx];
      real_type interp = test_spline.eval( x_grid[i_idx], y_grid[j_idx] );
      real_type error  = abs( interp - exact );

      max_interp_error = max( max_interp_error, error );
      avg_interp_error += error;
      ++count;
    }
  }

  avg_interp_error /= static_cast<real_type>( count );

  fmt::print(
    "\n"
    "Grid interpolation test (Biquintic spline):\n"
    "  Maximum error: {}\n"
    "  Average error: {}\n",
    format_error( max_interp_error ),
    format_error( avg_interp_error ) );

  if ( max_interp_error > 1e-10 )
  {
    fmt::print( fmt::fg( fmt::color::yellow ), "  Warning: Interpolation error larger than expected!\n" );
  }
  else
  {
    fmt::print( fmt::fg( fmt::color::green ), "  ✓ Interpolation is accurate at grid points\n" );
  }

  // Test derivative consistency using D() and DD() methods
  fmt::print(
    fmt::fg( fmt::color::light_blue ) | fmt::emphasis::bold,
    "\n\n"
    "╔══════════════════════════════════════════════════════════╗\n"
    "║          Derivative Method Consistency Test              ║\n"
    "╚══════════════════════════════════════════════════════════╝\n" );

  mt19937                              rng( 123 );
  uniform_real_distribution<real_type> dist_x( x_grid.front(), x_grid.back() );
  uniform_real_distribution<real_type> dist_y( y_grid.front(), y_grid.back() );

  real_type max_D_consistency_error  = 0.0;
  real_type max_DD_consistency_error = 0.0;

  for ( integer i = 0; i < 100; ++i )
  {
    real_type x = dist_x( rng );
    real_type y = dist_y( rng );

    // Skip points too close to boundaries
    real_type h_min = 1e-6;
    if (
      x <= x_grid.front() + h_min || x >= x_grid.back() - h_min || y <= y_grid.front() + h_min ||
      y >= y_grid.back() - h_min )
      continue;

    // Test D() method consistency
    real_type d[3];
    test_spline.D( x, y, d );

    real_type dx_indiv = test_spline.Dx( x, y );
    real_type dy_indiv = test_spline.Dy( x, y );

    real_type error_Dx = abs( d[1] - dx_indiv );
    real_type error_Dy = abs( d[2] - dy_indiv );

    max_D_consistency_error = max( max_D_consistency_error, max( error_Dx, error_Dy ) );

    // Test DD() method consistency
    real_type dd[6];
    test_spline.DD( x, y, dd );

    real_type dxx_indiv = test_spline.Dxx( x, y );
    real_type dyy_indiv = test_spline.Dyy( x, y );
    real_type dxy_indiv = test_spline.Dxy( x, y );

    real_type error_Dxx = abs( dd[3] - dxx_indiv );
    real_type error_Dyy = abs( dd[5] - dyy_indiv );
    real_type error_Dxy = abs( dd[4] - dxy_indiv );

    max_DD_consistency_error = max( max_DD_consistency_error, max( { error_Dxx, error_Dyy, error_Dxy } ) );
  }

  fmt::print(
    "\n"
    "Method consistency test (Biquintic spline):\n"
    "  D() method consistency error:  {}\n"
    "  DD() method consistency error: {}\n",
    format_error( max_D_consistency_error ),
    format_error( max_DD_consistency_error ) );

  if ( max_D_consistency_error < 1e-12 && max_DD_consistency_error < 1e-12 )
  {
    fmt::print( fmt::fg( fmt::color::green ), "  ✓ All derivative methods are consistent\n" );
  }

  fmt::print(
    fmt::fg( fmt::color::light_green ) | fmt::emphasis::bold,
    "\n\n"
    "════════════════════════════════════════════════════════════\n"
    "           All tests completed successfully!                \n"
    "════════════════════════════════════════════════════════════\n" );

  return 0;
}
