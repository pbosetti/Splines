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
#include "Utils_fmt.hh"

#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include <cmath>
#include <numeric>
#include <iomanip>
#include <sstream>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#endif

using namespace SplinesLoad;
using namespace std;
using Splines::integer;
using Splines::real_type;

static real_type x_data[] = { 0, 1, 2, 3, 4, 5 };
static real_type y_data[] = { 0, 1, 2, 3, 4, 5 };
static real_type z_data[] = { 10, 10, 10, 10, 10, 11, 10, 10, 10.5, 11, 9, 9, 11, 12, 13, 8, 7, 9,
                              11, 12, 13, 8,  7,  9,  11, 12, 13,   8,  7, 9, 11, 12, 13, 8, 7, 9 };

// Structure to hold 2D spline information
struct Spline2DInfo
{
  string    name;
  real_type x_min;
  real_type x_max;
  real_type y_min;
  real_type y_max;
  real_type z_min;
  real_type z_max;
  real_type z_avg;
};

// Structure to hold derivative errors statistics for interior and boundary points
struct DerivativeErrorStats
{
  real_type max_abs_error_interior;
  real_type avg_abs_error_interior;
  size_t    points_checked_interior;

  real_type max_abs_error_boundary;
  real_type avg_abs_error_boundary;
  size_t    points_checked_boundary;

  DerivativeErrorStats()
    : max_abs_error_interior( 0 )
    , avg_abs_error_interior( 0 )
    , points_checked_interior( 0 )
    , max_abs_error_boundary( 0 )
    , avg_abs_error_boundary( 0 )
    , points_checked_boundary( 0 )
  {
  }
};

// Structure to hold all derivative errors
struct DerivativeErrors
{
  DerivativeErrorStats Dx;
  DerivativeErrorStats Dy;
  DerivativeErrorStats Dxx;
  DerivativeErrorStats Dyy;
  DerivativeErrorStats Dxy;
};

// Helper function to format array as string
string format_array( const real_type * arr, size_t n )
{
  stringstream ss;
  for ( size_t i = 0; i < n; ++i )
  {
    if ( i > 0 ) ss << ", ";
    ss << fixed << setprecision( 2 ) << arr[i];
  }
  return ss.str();
}

// Function to print colored header
void print_header( const string & title )
{
  fmt::print(
    fg( fmt::color::cyan ) | fmt::emphasis::bold,
    "\n"
    "╔══════════════════════════════════════════════╗\n"
    "║{:^46}║\n"
    "╚══════════════════════════════════════════════╝\n",
    title );
}

// Function to print a table of results for 2D splines
void print_2D_spline_table( const vector<Spline2DInfo> & results )
{
  // Table header
  fmt::print(
    fg( fmt::color::yellow ) | fmt::emphasis::bold,
    "\n"
    "┌──────────────┬──────────┬──────────┬──────────┬──────────┬──────────┬──────────┬──────────┐\n"
    "│ {:^12} │ {:^8} │ {:^8} │ {:^8} │ {:^8} │ {:^8} │ {:^8} │ {:^8} │\n"
    "├──────────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┤\n",
    "Spline Type",
    "x_min",
    "x_max",
    "y_min",
    "y_max",
    "z_min",
    "z_max",
    "z_avg" );

  // Table rows with alternating colors
  for ( size_t i = 0; i < results.size(); ++i )
  {
    const auto & r         = results[i];
    auto         row_color = ( i % 2 == 0 ) ? fg( fmt::color::light_green ) : fg( fmt::color::light_blue );

    fmt::print( row_color, "│ {:<12} ", r.name );
    fmt::print( row_color, "│ {:>8.3f} ", r.x_min );
    fmt::print( row_color, "│ {:>8.3f} ", r.x_max );
    fmt::print( row_color, "│ {:>8.3f} ", r.y_min );
    fmt::print( row_color, "│ {:>8.3f} ", r.y_max );

    // Highlight z_min and z_max with different colors based on value
    auto z_min_color = ( r.z_min < 0 ) ? fg( fmt::color::red ) : fg( fmt::color::green );
    auto z_max_color = ( r.z_max < 0 ) ? fg( fmt::color::red ) : fg( fmt::color::green );
    auto z_avg_color = fg( fmt::color::orange );

    fmt::print( z_min_color, "│ {:>8.3f} ", r.z_min );
    fmt::print( z_max_color, "│ {:>8.3f} ", r.z_max );
    fmt::print( z_avg_color, "│ {:>8.3f} │\n", r.z_avg );
  }

  fmt::print(
    fg( fmt::color::yellow ) | fmt::emphasis::bold,
    "└──────────────┴──────────┴──────────┴──────────┴──────────┴──────────┴──────────┴──────────┘\n" );
}

// Function to format error with color and location (boundary vs interior)
string format_cell_with_stats(
  real_type      max_err,
  real_type      avg_err,
  size_t         points,
  const string & location,
  real_type      tol_max,
  real_type      tol_avg )
{
  stringstream ss;
  ss << location << ": ";

  // Determine color based on max error
  if ( max_err < tol_max && avg_err < tol_avg )
  {
    return fmt::format(
      fg( fmt::color::green ) | fmt::emphasis::bold,
      "{:52}",
      fmt::format( "{:<6} max:{:>8.1e} avg:{:>8.1e} ({:>3} pts)", location, max_err, avg_err, points ) );
  }
  else if ( max_err < 1e-4 && avg_err < 1e-5 )
  {
    return fmt::format(
      fg( fmt::color::light_green ),
      "{:52}",
      fmt::format( "{:<6} max:{:>8.1e} avg:{:>8.1e} ({:>3} pts)", location, max_err, avg_err, points ) );
  }
  else if ( max_err < 1e-2 && avg_err < 1e-3 )
  {
    return fmt::format(
      fg( fmt::color::yellow ),
      "{:52}",
      fmt::format( "{:<6} max:{:>8.1e} avg:{:>8.1e} ({:>3} pts)", location, max_err, avg_err, points ) );
  }
  else if ( max_err < 1.0 && avg_err < 0.1 )
  {
    return fmt::format(
      fg( fmt::color::orange ),
      "{:52}",
      fmt::format( "{:<6} max:{:>8.1e} avg:{:>8.1e} ({:>3} pts)", location, max_err, avg_err, points ) );
  }
  else
  {
    return fmt::format(
      fg( fmt::color::red ) | fmt::emphasis::bold,
      "{:52}",
      fmt::format( "{:<6} max:{:>8.1e} avg:{:>8.1e} ({:>3} pts)", location, max_err, avg_err, points ) );
  }
}

// Function to print a table for a specific derivative
void print_derivative_table(
  const string &                                 derivative_name,
  const vector<pair<string, DerivativeErrors>> & errors,
  real_type                                      tol_max_interior,
  real_type                                      tol_avg_interior,
  real_type                                      tol_max_boundary,
  real_type                                      tol_avg_boundary )
{
  auto color = fg( fmt::color::light_blue ) | fmt::emphasis::bold;
  fmt::print(
    color,
    "\n"
    "┌─{:─^122}─┐\n"
    "│ {:^122} │\n"
    "├─{:─^12}─┬─{:─^52}─┬─{:─^52}─┤\n"
    "│ {:^12} │ {:^52} │ {:^52} │\n"
    "├─{:─^12}─┼─{:─^52}─┼─{:─^52}─┤\n",
    "",
    fmt::format( "DERIVATIVE {} - FINITE DIFFERENCE ERRORS", derivative_name ),
    "",
    "",
    "",
    "Spline Type",
    "Interior Points (within patches)",
    "Boundary Points (between patches)",
    "",
    "",
    "" );

  // Function to get the appropriate error stats for the derivative
  auto get_error_stats = [&]( const DerivativeErrors & err, const string & deriv ) -> const DerivativeErrorStats &
  {
    if ( deriv == "Dx" ) return err.Dx;
    if ( deriv == "Dy" ) return err.Dy;
    if ( deriv == "Dxx" ) return err.Dxx;
    if ( deriv == "Dyy" ) return err.Dyy;
    if ( deriv == "Dxy" ) return err.Dxy;
    static DerivativeErrorStats empty{};
    return empty;
  };

  // Table rows
  for ( size_t i = 0; i < errors.size(); ++i )
  {
    const auto & [name, err] = errors[i];
    const auto & stats       = get_error_stats( err, derivative_name );
    auto         row_color   = ( i % 2 == 0 ) ? fg( fmt::color::light_green ) : fg( fmt::color::light_blue );

    // Interior errors
    string interior_cell = format_cell_with_stats(
      stats.max_abs_error_interior,
      stats.avg_abs_error_interior,
      stats.points_checked_interior,
      "Int",
      tol_max_interior,
      tol_avg_interior );

    // Boundary errors
    string boundary_cell = format_cell_with_stats(
      stats.max_abs_error_boundary,
      stats.avg_abs_error_boundary,
      stats.points_checked_boundary,
      "Bnd",
      tol_max_boundary,
      tol_avg_boundary );

    fmt::print( row_color, "│ {:<12} ", name );
    fmt::print( "│ {} ", interior_cell );
    fmt::print( "│ {} │\n", boundary_cell );

    if ( i < errors.size() - 1 ) { fmt::print( color, "├─{:─^12}─┼─{:─^52}─┼─{:─^52}─┤\n", "", "", "" ); }
  }

  fmt::print( color, "└─{:─^12}─┴─{:─^52}─┴─{:─^52}─┘\n", "", "", "" );
}

// Finite difference approximations

constexpr real_type eps = std::numeric_limits<real_type>::epsilon();

inline real_type h_first( real_type x )
{
  return std::cbrt( eps ) * std::max( real_type( 1 ), std::abs( x ) );
}

inline real_type h_second( real_type x )
{
  return std::sqrt( std::sqrt( eps ) ) * std::max( real_type( 1 ), std::abs( x ) );
}

namespace FiniteDifferences
{
  using real_type = double;

  namespace detail
  {
    constexpr real_type eps = std::numeric_limits<real_type>::epsilon();

    inline real_type h_first( real_type x, real_type y )
    {
      real_type scale = std::max( { real_type( 1 ), std::abs( x ), std::abs( y ) } );
      return std::cbrt( eps ) * scale;  // ~ 6e-6
    }

    inline real_type h_second( real_type x, real_type y )
    {
      real_type scale = std::max( { real_type( 1 ), std::abs( x ), std::abs( y ) } );
      return std::sqrt( std::sqrt( eps ) ) * scale;  // ~ 1e-4
    }
  }  // namespace detail

  // =========================
  // Derivate prime O(h^4)
  // =========================

  template <typename SplineType> real_type dx( const SplineType & spline, real_type x, real_type y )
  {
    real_type h = detail::h_first( x, y );
    return ( -spline( x + 2 * h, y ) + 8 * spline( x + h, y ) - 8 * spline( x - h, y ) + spline( x - 2 * h, y ) ) /
           ( 12 * h );
  }

  template <typename SplineType> real_type dy( const SplineType & spline, real_type x, real_type y )
  {
    real_type h = detail::h_first( x, y );
    return ( -spline( x, y + 2 * h ) + 8 * spline( x, y + h ) - 8 * spline( x, y - h ) + spline( x, y - 2 * h ) ) /
           ( 12 * h );
  }

  // =========================
  // Derivate seconde O(h^2)
  // =========================

  template <typename SplineType> real_type dxx( const SplineType & spline, real_type x, real_type y )
  {
    real_type h = detail::h_second( x, y );
    return ( spline( x + h, y ) - 2 * spline( x, y ) + spline( x - h, y ) ) / ( h * h );
  }

  template <typename SplineType> real_type dyy( const SplineType & spline, real_type x, real_type y )
  {
    real_type h = detail::h_second( x, y );
    return ( spline( x, y + h ) - 2 * spline( x, y ) + spline( x, y - h ) ) / ( h * h );
  }

  template <typename SplineType> real_type dxy( const SplineType & spline, real_type x, real_type y )
  {
    real_type h = detail::h_second( x, y );
    return ( spline( x + h, y + h ) - spline( x + h, y - h ) - spline( x - h, y + h ) + spline( x - h, y - h ) ) /
           ( 4 * h * h );
  }

}  // namespace FiniteDifferences


// Function to check if a point is on a patch boundary (including internal patch boundaries)
bool is_on_patch_boundary(
  real_type         x,
  real_type         y,
  const real_type * X,
  integer           nx,
  const real_type * Y,
  integer           ny,
  real_type         tolerance = 1e-6 )
{
  // Check if point is close to any knot in X direction
  for ( integer i = 0; i < nx; ++i )
  {
    if ( fabs( x - X[i] ) < tolerance ) { return true; }
  }

  // Check if point is close to any knot in Y direction
  for ( integer j = 0; j < ny; ++j )
  {
    if ( fabs( y - Y[j] ) < tolerance ) { return true; }
  }

  return false;
}

// Function to check derivatives using finite differences with patch boundary separation
template <typename SplineType> DerivativeErrors check_derivatives_with_patch_separation(
  SplineType &      spline,
  const string &    spline_name,
  const real_type * X_knots,
  integer           nx_knots,
  const real_type * Y_knots,
  integer           ny_knots,
  integer           grid_points_per_patch = 10,
  real_type         h                     = 1e-6 )
{
  DerivativeErrors errors{};

  // Initialize all statistics to zero
  auto init_stats = []() -> DerivativeErrorStats
  {
    DerivativeErrorStats stats;
    stats.max_abs_error_interior  = 0;
    stats.avg_abs_error_interior  = 0;
    stats.points_checked_interior = 0;
    stats.max_abs_error_boundary  = 0;
    stats.avg_abs_error_boundary  = 0;
    stats.points_checked_boundary = 0;
    return stats;
  };

  errors.Dx  = init_stats();
  errors.Dy  = init_stats();
  errors.Dxx = init_stats();
  errors.Dyy = init_stats();
  errors.Dxy = init_stats();

  real_type x_min = spline.x_min();
  real_type x_max = spline.x_max();
  real_type y_min = spline.y_min();
  real_type y_max = spline.y_max();

  // Generate test points with emphasis on patch boundaries
  vector<real_type> test_points_x, test_points_y;
  vector<bool>      is_boundary_point;

  // Generate points on patch boundaries (internal boundaries)
  for ( integer i = 0; i < nx_knots; ++i )
  {
    real_type x = X_knots[i];
    // Generate points along this x-boundary
    for ( integer j = 0; j <= grid_points_per_patch; ++j )
    {
      real_type y = y_min + ( y_max - y_min ) * j / grid_points_per_patch;
      if ( y >= y_min + h && y <= y_max - h && x >= x_min + h && x <= x_max - h )
      {
        test_points_x.push_back( x );
        test_points_y.push_back( y );
        is_boundary_point.push_back( true );
      }
    }
  }

  for ( integer j = 0; j < ny_knots; ++j )
  {
    real_type y = Y_knots[j];
    // Generate points along this y-boundary
    for ( integer i = 0; i <= grid_points_per_patch; ++i )
    {
      real_type x = x_min + ( x_max - x_min ) * i / grid_points_per_patch;
      if ( x >= x_min + h && x <= x_max - h && y >= y_min + h && y <= y_max - h )
      {
        test_points_x.push_back( x );
        test_points_y.push_back( y );
        is_boundary_point.push_back( true );
      }
    }
  }

  // Generate interior points (within patches, away from boundaries)
  integer patches_x = nx_knots - 1;
  integer patches_y = ny_knots - 1;

  for ( integer patch_i = 0; patch_i < patches_x; ++patch_i )
  {
    for ( integer patch_j = 0; patch_j < patches_y; ++patch_j )
    {
      real_type patch_x_min = X_knots[patch_i];
      real_type patch_x_max = X_knots[patch_i + 1];
      real_type patch_y_min = Y_knots[patch_j];
      real_type patch_y_max = Y_knots[patch_j + 1];

      // Generate interior points for this patch
      for ( integer i = 1; i < grid_points_per_patch; ++i )
      {
        for ( integer j = 1; j < grid_points_per_patch; ++j )
        {
          real_type x = patch_x_min + ( patch_x_max - patch_x_min ) * i / grid_points_per_patch;
          real_type y = patch_y_min + ( patch_y_max - patch_y_min ) * j / grid_points_per_patch;

          // Ensure we're not too close to patch boundaries
          real_type x_dist   = min( x - patch_x_min, patch_x_max - x );
          real_type y_dist   = min( y - patch_y_min, patch_y_max - y );
          real_type min_dist = min( x_dist, y_dist );

          if (
            min_dist > 0.1 * ( patch_x_max - patch_x_min ) && min_dist > 0.1 * ( patch_y_max - patch_y_min ) &&
            x >= x_min + h && x <= x_max - h && y >= y_min + h && y <= y_max - h )
          {
            test_points_x.push_back( x );
            test_points_y.push_back( y );
            is_boundary_point.push_back( false );
          }
        }
      }
    }
  }

  fmt::print(
    fg( fmt::color::blue ),
    "   Checking derivatives for {} ({} points total, with patch separation)...",
    spline_name,
    test_points_x.size() );

  size_t interior_points = 0;
  size_t boundary_points = 0;

  for ( size_t idx = 0; idx < test_points_x.size(); ++idx )
  {
    real_type x_val       = test_points_x[idx];
    real_type y_val       = test_points_y[idx];
    bool      is_boundary = is_boundary_point[idx];

    // Get analytical derivatives
    real_type Dx_analytic  = spline.Dx( x_val, y_val );
    real_type Dy_analytic  = spline.Dy( x_val, y_val );
    real_type Dxx_analytic = spline.Dxx( x_val, y_val );
    real_type Dyy_analytic = spline.Dyy( x_val, y_val );
    real_type Dxy_analytic = spline.Dxy( x_val, y_val );

    // Get finite difference approximations
    real_type Dx_fd  = FiniteDifferences::dx( spline, x_val, y_val );
    real_type Dy_fd  = FiniteDifferences::dy( spline, x_val, y_val );
    real_type Dxx_fd = FiniteDifferences::dxx( spline, x_val, y_val );
    real_type Dyy_fd = FiniteDifferences::dyy( spline, x_val, y_val );
    real_type Dxy_fd = FiniteDifferences::dxy( spline, x_val, y_val );

    // Calculate absolute errors
    real_type abs_error_Dx  = fabs( Dx_analytic - Dx_fd );
    real_type abs_error_Dy  = fabs( Dy_analytic - Dy_fd );
    real_type abs_error_Dxx = fabs( Dxx_analytic - Dxx_fd );
    real_type abs_error_Dyy = fabs( Dyy_analytic - Dyy_fd );
    real_type abs_error_Dxy = fabs( Dxy_analytic - Dxy_fd );

    // Update statistics based on boundary status
    if ( is_boundary )
    {
      boundary_points++;

      // Update Dx boundary stats
      errors.Dx.max_abs_error_boundary = max( errors.Dx.max_abs_error_boundary, abs_error_Dx );
      errors.Dx.avg_abs_error_boundary += abs_error_Dx;
      errors.Dx.points_checked_boundary++;

      // Update Dy boundary stats
      errors.Dy.max_abs_error_boundary = max( errors.Dy.max_abs_error_boundary, abs_error_Dy );
      errors.Dy.avg_abs_error_boundary += abs_error_Dy;
      errors.Dy.points_checked_boundary++;

      // Update Dxx boundary stats
      errors.Dxx.max_abs_error_boundary = max( errors.Dxx.max_abs_error_boundary, abs_error_Dxx );
      errors.Dxx.avg_abs_error_boundary += abs_error_Dxx;
      errors.Dxx.points_checked_boundary++;

      // Update Dyy boundary stats
      errors.Dyy.max_abs_error_boundary = max( errors.Dyy.max_abs_error_boundary, abs_error_Dyy );
      errors.Dyy.avg_abs_error_boundary += abs_error_Dyy;
      errors.Dyy.points_checked_boundary++;

      // Update Dxy boundary stats
      errors.Dxy.max_abs_error_boundary = max( errors.Dxy.max_abs_error_boundary, abs_error_Dxy );
      errors.Dxy.avg_abs_error_boundary += abs_error_Dxy;
      errors.Dxy.points_checked_boundary++;
    }
    else
    {
      interior_points++;

      // Update Dx interior stats
      errors.Dx.max_abs_error_interior = max( errors.Dx.max_abs_error_interior, abs_error_Dx );
      errors.Dx.avg_abs_error_interior += abs_error_Dx;
      errors.Dx.points_checked_interior++;

      // Update Dy interior stats
      errors.Dy.max_abs_error_interior = max( errors.Dy.max_abs_error_interior, abs_error_Dy );
      errors.Dy.avg_abs_error_interior += abs_error_Dy;
      errors.Dy.points_checked_interior++;

      // Update Dxx interior stats
      errors.Dxx.max_abs_error_interior = max( errors.Dxx.max_abs_error_interior, abs_error_Dxx );
      errors.Dxx.avg_abs_error_interior += abs_error_Dxx;
      errors.Dxx.points_checked_interior++;

      // Update Dyy interior stats
      errors.Dyy.max_abs_error_interior = max( errors.Dyy.max_abs_error_interior, abs_error_Dyy );
      errors.Dyy.avg_abs_error_interior += abs_error_Dyy;
      errors.Dyy.points_checked_interior++;

      // Update Dxy interior stats
      errors.Dxy.max_abs_error_interior = max( errors.Dxy.max_abs_error_interior, abs_error_Dxy );
      errors.Dxy.avg_abs_error_interior += abs_error_Dxy;
      errors.Dxy.points_checked_interior++;
    }
  }

  // Calculate averages
  auto calculate_averages = []( DerivativeErrorStats & stats )
  {
    if ( stats.points_checked_interior > 0 ) { stats.avg_abs_error_interior /= stats.points_checked_interior; }
    if ( stats.points_checked_boundary > 0 ) { stats.avg_abs_error_boundary /= stats.points_checked_boundary; }
  };

  calculate_averages( errors.Dx );
  calculate_averages( errors.Dy );
  calculate_averages( errors.Dxx );
  calculate_averages( errors.Dyy );
  calculate_averages( errors.Dxy );

  fmt::print( fg( fmt::color::green ), " ✓ ({} interior, {} boundary points)\n", interior_points, boundary_points );
  return errors;
}

// Function to generate grid data and save to files
template <typename SplineType> void generate_grid_data(
  SplineType &   spline,
  const string & base_filename,
  const string & spline_name,
  integer        grid_size = 100 )
{
  // Array of filenames
  array<string, 6> filenames = {
    fmt::format( "out/{}.txt", base_filename ),     fmt::format( "out/{}_Dx.txt", base_filename ),
    fmt::format( "out/{}_Dy.txt", base_filename ),  fmt::format( "out/{}_Dxx.txt", base_filename ),
    fmt::format( "out/{}_Dyy.txt", base_filename ), fmt::format( "out/{}_Dxy.txt", base_filename )
  };

  // Array of output streams
  array<ofstream, 6> files;

  // Open all files
  bool all_opened = true;
  for ( size_t idx = 0; idx < filenames.size(); ++idx )
  {
    files[idx].open( filenames[idx].data() );
    if ( !files[idx].is_open() )
    {
      fmt::print( fg( fmt::color::red ) | fmt::emphasis::bold, "❌ Error opening file: {}\n", filenames[idx] );
      all_opened = false;
    }
  }

  if ( !all_opened ) { return; }

  real_type x_min = spline.x_min();
  real_type x_max = spline.x_max();
  real_type y_min = spline.y_min();
  real_type y_max = spline.y_max();

  real_type dx = ( x_max - x_min ) / grid_size;
  real_type dy = ( y_max - y_min ) / grid_size;

  // Generate grid data
  fmt::print( fg( fmt::color::blue ), "   Generating {}x{} grid for {}...", grid_size + 1, grid_size + 1, spline_name );

  for ( integer i = 0; i <= grid_size; ++i )
  {
    real_type xxx = x_min + i * dx;
    for ( integer j = 0; j <= grid_size; ++j )
    {
      real_type yyy = y_min + j * dy;

      // Write to different files
      fmt::print( files[0], "{}\t", spline( xxx, yyy ) );
      fmt::print( files[1], "{}\t", spline.Dx( xxx, yyy ) );
      fmt::print( files[2], "{}\t", spline.Dy( xxx, yyy ) );
      fmt::print( files[3], "{}\t", spline.Dxx( xxx, yyy ) );
      fmt::print( files[4], "{}\t", spline.Dyy( xxx, yyy ) );
      fmt::print( files[5], "{}\t", spline.Dxy( xxx, yyy ) );
    }

    // New line in each file
    for ( auto & file : files ) { file << '\n'; }
  }

  // Close all files
  for ( auto & file : files ) { file.close(); }

  fmt::print( fg( fmt::color::green ), " ✓\n" );
}

// Function to analyze a 2D spline and collect statistics
template <typename SplineType>
Spline2DInfo analyze_2D_spline( SplineType & spline, const string & name, integer sample_points = 50 )
{
  Spline2DInfo info;
  info.name  = name;
  info.x_min = spline.x_min();
  info.x_max = spline.x_max();
  info.y_min = spline.y_min();
  info.y_max = spline.y_max();

  // Sample the spline to find statistics
  real_type z_min = numeric_limits<real_type>::max();
  real_type z_max = numeric_limits<real_type>::lowest();
  real_type z_sum = 0;

  real_type dx = ( info.x_max - info.x_min ) / sample_points;
  real_type dy = ( info.y_max - info.y_min ) / sample_points;

  for ( integer i = 0; i <= sample_points; ++i )
  {
    real_type xxx = info.x_min + i * dx;
    for ( integer j = 0; j <= sample_points; ++j )
    {
      real_type yyy = info.y_min + j * dy;
      real_type zzz = spline( xxx, yyy );

      z_sum += zzz;
      z_min = min( z_min, zzz );
      z_max = max( z_max, zzz );
    }
  }

  info.z_min = z_min;
  info.z_max = z_max;
  info.z_avg = z_sum / ( ( sample_points + 1 ) * ( sample_points + 1 ) );

  return info;
}

int main()
{
  print_header( "2D SPLINE INTERPOLATION TEST SUITE" );

  fmt::print( fg( fmt::color::magenta ) | fmt::emphasis::bold, "\n📊 Dataset Information:\n" );
  fmt::print( fg( fmt::color::gray ), "   Grid size: 6x6 points (5x5 patches)\n" );

  // Format knots arrays manually
  string x_knots_str = format_array( x_data, 6 );
  string y_knots_str = format_array( y_data, 6 );

  fmt::print( fg( fmt::color::gray ), "   x knots: {}\n", x_knots_str );
  fmt::print( fg( fmt::color::gray ), "   y knots: {}\n", y_knots_str );
  fmt::print(
    fg( fmt::color::gray ),
    "   Domain: x ∈ [{}, {}], y ∈ [{}, {}]\n",
    x_data[0],
    x_data[5],
    y_data[0],
    y_data[5] );

  BiQuinticSpline bq;
  BiCubicSpline   bc;
  Akima2Dspline   ak;
  BilinearSpline  bl;

  real_type X[6], Y[6], Z[6 * 6];

  std::copy_n( x_data, 6, X );
  std::copy_n( y_data, 6, Y );
  std::copy_n( z_data, 6 * 6, Z );

  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n🔨 Building 2D splines...\n" );

  // Build all splines
  bq.build( X, 1, Y, 1, Z, 6, 6, 6 );
  fmt::print( fg( fmt::color::green ), "   ✓ BiQuinticSpline built\n" );

  bc.build( X, 1, Y, 1, Z, 6, 6, 6 );
  fmt::print( fg( fmt::color::green ), "   ✓ BiCubicSpline built\n" );

  ak.build( X, 1, Y, 1, Z, 6, 6, 6 );
  fmt::print( fg( fmt::color::green ), "   ✓ Akima2Dspline built\n" );

  bl.build( X, 1, Y, 1, Z, 6, 6, 6 );
  fmt::print( fg( fmt::color::green ), "   ✓ BilinearSpline built\n" );

  // Analyze all splines
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n📈 Analyzing spline properties...\n" );

  vector<Spline2DInfo> spline_results;
  spline_results.push_back( analyze_2D_spline( bq, "BiQuintic" ) );
  spline_results.push_back( analyze_2D_spline( bc, "BiCubic" ) );
  spline_results.push_back( analyze_2D_spline( ak, "Akima2D" ) );
  spline_results.push_back( analyze_2D_spline( bl, "Bilinear" ) );

  // Print results table
  print_2D_spline_table( spline_results );

  // Check derivatives with finite differences and patch separation
  fmt::print(
    fg( fmt::color::cyan ) | fmt::emphasis::bold,
    "\n🔍 Checking derivatives with finite differences (patch separation)...\n" );

  vector<pair<string, DerivativeErrors>> derivative_errors;
  derivative_errors.emplace_back(
    "BiQuintic",
    check_derivatives_with_patch_separation( bq, "BiQuinticSpline", X, 6, Y, 6, 10 ) );
  derivative_errors.emplace_back(
    "BiCubic",
    check_derivatives_with_patch_separation( bc, "BiCubicSpline", X, 6, Y, 6, 10 ) );
  derivative_errors.emplace_back(
    "Akima2D",
    check_derivatives_with_patch_separation( ak, "Akima2Dspline", X, 6, Y, 6, 10 ) );
  derivative_errors.emplace_back(
    "Bilinear",
    check_derivatives_with_patch_separation( bl, "BilinearSpline", X, 6, Y, 6, 10 ) );

  // Print separate tables for each derivative
  fmt::print( fg( fmt::color::magenta ) | fmt::emphasis::bold, "\n📊 DERIVATIVE ERROR TABLES (Absolute Errors):\n" );

  // Define tolerances for interior and boundary points
  real_type first_deriv_tol_max_interior = 1e-6;
  real_type first_deriv_tol_avg_interior = 1e-7;
  real_type first_deriv_tol_max_boundary = 1e-5;  // More relaxed for boundaries
  real_type first_deriv_tol_avg_boundary = 1e-6;

  real_type second_deriv_tol_max_interior = 1e-3;
  real_type second_deriv_tol_avg_interior = 1e-4;
  real_type second_deriv_tol_max_boundary = 1e-2;  // More relaxed for boundaries
  real_type second_deriv_tol_avg_boundary = 1e-3;

  print_derivative_table(
    "Dx",
    derivative_errors,
    first_deriv_tol_max_interior,
    first_deriv_tol_avg_interior,
    first_deriv_tol_max_boundary,
    first_deriv_tol_avg_boundary );

  print_derivative_table(
    "Dy",
    derivative_errors,
    first_deriv_tol_max_interior,
    first_deriv_tol_avg_interior,
    first_deriv_tol_max_boundary,
    first_deriv_tol_avg_boundary );

  print_derivative_table(
    "Dxx",
    derivative_errors,
    second_deriv_tol_max_interior,
    second_deriv_tol_avg_interior,
    second_deriv_tol_max_boundary,
    second_deriv_tol_avg_boundary );

  print_derivative_table(
    "Dyy",
    derivative_errors,
    second_deriv_tol_max_interior,
    second_deriv_tol_avg_interior,
    second_deriv_tol_max_boundary,
    second_deriv_tol_avg_boundary );

  print_derivative_table(
    "Dxy",
    derivative_errors,
    second_deriv_tol_max_interior,
    second_deriv_tol_avg_interior,
    second_deriv_tol_max_boundary,
    second_deriv_tol_avg_boundary );

  // Generate output files
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n💾 Writing spline data to files...\n" );

  generate_grid_data( bq, "biquintic", "BiQuinticSpline" );
  generate_grid_data( bc, "bicubic", "BiCubicSpline" );
  generate_grid_data( ak, "akima2d", "Akima2Dspline" );
  generate_grid_data( bl, "bilinear", "BilinearSpline" );

  // Print summary of generated files
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n📁 Generated Files Summary:\n" );

  vector<string> spline_types = { "bicubic", "biquintic", "akima2d", "bilinear" };
  vector<string> derivatives  = { "", "_Dx", "_Dy", "_Dxx", "_Dyy", "_Dxy" };

  for ( const auto & spline_type : spline_types )
  {
    fmt::print( fg( fmt::color::gray ), "   {}:\n", spline_type );
    for ( const auto & deriv : derivatives )
    {
      fmt::print( fg( fmt::color::gray ), "     - out/{}{}.txt\n", spline_type, deriv );
    }
  }

  print_header( "2D SPLINE TEST COMPLETED" );
  fmt::print( fg( fmt::color::light_green ) | fmt::emphasis::bold, "\n🎉 All 2D spline analyses completed! 🎉\n" );
  fmt::print( fg( fmt::color::gray ), "   Grid data saved to: out/*.txt\n" );
  fmt::print( fg( fmt::color::gray ), "   Use visualization tools to plot the results.\n\n" );

  // Summary of derivative checks with color coding
  fmt::print( fg( fmt::color::magenta ) | fmt::emphasis::bold, "📋 DERIVATIVE CHECK SUMMARY (by patch location):\n" );

  for ( const auto & [name, err] : derivative_errors )
  {
    fmt::print( fg( fmt::color::white ) | fmt::emphasis::bold, "   {}:\n", name );

    // First derivatives
    bool first_deriv_interior_ok =
      ( err.Dx.max_abs_error_interior < first_deriv_tol_max_interior &&
        err.Dx.avg_abs_error_interior < first_deriv_tol_avg_interior &&
        err.Dy.max_abs_error_interior < first_deriv_tol_max_interior &&
        err.Dy.avg_abs_error_interior < first_deriv_tol_avg_interior );

    bool first_deriv_boundary_ok =
      ( err.Dx.max_abs_error_boundary < first_deriv_tol_max_boundary &&
        err.Dx.avg_abs_error_boundary < first_deriv_tol_avg_boundary &&
        err.Dy.max_abs_error_boundary < first_deriv_tol_max_boundary &&
        err.Dy.avg_abs_error_boundary < first_deriv_tol_avg_boundary );

    // Second derivatives
    bool second_deriv_interior_ok =
      ( err.Dxx.max_abs_error_interior < second_deriv_tol_max_interior &&
        err.Dxx.avg_abs_error_interior < second_deriv_tol_avg_interior &&
        err.Dyy.max_abs_error_interior < second_deriv_tol_max_interior &&
        err.Dyy.avg_abs_error_interior < second_deriv_tol_avg_interior &&
        err.Dxy.max_abs_error_interior < second_deriv_tol_max_interior &&
        err.Dxy.avg_abs_error_interior < second_deriv_tol_avg_interior );

    bool second_deriv_boundary_ok =
      ( err.Dxx.max_abs_error_boundary < second_deriv_tol_max_boundary &&
        err.Dxx.avg_abs_error_boundary < second_deriv_tol_avg_boundary &&
        err.Dyy.max_abs_error_boundary < second_deriv_tol_max_boundary &&
        err.Dyy.avg_abs_error_boundary < second_deriv_tol_avg_boundary &&
        err.Dxy.max_abs_error_boundary < second_deriv_tol_max_boundary &&
        err.Dxy.avg_abs_error_boundary < second_deriv_tol_avg_boundary );

    // Print results
    if ( first_deriv_interior_ok )
    {
      fmt::print( fg( fmt::color::green ), "     ✓ First derivatives (interior): OK\n" );
    }
    else
    {
      fmt::print( fg( fmt::color::yellow ), "     ⚠ First derivatives (interior): Some errors\n" );
    }

    if ( first_deriv_boundary_ok )
    {
      fmt::print( fg( fmt::color::green ), "     ✓ First derivatives (boundary): OK\n" );
    }
    else
    {
      fmt::print( fg( fmt::color::yellow ), "     ⚠ First derivatives (boundary): Some errors\n" );
    }

    if ( second_deriv_interior_ok )
    {
      fmt::print( fg( fmt::color::green ), "     ✓ Second derivatives (interior): OK\n" );
    }
    else
    {
      fmt::print( fg( fmt::color::yellow ), "     ⚠ Second derivatives (interior): Some errors\n" );
    }

    if ( second_deriv_boundary_ok )
    {
      fmt::print( fg( fmt::color::green ), "     ✓ Second derivatives (boundary): OK\n" );
    }
    else
    {
      fmt::print( fg( fmt::color::yellow ), "     ⚠ Second derivatives (boundary): Some errors\n" );
    }
  }

  fmt::print( fg( fmt::color::gray ), "\n   Tolerance levels (absolute errors):\n" );
  fmt::print( fg( fmt::color::gray ), "     - First derivatives (Dx, Dy):\n" );
  fmt::print( fg( fmt::color::gray ), "       • Interior: max < 1e-6, avg < 1e-7\n" );
  fmt::print( fg( fmt::color::gray ), "       • Boundary: max < 1e-5, avg < 1e-6\n" );
  fmt::print( fg( fmt::color::gray ), "     - Second derivatives (Dxx, Dyy, Dxy):\n" );
  fmt::print( fg( fmt::color::gray ), "       • Interior: max < 1e-3, avg < 1e-4\n" );
  fmt::print( fg( fmt::color::gray ), "       • Boundary: max < 1e-2, avg < 1e-3\n" );

  print_header( "ALL DONE FOLKS" );

  return 0;
}
