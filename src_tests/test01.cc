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

using namespace SplinesLoad;
using namespace std;
using Splines::integer;
using Splines::real_type;

// Test datasets
static real_type xx0[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
static real_type yy0[] = { 10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85 };

static real_type xx1[] = { 0, 1, 3, 4, 6, 7, 9, 10, 12, 13, 15 };
static real_type yy1[] = { 10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85 };

static real_type xx2[] = { 0, 2, 3, 5, 6, 8, 9, 11, 12, 14, 15 };
static real_type yy2[] = { 10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85 };

// RPN 14
static real_type xx3[] = { 7.99, 8.09, 8.19, 8.7, 9.2, 10, 12, 15, 20 };
static real_type yy3[] = { 0, 2.76429e-5, 4.37498e-2, 0.169183, 0.469428, 0.943740, 0.998636, 0.999919, 0.999994 };

// Titanium
static real_type xx4[] = { 595, 635, 695, 795, 855, 875, 895, 915, 935, 985, 1035, 1075 };
static real_type yy4[] = { 0.644, 0.652, 0.644, 0.694, 0.907, 1.336, 2.169, 1.598, 0.916, 0.607, 0.603, 0.608 };

// toolpath
static real_type xx5[] = { 0.11, 0.12, 0.15, 0.16 };
static real_type yy5[] = { 0.0003, 0.0003, 0.0004, 0.0004 };

static integer nn[] = { 11, 11, 11, 9, 12, 4 };

// ============================================================================
// STRUCTURES FOR DERIVATIVE ERROR ANALYSIS
// ============================================================================

// Structure to hold derivative error statistics
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

// Structure to hold all derivative errors for 1D spline
struct DerivativeErrors1D
{
  DerivativeErrorStats D;   // First derivative
  DerivativeErrorStats DD;  // Second derivative

  DerivativeErrors1D() : D(), DD() {}
};

// Structure to hold spline min/max results
struct SplineResult
{
  string    name;
  real_type x_min;
  real_type x_max;
  real_type y_min;
  real_type y_max;
  real_type x_at_min;
  real_type x_at_max;
};

// ============================================================================
// FINITE DIFFERENCE FUNCTIONS FOR 1D
// ============================================================================

namespace FiniteDifferences1D
{
  constexpr real_type eps = std::numeric_limits<real_type>::epsilon();

  // Adaptive step size for first derivative
  inline real_type h_first( real_type x )
  {
    real_type scale = std::max( real_type( 1 ), std::abs( x ) );
    return std::cbrt( eps ) * scale;  // ~ 6e-6 for double
  }

  // Adaptive step size for second derivative
  inline real_type h_second( real_type x )
  {
    real_type scale = std::max( real_type( 1 ), std::abs( x ) );
    return std::sqrt( std::sqrt( eps ) ) * scale;  // ~ 1e-4 for double
  }

  // Fourth-order central difference for first derivative
  template <typename SplineType> real_type dx( const SplineType & spline, real_type x )
  {
    real_type h = h_first( x );
    return ( -spline( x + 2 * h ) + 8 * spline( x + h ) - 8 * spline( x - h ) + spline( x - 2 * h ) ) / ( 12 * h );
  }

  // Second-order central difference for second derivative
  template <typename SplineType> real_type dxx( const SplineType & spline, real_type x )
  {
    real_type h = h_second( x );
    return ( spline( x + h ) - 2 * spline( x ) + spline( x - h ) ) / ( h * h );
  }
}  // namespace FiniteDifferences1D

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

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

// Helper function to format error with color and location
string format_cell_with_stats(
  real_type      max_err,
  real_type      avg_err,
  size_t         points,
  const string & location,
  real_type      tol_max,
  real_type      tol_avg )
{
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

// Function to print derivative error table
void print_derivative_table_1D(
  const string &                                   derivative_name,
  const vector<pair<string, DerivativeErrors1D>> & errors,
  integer                                          dataset )
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
    fmt::format( "DATASET {} - DERIVATIVE {} - FINITE DIFFERENCE ERRORS", dataset, derivative_name ),
    "",
    "",
    "",
    "Spline Type",
    "Interior Points (within segments)",
    "Boundary Points (at knots)",
    "",
    "",
    "" );

  // Tolerances
  real_type first_deriv_tol_max_interior = 1e-6;
  real_type first_deriv_tol_avg_interior = 1e-7;
  real_type first_deriv_tol_max_boundary = 1e-5;
  real_type first_deriv_tol_avg_boundary = 1e-6;

  real_type second_deriv_tol_max_interior = 1e-3;
  real_type second_deriv_tol_avg_interior = 1e-4;
  real_type second_deriv_tol_max_boundary = 1e-2;
  real_type second_deriv_tol_avg_boundary = 1e-3;

  // Table rows
  for ( size_t i = 0; i < errors.size(); ++i )
  {
    const auto & [name, err]           = errors[i];
    const DerivativeErrorStats & stats = ( derivative_name == "D" ) ? err.D : err.DD;

    auto row_color = ( i % 2 == 0 ) ? fg( fmt::color::light_green ) : fg( fmt::color::light_blue );

    // Choose tolerances based on derivative
    real_type tol_max_int = ( derivative_name == "D" ) ? first_deriv_tol_max_interior : second_deriv_tol_max_interior;
    real_type tol_avg_int = ( derivative_name == "D" ) ? first_deriv_tol_avg_interior : second_deriv_tol_avg_interior;
    real_type tol_max_bnd = ( derivative_name == "D" ) ? first_deriv_tol_max_boundary : second_deriv_tol_max_boundary;
    real_type tol_avg_bnd = ( derivative_name == "D" ) ? first_deriv_tol_avg_boundary : second_deriv_tol_avg_boundary;

    // Interior errors
    string interior_cell = format_cell_with_stats(
      stats.max_abs_error_interior,
      stats.avg_abs_error_interior,
      stats.points_checked_interior,
      "Int",
      tol_max_int,
      tol_avg_int );

    // Boundary errors
    string boundary_cell = format_cell_with_stats(
      stats.max_abs_error_boundary,
      stats.avg_abs_error_boundary,
      stats.points_checked_boundary,
      "Bnd",
      tol_max_bnd,
      tol_avg_bnd );

    fmt::print( row_color, "│ {:<12} ", name );
    fmt::print( "│ {} ", interior_cell );
    fmt::print( "│ {} │\n", boundary_cell );

    if ( i < errors.size() - 1 ) { fmt::print( color, "├─{:─^12}─┼─{:─^52}─┼─{:─^52}─┤\n", "", "", "" ); }
  }

  fmt::print( color, "└─{:─^12}─┴─{:─^52}─┴─{:─^52}─┘\n", "", "", "" );
}

// ============================================================================
// DERIVATIVE CHECK FUNCTION FOR 1D SPLINES
// ============================================================================

template <typename SplineType> DerivativeErrors1D check_derivatives_1D(
  SplineType &      spline,
  const string &    spline_name,
  const real_type * knots,
  integer           n_knots,
  integer           points_per_segment = 10 )
{
  DerivativeErrors1D errors;

  // Initialize statistics
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

  errors.D  = init_stats();
  errors.DD = init_stats();

  real_type x_min = spline.x_min();
  real_type x_max = spline.x_max();

  // Tolerance for boundary detection (relative to segment length)
  real_type boundary_tol = 1e-6 * ( x_max - x_min );

  fmt::print(
    fg( fmt::color::blue ),
    "   Checking derivatives for {} ({} knots, {} segments)...",
    spline_name,
    n_knots,
    n_knots - 1 );

  size_t total_points_checked = 0;

  // Test points in each segment
  for ( integer seg = 0; seg < n_knots - 1; ++seg )
  {
    real_type seg_start  = knots[seg];
    real_type seg_end    = knots[seg + 1];
    real_type seg_length = seg_end - seg_start;

    // Generate test points in this segment
    for ( integer p = 0; p <= points_per_segment; ++p )
    {
      real_type x;
      if ( p == 0 )
        x = seg_start + boundary_tol;  // Just after knot
      else if ( p == points_per_segment )
        x = seg_end - boundary_tol;  // Just before next knot
      else
        x = seg_start + p * seg_length / points_per_segment;

      // Check if point is interior or boundary
      bool is_boundary = false;
      for ( integer k = 0; k < n_knots; ++k )
      {
        if ( std::abs( x - knots[k] ) < boundary_tol )
        {
          is_boundary = true;
          break;
        }
      }

      // Ensure we're within spline domain
      if ( x < x_min || x > x_max ) continue;

      // Get analytical derivatives
      real_type D_analytic  = spline.D( x );
      real_type DD_analytic = spline.DD( x );

      // Get finite difference approximations
      real_type D_fd  = FiniteDifferences1D::dx( spline, x );
      real_type DD_fd = FiniteDifferences1D::dxx( spline, x );

      // Calculate absolute errors
      real_type abs_error_D  = std::abs( D_analytic - D_fd );
      real_type abs_error_DD = std::abs( DD_analytic - DD_fd );

      // Update statistics based on boundary status
      if ( is_boundary )
      {
        // First derivative boundary stats
        errors.D.max_abs_error_boundary = std::max( errors.D.max_abs_error_boundary, abs_error_D );
        errors.D.avg_abs_error_boundary += abs_error_D;
        errors.D.points_checked_boundary++;

        // Second derivative boundary stats
        errors.DD.max_abs_error_boundary = std::max( errors.DD.max_abs_error_boundary, abs_error_DD );
        errors.DD.avg_abs_error_boundary += abs_error_DD;
        errors.DD.points_checked_boundary++;
      }
      else
      {
        // First derivative interior stats
        errors.D.max_abs_error_interior = std::max( errors.D.max_abs_error_interior, abs_error_D );
        errors.D.avg_abs_error_interior += abs_error_D;
        errors.D.points_checked_interior++;

        // Second derivative interior stats
        errors.DD.max_abs_error_interior = std::max( errors.DD.max_abs_error_interior, abs_error_DD );
        errors.DD.avg_abs_error_interior += abs_error_DD;
        errors.DD.points_checked_interior++;
      }

      total_points_checked++;
    }
  }

  // Calculate averages
  auto calculate_averages = []( DerivativeErrorStats & stats )
  {
    if ( stats.points_checked_interior > 0 ) stats.avg_abs_error_interior /= stats.points_checked_interior;
    if ( stats.points_checked_boundary > 0 ) stats.avg_abs_error_boundary /= stats.points_checked_boundary;
  };

  calculate_averages( errors.D );
  calculate_averages( errors.DD );

  fmt::print(
    fg( fmt::color::green ),
    " ✓ ({} points total, {} interior, {} boundary)\n",
    total_points_checked,
    errors.D.points_checked_interior + errors.DD.points_checked_interior,
    errors.D.points_checked_boundary + errors.DD.points_checked_boundary );

  return errors;
}

// ============================================================================
// MAIN FUNCTION
// ============================================================================

// Function to print a table of results
void print_results_table( const vector<SplineResult> & results, integer dataset )
{
  // Table header
  fmt::print(
    fg( fmt::color::yellow ) | fmt::emphasis::bold,
    "\n"
    "┌──────────────┬──────────┬──────────┬──────────┬──────────┬──────────┬──────────┐\n"
    "│ {:^12} │ {:^8} │ {:^8} │ {:^8} │ {:^8} │ {:^8} │ {:^8} │\n"
    "├──────────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┤\n",
    "Spline Type",
    "x_min",
    "x_max",
    "y_min",
    "y_max",
    "x@min",
    "x@max" );

  // Table rows with alternating colors
  for ( size_t i = 0; i < results.size(); ++i )
  {
    const auto & r         = results[i];
    auto         row_color = ( i % 2 == 0 ) ? fg( fmt::color::light_green ) : fg( fmt::color::light_blue );

    fmt::print( row_color, "│ {:<12} ", r.name );
    fmt::print( row_color, "│ {:>8.3f} ", r.x_min );
    fmt::print( row_color, "│ {:>8.3f} ", r.x_max );

    // Highlight y_min and y_max with different colors based on value
    auto y_min_color = ( r.y_min < 0 ) ? fg( fmt::color::red ) : fg( fmt::color::green );
    auto y_max_color = ( r.y_max < 0 ) ? fg( fmt::color::red ) : fg( fmt::color::green );

    fmt::print( y_min_color, "│ {:>8.3f} ", r.y_min );
    fmt::print( y_max_color, "│ {:>8.3f} ", r.y_max );
    fmt::print( row_color, "│ {:>8.3f} ", r.x_at_min );
    fmt::print( row_color, "│ {:>8.3f} │\n", r.x_at_max );
  }

  fmt::print(
    fg( fmt::color::yellow ) | fmt::emphasis::bold,
    "└──────────────┴──────────┴──────────┴──────────┴──────────┴──────────┴──────────┘\n" );

  fmt::print(
    fg( fmt::color::gray ) | fmt::emphasis::italic,
    "Dataset {} - {} points analyzed\n",
    dataset,
    nn[dataset] );
}

int main()
{
  print_header( "SPLINE INTERPOLATION TEST SUITE" );

  // Create spline objects
  LinearSpline   li;
  ConstantSpline co;
  AkimaSpline    ak;
  CubicSpline    cs;
  BesselSpline   be;
  PchipSpline    pc;
  QuinticSpline  qs;

  // Open output files
  ofstream file_li, file_co, file_ak, file_cs, file_be, file_pc, file_qs;

  for ( integer k = 0; k < 6; ++k )
  {
    fmt::print( fg( fmt::color::magenta ) | fmt::emphasis::bold, "\n\n📊 DATASET {} ANALYSIS\n", k );

    real_type * xx{ nullptr };
    real_type * yy{ nullptr };

    switch ( k )
    {
      case 0:
        xx = xx0;
        yy = yy0;
        break;
      case 1:
        xx = xx1;
        yy = yy1;
        break;
      case 2:
        xx = xx2;
        yy = yy2;
        break;
      case 3:
        xx = xx3;
        yy = yy3;
        break;
      case 4:
        xx = xx4;
        yy = yy4;
        break;
      case 5:
        xx = xx5;
        yy = yy5;
        break;
    }

    // Open output files
    string fname;
    fname = fmt::format( "out/Linear{}.txt", k );
    file_li.open( fname.data() );
    fname = fmt::format( "out/Constant{}.txt", k );
    file_co.open( fname.data() );
    fname = fmt::format( "out/Akima{}.txt", k );
    file_ak.open( fname.data() );
    fname = fmt::format( "out/Cubic{}.txt", k );
    file_cs.open( fname.data() );
    fname = fmt::format( "out/Bessel{}.txt", k );
    file_be.open( fname.data() );
    fname = fmt::format( "out/Pchip{}.txt", k );
    file_pc.open( fname.data() );
    fname = fmt::format( "out/Quintic{}.txt", k );
    file_qs.open( fname.data() );

    real_type xmin{ xx[0] };
    real_type xmax{ xx[nn[k] - 1] };

    vector<SplineResult>                     results;
    vector<pair<string, DerivativeErrors1D>> derivative_errors;

    // Lambda function to process each spline type
    auto process_spline = [&]( auto & spline, const string & name, ofstream & file ) -> DerivativeErrors1D
    {
      fmt::print( fg( fmt::color::cyan ), "\n🔧 Processing {}...\n", name );

      // Build spline
      spline.clear();
      spline.reserve( nn[k] );
      for ( integer i = 0; i < integer( nn[k] ); ++i ) spline.push_back( xx[i], yy[i] );
      spline.build();

      // Get min/max information
      integer   i_min_pos, i_max_pos;
      real_type x_min_pos, x_max_pos, y_min, y_max;
      spline.y_min_max( i_min_pos, x_min_pos, y_min, i_max_pos, x_max_pos, y_max );

      // Save results
      results.push_back( { name, spline.x_min(), spline.x_max(), y_min, y_max, x_min_pos, x_max_pos } );

      // Write spline data to file
      file << "x\ty\tDy\tDDy\n";
      for ( real_type x = xmin; x <= xmax; x += ( xmax - xmin ) / 1000 )
        fmt::print( file, "{}\t{}\t{}\t{}\n", x, spline( x ), spline.D( x ), spline.DD( x ) );
      file.close();

      // Check derivatives
      auto errors = check_derivatives_1D( spline, name, xx, nn[k], 10 );
      derivative_errors.emplace_back( name, errors );

      fmt::print( fg( fmt::color::green ), "   ✓ {} completed\n", name );
      return errors;
    };

    // Process all spline types
    fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n📈 Building splines and checking derivatives...\n" );

    process_spline( li, "Linear", file_li );
    process_spline( co, "Constant", file_co );
    process_spline( ak, "Akima", file_ak );
    process_spline( cs, "Cubic", file_cs );
    process_spline( be, "Bessel", file_be );
    process_spline( pc, "Pchip", file_pc );
    process_spline( qs, "Quintic", file_qs );

    // Print min/max results table
    print_results_table( results, k );

    // Print derivative error tables
    fmt::print( fg( fmt::color::magenta ) | fmt::emphasis::bold, "\n📊 DERIVATIVE ERROR ANALYSIS - DATASET {}:\n", k );

    print_derivative_table_1D( "D", derivative_errors, k );
    print_derivative_table_1D( "DD", derivative_errors, k );

    // Summary of derivative checks
    fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n📋 DERIVATIVE CHECK SUMMARY - DATASET {}:\n", k );

    for ( const auto & [name, err] : derivative_errors )
    {
      fmt::print( fg( fmt::color::white ) | fmt::emphasis::bold, "   {}:\n", name );

      // First derivative summary
      bool first_deriv_interior_ok = ( err.D.max_abs_error_interior < 1e-6 && err.D.avg_abs_error_interior < 1e-7 );

      bool first_deriv_boundary_ok = ( err.D.max_abs_error_boundary < 1e-5 && err.D.avg_abs_error_boundary < 1e-6 );

      // Second derivative summary
      bool second_deriv_interior_ok = ( err.DD.max_abs_error_interior < 1e-3 && err.DD.avg_abs_error_interior < 1e-4 );

      bool second_deriv_boundary_ok = ( err.DD.max_abs_error_boundary < 1e-2 && err.DD.avg_abs_error_boundary < 1e-3 );

      // Print results with color coding
      auto print_check = []( bool ok, const string & text )
      {
        if ( ok )
          fmt::print( fg( fmt::color::green ), "     ✓ {}\n", text );
        else
          fmt::print( fg( fmt::color::yellow ), "     ⚠ {}\n", text );
      };

      print_check( first_deriv_interior_ok, "First derivative (interior)" );
      print_check( first_deriv_boundary_ok, "First derivative (boundary)" );
      print_check( second_deriv_interior_ok, "Second derivative (interior)" );
      print_check( second_deriv_boundary_ok, "Second derivative (boundary)" );
    }

    fmt::print( fg( fmt::color::gray ), "   Output files saved to: out/*{}.txt\n", k );
  }

  print_header( "TEST COMPLETED SUCCESSFULLY" );
  fmt::print( fg( fmt::color::light_green ) | fmt::emphasis::bold, "\n🎉 All spline analyses completed! 🎉\n" );
  fmt::print( fg( fmt::color::gray ), "Results saved in the 'out/' directory.\n\n" );

  return 0;
}
