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

#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

#include "Splines.hh"
#include "Utils_fmt.hh"
#include "Utils_string.hh"

#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cassert>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#endif

using namespace SplinesLoad;
using namespace std;
using Splines::integer;
using Splines::real_type;

// Test problem for Akima interpolation
// Ref. : Hiroshi Akima, Journal of the ACM, Vol. 17, No. 4, October 1970, pages 589-602.

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

// monotone
static real_type xx6[] = { 0, 1, 2, 3, 4 };
static real_type yy6[] = { 0, 1, 1.1, 2.0, 2.1 };

static integer nn[] = { 11, 11, 11, 9, 12, 4, 5 };

// ===========================================================================
// IMPROVED SUPPORT FUNCTIONS
// ===========================================================================

// Function for central finite differences
real_type finite_diff_central( SplineSet const & ss, real_type x, integer i, real_type h = 1e-6 )
{
  return ( ss( x + h, i ) - ss( x - h, i ) ) / ( 2 * h );
}

// Function for forward finite differences
real_type finite_diff_forward( SplineSet const & ss, real_type x, integer i, real_type h = 1e-6 )
{
  return ( ss( x + h, i ) - ss( x, i ) ) / h;
}

// Function for backward finite differences
real_type finite_diff_backward( SplineSet const & ss, real_type x, integer i, real_type h = 1e-6 )
{
  return ( ss( x, i ) - ss( x - h, i ) ) / h;
}

// Adaptive finite difference function (chooses optimal method)
real_type finite_diff_adaptive(
  SplineSet const & ss,
  real_type         x,
  integer           i,
  real_type         xmin,
  real_type         xmax,
  real_type         h = 1e-6 )
{
  if ( x <= xmin + h ) { return finite_diff_forward( ss, x, i, h ); }
  else if ( x >= xmax - h ) { return finite_diff_backward( ss, x, i, h ); }
  else
  {
    return finite_diff_central( ss, x, i, h );
  }
}

// Check if a point is near a knot (within tolerance)
bool is_near_knot( real_type x, real_type const * knots, integer n, real_type eps = 1e-10 )
{
  for ( integer i = 0; i < n; ++i )
  {
    if ( std::abs( x - knots[i] ) <= eps ) return true;
  }
  return false;
}

// Determine if a spline is differentiable at knots
bool is_spline_differentiable_at_knots( std::string_view spline_name )
{
  // Splines that are NOT differentiable at knots:
  // - CONSTANT: derivative is discontinuous (0 between knots, undefined at knots)
  // - LINEAR: derivative is discontinuous at knots
  return !( spline_name == "SPLINE_CONSTANT" || spline_name == "SPLINE_LINEAR" );
}

// Compute error robustly (avoid division by zero)
std::pair<real_type, real_type> compute_error( real_type exact, real_type approx, real_type abs_tol = 1e-12 )
{
  real_type abs_err = std::abs( exact - approx );
  real_type rel_err = 0.0;

  // If both values are very small, consider error zero
  if ( std::abs( exact ) < abs_tol && std::abs( approx ) < abs_tol ) { rel_err = 0.0; }
  // Otherwise compute relative error using maximum as denominator
  else
  {
    real_type denom = std::max( std::abs( exact ), std::abs( approx ) );
    if ( denom > abs_tol ) { rel_err = 100.0 * abs_err / denom; }
    else
    {
      rel_err = ( abs_err > abs_tol ) ? 100.0 : 0.0;
    }
  }

  return std::make_pair( abs_err, std::min( rel_err, 100.0 ) );  // Limit to 100%
}

// ===========================================================================
// TABLE FUNCTIONS
// ===========================================================================

// Print table header
void print_table_header( vector<string> const & headers )
{
  // Top line
  fmt::print( "┌{}", fmt::format( "{:─^{}}", "", 12 ) );
  for ( size_t i = 1; i < headers.size(); ++i ) { fmt::print( "┬{}", fmt::format( "{:─^{}}", "", 20 ) ); }
  fmt::print( "┐\n" );

  // Headers
  fmt::print( "│{:^12}", headers[0] );
  for ( size_t i = 1; i < headers.size(); ++i ) { fmt::print( "│{:^20}", headers[i] ); }
  fmt::print( "│\n" );

  // Separator
  fmt::print( "├{}", fmt::format( "{:─^{}}", "", 12 ) );
  for ( size_t i = 1; i < headers.size(); ++i ) { fmt::print( "┼{}", fmt::format( "{:─^{}}", "", 20 ) ); }
  fmt::print( "┤\n" );
}

// Print table row with color
template <typename COLOR> void print_table_row( vector<string> const & values, COLOR color )
{
  fmt::print( "│{:^12}", values[0] );
  for ( size_t i = 1; i < values.size() - 1; ++i ) { fmt::print( "│{:^20}", values[i] ); }
  fmt::print( "│" );
  fmt::print( fg( color ), "{:^20}", values.back() );
  fmt::print( "│\n" );
}

// Print table footer
void print_table_footer( size_t ncols )
{
  fmt::print( "└{}", fmt::format( "{:─^{}}", "", 12 ) );
  for ( size_t i = 1; i < ncols; ++i ) { fmt::print( "┴{}", fmt::format( "{:─^{}}", "", 20 ) ); }
  fmt::print( "┘\n" );
}

// ===========================================================================
// QUINTIC SPLINE TEST FUNCTIONS
// ===========================================================================

// Test quintic splines with different subtypes
void test_quintic_subtypes(
  SplineSet &                 ss,
  real_type                   x,
  const vector<string> &      quintic_subtype_names,
  vector<vector<real_type>> & quintic_results )
{
  // Get the base quintic spline index
  integer base_quintic_idx = -1;
  for ( integer i = 0; i < ss.num_splines(); ++i )
  {
    string name = string( ss.header( i ) );
    if ( name == "SPLINE_QUINTIC" )
    {
      base_quintic_idx = i;
      break;
    }
  }

  if ( base_quintic_idx < 0 ) return;

  // Test all quintic subtypes at this point
  for ( size_t sub_idx = 0; sub_idx < quintic_subtype_names.size(); ++sub_idx )
  {
    // In this test, we're simulating different quintic subtypes
    // In a real implementation, you would have separate splines for each subtype
    real_type value      = ss( x, base_quintic_idx );  // Using base quintic for all
    real_type derivative = ss.D( x, base_quintic_idx );

    quintic_results[sub_idx].push_back( value );
    quintic_results[sub_idx].push_back( derivative );
  }
}

// ===========================================================================
// MAIN FUNCTION
// ===========================================================================

int main()
{
  fmt::print(
    fg( fmt::color::cyan ) | fmt::emphasis::bold,
    "\n"
    "╔═══════════════════════════════════════════════════════════════════════════════════╗\n"
    "║                            TEST No.3 - SPLINE DERIVATIVES                         ║\n"
    "║                                                                                   ║\n"
    "║  This test verifies the correctness of spline derivatives by comparing            ║\n"
    "║  them with finite differences.                                                    ║\n"
    "╚═══════════════════════════════════════════════════════════════════════════════════╝\n\n" );

  // Print interpretation guidelines
  fmt::print(
    fg( fmt::color::yellow ) | fmt::emphasis::bold,
    "⚠️  INTERPRETATION GUIDELINES:\n"
    "   - CONSTANT and LINEAR splines are not differentiable at knots\n"
    "   - For these splines, knot tests are skipped or flagged as warnings\n"
    "   - Errors < 1e-4% are considered excellent for C1 splines (CUBIC, AKIMA, etc.)\n"
    "   - Errors < 1% are acceptable for less regular splines\n"
    "   - Colors: 🟢 Green = excellent, 🟡 Yellow = acceptable, 🔴 Red = problematic\n\n" );

  SplineSet ss;
  ofstream  file, file_D;

  // Configuration parameters
  const real_type h_fd           = 1e-6;   // Finite difference step
  const real_type rel_tol        = 1e-4;   // Relative tolerance for errors (0.01%)
  const real_type abs_tol        = 1e-8;   // Absolute tolerance for errors
  const real_type knot_tolerance = 1e-10;  // Tolerance for identifying knots

  // Quintic subtype names
  vector<string> quintic_subtype_names = { "QUINTIC_CUBIC", "QUINTIC_AKIMA", "QUINTIC_VANLEER", "QUINTIC_PCHIP" };

  for ( integer k = 0; k < 7; ++k )
  {
    real_type * xx{ nullptr };
    real_type * yy{ nullptr };
    string      dataset_name;

    switch ( k )
    {
      case 0:
        xx           = xx0;
        yy           = yy0;
        dataset_name = "Akima Test 0";
        break;
      case 1:
        xx           = xx1;
        yy           = yy1;
        dataset_name = "Akima Test 1";
        break;
      case 2:
        xx           = xx2;
        yy           = yy2;
        dataset_name = "Akima Test 2";
        break;
      case 3:
        xx           = xx3;
        yy           = yy3;
        dataset_name = "RPN 14";
        break;
      case 4:
        xx           = xx4;
        yy           = yy4;
        dataset_name = "Titanium";
        break;
      case 5:
        xx           = xx5;
        yy           = yy5;
        dataset_name = "Toolpath";
        break;
      case 6:
        xx           = xx6;
        yy           = yy6;
        dataset_name = "Monotone";
        break;
    }

    fmt::print(
      fg( fmt::color::yellow ) | fmt::emphasis::bold,
      "\n"
      "┌─────────────────────────────────────────────────────────────┐\n"
      "│ Dataset {:2}: {:47} │\n"
      "└─────────────────────────────────────────────────────────────┘\n",
      k,
      dataset_name );

    // Open files for output
    string fname = fmt::format( "out/SplineSet{}.txt", k );
    file.open( fname.data() );

    fname = fmt::format( "out/SplineSet{}_D.txt", k );
    file_D.open( fname.data() );

    real_type xmin = xx[0];
    real_type xmax = xx[nn[k] - 1];
    integer   npts = nn[k];

    // Number of splines: original 7 + 4 quintic subtypes
    integer nspl = 10;  // 7 original + 4 quintic subtypes

    // Define all splines including quintic subtypes
    vector<char const *> all_headers = { "SPLINE_CONSTANT", "SPLINE_LINEAR",

                                         "SPLINE_CUBIC",    "SPLINE_AKIMA",  "SPLINE_VANLEER",  "SPLINE_PCHIP",

                                         "QUINTIC_CUBIC",   "QUINTIC_AKIMA", "QUINTIC_VANLEER", "QUINTIC_PCHIP" };

    // For simplicity in this test, we'll use the same type for all quintic subtypes
    // In a real implementation, you would need different types or configurations
    vector<SplineType1D> all_stypes = {
      SplineType1D::CONSTANT,       SplineType1D::LINEAR,

      SplineType1D::CUBIC,          SplineType1D::AKIMA,  SplineType1D::VANLEER, SplineType1D::PCHIP,

      SplineType1D::QUINTIC_CUBIC,   // QUINTIC_CUBIC subtype
      SplineType1D::QUINTIC_AKIMA,   // QUINTIC_AKIMA subtype
      SplineType1D::QUINTIC_VANLEER, // QUINTIC_VANLEER subtype
      SplineType1D::QUINTIC_PCHIP    // QUINTIC_PCHIP subtype
    };

    // Prepare data for all splines (same Y values for all)
    vector<const real_type *> all_Y;
    for ( integer i = 0; i < nspl; ++i ) { all_Y.push_back( yy ); }

    // Build the spline set
    ss.build( nspl, npts, all_headers.data(), all_stypes.data(), xx, all_Y.data() );

    // Write files for compatibility (from original code)
    file << "x";
    file_D << "x";
    for ( integer i = 0; i < nspl; ++i )
    {
      file << '\t' << ss.header( i );
      file_D << '\t' << ss.header( i );
    }
    file << '\n';
    file_D << '\n';

    for ( real_type x = xmin; x <= xmax; x += ( xmax - xmin ) / 1000 )
    {
      file << x;
      file_D << x;
      for ( integer i = 0; i < nspl; ++i )
      {
        file << '\t' << ss( x, i );
        file_D << '\t' << ss.D( x, i );
      }
      file << '\n';
      file_D << '\n';
    }
    file.close();
    file_D.close();

    // ===========================================================================
    // DERIVATIVE TEST WITH FINITE DIFFERENCES
    // ===========================================================================

    // Prepare test points
    vector<real_type> test_points;

    // 1. Add all knots
    for ( integer i = 0; i < npts; ++i ) { test_points.push_back( xx[i] ); }

    // 2. Add internal points between knots
    for ( integer i = 0; i < npts - 1; ++i )
    {
      real_type x1 = xx[i];
      real_type x2 = xx[i + 1];
      real_type dx = x2 - x1;

      // Add 3 internal points per interval
      test_points.push_back( x1 + 0.25 * dx );
      test_points.push_back( x1 + 0.50 * dx );
      test_points.push_back( x1 + 0.75 * dx );
    }

    // Sort and remove duplicates
    sort( test_points.begin(), test_points.end() );
    test_points.erase( unique( test_points.begin(), test_points.end() ), test_points.end() );

    // Header for derivative check
    fmt::print( fg( fmt::color::green ), "\n✅ Derivative check - Dataset {}\n", k );

    vector<string> table_headers = { "x", "Spline", "D(x)", "FinDiff", "Err" };
    print_table_header( table_headers );

    // Global statistics
    integer total_tests   = 0;
    integer passed_tests  = 0;
    integer warning_tests = 0;
    integer failed_tests  = 0;
    integer skipped_tests = 0;

    // Results for quintic subtypes comparison
    vector<vector<real_type>> quintic_subtype_results( quintic_subtype_names.size() );

    for ( integer spline_idx = 0; spline_idx < nspl; ++spline_idx )
    {
      string spline_name    = string( ss.header( spline_idx ) );
      bool   differentiable = is_spline_differentiable_at_knots( spline_name );

      for ( size_t pt_idx = 0; pt_idx < test_points.size(); ++pt_idx )
      {
        real_type x       = test_points[pt_idx];
        bool      is_knot = is_near_knot( x, xx, npts, knot_tolerance );

        total_tests++;

        // ===================================================================
        // CASE 1: NON-DIFFERENTIABLE SPLINE AT KNOT
        // ===================================================================
        if ( !differentiable && is_knot )
        {
          skipped_tests++;

          // Only for debugging: print a warning
          if ( pt_idx == 0 )  // Only once per spline
          {
            fmt::print(
              fg( fmt::color::yellow ),
              "│{:^12}│{:^20}│{:^20}│{:^20}│{:^20}│\n",
              fmt::format( "{:.6f}", x ),
              spline_name,
              "N/A (discontinuous)",
              "N/A",
              "SKIPPED" );
          }
          continue;
        }

        // ===================================================================
        // CASE 2: DERIVATIVE CALCULATION
        // ===================================================================

        // Calculate derivative from spline
        real_type D_spline = ss.D( x, spline_idx );

        // Calculate derivative with finite differences
        real_type D_fd;
        try
        {
          D_fd = finite_diff_adaptive( ss, x, spline_idx, xmin, xmax, h_fd );
        }
        catch ( ... )
        {
          // Fallback: use appropriate method for position
          if ( x <= xmin + h_fd ) { D_fd = finite_diff_forward( ss, x, spline_idx, h_fd ); }
          else if ( x >= xmax - h_fd ) { D_fd = finite_diff_backward( ss, x, spline_idx, h_fd ); }
          else
          {
            D_fd = finite_diff_central( ss, x, spline_idx, h_fd );
          }
        }

        // Calculate error
        auto [abs_err, rel_err] = compute_error( D_spline, D_fd, abs_tol );

        // Determine if test passes
        bool test_passed = ( abs_err < abs_tol ) || ( rel_err < rel_tol );

        if ( test_passed ) { passed_tests++; }
        else if ( rel_err < 10.0 * rel_tol ) { warning_tests++; }
        else
        {
          failed_tests++;
        }

        // ===================================================================
        // CASE 3: PRINT RESULTS
        // ===================================================================
        // Print only for:
        // 1. All knots (for differentiable splines)
        // 2. Points with significant error
        // 3. First and last point of each spline
        if ( is_knot || !test_passed || pt_idx == 0 || pt_idx == test_points.size() - 1 )
        {
          // Format values
          vector<string> row_values = { fmt::format( "{:.6f}", x ),
                                        spline_name,
                                        fmt::format( "{:.12e}", D_spline ),
                                        fmt::format( "{:.12e}", D_fd ),
                                        fmt::format( "{:.2e}", abs_err ) };

          // Choose color based on error
          fmt::color color_code;
          if ( test_passed ) { color_code = fmt::color::green; }
          else if ( abs_err < 10.0 * abs_tol ) { color_code = fmt::color::yellow; }
          else
          {
            color_code = fmt::color::red;
          }

          print_table_row( row_values, color_code );
        }

        // Store results for quintic subtypes comparison
        if ( spline_name.find( "QUINTIC" ) != string::npos )
        {
          for ( size_t sub_idx = 0; sub_idx < quintic_subtype_names.size(); ++sub_idx )
          {
            if ( spline_name == quintic_subtype_names[sub_idx] )
            {
              quintic_subtype_results[sub_idx].push_back( D_spline );
            }
          }
        }
      }

      // Separator between different splines
      if ( spline_idx < nspl - 1 )
      {
        fmt::print(
          "├{}┼{}┼{}┼{}┼{}┤\n",
          fmt::format( "{:─^{}}", "", 12 ),
          fmt::format( "{:─^{}}", "", 20 ),
          fmt::format( "{:─^{}}", "", 20 ),
          fmt::format( "{:─^{}}", "", 20 ),
          fmt::format( "{:─^{}}", "", 20 ) );
      }
    }

    print_table_footer( table_headers.size() );

    // ===========================================================================
    // DATASET STATISTICS
    // ===========================================================================
    integer tested_points = total_tests - skipped_tests;

    fmt::print(
      fg( fmt::color::blue ),
      "\n📊 DATASET {} STATISTICS:\n"
      "   - Total test points: {}\n"
      "   - Actually tested points: {} ({} knots + {} internal points)\n"
      "   - Skipped tests: {} (knots for non-differentiable splines)\n"
      "   - Passed tests: {} ({:.1f}%)\n"
      "   - Warning tests: {} ({:.1f}%)\n"
      "   - Failed tests: {} ({:.1f}%)\n",
      k,
      total_tests,
      tested_points,
      npts,
      tested_points - npts,
      skipped_tests,
      passed_tests,
      ( tested_points > 0 ) ? 100.0 * passed_tests / tested_points : 0.0,
      warning_tests,
      ( tested_points > 0 ) ? 100.0 * warning_tests / tested_points : 0.0,
      failed_tests,
      ( tested_points > 0 ) ? 100.0 * failed_tests / tested_points : 0.0 );

    // ===========================================================================
    // QUINTIC SUBTYPE COMPARISON
    // ===========================================================================
    if ( !quintic_subtype_results.empty() && quintic_subtype_results[0].size() > 0 )
    {
      fmt::print(
        fg( fmt::color::magenta ),
        "\n🔬 QUINTIC SUBTYPE COMPARISON:\n"
        "   Comparing derivative values across quintic subtypes at {} test points\n",
        quintic_subtype_results[0].size() );

      // Compare derivatives between quintic subtypes
      fmt::print( "   Max differences between subtypes:\n" );

      for ( size_t i = 0; i < quintic_subtype_names.size(); ++i )
      {
        for ( size_t j = i + 1; j < quintic_subtype_names.size(); ++j )
        {
          if ( quintic_subtype_results[i].size() == quintic_subtype_results[j].size() )
          {
            real_type max_diff = 0.0;
            for ( size_t idx = 0; idx < quintic_subtype_results[i].size(); ++idx )
            {
              real_type diff = std::abs( quintic_subtype_results[i][idx] - quintic_subtype_results[j][idx] );
              max_diff       = std::max( max_diff, diff );
            }

            if ( max_diff > 1e-10 )
            {
              fmt::print( "     {} vs {}: {:.6e}\n", quintic_subtype_names[i], quintic_subtype_names[j], max_diff );
            }
          }
        }
      }
    }

    // Print specific recommendations for the dataset
    fmt::print( fg( fmt::color::cyan ), "\n💡 RECOMMENDATIONS FOR DATASET '{}':\n", dataset_name );

    if ( failed_tests > 0 )
    {
      fmt::print( fg( fmt::color::red ), "   ⚠️  Check implementations of splines with high errors\n" );
    }

    if ( warning_tests > 0 )
    {
      fmt::print( fg( fmt::color::yellow ), "   ⚠️  Some splines show moderate errors - check data regularity\n" );
    }

    if ( passed_tests == tested_points && tested_points > 0 )
    {
      fmt::print( fg( fmt::color::green ), "   ✅ All derivatives are correct within specified tolerances\n" );
    }

    fmt::print( "\n" );
  }

  // ===========================================================================
  // FINAL SUMMARY
  // ===========================================================================
  fmt::print(
    fg( fmt::color::cyan ) | fmt::emphasis::bold,
    "\n"
    "╔═════════════════════════════════════════════════════════════════════════════════════╗\n"
    "║                               FINAL TEST SUMMARY                                    ║\n"
    "║                                                                                     ║\n"
    "║  FINAL CONSIDERATIONS:                                                              ║\n"
    "║  1. CONSTANT and LINEAR splines are not tested at knots (discontinuous derivatives) ║\n"
    "║  2. For C1 splines (CUBIC, AKIMA, etc.) errors < 1e-4% are expected                 ║\n"
    "║  3. High errors may indicate:                                                       ║\n"
    "║     - Incorrect derivative implementation                                           ║\n"
    "║     - Low regularity data                                                           ║\n"
    "║     - Suboptimal finite difference step                                             ║\n"
    "║  4. For debugging, use:                                                             ║\n"
    "║     - Smaller/larger step h for finite differences                                  ║\n"
    "║     - Analytical derivative checks for simple cases                                 ║\n"
    "║     - C1 continuity verification at knots for splines that should have it           ║\n"
    "║  5. Quintic spline subtypes (CUBIC, AKIMA, VANLEER, PCHIP) are tested separately    ║\n"
    "╚═════════════════════════════════════════════════════════════════════════════════════╝\n\n" );

  return 0;
}
