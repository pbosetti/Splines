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

static integer nn[] = { 11, 11, 11, 9, 12, 4 };

// Structure to hold spline results
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

  LinearSpline   li;
  ConstantSpline co;
  AkimaSpline    ak;
  CubicSpline    cs;
  BesselSpline   be;
  PchipSpline    pc;
  QuinticSpline  qs;

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

    vector<SplineResult> results;

    // Lambda function to process each spline type
    auto process_spline = [&]( auto & spline, const string & name, ofstream & file )
    {
      fmt::print( fg( fmt::color::cyan ), "\n🔧 Processing {}...\n", name );

      spline.clear();
      spline.reserve( nn[k] );
      for ( integer i = 0; i < integer( nn[k] ); ++i ) spline.push_back( xx[i], yy[i] );

      spline.build();

      integer   i_min_pos, i_max_pos;
      real_type x_min_pos, x_max_pos, y_min, y_max;
      spline.y_min_max( i_min_pos, x_min_pos, y_min, i_max_pos, x_max_pos, y_max );

      // Save results
      results.push_back( { name, spline.x_min(), spline.x_max(), y_min, y_max, x_min_pos, x_max_pos } );

      // Write to file
      file << "x\ty\tDy\tDDy\n";
      for ( real_type x = xmin; x <= xmax; x += ( xmax - xmin ) / 1000 )
        fmt::print( file, "{}\t{}\t{}\t{}\n", x, spline( x ), spline.D( x ), spline.DD( x ) );

      file.close();

      fmt::print( fg( fmt::color::green ), "   ✓ {} completed\n", name );
    };

    // Process all spline types
    process_spline( li, "Linear", file_li );
    process_spline( co, "Constant", file_co );
    process_spline( ak, "Akima", file_ak );
    process_spline( cs, "Cubic", file_cs );
    process_spline( be, "Bessel", file_be );
    process_spline( pc, "Pchip", file_pc );
    process_spline( qs, "Quintic", file_qs );

    // Print results table
    print_results_table( results, k );

    fmt::print( fg( fmt::color::gray ), "   Output files saved to: out/*{}.txt\n", k );
  }

  print_header( "TEST COMPLETED SUCCESSFULLY" );
  fmt::print( fg( fmt::color::light_green ) | fmt::emphasis::bold, "\n🎉 All spline analyses completed! 🎉\n" );
  fmt::print( fg( fmt::color::gray ), "Results saved in the 'out/' directory.\n\n" );

  return 0;
}
