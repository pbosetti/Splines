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
#include "Utils_string.hh"

#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>

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

// Funzione per differenze finite centrali
real_type finite_diff_central( SplineSet const & ss, real_type x, integer i, real_type h = 1e-6 )
{
  return ( ss( x + h, i ) - ss( x - h, i ) ) / ( 2 * h );
}

// Funzione per differenze finite avanti/indietro
real_type finite_diff_forward( SplineSet const & ss, real_type x, integer i, real_type h = 1e-6 )
{
  return ( ss( x + h, i ) - ss( x, i ) ) / h;
}

real_type finite_diff_backward( SplineSet const & ss, real_type x, integer i, real_type h = 1e-6 )
{
  return ( ss( x, i ) - ss( x - h, i ) ) / h;
}

// Stampa tabella con bordi Unicode e colori
void print_table_header( vector<string> const & headers )
{
  // Linea superiore
  fmt::print( "┌{}", fmt::format( "{:─^{}}", "", 12 ) );
  for ( size_t i = 1; i < headers.size(); ++i ) { fmt::print( "┬{}", fmt::format( "{:─^{}}", "", 20 ) ); }
  fmt::print( "┐\n" );

  // Intestazioni
  fmt::print( "│{:^12}", headers[0] );
  for ( size_t i = 1; i < headers.size(); ++i ) { fmt::print( "│{:^20}", headers[i] ); }
  fmt::print( "│\n" );

  // Separatore
  fmt::print( "├{}", fmt::format( "{:─^{}}", "", 12 ) );
  for ( size_t i = 1; i < headers.size(); ++i ) { fmt::print( "┼{}", fmt::format( "{:─^{}}", "", 20 ) ); }
  fmt::print( "┤\n" );
}

template <typename COLOR> void print_table_row( vector<string> const & values, COLOR color )
{
  fmt::print( "│{:^12}", values[0] );
  for ( size_t i = 1; i < values.size() - 1; ++i ) { fmt::print( "│{:^20}", values[i] ); }
  fmt::print( "│" );
  fmt::print( fg( color ), "{:^20}", values.back() );
  fmt::print( "│\n" );
}

void print_table_footer( size_t ncols )
{
  fmt::print( "└{}", fmt::format( "{:─^{}}", "", 12 ) );
  for ( size_t i = 1; i < ncols; ++i ) { fmt::print( "┴{}", fmt::format( "{:─^{}}", "", 20 ) ); }
  fmt::print( "┘\n" );
}

int main()
{
  fmt::print(
    fg( fmt::color::cyan ) | fmt::emphasis::bold,
    "\n"
    "╔══════════════════════════════════════════════════════════╗\n"
    "║                         TEST N.3                         ║\n"
    "╚══════════════════════════════════════════════════════════╝\n\n" );

  SplineSet ss;
  ofstream  file, file_D;

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
      "┌────────────────────────────────────────────────────┐\n"
      "│ Dataset {:2}: {:38} │\n"
      "└────────────────────────────────────────────────────┘\n",
      k,
      dataset_name );

    string fname;
    fname = fmt::format( "out/SplineSet{}.txt", k );
    file.open( fname.data() );

    fname = fmt::format( "out/SplineSet{}_D.txt", k );
    file_D.open( fname.data() );

    real_type xmin{ xx[0] };
    real_type xmax{ xx[nn[k] - 1] };
    integer   nspl{ 7 };
    integer   npts{ nn[k] };

    char const * headers[] = { "SPLINE_CONSTANT", "SPLINE_LINEAR", "SPLINE_CUBIC",  "SPLINE_AKIMA",
                               "SPLINE_BESSEL",   "SPLINE_PCHIP",  "SPLINE_QUINTIC" };

    constexpr SplineType1D stype[] = { SplineType1D::CONSTANT, SplineType1D::LINEAR, SplineType1D::CUBIC,
                                       SplineType1D::AKIMA,    SplineType1D::BESSEL, SplineType1D::PCHIP,
                                       SplineType1D::QUINTIC };

    real_type const * Y[]{ yy, yy, yy, yy, yy, yy, yy, yy };

    ss.build( nspl, npts, headers, stype, xx, Y );
    // ss.info( cout );

    // Scrittura file (mantenuta per compatibilità)
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

    // ========== NUOVO: CONTROLLO DERIVATE CON DIFFERENZE FINITE ==========

    // Punti di test: punti di transizione (nodi) + punti interni
    vector<real_type> test_points;

    // Aggiungi tutti i nodi
    for ( integer i = 0; i < npts; ++i ) { test_points.push_back( xx[i] ); }

    // Aggiungi punti intermedi tra i nodi (evitando i nodi stessi)
    for ( integer i = 0; i < npts - 1; ++i )
    {
      real_type x1 = xx[i];
      real_type x2 = xx[i + 1];
      // Aggiungi 2 punti interni per ogni intervallo
      test_points.push_back( x1 + 0.25 * ( x2 - x1 ) );
      test_points.push_back( x1 + 0.75 * ( x2 - x1 ) );
    }

    // Ordina e rimuovi duplicati (nel caso i punti interni coincidano con i nodi)
    sort( test_points.begin(), test_points.end() );
    test_points.erase( unique( test_points.begin(), test_points.end() ), test_points.end() );

    // Tabella per il controllo delle derivate
    fmt::print( fg( fmt::color::green ), "\nControllo derivate - Dataset {}\n", k );

    vector<string> table_headers = { "x", "Spline", "D(x)", "FinDiff", "Err. Rel%" };
    print_table_header( table_headers );

    real_type const h   = 1e-6;
    real_type const tol = 1e-4;

    for ( integer spline_idx = 0; spline_idx < nspl; ++spline_idx )
    {
      string spline_name = string(
        ss.header( spline_idx ) );  // CORRETTO: conversione esplicita da string_view a string

      for ( size_t pt_idx = 0; pt_idx < test_points.size(); ++pt_idx )
      {
        real_type x = test_points[pt_idx];

        // Calcola derivata spline
        real_type D_spline = ss.D( x, spline_idx );

        // Calcola differenza finita appropriata
        real_type D_fd;
        if ( x <= xmin + h ) { D_fd = finite_diff_forward( ss, x, spline_idx, h ); }
        else if ( x >= xmax - h ) { D_fd = finite_diff_backward( ss, x, spline_idx, h ); }
        else
        {
          D_fd = finite_diff_central( ss, x, spline_idx, h );
        }

        // Calcola errore relativo
        real_type abs_err = abs( D_spline - D_fd );
        real_type rel_err = 0.0;
        if ( abs( D_spline ) > 1e-10 ) { rel_err = 100.0 * abs_err / abs( D_spline ); }
        else if ( abs( D_fd ) > 1e-10 ) { rel_err = 100.0 * abs_err / abs( D_fd ); }

        // Determina colore in base all'errore
        // Stampa solo se punto interessante (nodo o errore significativo)
        bool is_knot = ( find( xx, xx + npts, x ) != xx + npts );
        if ( is_knot || rel_err > tol )
        {
          vector<string> row_values = { fmt::format( "{:.6f}", x ),
                                        spline_name,
                                        fmt::format( "{:.6f}", D_spline ),
                                        fmt::format( "{:.6f}", D_fd ),
                                        fmt::format( "{:.2e}", rel_err ) };
          print_table_row(
            row_values,
            rel_err < tol ? fmt::color::green : ( rel_err < 10.0 * tol ? fmt::color::yellow : fmt::color::red ) );
        }
      }

      // Separatore tra spline diverse
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

    // Statistiche
    integer internal_points = static_cast<integer>( test_points.size() ) - npts;
    fmt::print(
      fg( fmt::color::blue ),
      "\nPunti testati: {} ({} nodi + {} punti interni)\n\n",
      test_points.size(),
      npts,
      internal_points );
  }

  fmt::print(
    fg( fmt::color::cyan ) | fmt::emphasis::bold,
    "\n"
    "╔══════════════════════════════════════════════════════════╗\n"
    "║                  TUTTI I TEST COMPLETATI                 ║\n"
    "╚══════════════════════════════════════════════════════════╝\n\n" );

  return 0;
}
