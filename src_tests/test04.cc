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

static real_type xx0[] = { 0, 1, 2, 3, 4, 5, 5, 6, 7, 8, 9, 10 };
static real_type yy0[] = { 10, 10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85 };

static real_type xx1[] = { 0, 1, 3, 4, 6, 7, 7, 9, 10, 12, 13, 15 };
static real_type yy1[] = { 10, 10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85 };

static real_type xx2[] = { 0, 2, 3, 5, 6, 8, 8, 9, 11, 12, 14, 15 };
static real_type yy2[] = { 10, 10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85 };

// RPN 14
static real_type xx3[] = { 7.99, 8.09, 8.19, 8.7, 8.7, 9.2, 10, 12, 15, 20 };
static real_type yy3[] = { 0,        2.76429e-5, 4.37498e-2, 0.169183, 0.169183,
                           0.469428, 0.943740,   0.998636,   0.999919, 0.999994 };

// Titanium
static real_type xx4[] = { 595, 635, 695, 795, 855, 875, 875, 895, 915, 935, 985, 1035, 1075 };
static real_type yy4[] = { 0.644, 0.652, 0.644, 0.694, 0.907, 1.336, 1.336, 2.169, 1.598, 0.916, 0.607, 0.603, 0.608 };

// toolpath
static real_type xx5[] = { 0.11, 0.12, 0.15, 0.16 };
static real_type yy5[] = { 0.0003, 0.0003, 0.0004, 0.0004 };

static integer nn[] = { 11 + 1, 11 + 1, 11 + 1, 9 + 1, 12 + 1, 4 };

// Funzione per differenze finite centrali
real_type finite_diff_central( SplineSet const & ss, real_type x, integer i, real_type h = 1e-6 ) {
  return ( ss(x + h, i) - ss(x - h, i) ) / (2 * h);
}

// Funzione per differenze finite avanti/indietro
real_type finite_diff_forward( SplineSet const & ss, real_type x, integer i, real_type h = 1e-6 ) {
  return ( ss(x + h, i) - ss(x, i) ) / h;
}

real_type finite_diff_backward( SplineSet const & ss, real_type x, integer i, real_type h = 1e-6 ) {
  return ( ss(x, i) - ss(x - h, i) ) / h;
}

// Stampa tabella con bordi Unicode e colori
void print_table_header( vector<string> const & headers ) {
  // Linea superiore
  fmt::print( "┌{}", fmt::format( "{:─^{}}", "", 12 ) );
  for ( size_t i = 1; i < headers.size()-1; ++i ) {
    fmt::print( "┬{}", fmt::format( "{:─^{}}", "", 16 ) );
  }
  fmt::print( "┬{:─^{}}┐\n", "", 25 );

  // Intestazioni
  fmt::print( "│{:^12}", headers[0] );
  for ( size_t i = 1; i < headers.size()-1; ++i ) {
    fmt::print( "│{:^16}", headers[i] );
  }
  fmt::print( "│{:^25}│\n", headers.back() );

  // Separatore
  fmt::print( "├{}", fmt::format( "{:─^{}}", "", 12 ) );
  for ( size_t i = 1; i < headers.size()-1; ++i ) {
    fmt::print( "┼{}", fmt::format( "{:─^{}}", "", 16 ) );
  }
  fmt::print( "┼{:─^{}}┤\n", "", 25 );
}

template <typename COLOR>
void print_table_row( vector<string> const & values, COLOR color ) {
  size_t i = 0;
  fmt::print( "│{:^12}", values[i] );
  for ( ++i; i < values.size()-2; ++i ) {
    fmt::print( "│{:^16}", values[i] );
  }
  fmt::print( "│");
  fmt::print( fg(color), "{:^16}", values[i] );
  fmt::print( "│{:^25}│\n", values.back() );
}

void print_table_footer( size_t ncols ) {
  fmt::print( "└{}", fmt::format( "{:─^{}}", "", 12 ) );
  for ( size_t i = 1; i < ncols-1; ++i ) {
    fmt::print( "┴{}", fmt::format( "{:─^{}}", "", 16 ) );
  }
  fmt::print( "┴{:─^{}}┘\n", "", 25 );
}

// Controlla se un punto è un nodo duplicato (punto di discontinuità)
bool is_duplicate_knot( real_type const * xx, integer npts, real_type x, integer & first_index, real_type eps = 1e-12 ) {
  first_index = -1;
  integer count = 0;
  for ( integer i = 0; i < npts; ++i ) {
    if ( abs(xx[i] - x) < eps ) {
      if ( first_index == -1 ) first_index = i;
      ++count;
      if ( count > 1 ) return true;
    }
  }
  return false;
}

// Verifica se un punto è troppo vicino a un nodo duplicato per evitare problemi nelle differenze finite
bool too_close_to_duplicate( real_type const * xx, integer npts, real_type x, real_type safety_margin = 1e-4 ) {
  for ( integer i = 0; i < npts; ++i ) {
    // Cerca nodi duplicati
    if ( i > 0 && abs(xx[i] - xx[i-1]) < 1e-12 ) {
      // Questo è un nodo duplicato
      if ( abs(x - xx[i]) < safety_margin ) {
        return true;
      }
    }
  }
  return false;
}

int main()
{
  fmt::print( fg(fmt::color::cyan) | fmt::emphasis::bold,
    "\n"
    "╔══════════════════════════════════════════════════════════╗\n"
    "║                         TEST N.4                         ║\n"
    "╚══════════════════════════════════════════════════════════╝\n\n" );

  SplineSet ss;
  ofstream  file, file_D;

  for ( integer k{ 0 }; k < 6; ++k )
  {
    real_type const * xx{ nullptr };
    real_type const * yy{ nullptr };
    string dataset_name;
    switch ( k )
    {
      case 0:
        xx = xx0; yy = yy0; dataset_name = "Akima Test 0 (con duplicato)"; break;
      case 1:
        xx = xx1; yy = yy1; dataset_name = "Akima Test 1 (con duplicato)"; break;
      case 2:
        xx = xx2; yy = yy2; dataset_name = "Akima Test 2 (con duplicato)"; break;
      case 3:
        xx = xx3; yy = yy3; dataset_name = "RPN 14 (con duplicato)"; break;
      case 4:
        xx = xx4; yy = yy4; dataset_name = "Titanium (con duplicato)"; break;
      case 5:
        xx = xx5; yy = yy5; dataset_name = "Toolpath"; break;
    }
    
    fmt::print( fg(fmt::color::yellow) | fmt::emphasis::bold,
      "\n"
      "┌──────────────────────────────────────────────────────────┐\n"
      "│ Dataset {}: {:42} │\n"
      "└──────────────────────────────────────────────────────────┘\n",
      k, dataset_name
    );

    string const fname{ fmt::format( "out/SplineSet{}.txt", k ) };
    file.open( fname.data() );
    
    string const fname_D{ fmt::format( "out/SplineSet{}_D.txt", k ) };
    file_D.open( fname_D.data() );

    real_type const xmin{ xx[0] };
    real_type const xmax{ xx[nn[k] - 1] };

    constexpr integer nspl      = 7;
    integer const     npts      = nn[k];
    char const *      headers[] = { "SPLINE_CONSTANT", "SPLINE_LINEAR", "SPLINE_AKIMA",  "SPLINE_BESSEL",
                                    "SPLINE_PCHIP",    "SPLINE_CUBIC",  "SPLINE_QUINTIC" };

    constexpr SplineType1D stype[]{ SplineType1D::CONSTANT, SplineType1D::LINEAR, SplineType1D::AKIMA,
                                    SplineType1D::BESSEL,   SplineType1D::PCHIP,  SplineType1D::CUBIC,
                                    SplineType1D::QUINTIC };

    real_type const * Y[] = { yy, yy, yy, yy, yy, yy, yy };

    ss.build( nspl, npts, headers, stype, xx, Y );

    // Scrittura file valori delle spline
    file << "x";
    for ( integer i{ 0 }; i < nspl; ++i ) file << '\t' << ss.header( i );
    file << '\n';
    for ( real_type x{ xmin }; x <= xmax; x += ( xmax - xmin ) / 1000 )
    {
      file << x;
      for ( integer i{ 0 }; i < nspl; ++i ) file << '\t' << ss( x, i );
      file << '\n';
    }
    file.close();
    
    // Scrittura file derivate delle spline
    file_D << "x";
    for ( integer i{ 0 }; i < nspl; ++i ) file_D << '\t' << ss.header( i );
    file_D << '\n';
    for ( real_type x{ xmin }; x <= xmax; x += ( xmax - xmin ) / 1000 )
    {
      file_D << x;
      for ( integer i{ 0 }; i < nspl; ++i ) file_D << '\t' << ss.D( x, i );
      file_D << '\n';
    }
    file_D.close();
    
    ss.info( cout );

    // ========== CONTROLLO DERIVATE CON DIFFERENZE FINITE ==========
    
    // Punti di test: punti interni lontani dai nodi duplicati
    vector<real_type> test_points;
    
    // Identifica intervalli validi (dove x[i] != x[i+1])
    vector<pair<real_type, real_type>> valid_intervals;
    for ( integer i = 0; i < npts-1; ++i ) {
      if ( abs(xx[i+1] - xx[i]) > 1e-12 ) {
        valid_intervals.emplace_back(xx[i], xx[i+1]);
      }
    }
    
    // Genera punti di test negli intervalli validi
    for ( const auto& interval : valid_intervals ) {
      real_type x1 = interval.first;
      real_type x2 = interval.second;
      real_type dx = (x2 - x1) / 10.0;
      
      // Aggiungi punti interni (evitando i bordi degli intervalli)
      for ( integer j = 1; j < 10; ++j ) {
        real_type x = x1 + j * dx;
        test_points.push_back(x);
      }
    }
    
    // Aggiungi anche i nodi singoli (non duplicati) se non sono troppo vicini a duplicati
    for ( integer i = 0; i < npts; ++i ) {
      integer first_index;
      if ( !is_duplicate_knot(xx, npts, xx[i], first_index) && 
           !too_close_to_duplicate(xx, npts, xx[i]) ) {
        test_points.push_back(xx[i]);
      }
    }
    
    // Ordina e rimuovi duplicati approssimati
    sort( test_points.begin(), test_points.end() );
    test_points.erase( unique( test_points.begin(), test_points.end(), 
                               [](real_type a, real_type b) { return abs(a - b) < 1e-12; } ), 
                       test_points.end() );

    // Tabella per il controllo delle derivate
    fmt::print( fg(fmt::color::green), "\nControllo derivate - Dataset {} ({} punti test validi)\n", 
                k, test_points.size() );
    
    vector<string> table_headers = {"x", "Spline", "D(x)", "FinDiff", "Err. Rel%", "Note"};
    print_table_header( table_headers );
    
    real_type const h = 1e-6;
    real_type const tol = 1e-4;
    integer total_points_tested = 0;
    integer points_with_error = 0;
    integer points_skipped = 0;
    
    for ( integer spline_idx = 0; spline_idx < nspl; ++spline_idx ) {
      string spline_name = string(ss.header( spline_idx ));
      
      for ( size_t pt_idx = 0; pt_idx < test_points.size(); ++pt_idx ) {
        real_type x = test_points[pt_idx];
        
        // Determina se è un nodo
        bool is_knot = ( find( xx, xx + npts, x ) != xx + npts );
        integer first_index;
        bool is_dup_knot = is_duplicate_knot(xx, npts, x, first_index);
        
        // Salta i nodi duplicati e punti troppo vicini ad essi
        if ( is_dup_knot || too_close_to_duplicate(xx, npts, x) ) {
          ++points_skipped;
          continue;
        }
        
        // Calcola derivata spline
        real_type D_spline = ss.D( x, spline_idx );
        
        // Calcola differenza finita appropriata
        real_type D_fd;
        string note = is_knot ? "nodo" : "interno";
        
        // Per i nodi (non duplicati), usa differenze finite unilaterali appropriate
        if ( is_knot ) {
          if ( x <= xmin + h ) {
            D_fd = finite_diff_forward( ss, x, spline_idx, h );
            note += " (estremo sinistro)";
          } else if ( x >= xmax - h ) {
            D_fd = finite_diff_backward( ss, x, spline_idx, h );
            note += " (estremo destro)";
          } else {
            // Per nodi interni, potremmo avere discontinuità, usa differenza finita dal lato appropriato
            // Controlla se è l'inizio o la fine di un intervallo
            bool is_start_of_interval = false;
            for ( integer i = 0; i < npts-1; ++i ) {
              if ( abs(xx[i] - x) < 1e-12 && abs(xx[i+1] - xx[i]) > 1e-12 ) {
                is_start_of_interval = true;
                break;
              }
            }
            if ( is_start_of_interval ) {
              D_fd = finite_diff_forward( ss, x, spline_idx, h );
              note += " (inizio intervallo)";
            } else {
              D_fd = finite_diff_backward( ss, x, spline_idx, h );
              note += " (fine intervallo)";
            }
          }
        } else {
          // Per punti interni, usa differenza finita centrale
          D_fd = finite_diff_central( ss, x, spline_idx, h );
        }
        
        // Calcola errore relativo
        real_type abs_err = abs( D_spline - D_fd );
        real_type rel_err = 0.0;
        if ( abs( D_spline ) > 1e-10 ) {
          rel_err = 100.0 * abs_err / abs( D_spline );
        } else if ( abs( D_fd ) > 1e-10 ) {
          rel_err = 100.0 * abs_err / abs( D_fd );
        } else {
          // Entrambi zero, errore zero
          rel_err = 0.0;
        }
        
        ++total_points_tested;
        if ( rel_err > tol ) ++points_with_error;
                
        // Stampa solo se: nodo, o errore significativo, o per campionamento ridotto
        bool should_print = ( is_knot || rel_err > tol || pt_idx % 5 == 0 );
        
        if ( should_print ) {
          vector<string> row_values = {
            fmt::format( "{:.6f}", x ),
            spline_name,
            fmt::format( "{:.6f}", D_spline ),
            fmt::format( "{:.6f}", D_fd ),
            fmt::format( "{:.2e}", rel_err ),
            note
          };
          print_table_row( row_values, rel_err < tol ? fmt::color::green : (rel_err < 10.0 * tol ? fmt::color::yellow : fmt::color::red )  );
        }
      }
      
      // Separatore tra spline diverse
      if ( spline_idx < nspl - 1 ) {
        fmt::print( "├{}┼{}┼{}┼{}┼{}┼{}┤\n",
          fmt::format( "{:─^{}}", "", 12 ),
          fmt::format( "{:─^{}}", "", 16 ),
          fmt::format( "{:─^{}}", "", 16 ),
          fmt::format( "{:─^{}}", "", 16 ),
          fmt::format( "{:─^{}}", "", 16 ),
          fmt::format( "{:─^{}}", "", 25 ) );
      }
    }
    
    print_table_footer( table_headers.size() );
    
    // Statistiche
    fmt::print( fg(fmt::color::blue), 
                "\n📊 Statistiche Dataset {}:\n"
                "   • Punti test validi: {}\n"
                "   • Punti saltati (vicini a discontinuità): {}\n"
                "   • Punti con errore > {:.0e}: {} ({:.1f}%)\n\n",
                k, total_points_tested, points_skipped, tol, points_with_error,
                total_points_tested > 0 ? 100.0 * points_with_error / total_points_tested : 0.0 );
  }

  fmt::print( fg(fmt::color::cyan) | fmt::emphasis::bold,
    "\n"
    "╔══════════════════════════════════════════════════════════╗\n"
    "║                  TUTTI I TEST COMPLETATI                 ║\n"
    "╚══════════════════════════════════════════════════════════╝\n\n"
  );

  return 0;
}
