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
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

#include "Splines.hh"
#include "Utils_fmt.hh"
#include "Utils_FD.hh"
#include "Utils_TridiagonalSolver.hh"
#include <set>

/*\
 |   ####  #    # #####  #  ####
 |  #    # #    # #    # # #    #
 |  #      #    # #####  # #
 |  #      #    # #    # # #
 |  #    # #    # #    # # #    #
 |   ####   ####  #####  #  ####
\*/

using namespace std;  // load standard namspace

namespace Splines
{

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*\
     Sistema lineare da risolvere

     D U UU
     L D U
       L D U
         L D U
           .....
              L D U
                L D U
                LL L D

  \*/
void CubicSpline_build(
  real_type const      X[],
  real_type const      Y[],
  real_type            Yp[],
  real_type            Ypp[],
  real_type            L[],
  real_type            D[],
  real_type            U[],
  integer const        npts,
  CubicSpline_BC const bc0,
  CubicSpline_BC const bcn )
{
  UTILS_ASSERT( npts >= 2, "CubicSpline_build, npts={} must be >= 2\n", npts );

  integer const n{ npts - 1 };
  real_type * Z{ Ypp };

  for ( integer i = 1; i < n; ++i )
  {
    real_type const HL{ X[i] - X[i - 1] };
    real_type const HR{ X[i + 1] - X[i] };
    real_type const HH{ HL + HR };
    L[i] = HL / HH;
    U[i] = HR / HH;
    D[i] = 2;
    Z[i] = 6 * ( ( Y[i + 1] - Y[i] ) / HR - ( Y[i] - Y[i - 1] ) / HL ) / HH;
  }

  real_type UU{ 0 }, LL{ 0 };

  switch ( bc0 )
  {
    case CubicSpline_BC::EXTRAPOLATE:
      L[0] = 0;
      D[0] = 1;
      U[0] = 0;
      if ( npts == 2 ) { 
        Z[0] = 0; 
      }
      else if ( npts == 3 )
      {
        Z[0] = Utils::second_derivative_3p(X[0], Y[0], X[1], Y[1], X[2], Y[2]);
      }
      else if ( npts == 4 )
      {
        Z[0] = Utils::second_derivative_4p(X[0], Y[0], X[1], Y[1], X[2], Y[2], X[3], Y[3]);
      }
      else
      {
        Z[0] = Utils::second_derivative_5p(X[0], Y[0], X[1], Y[1], X[2], Y[2], X[3], Y[3], X[4], Y[4]);
      }
      break;
    case CubicSpline_BC::NATURAL:
      L[0] = 0;
      D[0] = 1;
      U[0] = 0;
      Z[0] = 0;
      break;
    case CubicSpline_BC::PARABOLIC_RUNOUT:
      L[0] = 0;
      D[0] = 1;
      U[0] = -1;
      Z[0] = 0;
      break;
    case CubicSpline_BC::NOT_A_KNOT:
    {
      real_type const r = ( X[1] - X[0] ) / ( X[2] - X[1] );
      L[0] = 0;
      D[0] = 1;
      U[0] = -( 1 + r );
      UU    = r;  // elemento extra in posizione (0,2)
      Z[0] = 0;
    }
    break;
  }

  switch ( bcn )
  {
    case CubicSpline_BC::EXTRAPOLATE:
      L[n] = 0;
      D[n] = 1;
      U[n] = 0;
      if ( npts == 2 ) { 
        Z[n] = 0; 
      }
      else if ( npts == 3 )
      {
        Z[n] = Utils::second_derivative_3p(X[n], Y[n], X[n-1], Y[n-1], X[n-2], Y[n-2]);
      }
      else if ( npts == 4 )
      {
        Z[n] = Utils::second_derivative_4p(X[n], Y[n], X[n-1], Y[n-1], X[n-2], Y[n-2], X[n-3], Y[n-3]);
      }
      else
      {
        Z[n] = Utils::second_derivative_5p(X[n], Y[n], X[n-1], Y[n-1], X[n-2], Y[n-2], X[n-3], Y[n-3], X[n-4], Y[n-4]);
      }
      break;
    case CubicSpline_BC::NATURAL:
      L[n] = 0;
      D[n] = 1;
      U[n] = 0;
      Z[n] = 0;
      break;
    case CubicSpline_BC::PARABOLIC_RUNOUT:
      L[n] = -1;
      D[n] = 1;
      U[n] = 0;
      Z[n] = 0;
      break;
    case CubicSpline_BC::NOT_A_KNOT:
    {
      real_type const r = ( X[n - 1] - X[n - 2] ) / ( X[n] - X[n - 1] );
      U[n] = 0;
      D[n] = 1;
      L[n] = -( 1 + r );
      LL    = r;  // elemento extra in posizione (n, n-2)
      Z[n] = 0;
    }
    break;
  }

  // --- SEMPRE usa il solutore tridiagonale standard ---
  integer const m{ npts };
  Utils::TridiagonalSolver<real_type> solver;
  solver.resize( m );

  // Prepara i vettori per il sistema tridiagonale standard
  // a[i] = L[i+1] per i=0..m-2 (primo elemento a[0] corrisponde a L[1])
  // b[i] = D[i] per i=0..m-1
  // c[i] = U[i] per i=0..m-2 (c[m-1] non usato)
  Eigen::Matrix<real_type, Eigen::Dynamic, 1> a( m - 1 ), b( m ), c( m - 1 ), d( m );

  for ( integer i = 0; i < m - 1; ++i ) a( i ) = L[i + 1];
  for ( integer i = 0; i < m; ++i ) b( i ) = D[i];
  for ( integer i = 0; i < m - 1; ++i ) c( i ) = U[i];
  for ( integer i = 0; i < m; ++i ) d( i ) = Z[i];

  // Elimina separatamente UU e LL se presenti
  if ( UU != 0 ) {
    // UU è in posizione (0,2) → modifica la prima equazione
    // La prima equazione originale: D[0]*Z[0] + U[0]*Z[1] + UU*Z[2] = Z[0]
    // Portiamo UU*Z[2] a destra:
    d(0) -= UU * Z[2];  // Z[2] è ancora incognito, ma possiamo risolvere
    // Ora la matrice è tridiagonale standard
  }

  if ( LL != 0 ) {
    // LL è in posizione (n, n-2) → modifica l'ultima equazione
    // L'ultima equazione originale: LL*Z[n-2] + L[n]*Z[n-1] + D[n]*Z[n] = Z[n]
    // Portiamo LL*Z[n-2] a destra:
    d(n) -= LL * Z[n-2];
  }

  // Fattorizza e risolvi il sistema tridiagonale standard
  solver.factorize( a, b, c );
  Eigen::Matrix<real_type, Eigen::Dynamic, 1> x;
  solver.solve( a, b, d, x );

  // Copia la soluzione in Z
  for ( integer i = 0; i < m; ++i ) Z[i] = x( i );

  // Se UU o LL erano non zero, dobbiamo rifare una iterazione per aggiornare
  // i termini a destra che contenevano Z[2] e Z[n-2] ora noti
  if ( UU != 0 || LL != 0 ) {
    // Ricalcola i termini noti con i valori di Z appena trovati
    if ( UU != 0 ) d(0) = Z[0] - UU * Z[2];
    if ( LL != 0 ) d(n) = Z[n] - LL * Z[n-2];
    // Risolvi di nuovo (la fattorizzazione è già stata fatta)
    solver.solve( a, b, d, x );
    for ( integer i = 0; i < m; ++i ) Z[i] = x( i );
  }

  for ( integer i = 0; i < n; ++i )
  {
    real_type const DX = X[i + 1] - X[i];
    Yp[i]              = ( Y[i + 1] - Y[i] ) / DX - ( 2 * Z[i] + Z[i + 1] ) * ( DX / 6 );
  }
  real_type const DX2 = ( X[n] - X[n - 1] ) / 2;
  Yp[n]               = Yp[n - 1] + DX2 * ( Z[n - 1] + Z[n] );
}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void CubicSpline_build(
    real_type const      X[],
    real_type const      Y[],
    real_type            Yp[],
    integer const        npts,
    CubicSpline_BC const bc0,
    CubicSpline_BC const bcn )
  {
    Malloc_real mem( "CubicSpline_build" );
    mem.allocate( 4 * npts );
    real_type * L{ mem( npts ) };
    real_type * D{ mem( npts ) };
    real_type * U{ mem( npts ) };
    real_type * Z{ mem( npts ) };
    CubicSpline_build( X, Y, Yp, Z, L, D, U, npts, bc0, bcn );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void CubicSpline::build()
  {
    string msg{ fmt::format( "CubicSpline[{}]::build():", m_name ) };
    UTILS_ASSERT( m_npts > 1, "{} npts={} not enought points\n", msg, m_npts );
    Utils::check_NaN( m_X, msg + " X", m_npts, __LINE__, __FILE__ );
    Utils::check_NaN( m_Y, msg + " Y", m_npts, __LINE__, __FILE__ );
    integer ibegin{ 0 };
    integer iend{ 0 };
    do
    {
      // cerca intervallo monotono strettamente crescente
      for ( ++iend; iend < m_npts && m_X[iend - 1] < m_X[iend]; ++iend ) {}
      auto seg_bc0{ CubicSpline_BC::NOT_A_KNOT };
      auto seg_bcn{ CubicSpline_BC::NOT_A_KNOT };
      if ( ibegin == 0 ) seg_bc0 = m_bc0;
      if ( iend == m_npts ) seg_bcn = m_bcn;
      CubicSpline_build( m_X + ibegin, m_Y + ibegin, m_Yp + ibegin, iend - ibegin, seg_bc0, seg_bcn );
      ibegin = iend;
    } while ( iend < m_npts );

    Utils::check_NaN( m_Yp, msg + " Yp", m_npts, __LINE__, __FILE__ );
    m_search.must_reset();
  }

  using GC_namespace::GC_type;
  using GC_namespace::vec_real_type;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //!
  //!
  //! Setup a spline using a `GenericContainer`
  //!
  //! - gc("xdata") vector with the `x` coordinate of the data
  //! - gc("ydata") vector with the `y` coordinate of the data
  //!
  //! may contain
  //! - gc("bc_begin") and/or gc("bc_end")
  //!   - "extrapolate" extrapolate the boundary condition
  //!   - "natural"     make second derivative 0 at the border
  //!   - "parabolic"   make third derivative 0 at the border
  //!   - "not_a_knot"  not a knot condition of De Boor
  //!
  void CubicSpline::setup( GenericContainer const & gc )
  {
    /*
    // gc["xdata"]
    // gc["ydata"]
    //
    */
    string const where{ fmt::format( "CubicSpline[{}]::setup( gc ):", m_name ) };

    std::set<std::string> keywords;
    for ( auto const & pair : gc.get_map( where ) ) { keywords.insert( pair.first ); }
    keywords.erase( "spline_type" );

    GenericContainer const & gc_x{ gc( "xdata", where ) };
    keywords.erase( "xdata" );
    GenericContainer const & gc_y{ gc( "ydata", where ) };
    keywords.erase( "ydata" );

    vec_real_type x, y;
    {
      string const ff{ fmt::format( "{}, field `xdata'", where ) };
      gc_x.copyto_vec_real( x, ff );
    }
    {
      string const ff{ fmt::format( "{}, field `ydata'", where ) };
      gc_y.copyto_vec_real( y, ff );
    }
    if ( gc.exists( "bc_begin" ) )
    {
      keywords.erase( "bc_begin" );
      string_view bc{ gc.get_map_string( "bc_begin", where ) };
      if ( bc == "extrapolate" )
        m_bc0 = CubicSpline_BC::EXTRAPOLATE;
      else if ( bc == "natural" )
        m_bc0 = CubicSpline_BC::NATURAL;
      else if ( bc == "parabolic" )
        m_bc0 = CubicSpline_BC::PARABOLIC_RUNOUT;
      else if ( bc == "not_a_knot" )
        m_bc0 = CubicSpline_BC::NOT_A_KNOT;
      else
      {
        UTILS_ERROR( "{} unknow initial bc: {}\n", where, bc );
      }
    }
    else
    {
      UTILS_WARNING( false, "{}, missing field `bc_begin` using `extrapolate` as default value\n", where );
    }

    if ( gc.exists( "bc_end" ) )
    {
      keywords.erase( "bc_end" );
      string_view bc{ gc.get_map_string( "bc_end", where ) };
      if ( bc == "extrapolate" )
        m_bcn = CubicSpline_BC::EXTRAPOLATE;
      else if ( bc == "natural" )
        m_bcn = CubicSpline_BC::NATURAL;
      else if ( bc == "parabolic" )
        m_bcn = CubicSpline_BC::PARABOLIC_RUNOUT;
      else if ( bc == "not_a_knot" )
        m_bcn = CubicSpline_BC::NOT_A_KNOT;
      else
      {
        UTILS_ERROR( "{} unknow final bc: {}\n", where, bc );
      }
    }
    else
    {
      UTILS_WARNING( false, "{}, missing field `bc_end` using `extrapolate` as default value\n", where );
    }

    UTILS_WARNING(
      keywords.empty(),
      "{}: unused keys\n{}\n",
      where,
      [&keywords]() -> string
      {
        string res;
        for ( auto const & it : keywords )
        {
          res += it;
          res += ' ';
        };
        return res;
      }() );
    this->build( x, y );
  }

}  // namespace Splines
