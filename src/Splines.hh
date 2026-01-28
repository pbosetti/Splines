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

#pragma once

#ifndef SPLINES_HH
#define SPLINES_HH

#ifdef __GNUC__
#pragma GCC diagnostic push
#endif

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

#include "SplinesConfig.hh"

#include "Utils_fmt.hh"
#include "Utils_FD.hh"
#include "Utils_search_intervals2.hh"
#include "Utils_TridiagonalSolver.hh"
#include "Utils_AlgoHNewton.hh"

//!
//! Namespace of Splines library
//!
namespace Splines
{

  using std::basic_istream;
  using std::basic_ostream;
  using std::cerr;
  using std::cin;
  using std::copy_n;
  using std::cout;
  using std::exception;
  using std::lower_bound;
  using std::pair;
  using std::runtime_error;
  using std::string;
  using std::string_view;
  using std::vector;

  typedef double real_type;  //!< Floating point type for splines
  typedef int    integer;    //!< Signed integer type for splines

  using Malloc_real  = Utils::Malloc<real_type>;
  using ostream_type = basic_ostream<char>;
  using istream_type = basic_istream<char>;

  using Vec  = Eigen::Array<real_type, Eigen::Dynamic, 1>;
  using Mat  = Eigen::Array<real_type, Eigen::Dynamic, Eigen::Dynamic>;
  using MatC = Eigen::Array<real_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

  using Mat2x2 = Eigen::Matrix<real_type, 2, 2>;
  using Mat4x4 = Eigen::Matrix<real_type, 4, 4>;
  using Mat6x6 = Eigen::Matrix<real_type, 6, 6>;

  using GC_namespace::GC_type;
  using GC_namespace::map_type;
  using GC_namespace::mat_real_type;
  using GC_namespace::vec_int_type;
  using GC_namespace::vec_real_type;
  using GC_namespace::vec_string_type;
  using GC_namespace::vector_type;

  void backtrace( ostream_type & );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! Associate a number for each type of splines implemented
  using SplineType1D = enum class SplineType1D : integer {
    CONSTANT       = 0,
    LINEAR         = 1,
    CUBIC          = 2,
    AKIMA          = 3,
    BESSEL         = 4,
    PCHIP          = 5,
    QUINTIC_CUBIC  = 6,
    QUINTIC_AKIMA  = 7,
    QUINTIC_BESSEL = 8,
    QUINTIC_PCHIP  = 9,
    HERMITE        = 10,
    SPLINE_SET     = 11,
    SPLINE_VEC     = 12
  };

  using Spline_sub_type = enum class Spline_sub_type : integer { CUBIC = 0, PCHIP = 1, AKIMA = 2, BESSEL = 3 };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! Associate a number for each type of splines implemented
  using SplineType2D = enum class SplineType2D : integer {
    BILINEAR         = 0,
    BICUBIC_CUBIC    = 1,
    BICUBIC_AKIMA    = 2,
    BICUBIC_BESSEL   = 3,
    BICUBIC_PCHIP    = 4,
    BIQUINTIC_CUBIC  = 5,
    BIQUINTIC_AKIMA  = 6,
    BIQUINTIC_BESSEL = 7,
    BIQUINTIC_PCHIP  = 8,
  };

  inline SplineType1D string_to_splineType1D( string_view nin )
  {
    string n{ nin };
    std::transform( n.begin(), n.end(), n.begin(), ::tolower );
    if ( n == "constant" ) return SplineType1D::CONSTANT;
    if ( n == "linear" ) return SplineType1D::LINEAR;

    if ( n == "cubic" ) return SplineType1D::CUBIC;
    if ( n == "akima" ) return SplineType1D::AKIMA;
    if ( n == "bessel" ) return SplineType1D::BESSEL;
    if ( n == "pchip" ) return SplineType1D::PCHIP;

    if ( n == "quintic" ) return SplineType1D::QUINTIC_CUBIC;
    if ( n == "quintic_akima" ) return SplineType1D::QUINTIC_AKIMA;
    if ( n == "quintic_bessel" ) return SplineType1D::QUINTIC_BESSEL;
    if ( n == "quintic_pchip" ) return SplineType1D::QUINTIC_PCHIP;

    if ( n == "hermite" ) return SplineType1D::HERMITE;
    if ( n == "spline_set" ) return SplineType1D::SPLINE_SET;
    if ( n == "spline_vec" ) return SplineType1D::SPLINE_VEC;
    throw std::runtime_error( fmt::format( "string_to_splineType1D({}) unknown type\n", n ) );
  }

  inline std::string to_string( Spline_sub_type t )
  {
    switch ( t )
    {
      case Spline_sub_type::CUBIC: return "CUBIC";
      case Spline_sub_type::PCHIP: return "PCHIP";
      case Spline_sub_type::AKIMA: return "AKIMA";
      case Spline_sub_type::BESSEL: return "BESSEL";
    }
    return "NOTYPE";
  }

  inline SplineType2D string_to_splineType2D( string_view nin )
  {
    string n{ nin };
    std::transform( n.begin(), n.end(), n.begin(), ::tolower );
    if ( n == "bilinear" ) return SplineType2D::BILINEAR;

    if ( n == "bicubic" ) return SplineType2D::BICUBIC_CUBIC;
    if ( n == "bicubic_akima" ) return SplineType2D::BICUBIC_AKIMA;
    if ( n == "bicubic_bessel" ) return SplineType2D::BICUBIC_BESSEL;
    if ( n == "bicubic_pchip" ) return SplineType2D::BICUBIC_PCHIP;

    if ( n == "biquintic" ) return SplineType2D::BIQUINTIC_CUBIC;
    if ( n == "biquintic_akima" ) return SplineType2D::BIQUINTIC_AKIMA;
    if ( n == "biquintic_bessel" ) return SplineType2D::BIQUINTIC_BESSEL;
    if ( n == "biquintic_pchip" ) return SplineType2D::BIQUINTIC_PCHIP;

    throw std::runtime_error( fmt::format( "string_to_splineType2D({}) unknown type\n", n ) );
  }

  inline char const * to_string( SplineType1D const t )
  {
    switch ( t )
    {
      case SplineType1D::CONSTANT: return "SPLINE_CONSTANT";
      case SplineType1D::LINEAR: return "SPLINE_LINEAR";
      case SplineType1D::CUBIC: return "SPLINE_CUBIC";
      case SplineType1D::AKIMA: return "SPLINE_AKIMA";
      case SplineType1D::BESSEL: return "SPLINE_BESSEL";
      case SplineType1D::PCHIP: return "SPLINE_PCHIP";
      case SplineType1D::QUINTIC_CUBIC: return "SPLINE_QUINTIC";
      case SplineType1D::QUINTIC_AKIMA: return "SPLINE_QUINTIC_AKIMA";
      case SplineType1D::QUINTIC_BESSEL: return "SPLINE_QUINTIC_BESSEL";
      case SplineType1D::QUINTIC_PCHIP: return "SPLINE_QUINTIC_PCHIP";
      case SplineType1D::HERMITE: return "SPLINE_HERMITE";
      case SplineType1D::SPLINE_SET: return "SPLINE_SPLINE_SET";
      case SplineType1D::SPLINE_VEC: return "SPLINE_SPLINE_VEC";
    }
    return "NO_TYPE";
  }

  inline char const * to_string( SplineType2D const t )
  {
    switch ( t )
    {
      case SplineType2D::BILINEAR: return "SPLINE2D_BILINEAR";
      case SplineType2D::BICUBIC_CUBIC: return "SPLINE2D_BICUBIC_CUBIC";
      case SplineType2D::BICUBIC_AKIMA: return "SPLINE2D_BICUBIC_AKIMA";
      case SplineType2D::BICUBIC_BESSEL: return "SPLINE2D_BICUBIC_BESSEL";
      case SplineType2D::BICUBIC_PCHIP: return "SPLINE2D_BICUBIC_PCHIP";
      case SplineType2D::BIQUINTIC_CUBIC: return "SPLINE2D_BIQUINTIC_CUBIC";
      case SplineType2D::BIQUINTIC_AKIMA: return "SPLINE2D_BIQUINTIC_AKIMA";
      case SplineType2D::BIQUINTIC_BESSEL: return "SPLINE2D_BIQUINTIC_BESSEL";
      case SplineType2D::BIQUINTIC_PCHIP: return "SPLINE2D_BIQUINTIC_PCHIP";
    }
    return "NO_TYPE";
  }

  using GC_namespace::GC_type;
  using GC_namespace::GenericContainer;
  using GC_namespace::map_type;
  using GC_namespace::vec_real_type;
  using GC_namespace::vec_string_type;
  using GC_namespace::vector_type;

  /*
  //   _   _                     _ _
  //  | | | | ___ _ __ _ __ ___ (_) |_ ___
  //  | |_| |/ _ \ '__| '_ ` _ \| | __/ _ \
  //  |  _  |  __/ |  | | | | | | | ||  __/
  //  |_| |_|\___|_|  |_| |_| |_|_|\__\___|
  */

#ifdef AUTODIFF_SUPPORT
  template <typename T> inline void Hermite3( T const & x, real_type const H, T base[4] )
  {
    T const X = x / H;
    base[1]   = X * X * ( 3 - 2 * X );
    base[0]   = 1 - base[1];
    base[2]   = x * ( X * ( X - 2 ) + 1 );
    base[3]   = x * X * ( X - 1 );
  }

  template <typename T> inline void Hermite5( T const & x, real_type H, T base[6] )
  {
    const real_type invH  = real_type( 1 ) / H;
    const real_type invH2 = invH * invH;
    const real_type invH4 = invH2 * invH2;

    const T x2 = x * x;
    const T x3 = x2 * x;
    const T x4 = x2 * x2;
    const T x5 = x4 * x;

    const T xm  = H - x;
    const T xm2 = xm * xm;
    const T xm3 = xm2 * xm;

    base[0] = invH4 * xm3 * ( 3 * x * H + H * H + 6 * x2 );
    base[1] = invH4 * ( -15 * H * x4 + 6 * x5 + 10 * H * H * x3 );
    base[2] = invH4 * xm3 * x * ( H + 3 * x );
    base[3] = invH4 * ( 3 * x - 4 * H ) * xm * x3;

    const real_type c = real_type( 0.5 ) * invH2 * invH;

    base[4] = c * xm3 * x2;
    base[5] = c * xm2 * x3;
  }
#endif

  inline void Hermite3( real_type const x, real_type const H, real_type base[4] )
  {
    real_type const X = x / H;
    base[1]           = X * X * ( 3 - 2 * X );
    base[0]           = 1 - base[1];
    base[2]           = x * ( X * ( X - 2 ) + 1 );
    base[3]           = x * X * ( X - 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline void Hermite3_D( real_type const x, real_type const H, real_type base_D[4] )
  {
    real_type const X = x / H;
    base_D[0]         = 6.0 * X * ( X - 1.0 ) / H;
    base_D[1]         = -base_D[0];
    base_D[2]         = ( ( 3 * X - 4 ) * X + 1 );
    base_D[3]         = X * ( 3 * X - 2 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline void Hermite3_DD( real_type const x, real_type const H, real_type base_DD[4] )
  {
    real_type const X = x / H;
    base_DD[0]        = ( 12 * X - 6 ) / ( H * H );
    base_DD[1]        = -base_DD[0];
    base_DD[2]        = ( 6 * X - 4 ) / H;
    base_DD[3]        = ( 6 * X - 2 ) / H;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline void Hermite3_DDD( real_type, real_type const H, real_type base_DDD[4] )
  {
    base_DDD[0] = 12 / ( H * H * H );
    base_DDD[1] = -base_DDD[0];
    base_DDD[2] = 6 / ( H * H );
    base_DDD[3] = base_DDD[2];
  }

  // --------------------------------------------------------------------------

  inline void Hermite5( real_type x, real_type H, real_type base[6] )
  {
    // 1. Normalizzazione: t va da 0 a 1
    // Sostituisce le costose divisioni ripetute (invH, invH4, invH5)
    const real_type s = 1 / H;
    const real_type t = x * s;
    const real_type u = 1 - t;  // u è il complemento (H-x)/H

    // 2. Precalcolo delle potenze di t e u
    // Usiamo t*t invece di pow() per velocità
    const real_type t2 = t * t;
    const real_type t3 = t2 * t;

    const real_type u2 = u * u;
    const real_type u3 = u2 * u;

    // 3. Calcolo delle funzioni di base
    // Polinomio base: h00 = (1 + 3t + 6t^2) * u^3
    const real_type poly_t = 1 + t * ( 3 + 6 * t );
    const real_type poly_u = 1 + u * ( 3 + 6 * u );

    // base[0]: Valore al nodo sinistro (x=0)
    base[0] = u3 * poly_t;

    // base[1]: Valore al nodo destro (x=H)
    // Sfrutta la simmetria: scambia t con u
    base[1] = t3 * poly_u;

    // base[2]: Derivata prima al nodo sinistro
    // Scala originale: H * t * (1+3t) * u^3
    base[2] = H * t * u3 * ( 1 + 3 * t );

    // base[3]: Derivata prima al nodo destro
    // Scala originale: -H * u * (1+3u) * t^3
    base[3] = -H * u * t3 * ( 1 + 3 * u );

    // base[4]: Derivata seconda al nodo sinistro
    // Scala originale: 0.5 * H^2 * t^2 * u^3
    const real_type H2_half = ( H * H ) / 2;
    base[4]                 = H2_half * t2 * u3;

    // base[5]: Derivata seconda al nodo destro
    // Scala originale: 0.5 * H^2 * u^2 * t^3
    base[5] = H2_half * u2 * t3;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline void Hermite5_D( real_type x, real_type H, real_type base_D[6] )
  {
    // 1. Normalizzazione
    // t = x/H, u = (H-x)/H = 1-t
    const real_type s = 1 / H;
    const real_type t = x * s;
    const real_type u = 1 - t;

    // 2. Precalcolo potenze
    const real_type t2 = t * t;
    const real_type u2 = u * u;

    // 3. Calcolo Derivate
    // Nota: La derivata rispetto a x scala di (1/H) rispetto alla derivata su t.

    // --- Basi 0 e 1 (Valori ai nodi) ---
    // Derivata del Quintic Step: 30 * t^2 * u^2 / H
    // common = 30 * t^2 * u^2 * s
    const real_type common = ( 30 * s ) * t2 * u2;

    base_D[0] = -common;
    base_D[1] = common;

    // --- Basi 2 e 3 (Derivate Prime ai nodi) ---
    // Polinomio derivato base: p(t) = (1 - 3t) * (1 + 5t)
    // base_D[2] = u^2 * p(t)
    // base_D[3] = t^2 * p(u)  <-- Simmetria perfetta

    // Espressione fattorizzata per stabilità e velocità
    base_D[2] = u2 * ( 1 - 3 * t ) * ( 1 + 5 * t );

    base_D[3] = t2 * ( 1 - 3 * u ) * ( 1 + 5 * u );

    // --- Basi 4 e 5 (Derivate Seconde ai nodi) ---
    // Polinomio derivato base: q(t) = t * (2 - 5t)
    // Scala fisica: 0.5 * H

    const real_type h_half = H / 2;

    base_D[4] = h_half * u2 * t * ( 2 - 5 * t );

    // Nota il segno meno per simmetria speculare sulla derivata
    base_D[5] = -h_half * t2 * u * ( 2 - 5 * u );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline void Hermite5_DD( real_type x, real_type H, real_type base_DD[6] )
  {
    // 1. Normalizzazione
    // s = 1/H, s2 = 1/H^2
    const real_type s  = 1 / H;
    const real_type s2 = s * s;

    // t = x/H (da 0 a 1), u = 1 - t
    const real_type t = x * s;
    const real_type u = 1 - t;

    // 2. Termini comuni
    // Molte basi condividono il termine t * u
    const real_type tu = t * u;

    // --- Basi 0 e 1 (Valori ai nodi) ---
    // Formula originale normalizzata: 60/H^2 * t * u * (1 - 2t)
    // Nota: (1 - 2t) è equivalente a (u - t)

    const real_type common = ( 60 * s2 ) * tu * ( u - t );

    base_DD[0] = -common;
    base_DD[1] = common;

    // --- Basi 2 e 3 (Derivate Prime ai nodi) ---
    // Scala fisica: 1/H
    // Base 2: 12/H * t * u * (5t - 3)
    // Base 3: 12/H * t * u * (5t - 2)

    const real_type k_tan = ( 12 * s ) * tu;

    base_DD[2] = k_tan * ( 5 * t - 3 );
    base_DD[3] = k_tan * ( 5 * t - 2 );

    // --- Basi 4 e 5 (Derivate Seconde ai nodi) ---
    // Scala fisica: 1 (Adimensionale rispetto a H nella derivata seconda)
    // Base 4: u * (1 - 8t + 10t^2)
    // Base 5: t * (3 - 12t + 10t^2)

    // Sfruttiamo la simmetria: Base5(t) = Base4(u)
    // Calcoliamo il polinomio 'p' usando t per base 4, e u per base 5.
    // p(v) = 1 - 8v + 10v^2

    // Base 4 (lato sinistro, usa t per il polinomio interno)
    // Nota: L'espressione u * (...) deriva dalla fattorizzazione
    base_DD[4] = u * ( 1 + t * ( 10 * t - 8 ) );

    // Base 5 (lato destro, usa u per il polinomio interno per simmetria)
    // Equivale a: t * (3 - 12t + 10t^2)
    base_DD[5] = t * ( 1 + u * ( 10 * u - 8 ) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  inline void Hermite5_DDD( real_type x, real_type H, real_type base_DDD[6] )
  {
    // 1. Normalizzazione
    const real_type s = 1 / H;
    const real_type t = x * s;

    // t^2 per valutare le quadratiche
    const real_type t2 = t * t;

    // 2. Fattori di scala precalcolati
    // Base 0,1 scalano con 1/H^3
    // Base 2,3 scalano con 1/H^2
    // Base 4,5 scalano con 1/H
    const real_type s2 = s * s;
    const real_type s3 = s2 * s;

    // --- Basi 0 e 1 ---
    // Polinomio: 60 * (6t - 6t^2 - 1)
    const real_type k0     = 60 * s3;
    const real_type common = k0 * ( 6 * ( t - t2 ) - 1 );

    base_DDD[0] = common;
    base_DDD[1] = -common;

    // --- Basi 2 e 3 ---
    // Base 2: 12 * (16t - 15t^2 - 3)
    // Base 3: 12 * (14t - 15t^2 - 2)
    const real_type k1 = real_type( 12 ) * s2;

    // Fattorizzazione di Horner per risparmiare 1 mul: t*(16 - 15t) - 3
    base_DDD[2] = k1 * ( t * ( 16 - 15 * t ) - real_type( 3 ) );

    // Fattorizzazione: t*(14 - 15t) - 2
    base_DDD[3] = k1 * ( t * ( 14 - 15 * t ) - real_type( 2 ) );

    // --- Basi 4 e 5 ---
    // Base 4: 3 * (12t - 10t^2 - 3)
    // Base 5: 3 * (10t^2 - 8t + 1)
    const real_type k2 = 3 * s;

    base_DDD[4] = k2 * ( t * ( 12 - 10 * t ) - 3 );

    base_DDD[5] = k2 * ( t * ( 10 * t - 8 ) + 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline void Hermite5_DDDD( real_type x, real_type H, real_type base_DDDD[6] )
  {
    // 1. Normalizzazione
    // s = 1/H, t = x/H
    const real_type s = 1 / H;
    const real_type t = x * s;

    // 2. Fattori di scala
    // Base 0,1 scalano con 1/H^4
    // Base 2,3 scalano con 1/H^3
    // Base 4,5 scalano con 1/H^2
    const real_type s2 = s * s;
    const real_type s3 = s2 * s;
    const real_type s4 = s2 * s2;

    // --- Basi 0 e 1 ---
    // Originale: (360*H - 720*x) / H^5
    // Normalizzato: 360/H^4 * (1 - 2t)
    const real_type k0    = 360 * s4;
    const real_type term0 = k0 * ( 1 - 2 * t );

    base_DDDD[0] = term0;
    base_DDDD[1] = -term0;

    // --- Basi 2 e 3 ---
    // Costante comune: 24 / H^3
    const real_type k1 = 24 * s3;

    // Base 2: k1 * (8 - 15t)
    base_DDDD[2] = k1 * ( 8 - 15 * t );

    // Base 3: k1 * (7 - 15t)
    base_DDDD[3] = k1 * ( 7 - 15 * t );

    // --- Basi 4 e 5 ---
    // Costante comune: 12 / H^2
    const real_type k2 = 12 * s2;

    // Base 4: k2 * (3 - 5t)
    base_DDDD[4] = k2 * ( 3 - 5 * t );

    // Base 5: k2 * (5t - 2)
    base_DDDD[5] = k2 * ( 5 * t - 2 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline void Hermite5_DDDDD( real_type /*x*/, real_type H, real_type base_DDDDD[6] )
  {
    // 1. Unica divisione necessaria
    const real_type s = 1 / H;

    // 2. Calcolo delle potenze di s (moltiplicazioni)
    // s3 = 1/H^3, s4 = 1/H^4, s5 = 1/H^5
    const real_type s2 = s * s;
    const real_type s3 = s2 * s;
    const real_type s4 = s2 * s2;  // Oppure s3 * s
    const real_type s5 = s4 * s;

    // 3. Calcolo delle costanti
    // Nota: Le derivate sono costanti, quindi calcoliamo solo i valori scalari.

    // --- Basi 0 e 1 (Valore) ---
    // Derivata 5a scala con 1/H^5
    const real_type c0 = 720 * s5;

    base_DDDDD[0] = -c0;
    base_DDDDD[1] = c0;

    // --- Basi 2 e 3 (Derivata) ---
    // Derivata 5a scala con 1/H^4 (poiché H^1 * H^-5 = H^-4)
    const real_type c1 = 360 * s4;

    // Entrambe negative perché il coefficiente di x^5 nelle basi originali
    // ha lo stesso segno relativo per le derivate (antisimmetria nello sviluppo)
    base_DDDDD[2] = -c1;
    base_DDDDD[3] = -c1;

    // --- Basi 4 e 5 (Curvatura) ---
    // Derivata 5a scala con 1/H^3 (poiché H^2 * H^-5 = H^-3)
    const real_type c2 = 60 * s3;

    base_DDDDD[4] = -c2;
    base_DDDDD[5] = c2;
  }

  //!
  //! Convert polynomial defined using Hermite base
  //!
  //! \f[ p(x) = p_0 h_0(t) + p_1 h_1(t) +  p'_0 h_2(t) + p'_1 h_3(t) \f]
  //!
  //! to standard form
  //!
  //! \f[ p(x) = A t^3 + B t^2 + C t + D \f]
  //!
  //!
  //!
  static inline void Hermite3_to_poly(
    real_type   H,
    real_type   P0,
    real_type   P1,
    real_type   DP0,
    real_type   DP1,
    real_type & A,
    real_type & B,
    real_type & C,
    real_type & D )
  {
    const real_type invH  = real_type( 1 ) / H;
    const real_type slope = ( P1 - P0 ) * invH;  // Pendenza media (Secante)

    // Precalcolo inverso quadrato
    const real_type invH2 = invH * invH;

    // A = (Sum(derivs) - 2*slope) / H^2
    // Corrisponde a: (DP0 + DP1 - 2*(P1-P0)/H) / H^2
    A = ( ( DP0 + DP1 ) - 2 * slope ) * invH2;

    // B = (3*slope - (2*DP0 + DP1)) / H
    // Corrisponde a: (3*(P1-P0)/H - 2*DP0 - DP1) / H
    B = ( 3 * slope - ( 2 * DP0 + DP1 ) ) * invH;

    // C = Derivata iniziale
    C = DP0;

    // D = Posizione iniziale
    D = P0;
  }

  //!
  //! Convert polynomial defined using Hermite base
  //!
  //! \f[
  //! p(x) = p_0 h_0(t) + p_1 h_1(t) +
  //! p'_0 h_2(t) + p'_1 h_3(t) +
  //! p''_0 h_4(t) + p''_1 h_5(t)
  //! \f]
  //!
  //! to standard form
  //!
  //! \f[ p(x) = A t^5 + B t^4 + C t^3 + D t^2 + E t + F \f]
  //!
  //!
  static inline void Hermite5_to_poly(
    real_type   h,
    real_type   P0,
    real_type   P1,
    real_type   DP0,
    real_type   DP1,
    real_type   DDP0,
    real_type   DDP1,
    real_type & A,
    real_type & B,
    real_type & C,
    real_type & D,
    real_type & E,
    real_type & F )
  {
    // 1. Inverso di h (type-safe)
    const real_type invH = 1 / h;

    // 2. Precalcoli
    // Pendenza media (Velocità secante)
    const real_type slope = ( P1 - P0 ) * invH;

    // Accelerazioni dimezzate (moltiplicazione vs divisione)
    const real_type half_DDP0 = DDP0 / 2;
    const real_type half_DDP1 = DDP1 / 2;

    // 3. Potenze dell'inverso
    const real_type invH2 = invH * invH;
    const real_type invH3 = invH2 * invH;

    // 4. Calcolo Coefficienti

    // Grado 0, 1, 2
    F = P0;
    E = DP0;
    D = half_DDP0;

    // Grado 3 (C)
    // Correzione: Parentesi aggiunte per moltiplicare TUTTO per invH
    // Struttura: ( [Accelerazione] + [Velocità]*invH ) * invH
    // Risultato dimensionale: Accelerazione/H -> P/H^3
    C = ( ( half_DDP1 - 3 * half_DDP0 ) + invH * ( 10 * slope - ( 6 * DP0 + 4 * DP1 ) ) ) * invH;

    // Grado 4 (B)
    // Struttura: ( [Accelerazione] + [Velocità]*invH ) * invH^2
    // Risultato dimensionale: Accelerazione/H^2 -> P/H^4
    B = ( ( 3 * half_DDP0 - 2 * half_DDP1 ) + invH * ( ( 8 * DP0 + 7 * DP1 ) - 15 * slope ) ) * invH2;

    // Grado 5 (A)
    // Struttura: ( [Accelerazione] + [Velocità]*invH ) * invH^3
    // Risultato dimensionale: Accelerazione/H^3 -> P/H^5
    A = ( ( half_DDP1 - half_DDP0 ) + invH * ( 6 * slope - 3 * ( DP0 + DP1 ) ) ) * invH3;
  }

  // --------------------------------------------------------------------------

  static inline real_type bilinear3( real_type const p[4], real_type const M[4][4], real_type const q[4] )
  {
    return p[0] * ( M[0][0] * q[0] + M[0][1] * q[1] + M[0][2] * q[2] + M[0][3] * q[3] ) +
           p[1] * ( M[1][0] * q[0] + M[1][1] * q[1] + M[1][2] * q[2] + M[1][3] * q[3] ) +
           p[2] * ( M[2][0] * q[0] + M[2][1] * q[1] + M[2][2] * q[2] + M[2][3] * q[3] ) +
           p[3] * ( M[3][0] * q[0] + M[3][1] * q[1] + M[3][2] * q[2] + M[3][3] * q[3] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static inline real_type bilinear5( real_type const p[6], real_type const M[6][6], real_type const q[6] )
  {
    return p[0] *
             ( M[0][0] * q[0] + M[0][1] * q[1] + M[0][2] * q[2] + M[0][3] * q[3] + M[0][4] * q[4] + M[0][5] * q[5] ) +
           p[1] *
             ( M[1][0] * q[0] + M[1][1] * q[1] + M[1][2] * q[2] + M[1][3] * q[3] + M[1][4] * q[4] + M[1][5] * q[5] ) +
           p[2] *
             ( M[2][0] * q[0] + M[2][1] * q[1] + M[2][2] * q[2] + M[2][3] * q[3] + M[2][4] * q[4] + M[2][5] * q[5] ) +
           p[3] *
             ( M[3][0] * q[0] + M[3][1] * q[1] + M[3][2] * q[2] + M[3][3] * q[3] + M[3][4] * q[4] + M[3][5] * q[5] ) +
           p[4] *
             ( M[4][0] * q[0] + M[4][1] * q[1] + M[4][2] * q[2] + M[4][3] * q[3] + M[4][4] * q[4] + M[4][5] * q[5] ) +
           p[5] *
             ( M[5][0] * q[0] + M[5][1] * q[1] + M[5][2] * q[2] + M[5][3] * q[3] + M[5][4] * q[4] + M[5][5] * q[5] );
  }

}  // namespace Splines

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#include "Splines/SplineUtils.hxx"
#include "Splines/SplineParametrization.hxx"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

namespace Splines
{

  /*\
   |   ____        _ _
   |  / ___| _ __ | (_)_ __   ___
   |  \___ \| '_ \| | | '_ \ / _ \
   |   ___) | |_) | | | | | |  __/
   |  |____/| .__/|_|_|_| |_|\___|
   |        |_|
  \*/
  //!
  //! Spline Management Class
  //!
  class Spline
  {
    friend class SplineSet;

  protected:
    string const m_name;
    bool         m_curve_is_closed         = false;
    bool         m_curve_can_extend        = true;
    bool         m_curve_extended_constant = false;

    integer     m_npts          = 0;
    integer     m_npts_reserved = 0;
    real_type * m_X             = nullptr;  // allocated in the derived class!
    real_type * m_Y             = nullptr;  // allocated in the derived class!

    Utils::SearchInterval<real_type, integer> m_search;

  protected:
    void copy_flags( Spline const & S )
    {
      m_curve_is_closed         = S.m_curve_is_closed;
      m_curve_can_extend        = S.m_curve_can_extend;
      m_curve_extended_constant = S.m_curve_extended_constant;
    }

  public:
    Spline( Spline const & )                   = delete;
    Spline const & operator=( Spline const & ) = delete;

    //! \name Constructors
    ///@{

    //!
    //! spline constructor
    //!
    explicit Spline( string_view const name = "Spline" ) : m_name( name )
    {
      m_search.setup( &m_name, &m_npts, &m_X, &m_curve_is_closed, &m_curve_can_extend );
    }

    //!
    //! spline destructor
    //!
    virtual ~Spline() = default;

    ///@}

    //!
    //! \name Open/Close
    //!
    ///@{

    //!
    //! \return string with the name of the spline
    //!
    string_view name() const { return m_name; }

    //! \return `true` if spline is a closed spline
    bool is_closed() const { return m_curve_is_closed; }
    //!
    //! Set spline as a closed spline.
    //! When evaluated if parameter is outside the domain
    //! is wrapped cyclically before evalation.
    //!
    void make_closed() { m_curve_is_closed = true; }
    //!
    //! Set spline as an opened spline.
    //! When evaluated if parameter is outside the domain
    //! an error is produced.
    //!
    void make_opened() { m_curve_is_closed = false; }

    //!
    //! \return `true` if spline cannot extend outside interval of definition
    //!
    bool is_bounded() const { return !m_curve_can_extend; }
    //!
    //! Set spline as unbounded.
    //! When evaluated if parameter is outside the domain
    //! an extrapolated value is used.
    //!
    void make_unbounded() { m_curve_can_extend = true; }
    //!
    //! Set spline as bounded.
    //! When evaluated if parameter is outside the domain
    //! an error is issued.
    //!
    void make_bounded() { m_curve_can_extend = false; }

    //!
    //! \return `true` if the spline extend with a constant value
    //!
    bool is_extended_constant() const { return m_curve_extended_constant; }
    //!
    //! Set spline to extend constant.
    //! When evaluated if parameter is outside the domain
    //! the value returned is the value of the closed border.
    //!
    void make_extended_constant() { m_curve_extended_constant = true; }
    //!
    //! Set spline to extend NOT constant.
    //! When evaluated if parameter is outside the domain
    //! teh value returned is extrapolated using the last spline polynomial.
    //!
    void make_extended_not_constant() { m_curve_extended_constant = false; }

    ///@}

    //! \name Spline Data Info
    ///@{

    //!
    //! return true if spline is empty (no points)
    //!
    bool empty() const { return m_npts == 0; }

    //!
    //! the number of support points of the spline.
    //!
    integer num_points() const { return m_npts; }

    //!
    //! Return the pointer of values of x-nodes.
    //!
    real_type const * x_nodes() const { return m_X; }

    //!
    //! Return the pointer of values of y-nodes.
    //!
    real_type const * y_nodes() const { return m_Y; }

    //!
    //! the i-th node of the spline (x component).
    //!
    real_type x_node( integer const i ) const { return m_X[i]; }

    //!
    //! the i-th node of the spline (y component).
    //!
    real_type y_node( integer const i ) const { return m_Y[i]; }

    //!
    //! first node of the spline (x component).
    //!
    real_type x_begin() const { return m_X[0]; }

    //!
    //! first node of the spline (y component).
    //!
    real_type y_begin() const { return m_Y[0]; }

    //!
    //! last node of the spline (x component).
    //!
    real_type x_end() const { return m_X[m_npts - 1]; }

    //!
    //! last node of the spline (y component).
    //!
    real_type y_end() const { return m_Y[m_npts - 1]; }

    //!
    //! x-minumum spline value
    //!
    real_type x_min() const { return m_X[0]; }

    //!
    //! x-maximum spline value
    //!
    real_type x_max() const { return m_X[m_npts - 1]; }

    //!
    //! y-minumum spline value
    //!
    real_type y_min() const
    {
      integer N = ( type() == SplineType1D::CONSTANT ) ? m_npts - 1 : m_npts;
      return *std::min_element( m_Y, m_Y + N );
    }

    //!
    //! return y-maximum spline value
    //!
    real_type y_max() const
    {
      integer N = ( type() == SplineType1D::CONSTANT ) ? m_npts - 1 : m_npts;
      return *std::max_element( m_Y, m_Y + N );
    }

    //!
    //! Search the max and min values of `y` along the spline
    //! with the corresponding `x` position
    //!
    //! \param[out] i_min_pos interval where is the minimum
    //! \param[out] x_min_pos where is the minimum
    //! \param[out] y_min     the minimum value
    //! \param[out] i_max_pos interval where is the maximum
    //! \param[out] x_max_pos where is the maximum
    //! \param[out] y_max     the maximum value
    //!
    virtual void y_min_max(
      integer &   i_min_pos,
      real_type & x_min_pos,
      real_type & y_min,
      integer &   i_max_pos,
      real_type & x_max_pos,
      real_type & y_max ) const
    {
      i_min_pos = i_max_pos = 0;
      x_min_pos = y_min = x_max_pos = y_max = 0;
      UTILS_ERROR( "In spline: {} y_min_max not implemented\n", info() );
    }

    //!
    //! Search the max and min values of `y` along the spline
    //! with the corresponding `x` position
    //!
    //! \param[out] i_min_pos interval where is the minimum
    //! \param[out] x_min_pos where is the minimum
    //! \param[out] y_min     the minimum value
    //! \param[out] i_max_pos interval where is the maximum
    //! \param[out] x_max_pos where is the maximum
    //! \param[out] y_max     the maximum value
    //!
    virtual void y_min_max(
      vector<integer> &   i_min_pos,
      vector<real_type> & x_min_pos,
      vector<real_type> & y_min,
      vector<integer> &   i_max_pos,
      vector<real_type> & x_max_pos,
      vector<real_type> & y_max ) const
    {
      i_min_pos.clear();
      i_max_pos.clear();
      x_min_pos.clear();
      x_max_pos.clear();
      y_min.clear();
      y_max.clear();
      UTILS_ERROR( "In spline: {} y_min_max not implemented\n", info() );
    }
    ///@}

    //! \name Build
    ///@{

    //!
    //! Build a spline using data in `GenericContainer`
    //!
    void build( GenericContainer const & gc ) { setup( gc ); }

    //!
    //! Build a spline using data in a file `file_name`
    //!
    void build( string const & file_name ) { setup( file_name ); }

    //!
    //! Build a spline.
    //!
    //! \param x    vector of x-coordinates
    //! \param incx access elements as `x[0]`, `x[incx]`, `x[2*incx]`,...
    //! \param y    vector of y-coordinates
    //! \param incy access elements as `y[0]`, `y[incx]`, `y[2*incx]`,...
    //! \param n    total number of points
    //!
    virtual void build( real_type const x[], integer incx, real_type const y[], integer incy, integer const n )
    {
      reserve( n );
      for ( integer i = 0; i < n; ++i ) m_X[i] = x[i * incx];
      for ( integer i = 0; i < n; ++i ) m_Y[i] = y[i * incy];
      m_npts = n;
      build();
    }

    //!
    //! Build a spline.
    //!
    //! \param x vector of x-coordinates
    //! \param y vector of y-coordinates
    //! \param n total number of points
    //!
    inline void build( real_type const x[], real_type const y[], integer const n ) { this->build( x, 1, y, 1, n ); }

    //!
    //! Build a spline.
    //!
    //! \param x vector of x-coordinates
    //! \param y vector of y-coordinates
    //!
    void build( vector<real_type> const & x, vector<real_type> const & y )
    {
      integer N = std::min<integer>( x.size(), y.size() );
      this->build( x.data(), 1, y.data(), 1, N );
    }

    //!
    //! Build a spline using internal stored data
    //!
    virtual void build() = 0;

    //!
    //! Setup a spline using a `GenericContainer`
    //!
    //! - gc("xdata") vector with the `x` coordinate of the data
    //! - gc("ydata") vector with the `y` coordinate of the data
    //!
    virtual void setup( GenericContainer const & gc )
    {
      /*
      // gc["xdata"]
      // gc["ydata"]
      //
      */
      string const where = fmt::format( "Spline[{}]::setup( gc ):", m_name );

      std::set<std::string> keywords;
      for ( auto const & pair : gc.get_map( where ) ) { keywords.insert( pair.first ); }

      GenericContainer const & gc_x = gc( "xdata", where );
      keywords.erase( "xdata" );
      GenericContainer const & gc_y = gc( "ydata", where );
      keywords.erase( "ydata" );
      keywords.erase( "spline_type" );

      vec_real_type x, y;
      {
        string const ff = fmt::format( "{}, field `xdata'", where );
        gc_x.copyto_vec_real( x, ff );
      }
      {
        string const ff = fmt::format( "{}, field `ydata'", where );
        gc_y.copyto_vec_real( y, ff );
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
      build( x, y );
    }
    //!
    //! Setup a spline using a `GenericContainer` readed from file
    //!
    //! - file_name file name of the file with data
    //!
    //! - `.json`
    //! - `.yaml` or `.yml`
    //! - `.toml`
    //!
    void setup( string const & file_name )
    {
      GenericContainer gc;
      UTILS_ASSERT( gc.from_file( file_name ), "Spline::setup( '{}' ) failed to read\n", file_name );
      setup( gc );
    }

    ///@}

    //! \name Incremental Build
    ///@{

    //!
    //! Allocate memory for `npts` points
    //!
    virtual void reserve( integer const npts ) = 0;

    //!
    //! Add a support point (x,y) to the spline.
    //!
    void push_back( real_type const x, real_type const y )
    {
      if ( m_npts > 0 )
      {
        UTILS_ASSERT(
          x >= m_X[m_npts - 1],  // ammetto punti doppi
          "Spline[{}]::push_back, non monotone insert at insert N.{}"
          "\nX[{}] = {:.5}\nX[{}] = {:>5}\n",
          m_name,
          m_npts,
          m_npts - 1,
          m_X[m_npts - 1],
          m_npts,
          x );
      }
      if ( m_npts_reserved == 0 ) { reserve( 2 ); }
      else if ( m_npts >= m_npts_reserved )
      {
        // Versione "Fall-back" se devi per forza usare il buffer temporaneo
        integer const saved_npts = m_npts;
        Vec Xsaved( m_npts );
        Vec Ysaved( m_npts );

        // memcpy è più veloce di copy_n per array raw
        std::memcpy( Xsaved.data(), m_X, m_npts * sizeof( real_type ) );
        std::memcpy( Ysaved.data(), m_Y, m_npts * sizeof( real_type ) );

        reserve( ( m_npts + 1 ) * 2 );  // La tua funzione distruttiva

        m_npts = saved_npts;
        std::memcpy( m_X, Xsaved.data(), m_npts * sizeof( real_type ) );
        std::memcpy( m_Y, Ysaved.data(), m_npts * sizeof( real_type ) );
      }
      m_X[m_npts] = x;
      m_Y[m_npts] = y;
      ++m_npts;
      m_search.must_reset();
    }

    //!
    //! Drop last inserted point of the spline.
    //!
    void drop_back()
    {
      if ( m_npts > 0 ) --m_npts;
      m_search.must_reset();
    }

    //!
    //! Delete the support points, empty the spline.
    //!
    virtual void clear() = 0;

    ///@}

    //! \name Manipulate
    ///@{

    //!
    //! change X-origin of the spline
    //!
    void set_origin( real_type const x0 ) const
    {
      real_type const Tx = x0 - m_X[0];
      real_type *     ix = m_X;
      while ( ix < m_X + m_npts ) *ix++ += Tx;
    }

    //!
    //! change X-range of the spline
    //!
    void set_range( real_type xmin, real_type xmax )
    {
      UTILS_ASSERT( xmax > xmin, "Spline[{}]::set_range({},{}) bad range ", m_name, xmin, xmax );

      // Calcolo range attuale (delta)
      real_type const dx_old = m_X[m_npts - 1] - m_X[0];

      // Opzionale: Protezione divisione per zero se la spline è collassata su un punto
      // UTILS_ASSERT( std::abs(dx_old) > epsilon, ... );

      real_type const S  = ( xmax - xmin ) / dx_old;
      real_type const Tx = xmin - S * m_X[0];

      // EIGEN IMPLEMENTATION
      // Usiamo Eigen::Array invece di Matrix perché l'operazione è element-wise (scalare)
      // Map crea una vista sui dati esistenti senza copiare memoria.
      Eigen::Map<Vec> map_X{ m_X, m_npts };

      // Operazione in-place: x[i] = x[i] * S + Tx
      // Eigen utilizzerà istruzioni vettoriali (es. vfmadd su AVX2)
      map_X = map_X * S + Tx;
    }
    ///@}

    //! \name Dump Data
    ///@{

    //!
    //! dump a sample of the spline
    //!
    void dump( ostream_type & s, integer const nintervals, string_view header = "x\ty" ) const
    {
      s << header << '\n';
      real_type const dx = ( x_max() - x_min() ) / nintervals;
      for ( integer i = 0; i <= nintervals; ++i )
      {
        real_type x = x_min() + i * dx;
        fmt::print( s, "{}\t{}\n", x, this->eval( x ) );
      }
    }

    //!
    //! dump a sample of the spline
    //!
    void dump( string_view fname, integer const nintervals, string_view header = "x\ty" ) const
    {
      std::ofstream file( fname.data() );
      this->dump( file, nintervals, header );
      file.close();
    }

    //!
    //! Print spline coefficients
    //!
    virtual void write_to_stream( ostream_type & s ) const = 0;

    ///@}

    //! \name Evaluation
    ///@{

#ifdef AUTODIFF_SUPPORT
    //!
    //! Evaluate spline value
    //!
    virtual real_type         eval( real_type const x ) const           = 0;
    virtual autodiff::dual1st eval( autodiff::dual1st const & x ) const = 0;
    virtual autodiff::dual2nd eval( autodiff::dual2nd const & x ) const = 0;

    // Template unificato per tutti i tipi
    template <typename T> auto eval( T const & x ) const
    {
      if constexpr ( std::is_arithmetic<T>::value )
      {
        // Se T è un tipo numerico (int, float, double, etc.), promuovi a real_type
        return eval( static_cast<real_type>( x ) );
      }
      else
      {
        // Altrimenti deduce automaticamente il tipo duale appropriato
        return eval( autodiff::detail::to_dual( x ) );
      }
    }

    template <typename T> auto operator()( T const & x ) const -> decltype( eval( x ) ) { return eval( x ); }
#endif

    //!
    //! First derivative
    //!
    virtual real_type D( real_type const x ) const = 0;

    //!
    //! Second derivative
    //!
    virtual real_type DD( real_type const x ) const = 0;

    //!
    //! Third derivative
    //!
    virtual real_type DDD( real_type const x ) const = 0;

    //!
    //! 4th derivative
    //!
    virtual real_type DDDD( real_type const ) const { return real_type( 0 ); }

    //!
    //! 5th derivative
    //!
    virtual real_type DDDDD( real_type const ) const { return real_type( 0 ); }

    virtual void D( real_type const x, real_type dd[2] ) const  = 0;
    virtual void DD( real_type const x, real_type dd[3] ) const = 0;

    ///@}

    //!
    //! \name Evaluation Aliases
    //!
    ///@{
    //!
    //! Alias for `real_type eval( real_type x )`
    //!
    real_type operator()( real_type const x ) const { return this->eval( x ); }
    //!
    //! Alias for `real_type D( real_type x )`
    //!
    real_type eval_D( real_type const x ) const { return this->D( x ); }
    //!
    //! Alias for `real_type DD( real_type x )`
    //!
    real_type eval_DD( real_type const x ) const { return this->DD( x ); }
    //!
    //! Alias for `real_type DDD( real_type x )`
    //!
    real_type eval_DDD( real_type const x ) const { return this->DDD( x ); }
    //!
    //! Alias for `real_type DDDD( real_type x )`
    //!
    real_type eval_DDDD( real_type const x ) const { return this->DDDD( x ); }
    //!
    //! Alias for `real_type DDDD( real_type x )`
    //!
    real_type eval_DDDDD( real_type const x ) const { return this->DDDDD( x ); }
    ///@}

    //!
    //! \name Evaluation when segment is known
    //!
    ///@{

    //!
    //! Evaluate spline value
    //!
    virtual real_type id_eval( integer const ni, real_type const x ) const = 0;

    //!
    //! First derivative
    //!
    virtual real_type id_D( integer const ni, real_type const x ) const = 0;

    //!
    //! Second derivative
    //!
    virtual real_type id_DD( integer const ni, real_type const x ) const = 0;

    //!
    //! Third derivative
    //!
    virtual real_type id_DDD( integer const ni, real_type const x ) const = 0;

    //!
    //! 4th derivative
    //!
    virtual real_type id_DDDD( integer const, real_type const ) const { return real_type( 0 ); }

    //!
    //! 5th derivative
    //!
    virtual real_type id_DDDDD( integer const, real_type const ) const { return real_type( 0 ); }

    ///@}

    //! \name Get Info
    ///@{

    //!
    //! get the piecewise polinomials of the spline
    //!
    virtual integer  // order
    coeffs( real_type cfs[], real_type nodes[], bool transpose = false ) const = 0;

    //!
    //! Spline order of the piecewise polynomial
    //!
    virtual integer order() const = 0;

    //!
    //! spline type returned as a string
    //!
    char const * type_name() const { return to_string( type() ); }

    //!
    //! spline type returned as integer
    //!
    virtual SplineType1D type() const = 0;

    //!
    //! String information of the kind and order of the spline
    //!
    string info() const
    {
      string res = fmt::format( "Spline `{}` of type: {} of order: {}", m_name, type_name(), order() );
      if ( m_npts > 0 )
        res += fmt::format( "\nx_min={:.5} x_max={:.5} y_min={:.5} y_max={:.5}", x_min(), x_max(), y_min(), y_max() );
      return res;
    }

    //!
    //! Print information of the kind and order of the spline
    //!
    void info( ostream_type & stream ) const { stream << this->info() << '\n'; }

    ///@}

#ifdef SPLINES_BACK_COMPATIBILITY
    void pushBack( real_type x, real_type y ) { push_back( x, y ); }
    void dropBack() { drop_back(); }
    void setOrigin( real_type x0 ) { set_origin( x0 ); }
    void setRange( real_type xmin, real_type xmax ) { set_range( xmin, xmax ); }
    void writeToStream( ostream_type & s ) const { write_to_stream( s ); }
#endif
  };

  //!
  //! Compute curvature of a planar curve
  //! Formula: k = (v x a) / |v|^3
  //!
  static inline real_type curvature( real_type s, Spline const & X, Spline const & Y )
  {
    // Definizione di un vettore 2D a dimensione fissa (allocato sullo stack, zero overhead)
    using Vec2 = Eigen::Matrix<real_type, 2, 1>;

    // Helper function: Determinante 2D (analogo al prodotto vettoriale "cross" in 3D)
    // Restituisce v1.x * v2.y - v1.y * v2.x
    auto kross = []( const Vec2 & a, const Vec2 & b ) -> real_type { return a.x() * b.y() - a.y() * b.x(); };

    // Eigen caricherà questi valori direttamente nei registri CPU
    Vec2 v( X.D( s ), Y.D( s ) );    // Velocità
    Vec2 a( X.DD( s ), Y.DD( s ) );  // Accelerazione

    const real_type speed2 = v.squaredNorm();  // |v|^2

    // Protezione divisione per zero
    if ( speed2 <= 1e-12 ) return 0;

    // Numeratore: v x a
    const real_type num = kross( v, a );

    // Denominatore: |v|^3 = speed2 * sqrt(speed2)
    return num / ( speed2 * std::sqrt( speed2 ) );
  }

  //!
  //! Compute curvature derivative
  //! Formula vettoriale compatta: k' = ( |v|^2 (v x j) - 3 (v x a) (v . a) ) / |v|^5
  //!
  static inline real_type curvature_D( real_type s, Spline const & X, Spline const & Y )
  {
    // Definizione di un vettore 2D a dimensione fissa (allocato sullo stack, zero overhead)
    using Vec2 = Eigen::Matrix<real_type, 2, 1>;

    // Helper function: Determinante 2D (analogo al prodotto vettoriale "cross" in 3D)
    // Restituisce v1.x * v2.y - v1.y * v2.x
    auto kross = []( const Vec2 & a, const Vec2 & b ) -> real_type { return a.x() * b.y() - a.y() * b.x(); };

    Vec2 v( X.D( s ), Y.D( s ) );
    Vec2 a( X.DD( s ), Y.DD( s ) );
    Vec2 j( X.DDD( s ), Y.DDD( s ) );  // Jerk (Derivata terza)

    const real_type v2 = v.squaredNorm();  // |v|^2
    if ( v2 <= 1e-12 ) return 0;

    // Termini vettoriali comuni
    const real_type v_cross_a = kross( v, a );
    const real_type v_dot_a   = v.dot( a );
    const real_type v_cross_j = kross( v, j );

    // Applicazione formula derivata del quoziente vettoriale
    const real_type num = v2 * v_cross_j - 3 * v_cross_a * v_dot_a;

    // Denominatore: |v|^5
    return num / ( v2 * v2 * std::sqrt( v2 ) );
  }

  //!
  //! Compute curvature second derivative
  //! Sostituisce l'espansione polinomiale manuale con algebra vettoriale esatta
  //! Riduce drasticamente il numero di moltiplicazioni e migliora la precisione.
  //!
  static inline real_type curvature_DD( real_type s, Spline const & X, Spline const & Y )
  {
    // Definizione di un vettore 2D a dimensione fissa (allocato sullo stack, zero overhead)
    using Vec2 = Eigen::Matrix<real_type, 2, 1>;

    // Helper function: Determinante 2D (analogo al prodotto vettoriale "cross" in 3D)
    // Restituisce v1.x * v2.y - v1.y * v2.x
    auto kross = []( const Vec2 & a, const Vec2 & b ) -> real_type { return a.x() * b.y() - a.y() * b.x(); };

    // Recupero dati (idealmente la classe Spline dovrebbe avere una funzione eval_all(s, ...) per farlo in una
    // chiamata)
    Vec2 v( X.D( s ), Y.D( s ) );
    Vec2 a( X.DD( s ), Y.DD( s ) );
    Vec2 j( X.DDD( s ), Y.DDD( s ) );
    Vec2 s_snap( X.DDDD( s ), Y.DDDD( s ) );  // Snap/Jounce (Derivata quarta)

    const real_type v2 = v.squaredNorm();
    if ( v2 <= 1e-12 ) return 0;

    // --- Precalcoli Termini vettoriali (Cache nei registri) ---
    const real_type v_dot_a   = v.dot( a );
    const real_type v_cross_a = kross( v, a );
    const real_type v_cross_j = kross( v, j );

    // Termini misti derivata seconda
    // T2 = (a x j) + (v x s)
    const real_type T2 = kross( a, j ) + kross( v, s_snap );

    // D2 = |a|^2 + v . j
    const real_type D2 = a.squaredNorm() + v.dot( j );

    // --- Calcolo Numeratore ---
    // Formula derivata analiticamente:
    // N = v^2 [ v^2 * T2  -  6 * (v.a) * (v x j)  -  3 * (v x a) * D2 ]  +  15 * (v x a) * (v.a)^2

    const real_type term_bracket = v2 * T2 - 6 * v_dot_a * v_cross_j - 3 * v_cross_a * D2;
    const real_type num          = v2 * term_bracket + 15 * v_cross_a * ( v_dot_a * v_dot_a );

    // --- Risultato Finale ---
    // Denominatore: |v|^7
    return num / ( v2 * v2 * v2 * std::sqrt( v2 ) );
  }

}  // namespace Splines

#include "Splines/SplineCubicBase.hxx"
#include "Splines/SplineCubic.hxx"

#include "Splines/SplineBuild.hxx"

#include "Splines/SplineHermite.hxx"
#include "Splines/SplineAkima.hxx"
#include "Splines/SplineBessel.hxx"
#include "Splines/SplineConstant.hxx"
#include "Splines/SplineLinear.hxx"
#include "Splines/SplinePchip.hxx"
#include "Splines/SplineQuinticBase.hxx"
#include "Splines/SplineQuintic.hxx"

#include "Splines/SplineSurf.hxx"

#include "Splines/SplineBilinear.hxx"

#include "Splines/SplineBiCubicBase.hxx"
#include "Splines/SplineBiCubic.hxx"

#include "Splines/SplineBiQuinticBase.hxx"
#include "Splines/SplineBiQuintic.hxx"

#include "Splines/SplineVec.hxx"
#include "Splines/SplineSet.hxx"
#include "Splines/Splines1D.hxx"
#include "Splines/Splines2D.hxx"

#include "Splines/Splines1Dblend.hxx"
#include "Splines/Splines2Dblend.hxx"

namespace SplinesLoad
{

  using Splines::AkimaSpline;
  using Splines::BesselSpline;
  using Splines::ConstantSpline;
  using Splines::CubicSpline;
  using Splines::CubicSplineBase;
  using Splines::LinearSpline;
  using Splines::PchipSpline;
  using Splines::QuinticSpline;
  using Splines::Spline;
  using Splines::Spline1D;

  using Splines::BiCubicSpline;
  using Splines::BilinearSpline;
  using Splines::BiQuinticSpline;
  using Splines::Spline2D;

  using Splines::SplineSet;
  using Splines::SplineVec;

  using Splines::SplineType1D;
  using Splines::SplineType2D;
}  // namespace SplinesLoad

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif

#endif
