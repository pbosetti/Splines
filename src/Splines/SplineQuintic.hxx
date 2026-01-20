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

/*\
 |    ___        _       _   _      ____        _ _
 |   / _ \ _   _(_)_ __ | |_(_) ___/ ___| _ __ | (_)_ __   ___
 |  | | | | | | | | '_ \| __| |/ __\___ \| '_ \| | | '_ \ / _ \
 |  | |_| | |_| | | | | | |_| | (__ ___) | |_) | | | | | |  __/
 |   \__\_\\__,_|_|_| |_|\__|_|\___|____/| .__/|_|_|_| |_|\___|
 |                                       |_|
 |
\*/

namespace Splines
{

  using QuinticSpline_sub_type = enum class QuinticSpline_sub_type : integer {
    CUBIC  = 0,
    PCHIP  = 1,
    AKIMA  = 2,
    BESSEL = 3
  };

  inline std::string to_string( QuinticSpline_sub_type t )
  {
    switch ( t )
    {
      case QuinticSpline_sub_type::CUBIC: return "CUBIC";
      case QuinticSpline_sub_type::PCHIP: return "PCHIP";
      case QuinticSpline_sub_type::AKIMA: return "AKIMA";
      case QuinticSpline_sub_type::BESSEL: return "BESSEL";
    }
    return "NOTYPE";
  }

  // ============================================================================
  // OTTIMIZZAZIONE 1: Funzione helper inline per calcoli ripetuti
  // ============================================================================
  namespace detail {
    // Calcola termini comuni per evitare ricalcoli
    struct IntervalData {
      real_type h, h2, h3, inv_h, inv_h2, inv_h3;
      
      inline void compute(real_type h_val) {
        h = h_val;
        inv_h = 1.0 / h;
        h2 = h * h;
        inv_h2 = inv_h * inv_h;
        h3 = h2 * h;
        inv_h3 = inv_h2 * inv_h;
      }
    };
  }

  /*!
   * \brief Construct a complete quintic spline from cubic Hermite data
   */
  inline void Quintic_build(
    QuinticSpline_sub_type const q_sub_type,
    real_type const              X[],
    real_type const              Y[],
    real_type                    Yp[],
    real_type                    Ypp[],
    integer const                npts )
  {
    UTILS_ASSERT( npts >= 2, "Quintic_build, npts={} must be >= 2\n", npts );

    // ========================================================================
    // FASE 1: Stima derivate prime
    // ========================================================================
    switch ( q_sub_type )
    {
      case QuinticSpline_sub_type::CUBIC:
        CubicSpline_build( X, Y, Yp, Ypp, npts, CubicSpline_BC::EXTRAPOLATE, CubicSpline_BC::EXTRAPOLATE );
        break;
      case QuinticSpline_sub_type::PCHIP: 
        Pchip_build( X, Y, Yp, npts ); 
        break;
      case QuinticSpline_sub_type::AKIMA:
      {
        Malloc_real mem( "Quintic_build::work memory" );
        Akima_build( X, Y, Yp, mem.malloc( npts ), npts );
      }
      break;
      case QuinticSpline_sub_type::BESSEL: 
        Bessel_build( X, Y, Yp, npts ); 
        break;
      default: 
        UTILS_ERROR( "Unknown QuinticSpline_sub_type value\n" );
    }

    // ========================================================================
    // CASO BANALE: 2 punti
    // ========================================================================
    if ( npts == 2 )
    {
      Ypp[0] = Ypp[1] = 0;
      return;
    }

    integer const n = npts - 1;

    // ========================================================================
    // OTTIMIZZAZIONE 2: Allocazione unica per tutti i dati degli intervalli
    // Invece di allocare solo h, allochiamo tutto in un blocco contiguo
    // ========================================================================
    Malloc_real mem_intervals( "QuinticSpline::intervals" );
    real_type * interval_data = mem_intervals.malloc( 4 * n ); // h, h2, h3, inv_h
    
    real_type * h      = interval_data;
    real_type * h2     = interval_data + n;
    real_type * h3     = interval_data + 2 * n;
    real_type * inv_h  = interval_data + 3 * n;

    // Calcola tutti i termini degli intervalli in un'unica passata
    // Migliora cache locality
    for ( integer i = 0; i < n; ++i )
    {
      h[i] = X[i + 1] - X[i];
      UTILS_ASSERT( h[i] > 0, "X must be strictly increasing" );
      
      inv_h[i] = 1.0 / h[i];
      h2[i] = h[i] * h[i];
      h3[i] = h2[i] * h[i];
    }

    // ========================================================================
    // OTTIMIZZAZIONE 3: Sistema tridiagonale con allocazione ottimizzata
    // ========================================================================
    Utils::TridiagonalSolver<real_type> solver;
    solver.resize( npts );

    // Alloca vettori contigui per migliorare cache performance
    Malloc_real mem_system( "QuinticSpline::system" );
    real_type * system_data = mem_system.malloc( 3 * npts - 2 );
    
    // Map diretto su array allocati (zero overhead)
    Utils::TridiagonalSolver<real_type>::VecS a = 
      Eigen::Map<Eigen::Matrix<real_type, Eigen::Dynamic, 1>>( system_data, npts - 1 );
    Utils::TridiagonalSolver<real_type>::VecS b = 
      Eigen::Map<Eigen::Matrix<real_type, Eigen::Dynamic, 1>>( system_data + npts - 1, npts );
    Utils::TridiagonalSolver<real_type>::VecS c = 
      Eigen::Map<Eigen::Matrix<real_type, Eigen::Dynamic, 1>>( system_data + 2 * npts - 1, npts - 1 );
    
    // d mappa direttamente su Ypp (zero-copy!)
    Eigen::Map<Eigen::Matrix<real_type, Eigen::Dynamic, 1>> d( Ypp, npts );

    // ========================================================================
    // BORDO SINISTRO (i=0)
    // ========================================================================
    {
      real_type const h0_inv = inv_h[0];
      real_type const dY = Y[1] - Y[0];
      
      d( 0 ) = -10.0 * dY * h0_inv + 6.0 * Yp[0] + 4.0 * Yp[1];
      b( 0 ) = 2.0 * h[0];
      c( 0 ) = h[0];
    }

    // ========================================================================
    // EQUAZIONI INTERNE (i=1..n-1)
    // OTTIMIZZAZIONE 4: Riduzione operazioni e riuso variabili
    // ========================================================================
    for ( integer i = 1; i < n; ++i )
    {
      integer const iL = i - 1;
      
      real_type const hL_inv = inv_h[iL];
      real_type const hR_inv = inv_h[i];
      real_type const hL2_inv = hL_inv * hL_inv;
      real_type const hR2_inv = hR_inv * hR_inv;
      real_type const hL3_inv = hL2_inv * hL_inv;
      real_type const hR3_inv = hR2_inv * hR_inv;

      // Coefficienti matrice
      a( iL ) = -hL_inv;
      b( i )  = 3.0 * ( hL_inv + hR_inv );
      c( i )  = -hR_inv;

      // RHS: fattorizza calcoli comuni
      real_type const dY_L = Y[i] - Y[iL];
      real_type const dY_R = Y[i + 1] - Y[i];
      
      real_type const termY = 20.0 * ( dY_R * hR3_inv - dY_L * hL3_inv );
      
      // Calcola termini derivate con meno operazioni
      real_type const termYp_R = ( 3.0 * Yp[i] + 2.0 * Yp[i + 1] ) * hR2_inv;
      real_type const termYp_L = ( 2.0 * Yp[iL] + 3.0 * Yp[i] ) * hL2_inv;

      d( i ) = termY - 4.0 * ( termYp_R - termYp_L );
    }

    // ========================================================================
    // BORDO DESTRO (i=n)
    // ========================================================================
    {
      real_type const hn_inv = inv_h[n - 1];
      real_type const dY = Y[n] - Y[n - 1];
      
      a( n - 1 ) = h[n - 1];
      b( n )     = 2.0 * h[n - 1];
      d( n )     = 10.0 * dY * hn_inv - 4.0 * Yp[n - 1] - 6.0 * Yp[n];
    }

    // ========================================================================
    // RISOLUZIONE: solve_inplace modifica d che è già mappato su Ypp
    // ========================================================================
    solver.factorize( a, b, c );
    solver.solve_inplace( a, b, d );
  }

  //! Quintic spline class
  class QuinticSpline : public QuinticSplineBase
  {
    QuinticSpline_sub_type m_q_sub_type{ QuinticSpline_sub_type::CUBIC };

  public:
    //!
    //! \name Constructors
    //!
    ///@{

    using QuinticSplineBase::build;
    using QuinticSplineBase::reserve;

    //!
    //! Build an empty spline of `QuinticSpline` type
    //!
    //! \param name the name of the spline
    //!
    explicit QuinticSpline( string_view name = "Spline" ) : QuinticSplineBase( name ) {}

    //!
    //! Spline destructor.
    //!
    ~QuinticSpline() override = default;

    ///@}

    //!
    //! Set spline type
    //!
    //! - CUBIC
    //! - PCHIP
    //! - AKIMA
    //! - BESSEL
    //!
    void set_quintic_type( QuinticSpline_sub_type qt ) { m_q_sub_type = qt; }

    void setQuinticType( QuinticSpline_sub_type qt ) { m_q_sub_type = qt; }

    // --------------------------- VIRTUALS -----------------------------------
    //! Build a Monotone quintic spline from previously inserted points
    void build() override
    {
      string msg{ fmt::format( "QuinticSpline[{}]::build():", m_name ) };
      UTILS_ASSERT( m_npts > 1, "{} npts = {} not enought points\n", msg, m_npts );
      Utils::check_NaN( m_X, msg + " X", m_npts, __LINE__, __FILE__ );
      Utils::check_NaN( m_Y, msg + " Y", m_npts, __LINE__, __FILE__ );
      
      integer ibegin{ 0 };
      integer iend{ 0 };
      do
      {
        // cerca intervallo monotono strettamente crescente
        for ( ++iend; iend < m_npts && m_X[iend - 1] < m_X[iend]; ++iend ) {}
        Quintic_build( m_q_sub_type, m_X + ibegin, m_Y + ibegin, m_Yp + ibegin, m_Ypp + ibegin, iend - ibegin );
        ibegin = iend;
      } while ( iend < m_npts );

      Utils::check_NaN( m_Yp, msg + " Yp", m_npts, __LINE__, __FILE__ );
      Utils::check_NaN( m_Ypp, msg + " Ypp", m_npts, __LINE__, __FILE__ );
      m_search.must_reset();
    }

    //! Build a Monotone quintic spline from data from `gc`
    void setup( GenericContainer const & gc ) override
    {
      string const where{ fmt::format( "QuinticSpline[{}]::setup( gc ):", m_name ) };

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
      
      if ( gc.exists( "spline_sub_type" ) )
      {
        string_view st{ gc.get_map_string( "spline_sub_type", where ) };
        keywords.erase( "spline_sub_type" );
        
        // OTTIMIZZAZIONE 5: Switch più efficiente del if-else chain
        switch ( st[0] ) // First character hint
        {
          case 'c': // cubic
            if ( st == "cubic" ) m_q_sub_type = QuinticSpline_sub_type::CUBIC;
            else goto unknown_type;
            break;
          case 'p': // pchip
            if ( st == "pchip" ) m_q_sub_type = QuinticSpline_sub_type::PCHIP;
            else goto unknown_type;
            break;
          case 'a': // akima
            if ( st == "akima" ) m_q_sub_type = QuinticSpline_sub_type::AKIMA;
            else goto unknown_type;
            break;
          case 'b': // bessel
            if ( st == "bessel" ) m_q_sub_type = QuinticSpline_sub_type::BESSEL;
            else goto unknown_type;
            break;
          default:
          unknown_type:
            UTILS_ERROR( "{} unknow sub type: {}\n", where, st );
        }
      }
      else
      {
        UTILS_WARNING( false, "{}, missing field `spline_sub_type` using `cubic` as default value\n", where );
      }

      UTILS_WARNING(
        keywords.empty(),
        "{}: unused keys\n{}\n",
        where,
        [&keywords]() -> string
        {
          string res;
          res.reserve( keywords.size() * 10 ); // OTTIMIZZAZIONE 6: Pre-allocazione
          for ( auto const & it : keywords )
          {
            res += it;
            res += ' ';
          }
          return res;
        }() );

      this->build( x, y );
    }
  };

}  // namespace Splines

// EOF: SplineQuintic.hxx
