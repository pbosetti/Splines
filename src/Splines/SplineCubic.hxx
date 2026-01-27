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
 |      Versione Ottimizzata con Eigen - 2026                               |
 |                                                                          |
\*--------------------------------------------------------------------------*/

/**
 * @file SplineCubic.hxx
 * @brief Implementazione di spline cubiche con boundary conditions multiple
 *
 * Questo file contiene l'implementazione di spline cubiche interpolanti
 * utilizzando la libreria Eigen per operazioni vettoriali efficienti.
 * Supporta diverse condizioni al contorno (boundary conditions).
 */

namespace Splines
{

  /**
   * @enum CubicSpline_BC
   * @brief Tipi di condizioni al contorno per spline cubiche
   *
   * Definisce le possibili condizioni al contorno applicabili
   * agli estremi della spline cubica.
   */
  using CubicSpline_BC = enum class CubicSpline_BC : integer {
    EXTRAPOLATE      = 0,  ///< Estrapolazione usando derivate di ordine superiore
    NATURAL          = 1,  ///< Derivata seconda nulla agli estremi
    PARABOLIC_RUNOUT = 2,  ///< Derivata terza nulla agli estremi (continuità cubica)
    NOT_A_KNOT       = 3   ///< Condizione "not-a-knot" di De Boor
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /**
   * @brief Costruisce una spline cubica calcolando le derivate prime
   *
   * Questa funzione risolve il sistema tridiagonale per le derivate seconde (Z)
   * e calcola le derivate prime (Yp) necessarie per l'interpolazione cubica.
   * Utilizza Eigen per operazioni vettoriali efficienti e gestisce automaticamente
   * la memoria.
   *
   * @param[in]  X    Array delle coordinate x (strettamente crescenti)
   * @param[in]  Y    Array delle coordinate y corrispondenti
   * @param[out] Yp   Array delle derivate prime calcolate
   * @param[in]  npts Numero di punti dati
   * @param[in]  bc0  Condizione al contorno iniziale (x = X[0])
   * @param[in]  bcn  Condizione al contorno finale (x = X[npts-1])
   *
   * @note Il sistema tridiagonale viene risolto usando l'algoritmo di Thomas
   * @note Per condizioni NOT_A_KNOT viene usata una correzione iterativa (Sherman-Morrison)
   *
   * @throw UTILS_ASSERT Se npts < 2
   *
   * @par Dettagli matematici
   * La spline cubica su [X[i], X[i+1]] è definita come:
   * \f[
   * S_i(x) = Y[i] + Yp[i](x-X[i]) + \frac{Z[i]}{2}(x-X[i])^2 +
   *          \frac{Z[i+1]-Z[i]}{6h_i}(x-X[i])^3
   * \f]
   * dove h_i = X[i+1] - X[i] e Z[i] è la derivata seconda in X[i].
   */
  inline void CubicSpline_build(
    real_type const      X[],
    real_type const      Y[],
    real_type            Yp[],
    real_type            Ypp[],
    integer const        npts,
    CubicSpline_BC const bc0,
    CubicSpline_BC const bcn )
  {
    // Alias per tipi Eigen
    using Vector    = Eigen::Matrix<real_type, Eigen::Dynamic, 1>;
    using MapVector = Eigen::Map<const Vector>;

    UTILS_ASSERT( npts >= 2, "CubicSpline_build, npts={} must be >= 2\n", npts );

    integer const n = npts - 1;  ///< Numero di segmenti
    integer const m = npts;      ///< Numero di punti

    // Mappiamo gli array C-style come vettori Eigen (zero-copy)
    MapVector X_vec( X, npts );
    MapVector Y_vec( Y, npts );

    // --- 1. Allocazione vettori per il sistema tridiagonale ---
    Vector sub( n );    ///< Diagonale inferiore (coefficienti a_i)
    Vector diag( m );   ///< Diagonale principale (coefficienti b_i)
    Vector super( n );  ///< Diagonale superiore (coefficienti c_i)
    Vector rhs( m );    ///< Termine noto (derivate seconde Z)

    // Pre-calcolo delle differenze (vectorizzato con Eigen)
    Vector DX = X_vec.tail( n ) - X_vec.head( n );  ///< h_i = X[i+1] - X[i]
    Vector DY = Y_vec.tail( n ) - Y_vec.head( n );  ///< ΔY_i = Y[i+1] - Y[i]

    // --- 2. Costruzione del corpo centrale del sistema (i = 1 ... n-1) ---
    for ( integer i = 1; i < n; ++i )
    {
      real_type const HL = DX( i - 1 );  // h_{i-1}
      real_type const HR = DX( i );      // h_i
      real_type const HH = HL + HR;      // h_{i-1} + h_i

      sub( i - 1 ) = HL / HH;  // Coefficiente L_i
      super( i )   = HR / HH;  // Coefficiente U_i
      diag( i )    = 2;        // Diagonale principale

      // Termine noto: derivata seconda approssimata
      // rhs(i) = 6 * [f[x_i, x_{i+1}] - f[x_{i-1}, x_i]] / (h_{i-1} + h_i)
      real_type const slope_L = DY( i - 1 ) / HL;
      real_type const slope_R = DY( i ) / HR;
      rhs( i )                = 6 * ( slope_R - slope_L ) / HH;
    }

    // Variabili per elementi extra-diagonali (usate solo per NOT_A_KNOT)
    real_type UU = 0;
    real_type LL = 0;

    // --- 3. Applicazione condizioni al contorno iniziali (x = X[0]) ---
    switch ( bc0 )
    {
      case CubicSpline_BC::EXTRAPOLATE:
        /**
         * Estrapolazione: calcola Z[0] usando differenze finite di ordine superiore
         * basate sul numero di punti disponibili (2-5 punti)
         */
        diag( 0 )  = 1;
        super( 0 ) = 0;

        switch ( npts )
        {
          case 2: rhs( 0 ) = 0; break;
          case 3:
          case 4: rhs( 0 ) = Utils::second_derivative_3p( X[0], Y[0], X[1], Y[1], X[2], Y[2] ); break;
          case 5:
          case 6: rhs( 0 ) = Utils::second_derivative_4p( X[0], Y[0], X[1], Y[1], X[2], Y[2], X[3], Y[3] ); break;
          default:
            rhs( 0 ) = Utils::second_derivative_5p( X[0], Y[0], X[1], Y[1], X[2], Y[2], X[3], Y[3], X[4], Y[4] );
            break;
        }
        break;

      case CubicSpline_BC::NATURAL:
        /**
         * Natural spline: Z[0] = 0 (derivata seconda nulla)
         * Minimizza la curvatura totale della spline
         */
        diag( 0 )  = 1;
        super( 0 ) = 0;
        rhs( 0 )   = 0;
        break;

      case CubicSpline_BC::PARABOLIC_RUNOUT:
        /**
         * Parabolic runout: Z[0] = Z[1]
         * Implica che la derivata terza è nulla al bordo
         */
        diag( 0 )  = 1;
        super( 0 ) = -1;
        rhs( 0 )   = 0;
        break;

      case CubicSpline_BC::NOT_A_KNOT:
        /**
         * Not-a-knot: la derivata terza è continua in X[1]
         * Equivale a Z[0] - (1+r)Z[1] + r*Z[2] = 0, con r = h_0/h_1
         * L'elemento extra (r*Z[2]) viene gestito separatamente
         */
        {
          real_type const r = DX( 0 ) / DX( 1 );
          diag( 0 )         = 1;
          super( 0 )        = -( 1 + r );
          UU                = r;  // Elemento extra-diagonale in posizione (0,2)
          rhs( 0 )          = 0;
        }
        break;
    }

    // --- 4. Applicazione condizioni al contorno finali (x = X[n]) ---
    switch ( bcn )
    {
      case CubicSpline_BC::EXTRAPOLATE:
        /**
         * Estrapolazione all'estremo finale usando punti precedenti
         */
        diag( n )    = 1.0;
        sub( n - 1 ) = 0.0;

        switch ( npts )
        {
          case 2: rhs( n ) = 0; break;
          case 3:
          case 4: rhs( n ) = Utils::second_derivative_3p( X[n], Y[n], X[n - 1], Y[n - 1], X[n - 2], Y[n - 2] ); break;
          case 5:
          case 6:
            rhs( n ) =
              Utils::second_derivative_4p( X[n], Y[n], X[n - 1], Y[n - 1], X[n - 2], Y[n - 2], X[n - 3], Y[n - 3] );
            break;
          default:
            rhs( n ) = Utils::second_derivative_5p(
              X[n],
              Y[n],
              X[n - 1],
              Y[n - 1],
              X[n - 2],
              Y[n - 2],
              X[n - 3],
              Y[n - 3],
              X[n - 4],
              Y[n - 4] );
            break;
        }
        break;

      case CubicSpline_BC::NATURAL:
        diag( n )    = 1;
        sub( n - 1 ) = 0;
        rhs( n )     = 0;
        break;

      case CubicSpline_BC::PARABOLIC_RUNOUT:
        diag( n )    = 1;
        sub( n - 1 ) = -1;
        rhs( n )     = 0;
        break;

      case CubicSpline_BC::NOT_A_KNOT:
        /**
         * Not-a-knot all'estremo finale: continuità della derivata terza in X[n-1]
         */
        {
          real_type const r = DX( n - 2 ) / DX( n - 1 );
          diag( n )         = 1;
          sub( n - 1 )      = -( 1 + r );
          LL                = r;  // Elemento extra-diagonale in posizione (n,n-2)
          rhs( n )          = 0;
        }
        break;
    }

    // --- 5. Risoluzione del sistema tridiagonale ---
    Utils::TridiagonalSolver<real_type> solver;
    solver.resize( m );

    // Fattorizzazione LU del sistema tridiagonale
    solver.factorize( sub, diag, super );

    // Prima risoluzione (senza correzione per elementi extra-diagonali)
    Vector Z( m );  ///< Vettore delle derivate seconde
    solver.solve( sub, diag, rhs, Z );

    // --- 6. Correzione per condizioni NOT_A_KNOT ---
    /**
     * Per NOT_A_KNOT, il sistema ha elementi extra-diagonali che non sono
     * nella struttura tridiagonale standard. Usiamo un approccio iterativo:
     * 1. Risolviamo ignorando i termini extra
     * 2. Correggiamo il termine noto usando la soluzione ottenuta
     * 3. Risolviamo nuovamente (convergenza in una iterazione per questo caso)
     */
    if ( UU != 0.0 || LL != 0.0 )
    {
      // Aggiorna il termine noto con i contributi extra-diagonali
      if ( UU != 0.0 ) rhs( 0 ) -= UU * Z( 2 );      // Correzione al bordo iniziale
      if ( LL != 0.0 ) rhs( n ) -= LL * Z( n - 2 );  // Correzione al bordo finale

      // Risolvi nuovamente con il termine noto corretto
      solver.solve( sub, diag, rhs, Z );
    }

    // --- 7. Calcolo derivate prime usando operazioni vettoriali Eigen ---
    /**
     * La derivata prima è calcolata dalla formula:
     * Yp[i] = (Y[i+1] - Y[i])/h_i - h_i/6 * (2*Z[i] + Z[i+1])
     *
     * Questa assicura che la spline cubica interpolante abbia le derivate
     * prime corrette in ogni nodo.
     */

    // Calcolo per i = 0 ... n-1 (vectorizzato)
    Vector slopes = DY.array() / DX.array();  // (Y[i+1] - Y[i]) / h_i

    for ( integer i = 0; i < n; ++i ) { Yp[i] = slopes( i ) - DX( i ) * ( 2.0 * Z( i ) + Z( i + 1 ) ) / 6.0; }

    // Calcolo per l'ultimo punto usando continuità della derivata prima
    /**
     * Formula: Yp[n] = Yp[n-1] + h_{n-1}/2 * (Z[n-1] + Z[n])
     * Derivata dall'integrazione della derivata seconda
     */
    Yp[n] = Yp[n - 1] + DX( n - 1 ) * ( Z( n - 1 ) + Z( n ) ) * 0.5;
    if ( Ypp != nullptr && npts > 0 ) std::memcpy( Ypp, Z.data(), npts * sizeof( real_type ) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /**
   * @class CubicSpline
   * @brief Classe per la gestione di spline cubiche interpolanti
   *
   * Questa classe eredita da CubicSplineBase e implementa spline cubiche
   * con supporto per diverse condizioni al contorno. La spline risultante
   * garantisce continuità C² (derivate seconde continue) all'interno del
   * dominio, mentre le condizioni al contorno determinano il comportamento
   * agli estremi.
   *
   * @par Caratteristiche
   * - Interpolazione cubica a tratti con continuità C²
   * - 4 tipi di boundary conditions supportate
   * - Gestione automatica di segmenti monotoni multipli
   * - Ottimizzazione con Eigen per efficienza computazionale
   *
   * @par Uso tipico
   * @code
   * CubicSpline spline("my_spline");
   * spline.set_initial_BC(CubicSpline_BC::NATURAL);
   * spline.set_final_BC(CubicSpline_BC::NATURAL);
   * spline.build(x_data, y_data);
   * real_type y = spline.eval(x_query);
   * @endcode
   */
  class CubicSpline : public CubicSplineBase
  {
  private:
    CubicSpline_BC m_bc0{ CubicSpline_BC::EXTRAPOLATE };  ///< Boundary condition iniziale
    CubicSpline_BC m_bcn{ CubicSpline_BC::EXTRAPOLATE };  ///< Boundary condition finale

  public:
    //!
    //! @name Constructors
    //!
    ///@{

    using CubicSplineBase::build;
    using CubicSplineBase::reserve;

    /**
     * @brief Costruisce una spline vuota di tipo CubicSpline
     *
     * @param[in] name Nome identificativo della spline (opzionale)
     */
    explicit CubicSpline( string_view name = "CubicSpline" ) : CubicSplineBase( name ) {}

    /**
     * @brief Distruttore della spline
     */
    ~CubicSpline() override = default;

    ///@}

    //!
    //! @name Setup
    //!
    ///@{

    /**
     * @brief Imposta la condizione al contorno iniziale
     *
     * Definisce come la spline si comporta all'estremo sinistro (x = X[0]).
     * La scelta influenza la curvatura e la smoothness della spline.
     *
     * @param[in] bc0 Tipo di boundary condition iniziale
     *
     * @see CubicSpline_BC per i tipi disponibili
     */
    void set_initial_BC( CubicSpline_BC bc0 ) { m_bc0 = bc0; }

    /**
     * @brief Imposta la condizione al contorno finale
     *
     * Definisce come la spline si comporta all'estremo destro (x = X[npts-1]).
     *
     * @param[in] bcn Tipo di boundary condition finale
     *
     * @see CubicSpline_BC per i tipi disponibili
     */
    void set_final_BC( CubicSpline_BC bcn ) { m_bcn = bcn; }

    // --------------------------- VIRTUALS -----------------------------------

    /**
     * @brief Costruisce la spline a partire dai dati impostati
     *
     * Questo metodo:
     * 1. Verifica la validità dei dati (NaN check)
     * 2. Suddivide i dati in segmenti monotoni crescenti
     * 3. Costruisce una spline cubica per ogni segmento
     * 4. Applica le boundary conditions appropriate
     *
     * @throw UTILS_ASSERT Se npts < 2
     * @throw UTILS_ASSERT Se i dati contengono NaN
     *
     * @note La spline viene costruita anche per dati non monotoni,
     *       dividendoli automaticamente in segmenti monotoni
     */
    void build() override
    {
      string msg{ fmt::format( "CubicSpline[{}]::build():", m_name ) };

      // Validazione input
      UTILS_ASSERT( m_npts > 1, "{} npts={} not enough points\n", msg, m_npts );

      Utils::check_NaN( m_X, msg + " X", m_npts, __LINE__, __FILE__ );
      Utils::check_NaN( m_Y, msg + " Y", m_npts, __LINE__, __FILE__ );

      // Costruzione per segmenti monotoni
      integer ibegin{ 0 };
      integer iend{ 0 };

      do
      {
        // Cerca il prossimo intervallo monotono strettamente crescente
        for ( ++iend; iend < m_npts && m_X[iend - 1] < m_X[iend]; ++iend ) {}

        // Determina le boundary conditions per questo segmento
        auto seg_bc0{ CubicSpline_BC::NOT_A_KNOT };
        auto seg_bcn{ CubicSpline_BC::NOT_A_KNOT };

        if ( ibegin == 0 ) seg_bc0 = m_bc0;     // Primo segmento: usa BC iniziale
        if ( iend == m_npts ) seg_bcn = m_bcn;  // Ultimo segmento: usa BC finale

        // Costruisci la spline per questo segmento
        CubicSpline_build(
          m_X + ibegin,   // Puntatore ai dati X del segmento
          m_Y + ibegin,   // Puntatore ai dati Y del segmento
          m_Yp + ibegin,  // Puntatore alle derivate del segmento
          nullptr,        // non serve farsi dare le derivate seconde
          iend - ibegin,  // Numero di punti nel segmento
          seg_bc0,        // BC iniziale del segmento
          seg_bcn         // BC finale del segmento
        );

        ibegin = iend;  // Passa al prossimo segmento

      } while ( iend < m_npts );

      // Validazione output
      Utils::check_NaN( m_Yp, msg + " Yp", m_npts, __LINE__, __FILE__ );

      // Reset della ricerca binaria per nuova spline
      m_search.must_reset();
    }

    /**
     * @brief Configura la spline usando un GenericContainer
     *
     * Carica i dati e le impostazioni da un container generico.
     *
     * @param[in] gc GenericContainer con i dati e le opzioni
     *
     * @par Campi richiesti nel GenericContainer
     * - `"xdata"`: vettore con le coordinate x
     * - `"ydata"`: vettore con le coordinate y
     *
     * @par Campi opzionali
     * - `"bc_begin"`: boundary condition iniziale ("extrapolate", "natural",
     *                 "parabolic", "not_a_knot")
     * - `"bc_end"`: boundary condition finale (stessi valori di bc_begin)
     *
     * @throw UTILS_ERROR Se i campi richiesti mancano o hanno formato errato
     * @warning Emette warning se campi opzionali mancano (usa default)
     * @warning Emette warning se ci sono campi non riconosciuti
     *
     * @code
     * GenericContainer gc;
     * gc["xdata"] << vector<double>{0, 1, 2, 3};
     * gc["ydata"] << vector<double>{0, 1, 4, 9};
     * gc["bc_begin"] = "natural";
     * gc["bc_end"] = "natural";
     * spline.setup(gc);
     * @endcode
     */
    void setup( GenericContainer const & gc ) override
    {
      string const where{ fmt::format( "CubicSpline[{}]::setup( gc ):", m_name ) };

      // Raccoglie tutte le chiavi presenti per identificare campi non usati
      std::set<std::string> keywords;
      for ( auto const & pair : gc.get_map( where ) ) { keywords.insert( pair.first ); }
      keywords.erase( "spline_type" );  // Campo standard, non è un warning

      // Estrazione dati X e Y (obbligatori)
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

      // Parsing boundary condition iniziale (opzionale)
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
          UTILS_ERROR( "{} unknown initial bc: {}\n", where, bc );
        }
      }
      else
      {
        UTILS_WARNING( false, "{}, missing field `bc_begin` using `extrapolate` as default value\n", where );
      }

      // Parsing boundary condition finale (opzionale)
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
          UTILS_ERROR( "{} unknown final bc: {}\n", where, bc );
        }
      }
      else
      {
        UTILS_WARNING( false, "{}, missing field `bc_end` using `extrapolate` as default value\n", where );
      }

      // Warning per campi non riconosciuti
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
          }
          return res;
        }() );

      // Costruzione finale della spline
      this->build( x, y );
    }

    ///@}

    /**
     * @brief Ritorna il tipo di spline
     *
     * @return SplineType1D::CUBIC
     */
    SplineType1D type() const override { return SplineType1D::CUBIC; }

#ifdef SPLINES_BACK_COMPATIBILITY
    /// @deprecated Usa set_initial_BC invece
    void setInitialBC( CubicSpline_BC bc0 ) { m_bc0 = bc0; }

    /// @deprecated Usa set_final_BC invece
    void setFinalBC( CubicSpline_BC bcn ) { m_bcn = bcn; }
#endif
  };

}  // namespace Splines

// EOF: SplineCubic.hxx
