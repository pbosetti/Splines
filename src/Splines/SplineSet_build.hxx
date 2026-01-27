///////////////////////////////////////////////////////////////////////////
/**
 * \brief Build the spline set from raw data arrays
 *
 * Constructs multiple splines sharing the same x-nodes. Each spline can be
 * of a different type (linear, cubic, etc.) and may have optional derivative data.
 *
 * \param[in] nspl    Number of splines
 * \param[in] npts    Number of points per spline
 * \param[in] headers Array of spline names (size = nspl)
 * \param[in] stype   Array of spline types (size = nspl)
 * \param[in] X       Array of x-node values (size = npts)
 * \param[in] Y       Array of pointers to y-value arrays (size = nspl)
 * \param[in] Yp      Optional array of pointers to derivative arrays (size = nspl, default: nullptr)
 *
 * \par Memory Management:
 * Allocates internal storage and copies input data. Input arrays can be freed
 * after calling build().
 *
 * \par Supported Spline Types:
 * - CONSTANT: Piecewise constant
 * - LINEAR: Linear interpolation
 * - CUBIC: Cubic spline (natural)
 * - AKIMA: Akima spline
 * - BESSEL: Bessel spline
 * - PCHIP: Piecewise Cubic Hermite Interpolating Polynomial
 * - HERMITE: Hermite spline (requires Yp)
 * - QUINTIC: Quintic spline
 *
 * \throw UTILS_ASSERT if nspl ≤ 0 or npts ≤ 1
 * \throw UTILS_ERROR for unsupported spline types
 *
 * \note For HERMITE splines, Yp must be provided.
 * \note Automatically checks monotonicity for applicable spline types.
 */
void build(
  integer const           nspl,
  integer const           npts,
  char const * const      headers[],
  SplineType1D const      stype[],
  real_type const         X[],
  real_type const * const Y[],
  real_type const * const Yp[] = nullptr )
{
  string const msg = fmt::format( "SplineSet[{}]::build(...):", m_name );

  // Validazione input
  UTILS_ASSERT( nspl > 0, "{} expected positive nspl = {}\n", msg, nspl );
  UTILS_ASSERT( npts > 1, "{} expected npts = {} greater than 1", msg, npts );
  UTILS_ASSERT( headers != nullptr, "{} headers array is null\n", msg );
  UTILS_ASSERT( stype != nullptr, "{} stype array is null\n", msg );
  UTILS_ASSERT( X != nullptr, "{} X array is null\n", msg );
  UTILS_ASSERT( Y != nullptr, "{} Y array is null\n", msg );

  m_nspl = nspl;
  m_npts = npts;

  // Allocazione memoria per le strutture principali
  m_splines.resize( m_nspl );
  m_is_monotone = m_mem_int.realloc( m_nspl );
  m_header_to_position.clear();

  // ===================================================================
  // FASE 1: Calcolo della memoria totale necessaria
  // ===================================================================
  integer mem = npts;  // Spazio base per X

  for ( integer spl = 0; spl < nspl; ++spl )
  {
    UTILS_ASSERT( headers[spl] != nullptr, "{} headers[{}] array is null\n", msg, spl );

    switch ( stype[spl] )
    {
      // Spline quintiche: necessitano Y, Yp, Ypp (3 * npts)
      case SplineType1D::QUINTIC_CUBIC:
      case SplineType1D::QUINTIC_AKIMA:
      case SplineType1D::QUINTIC_BESSEL:
      case SplineType1D::QUINTIC_PCHIP:
        mem += npts;  // Spazio per Ypp
        [[fallthrough]];

      // Spline cubiche/hermitiane: necessitano Y, Yp (2 * npts)
      case SplineType1D::CUBIC:
      case SplineType1D::AKIMA:
      case SplineType1D::BESSEL:
      case SplineType1D::PCHIP:
      case SplineType1D::HERMITE:
        mem += npts;  // Spazio per Yp
        [[fallthrough]];

      // Spline lineari/costanti: necessitano solo Y (npts)
      case SplineType1D::CONSTANT:
      case SplineType1D::LINEAR:
        mem += npts;  // Spazio per Y
        break;

      // Tipi non supportati
      case SplineType1D::SPLINE_SET:
      case SplineType1D::SPLINE_VEC:
        UTILS_ERROR(
          "{} At spline n.{} named {} cannot be done for type = {}\n",
          msg,
          spl,
          headers[spl],
          to_string( stype[spl] ) );
    }
  }

  // ===================================================================
  // FASE 2: Allocazione memoria
  // ===================================================================
  m_mem.reallocate( mem + 2 * nspl );
  m_mem_p.reallocate( 3 * nspl );  // Per Y, Yp, Ypp pointers

  // Assegnazione dei puntatori alle sezioni di memoria
  m_Y    = m_mem_p( m_nspl );
  m_Yp   = m_mem_p( m_nspl );
  m_Ypp  = m_mem_p( m_nspl );
  m_X    = m_mem( m_npts );
  m_Ymin = m_mem( m_nspl );
  m_Ymax = m_mem( m_nspl );

  // Copia dei nodi X (condivisi da tutte le spline)
  if ( npts > 0 ) std::memcpy( m_X, X, npts * sizeof( *m_X ) );

  // ===================================================================
  // FASE 3: Costruzione delle singole spline
  // ===================================================================
  for ( integer spl{ 0 }; spl < nspl; ++spl )
  {
    UTILS_ASSERT( Y[spl] != nullptr, "{} Y[{}] array is null\n", msg, spl );

    // Riferimenti ai puntatori per questa spline
    real_type *& pY   = m_Y[spl];
    real_type *& pYp  = m_Yp[spl];
    real_type *& pYpp = m_Ypp[spl];

    // Allocazione e copia dei valori Y
    pY = m_mem( m_npts );
    if ( npts > 0 ) std::memcpy( pY, Y[spl], npts * sizeof( *pY ) );

    // ---------------------------------------------------------------
    // Calcolo min/max per questa spline
    // ---------------------------------------------------------------
    // NOTA: Per CONSTANT si esclude l'ultimo punto (comportamento corretto?)
    if ( stype[spl] == SplineType1D::CONSTANT )
    {
      m_Ymin[spl] = *std::min_element( pY, pY + npts - 1 );
      m_Ymax[spl] = *std::max_element( pY, pY + npts - 1 );
    }
    else
    {
      m_Ymin[spl] = *std::min_element( pY, pY + npts );
      m_Ymax[spl] = *std::max_element( pY, pY + npts );
    }

    // Inizializzazione puntatori derivate
    pYpp = pYp = nullptr;

    // ---------------------------------------------------------------
    // Allocazione memoria per derivate se necessario
    // ---------------------------------------------------------------
    switch ( stype[spl] )
    {
      case SplineType1D::QUINTIC_CUBIC:
      case SplineType1D::QUINTIC_AKIMA:
      case SplineType1D::QUINTIC_BESSEL:
      case SplineType1D::QUINTIC_PCHIP:
        pYpp = m_mem( m_npts );  // Alloca seconda derivata
        [[fallthrough]];

      case SplineType1D::CUBIC:
      case SplineType1D::AKIMA:
      case SplineType1D::BESSEL:
      case SplineType1D::PCHIP:
      case SplineType1D::HERMITE:
        pYp = m_mem( m_npts );  // Alloca prima derivata

        // Per HERMITE le derivate devono essere fornite dall'utente
        if ( stype[spl] == SplineType1D::HERMITE )
        {
          UTILS_ASSERT(
            Yp != nullptr && Yp[spl] != nullptr,
            "{} At spline n.{} named {}\n"
            "expect to find derivative values",
            msg,
            spl,
            headers[spl] );
          if ( npts > 0 ) std::memcpy( pYp, Yp[spl], npts * sizeof( *pYp ) );
        }
        [[fallthrough]];

      case SplineType1D::CONSTANT:
      case SplineType1D::LINEAR:
      case SplineType1D::SPLINE_SET:
      case SplineType1D::SPLINE_VEC: break;
    }

    // ---------------------------------------------------------------
    // Creazione e costruzione della spline specifica
    // ---------------------------------------------------------------
    string_view               h = headers[spl];
    std::unique_ptr<Spline> & s = m_splines[spl];

    m_is_monotone[spl] = -1;  // Valore di default (sconosciuto)

    // MIGLIORAMENTO: Questo switch è molto ripetitivo
    // Potrebbe essere refactorizzato con template o funzioni helper
    switch ( stype[spl] )
    {
      case SplineType1D::CONSTANT:
      {
        auto S = std::make_unique<ConstantSpline>( h );
        S->reserve_external( m_npts, m_X, pY );
        S->m_npts = m_npts;
        S->build();
        s = std::move( S );
      }
      break;

      case SplineType1D::LINEAR:
      {
        auto S = std::make_unique<LinearSpline>( h );
        S->reserve_external( m_npts, m_X, pY );
        S->m_npts = m_npts;
        S->build();

        // Check manuale monotonia per spline lineari
        // NOTA: Potrebbe essere spostato in LinearSpline::is_monotone()
        integer flag{ 1 };  // 1 = strettamente monotona crescente
        for ( integer j = 1; j < m_npts; ++j )
        {
          if ( pY[j - 1] > pY[j] )
          {
            flag = -1;  // Non monotona
            break;
          }
          // MIGLIORAMENTO: La condizione è complicata, meglio separare
          if ( Utils::is_zero( pY[j - 1] - pY[j] ) && m_X[j - 1] < m_X[j] ) flag = 0;  // Monotona non stretta (plateau)
        }
        m_is_monotone[spl] = flag;
        s                  = std::move( S );
      }
      break;

      // NOTA: Tutti i casi seguenti sono quasi identici
      // Solo il tipo di spline cambia
      case SplineType1D::CUBIC:
      {
        auto S = std::make_unique<CubicSpline>( h );
        S->reserve_external( m_npts, m_X, pY, pYp );
        S->m_npts = m_npts;
        S->build();
        m_is_monotone[spl] = S->is_monotone();
        s                  = std::move( S );
      }
      break;

      case SplineType1D::AKIMA:
      {
        auto S = std::make_unique<AkimaSpline>( h );
        S->reserve_external( m_npts, m_X, pY, pYp );
        S->m_npts = m_npts;
        S->build();
        m_is_monotone[spl] = S->is_monotone();
        s                  = std::move( S );
      }
      break;

      case SplineType1D::BESSEL:
      {
        auto S = std::make_unique<BesselSpline>( h );
        S->reserve_external( m_npts, m_X, pY, pYp );
        S->m_npts = m_npts;
        S->build();
        m_is_monotone[spl] = S->is_monotone();
        s                  = std::move( S );
      }
      break;

      case SplineType1D::PCHIP:
      {
        auto S = std::make_unique<PchipSpline>( h );
        S->reserve_external( m_npts, m_X, pY, pYp );
        S->m_npts = m_npts;
        S->build();
        m_is_monotone[spl] = S->is_monotone();
        s                  = std::move( S );
      }
      break;

      case SplineType1D::HERMITE:
      {
        auto S = std::make_unique<HermiteSpline>( h );
        S->reserve_external( m_npts, m_X, pY, pYp );
        S->m_npts = m_npts;
        S->build();
        m_is_monotone[spl] = S->is_monotone();
        s                  = std::move( S );
      }
      break;

      // Casi quintic con seconda derivata
      case SplineType1D::QUINTIC_CUBIC:
      {
        auto S = std::make_unique<QuinticSpline>( Spline_sub_type::CUBIC, h );
        S->reserve_external( m_npts, m_X, pY, pYp, pYpp );
        S->m_npts = m_npts;
        S->build();
        m_is_monotone[spl] = S->is_monotone();
        s                  = std::move( S );
      }
      break;

      case SplineType1D::QUINTIC_AKIMA:
      {
        auto S = std::make_unique<QuinticSpline>( Spline_sub_type::AKIMA, h );
        S->reserve_external( m_npts, m_X, pY, pYp, pYpp );
        S->m_npts = m_npts;
        S->build();
        m_is_monotone[spl] = S->is_monotone();
        s                  = std::move( S );
      }
      break;

      case SplineType1D::QUINTIC_BESSEL:
      {
        auto S = std::make_unique<QuinticSpline>( Spline_sub_type::BESSEL, h );
        S->reserve_external( m_npts, m_X, pY, pYp, pYpp );
        S->m_npts = m_npts;
        S->build();
        m_is_monotone[spl] = S->is_monotone();
        s                  = std::move( S );
      }
      break;

      case SplineType1D::QUINTIC_PCHIP:
      {
        auto S = std::make_unique<QuinticSpline>( Spline_sub_type::PCHIP, h );
        S->reserve_external( m_npts, m_X, pY, pYp, pYpp );
        S->m_npts = m_npts;
        S->build();
        m_is_monotone[spl] = S->is_monotone();
        s                  = std::move( S );
      }
      break;

      case SplineType1D::SPLINE_SET:
      case SplineType1D::SPLINE_VEC:
        UTILS_ERROR(
          "{} At spline n.{} named {}\n"
          "{} not allowed as spline type\n"
          "in SplineSet::build for {}-th spline\n",
          msg,
          spl,
          headers[spl],
          to_string( stype[spl] ),
          spl );
    }

    // ---------------------------------------------------------------
    // Registrazione del mapping nome -> indice
    // ---------------------------------------------------------------
    // ERRORE POTENZIALE: .data() potrebbe non essere sicuro
    // MIGLIORAMENTO: Usare s->name() direttamente se possibile
    // ATTENZIONE: Possibile duplicazione di nomi non gestita
    m_header_to_position.insert( { s->name().data(), static_cast<integer>( spl ) } );
  }

  // ===================================================================
  // FASE 4: Verifica finale allocazione memoria
  // ===================================================================
  // Controllo che tutta la memoria allocata sia stata effettivamente usata
  m_mem.must_be_empty( "SplineSet::build, baseValue" );
  m_mem_p.must_be_empty( "SplineSet::build, basePointer" );
}
