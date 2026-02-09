/**
 * \brief Build spline set from GenericContainer
 *
 * Constructs a set of splines from a GenericContainer with validated input data.
 * Supports multiple input formats for flexibility and provides comprehensive error checking.
 *
 * \param[in] gc GenericContainer containing spline set configuration
 *
 * \par Required Fields:
 * - "spline_type": vec_string_type - Array of spline type names (e.g., "linear", "cubic")
 * - "xdata": vec_real_type - X-coordinates shared by all splines (must be strictly increasing)
 * - "ydata": Multiple formats supported:
 *   * vec_real_type: Single spline data
 *   * mat_real_type: Matrix where each column is a spline
 *   * vector_type: Vector of vectors, one per spline
 *   * map_type: Map of {spline_name: y_values}
 *
 * \par Optional Fields:
 * - "headers": vec_string_type - Names for splines (required unless ydata is map_type)
 * - "ypdata": map_type - Derivative values for Hermite splines {spline_name: yp_values}
 * - "boundary": vector_type - Boundary conditions per spline (closed, extend, etc.)
 *
 * \par Boundary Configuration:
 * Each boundary element can contain:
 * - "closed": bool - If true, spline wraps around
 * - "extend"/"can_extend": bool - Allow extrapolation beyond data range
 * - "extend_constant": bool - Use constant extrapolation (vs linear)
 *
 * \throw UTILS_ASSERT if required fields missing or data dimensions inconsistent
 * \throw UTILS_ERROR if unsupported data types provided
 *
 * \note CONSTANT splines use n-1 points (interval-based values)
 * \note Warns about unused configuration keys to catch typos
 *
 * \see build() for the underlying construction mechanism
 */
void setup( GenericContainer const & gc )
{
  string const where = fmt::format( "SplineSet[{}]::setup( gc ): ", m_name );

  // ===================================================================
  // FASE 1: Raccolta e tracciamento delle chiavi disponibili
  // ===================================================================
  std::set<std::string> keywords;
  for ( auto const & pair : gc.get_map( where ) ) { keywords.insert( pair.first ); }

  // Strutture dati temporanee per la costruzione
  vec_string_type       spline_type_vec;  // Tipi di spline come stringhe
  vec_real_type         X;                // Nodi x condivisi
  vector<SplineType1D>  stype;            // Tipi di spline convertiti (enum)
  vec_string_type       headers;          // Nomi delle spline
  vector<vec_real_type> Y, Yp;            // Valori y e derivate

  // ===================================================================
  // FASE 2: Lettura e validazione dei campi obbligatori
  // ===================================================================

  // --- Lettura "spline_type" (obbligatorio) ---
  GenericContainer const & gc_stype = gc( "spline_type", where );
  keywords.erase( "spline_type" );

  // --- Lettura "xdata" (obbligatorio) ---
  GenericContainer const & gc_xdata = gc( "xdata", where );
  keywords.erase( "xdata" );

  // --- Lettura "ydata" (obbligatorio) ---
  GenericContainer const & gc_ydata = gc( "ydata", where );
  keywords.erase( "ydata" );

  // ===================================================================
  // FASE 3: Parsing e conversione dei tipi di spline
  // ===================================================================
  gc_stype.copyto_vec_string( spline_type_vec, ( where + "in reading 'spline_type'\n" ) );
  m_nspl = static_cast<integer>( spline_type_vec.size() );

  // Validazione: almeno una spline richiesta
  UTILS_ASSERT( m_nspl > 0, "{} expected at least one spline, got nspl = {}\n", where, m_nspl );

  // Conversione da stringhe a enum SplineType1D
  stype.resize( m_nspl );
  for ( integer spl = 0; spl < m_nspl; ++spl )
  {
    // NOTA: string_to_splineType1D lancia eccezione se tipo non riconosciuto
    stype[spl] = string_to_splineType1D( spline_type_vec[spl] );
  }

  // ===================================================================
  // FASE 4: Lettura headers (se ydata non è MAP)
  // ===================================================================
  // Quando ydata è una mappa, le chiavi diventano i nomi delle spline.
  // Altrimenti, gli headers devono essere forniti esplicitamente.
  if ( GC_type::MAP != gc_ydata.get_type() )
  {
    GenericContainer const & gc_headers = gc( "headers", where );
    keywords.erase( "headers" );
    gc_headers.copyto_vec_string( headers, ( where + ", reading 'headers'" ) );

    // Validazione: numero di headers deve corrispondere al numero di spline
    UTILS_ASSERT(
      headers.size() == static_cast<size_t>( m_nspl ),
      "{} field 'headers' expected to be of size {} found of size {}\n",
      where,
      m_nspl,
      headers.size() );

    // Validazione: nessun header può essere vuoto o null
    for ( integer spl = 0; spl < m_nspl; ++spl )
    {
      UTILS_ASSERT( !headers[spl].empty(), "{} header[{}] cannot be empty\n", where, spl );
    }
  }

  // ===================================================================
  // FASE 5: Lettura e validazione dei dati X (condivisi)
  // ===================================================================
  gc_xdata.copyto_vec_real( X, ( where + "reading 'xdata'" ) );
  m_npts = static_cast<integer>( X.size() );

  // Validazione: almeno 2 punti necessari per una spline
  UTILS_ASSERT( m_npts > 1, "{} expected at least 2 points, got npts = {}\n", where, m_npts );

  // Validazione: X deve essere strettamente crescente
  for ( integer i = 1; i < m_npts; ++i )
  {
    UTILS_ASSERT(
      X[i] > X[i - 1],
      "{} xdata must be strictly increasing: X[{}] = {} <= X[{}] = {}\n",
      where,
      i - 1,
      X[i - 1],
      i,
      X[i] );
  }

  // Pre-allocazione per le Y e Yp di tutte le spline
  Y.resize( m_nspl );
  Yp.resize( m_nspl );

  // ===================================================================
  // FASE 6: Lettura dei dati Y (supporta formati multipli)
  // ===================================================================
  switch ( gc_ydata.get_type() )
  {
    // ---------------------------------------------------------------
    // CASO 1: Vettore singolo → Una sola spline
    // ---------------------------------------------------------------
    case GC_type::VEC_BOOL:
    case GC_type::VEC_INTEGER:
    case GC_type::VEC_LONG:
    case GC_type::VEC_REAL:
    case GC_type::VEC_COMPLEX:
    {
      vec_real_type data;
      gc_ydata.copyto_vec_real( data, ( where + "reading 'ydata'" ) );

      // Deve esserci esattamente una spline
      UTILS_ASSERT(
        m_nspl == 1,
        "{} number of splines [{}] is incompatible with ydata as vector (expected 1 spline)\n",
        where,
        m_nspl );

      // Numero di punti deve corrispondere (o n-1 per CONSTANT)
      integer expected_npts = ( stype[0] == SplineType1D::CONSTANT ) ? m_npts - 1 : m_npts;
      UTILS_ASSERT(
        static_cast<size_t>( expected_npts ) == data.size(),
        "{} number of points [{}] differs from ydata size [{}] for spline type '{}'\n",
        where,
        expected_npts,
        data.size(),
        spline_type_vec[0] );

      // Usa move semantics invece di copy
      Y[0] = std::move( data );
    }
    break;

    // ---------------------------------------------------------------
    // CASO 2: Matrice → Ogni colonna è una spline
    // ---------------------------------------------------------------
    case GC_type::MAT_INTEGER:
    case GC_type::MAT_LONG:
    case GC_type::MAT_REAL:
    {
      mat_real_type data;
      gc_ydata.copyto_mat_real( data, ( where + "reading 'ydata'" ) );

      // Numero colonne deve corrispondere al numero di spline
      UTILS_ASSERT(
        static_cast<size_t>( m_nspl ) == data.num_cols(),
        "{} number of splines [{}] differs from ydata columns [{}]\n",
        where,
        m_nspl,
        data.num_cols() );

      // Numero righe deve corrispondere al numero di punti
      UTILS_ASSERT(
        static_cast<size_t>( m_npts ) == data.num_rows(),
        "{} number of points [{}] differs from ydata rows [{}]\n",
        where,
        m_npts,
        data.num_rows() );

      // Estrazione delle colonne
      for ( integer i = 0; i < m_nspl; ++i ) { data.get_column( static_cast<unsigned>( i ), Y[i] ); }
    }
    break;

    // ---------------------------------------------------------------
    // CASO 3: Vettore di vettori → Array di spline
    // ---------------------------------------------------------------
    case GC_type::VECTOR:
    {
      vector_type const & data = gc_ydata.get_vector();

      // Dimensione del vettore deve corrispondere al numero di spline
      UTILS_ASSERT(
        static_cast<size_t>( m_nspl ) == data.size(),
        "{} field 'ydata' expected of size {} found of size {}\n",
        where,
        m_nspl,
        data.size() );

      string msg1 = where + " reading 'ydata' columns";
      for ( integer spl = 0; spl < m_nspl; ++spl )
      {
        GenericContainer const & datai = data[spl];

        // Calcola numero atteso di righe (n-1 per CONSTANT, n altrimenti)
        integer expected_npts = m_npts;
        if ( stype[spl] == SplineType1D::CONSTANT ) { --expected_npts; }

        datai.copyto_vec_real( Y[spl], msg1 );

        // Validazione dimensione per questa spline
        UTILS_ASSERT(
          static_cast<size_t>( expected_npts ) == Y[spl].size(),
          "{} column {} of 'ydata' of type '{}' expected of size {} found of size {}\n",
          where,
          spl,
          spline_type_vec[spl],
          expected_npts,
          Y[spl].size() );
      }
    }
    break;

    // ---------------------------------------------------------------
    // CASO 4: Mappa → Chiavi = nomi spline, valori = dati
    // ---------------------------------------------------------------
    case GC_type::MAP:
    {
      map_type const & data = gc_ydata.get_map();

      // Dimensione mappa deve corrispondere al numero di spline
      UTILS_ASSERT(
        data.size() == static_cast<size_t>( m_nspl ),
        "{} field 'ydata' expected of size {} found of size {}\n",
        where,
        m_nspl,
        data.size() );

      // In modalità MAP, gli headers vengono estratti dalle chiavi della mappa
      headers.clear();
      headers.reserve( data.size() );

      string  msg1 = where + " reading 'ydata' columns";
      integer spl  = 0;

      for ( auto const & [name, container] : data )
      {
        // Estrai nome dalla chiave
        headers.emplace_back( name );

        // Validazione: nome non vuoto
        UTILS_ASSERT( !name.empty(), "{} spline name cannot be empty in ydata map\n", where );

        // Calcola numero atteso di righe
        integer expected_npts = m_npts;
        if ( stype[spl] == SplineType1D::CONSTANT ) { --expected_npts; }

        container.copyto_vec_real( Y[spl], msg1 );

        // Validazione dimensione (FIXED: era "ot type")
        UTILS_ASSERT(
          static_cast<size_t>( expected_npts ) == Y[spl].size(),
          "{} column '{}' of 'ydata' of type '{}' expected of size {} found of size {}\n",
          where,
          name,
          spline_type_vec[spl],
          expected_npts,
          Y[spl].size() );

        ++spl;
      }
    }
    break;

    // ---------------------------------------------------------------
    // CASO DEFAULT: Tipo non supportato
    // ---------------------------------------------------------------
    default:
      UTILS_ERROR(
        "{} field 'ydata' expected to be of type:\n"
        "  - vec_[int/long/real]_type (single spline)\n"
        "  - mat_[int/long/real]_type (matrix of splines)\n"
        "  - vector_type (vector of spline data)\n"
        "  - map_type (named splines)\n"
        "Found: '{}'\n",
        where,
        gc_ydata.get_type_name() );
      break;
  }

  // ===================================================================
  // FASE 7: Lettura derivate (opzionale, per spline Hermite)
  // ===================================================================
  if ( gc.exists( "ypdata" ) )
  {
    GenericContainer const & gc_ypdata = gc( "ypdata", where );
    keywords.erase( "ypdata" );

    // Le derivate devono essere fornite come MAP (nome → valori)
    UTILS_ASSERT(
      GC_type::MAP == gc_ypdata.get_type(),
      "{} field 'ypdata' expected to be of type 'map_type' found: '{}'\n",
      where,
      gc_ypdata.get_type_name() );

    // ---------------------------------------------------------------
    // Costruzione della mappa nome → indice (usa unordered_map)
    // ---------------------------------------------------------------
    std::unordered_map<string, integer> h_to_pos;
    h_to_pos.reserve( headers.size() );

    for ( integer idx = 0; idx < static_cast<integer>( headers.size() ); ++idx ) { h_to_pos[headers[idx]] = idx; }

    string           msg1 = where + " reading 'ypdata' columns";
    map_type const & data = gc_ypdata.get_map();

    // Per ogni derivata fornita nella mappa
    for ( auto const & [name, container] : data )
    {
      // Cerca la posizione corrispondente negli headers
      auto it_pos = h_to_pos.find( name );
      UTILS_ASSERT(
        it_pos != h_to_pos.end(),
        "{} column '{}' of 'ypdata' does not match any spline name in headers\n",
        where,
        name );

      integer spl = it_pos->second;

      // Calcola numero atteso di punti
      integer expected_npts = m_npts;
      if ( stype[spl] == SplineType1D::CONSTANT ) { --expected_npts; }

      container.copyto_vec_real( Yp[spl], msg1 );

      // FIXED: Confronta con Yp[spl].size() invece di Y[spl].size()
      // FIXED: Typo "or type" → "of type"
      UTILS_ASSERT(
        static_cast<size_t>( expected_npts ) == Yp[spl].size(),
        "{} column '{}' of 'ypdata' of type '{}' expected of size {} found of size {}\n",
        where,
        name,
        spline_type_vec[spl],
        expected_npts,
        Yp[spl].size() );
    }
  }

  // ===================================================================
  // FASE 8: Preparazione puntatori per chiamata a build()
  // ===================================================================
  // NOTA: build() richiede array di puntatori C-style (char**, real_type**)
  // Questa sezione prepara questi array dai vettori C++ moderni

  Utils::Malloc<void *> mem( where );
  mem.allocate( 3 * m_nspl );

  void ** pp__headers = mem( m_nspl );
  void ** pp__Y       = mem( m_nspl );
  void ** pp__Yp      = mem( m_nspl );

  // Conversione da vector<string>/vector<vector<real_type>> a puntatori C
  // ATTENZIONE: headers deve rimanere valido fino alla fine di build()
  for ( integer spl = 0; spl < m_nspl; ++spl )
  {
    // Puntatore al nome della spline (const char*)
    // SICUREZZA: headers[spl] deve rimanere valido durante build()
    pp__headers[spl] = const_cast<void *>( reinterpret_cast<void const *>( headers[spl].c_str() ) );

    // Puntatore ai dati Y
    pp__Y[spl] = reinterpret_cast<void *>( Y[spl].data() );

    // Puntatore ai dati Yp (nullptr se non forniti)
    pp__Yp[spl] = reinterpret_cast<void *>( Yp[spl].empty() ? nullptr : Yp[spl].data() );
  }

  // ===================================================================
  // FASE 9: Chiamata a build() con i dati preparati
  // ===================================================================
  // NOTA: I cast multipli sono necessari per l'interfaccia C-style di build()
  // ma sono sicuri perché i dati sorgente rimangono validi per tutta la chiamata
  this->build(
    m_nspl,
    m_npts,
    reinterpret_cast<char const **>( const_cast<void const **>( pp__headers ) ),
    stype.data(),
    X.data(),
    reinterpret_cast<real_type const **>( const_cast<void const **>( pp__Y ) ),
    reinterpret_cast<real_type const **>( const_cast<void const **>( pp__Yp ) ) );

  // mem viene deallocato automaticamente qui (RAII)

  // ===================================================================
  // FASE 10: Configurazione condizioni al contorno (opzionale)
  // ===================================================================
  if ( gc.exists( "boundary" ) )
  {
    GenericContainer const & gc_boundary = gc( "boundary" );
    keywords.erase( "boundary" );

    integer ne = static_cast<integer>( gc_boundary.get_num_elements() );

    // Validazione: una configurazione boundary per ogni spline
    // FIXED: Ordine corretto dei parametri nel messaggio
    UTILS_ASSERT(
      ne == m_nspl,
      "{} field 'boundary' expected a generic vector of size: {} but is of size: {}\n",
      where,
      m_nspl,
      ne );

    // Configurazione boundary per ogni spline
    for ( integer ispl = 0; ispl < ne; ++ispl )
    {
      Spline *                 S    = m_splines[ispl].get();
      GenericContainer const & item = gc_boundary( ispl, "SplineSet boundary data" );

      // Raccolta keywords per questo elemento boundary
      std::set<std::string> keywords2;
      for ( auto const & pair : item.get_map( where ) ) { keywords2.insert( pair.first ); }

      // ---------------------------------------------------------------
      // Configurazione closed/open
      // ---------------------------------------------------------------
      bool is_closed = false;
      item.get_if_exists( "closed", is_closed );
      keywords2.erase( "closed" );

      if ( is_closed )
      {
        // Spline periodica (endpoint = startpoint)
        S->make_closed();
      }
      else
      {
        // Spline aperta (gestione extrapolation)
        keywords2.erase( "extend" );
        keywords2.erase( "can_extend" );
        S->make_opened();

        bool can_extend = false;
        // Supporta sia "extend" che "can_extend" come alias
        if ( !item.get_if_exists( "extend", can_extend ) ) { item.get_if_exists( "can_extend", can_extend ); }

        if ( can_extend )
        {
          // Permette valutazione oltre il range dei dati
          S->make_unbounded();

          bool extend_constant = false;
          keywords2.erase( "extend_constant" );
          item.get_if_exists( "extend_constant", extend_constant );

          if ( extend_constant )
          {
            // Extrapolazione costante (usa valore endpoint)
            S->make_extended_constant();
          }
          else
          {
            // Extrapolazione con continuità della derivata
            S->make_extended_not_constant();
          }
        }
        else
        {
          // Non permette valutazione oltre il range
          S->make_bounded();
        }
      }

      // ---------------------------------------------------------------
      // Warning per chiavi non utilizzate in questo boundary
      // ---------------------------------------------------------------
      if ( !keywords2.empty() )
      {
        // Usa std::accumulate per costruire la stringa in modo efficiente
        string unused_keys{ std::accumulate(
          keywords2.begin(),
          keywords2.end(),
          string{},
          []( string acc, string const & key ) { return acc.empty() ? key : acc + " " + key; } ) };

        UTILS_WARNING( false, "{} spline N.{} of {} has unused boundary keys: {}\n", where, ispl, ne, unused_keys );
      }
    }
  }

  // ===================================================================
  // FASE 11: Warning per chiavi non utilizzate nel GC principale
  // ===================================================================
  if ( !keywords.empty() )
  {
    // Usa std::accumulate per costruire la stringa in modo efficiente
    string unused_keys{ std::accumulate(
      keywords.begin(),
      keywords.end(),
      string{},
      []( string acc, string const & key ) { return acc.empty() ? key : acc + " " + key; } ) };

    UTILS_WARNING( false, "{} unused configuration keys: {}\n", where, unused_keys );
  }
}
