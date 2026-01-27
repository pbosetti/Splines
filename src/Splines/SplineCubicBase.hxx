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
 |    ____      _     _        ____        _ _              ____
 |   / ___|   _| |__ (_) ___  / ___| _ __ | (_)_ __   ___  | __ )  __ _ ___  ___
 |  | |  | | | | '_ \| |/ __| \___ \| '_ \| | | '_ \ / _ \ |  _ \ / _` / __|/ _ \
 |  | |__| |_| | |_) | | (__   ___) | |_) | | | | | |  __/ | |_) | (_| \__ \  __/
 |   \____\__,_|_.__/|_|\___| |____/| .__/|_|_|_| |_|\___| |____/ \__,_|___/\___|
 |                                  |_|
\*/

namespace Splines
{

  //!
  //! cubic spline base class
  //!
  class CubicSplineBase : public Spline
  {
  protected:
    Malloc_real m_mem_cubic;
    real_type * m_Yp{ nullptr };
    bool        m_external_alloc{ false };

  public:
    using Spline::build;

    //!
    //! \name Contructors/Destructors
    ///@{
    //!
    //! Spline constructor.
    //!
    explicit CubicSplineBase( string_view name = "CubicSplineBase" )
      : Spline( name ), m_mem_cubic( fmt::format( "CubicSplineBase[{}]::m_mem_cubic", name ) )
    {
    }

    ~CubicSplineBase() override {}
    ///@}

    //!
    //! Build a copy of spline `S`
    //!
    void copy_spline( CubicSplineBase const & S )
    {
      // 1. Protezione contro l'auto-assegnazione (Cruciale!)
      if ( this == &S ) return;

      // 2. Alloca memoria
      CubicSplineBase::reserve( S.m_npts );
      m_npts = S.m_npts;

      // 3. Copia in blocco (Bulk Copy)
      if ( m_npts > 0 )
      {
        // Calcolo size in bytes una sola volta
        size_t const byte_size = m_npts * sizeof( real_type );

        std::memcpy( m_X, S.m_X, byte_size );
        std::memcpy( m_Y, S.m_Y, byte_size );
        std::memcpy( m_Yp, S.m_Yp, byte_size );
      }

      copy_flags( S );
    }

    //!
    //! Return the pointer of values of yp-nodes.
    //!
    real_type const * yp_nodes() const { return m_Yp; }

    //!
    //! Return the i-th node of the spline (y' component).
    //!
    real_type yp_node( integer i ) const { return m_Yp[i]; }

    //!
    //! Change X-range of the spline.
    //!
    void set_range( real_type xmin, real_type xmax )
    {
      // 1. Calcola il fattore di scala PRIMA che m_X venga modificato.
      //    Se stiamo scalando le derivate (y' = dy/dx), il fattore è OldRange / NewRange.
      real_type const dx_old = m_X[m_npts - 1] - m_X[0];
      real_type const dx_new = xmax - xmin;

      // Protezione divisione per zero e calcolo rapporto
      // Se dx_new è 0, evitiamo NaN.
      real_type const recS = ( std::abs( dx_new ) > 1e-12 ) ? dx_new / dx_old : 1.0;

      // 2. Chiama la classe base per scalare m_X
      Spline::set_range( xmin, xmax );

      // 3. Check ottimizzazione: Se il fattore è 1 (nessun cambio scala), esci.
      if ( std::abs( recS - 1.0 ) < 1e-12 ) return;

      // 4. Applica il fattore di scala usando Eigen (SIMD)
      Eigen::Map<Eigen::Array<real_type, Eigen::Dynamic, 1>> map_vals( m_X, m_npts );
      map_vals -= m_X[0];
      map_vals *= recS;
      map_vals += xmin;
      m_X[m_npts - 1] = xmax;
    }

    //!
    //! Use externally allocated memory for `npts` points.
    //!
    void reserve_external( integer n, real_type *& p_x, real_type *& p_y, real_type *& p_dy )
    {
      m_npts_reserved  = n;
      m_X              = p_x;
      m_Y              = p_y;
      m_Yp             = p_dy;
      m_external_alloc = true;
      m_npts           = 0;
    }

    void y_min_max(
      integer &   i_min_pos,
      real_type & x_min_pos,
      real_type & y_min,
      integer &   i_max_pos,
      real_type & x_max_pos,
      real_type & y_max ) const override
    {
      UTILS_ASSERT( m_npts > 0, "CubicSplineBase[{}]::y_min_max() empty spline!", m_name );
      // find max min alongh the nodes
      i_min_pos = i_max_pos = 0;
      x_min_pos = x_max_pos = m_X[0];
      y_min = y_max = m_Y[0];
      PolynomialRoots::Quadratic q;
      for ( integer i = 1; i < m_npts; ++i )
      {
        real_type const & X0  = m_X[i - 1];
        real_type const & X1  = m_X[i];
        real_type const & P0  = m_Y[i - 1];
        real_type const & P1  = m_Y[i];
        real_type const & DP0 = m_Yp[i - 1];
        real_type const & DP1 = m_Yp[i];
        real_type const   H   = X1 - X0;
        real_type         A, B, C, D;
        Hermite3_to_poly( H, P0, P1, DP0, DP1, A, B, C, D );
        q.setup( 3 * A, 2 * B, C );
        real_type     r[2];
        integer const nr = q.getRootsInOpenRange( 0, H, r );
        for ( integer j = 0; j < nr; ++j )
        {
          real_type const rr = r[j];
          real_type const yy = ( ( ( A * rr ) + B ) * rr + C ) * rr + D;
          if ( yy > y_max )
          {
            y_max     = yy;
            x_max_pos = X0 + rr;
            i_max_pos = i;
          }
          else if ( yy < y_min )
          {
            y_min     = yy;
            x_min_pos = X0 + rr;
            i_min_pos = i;
          }
        }
        if ( P1 > y_max )
        {
          y_max     = P1;
          x_max_pos = X1;
          i_max_pos = i;
        }
        else if ( P1 < y_min )
        {
          y_min     = P1;
          x_min_pos = X1;
          i_min_pos = i;
        }
      }
    }

    void y_min_max(
      vector<integer> &   i_min_pos,
      vector<real_type> & x_min_pos,
      vector<real_type> & y_min,
      vector<integer> &   i_max_pos,
      vector<real_type> & x_max_pos,
      vector<real_type> & y_max ) const override
    {
      constexpr real_type epsi{ 1e-8 };
      i_min_pos.clear();
      i_max_pos.clear();
      x_min_pos.clear();
      x_max_pos.clear();
      y_min.clear();
      y_max.clear();
      UTILS_ASSERT( m_npts > 0, "CubicSplineBase[{}]::y_min_max() empty spline!", m_name );
      // find max min along the nodes
      if ( m_Yp[0] >= 0 )
      {
        y_min.emplace_back( m_Y[0] );
        x_min_pos.emplace_back( m_X[0] );
        i_min_pos.emplace_back( 0 );
      }
      if ( m_Yp[0] <= 0 )
      {
        y_max.emplace_back( m_Y[0] );
        x_max_pos.emplace_back( m_X[0] );
        i_max_pos.emplace_back( 0 );
      }
      PolynomialRoots::Quadratic q;
      for ( integer i = 1; i < m_npts; ++i )
      {
        real_type const & X0  = m_X[i - 1];
        real_type const & X1  = m_X[i];
        real_type const & P0  = m_Y[i - 1];
        real_type const & P1  = m_Y[i];
        real_type const & DP0 = m_Yp[i - 1];
        real_type const & DP1 = m_Yp[i];
        real_type const   H   = X1 - X0;
        real_type         A, B, C, D;
        Hermite3_to_poly( H, P0, P1, DP0, DP1, A, B, C, D );
        q.setup( 3 * A, 2 * B, C );
        real_type     r[2];
        integer const nr = q.getRootsInOpenRange( 0, H, r );
        for ( integer j = 0; j < nr; ++j )
        {
          real_type const rr = r[j];
          real_type const yy = ( ( ( A * rr ) + B ) * rr + C ) * rr + D;
          real_type const ddy{ 3 * A * rr + B };
          if ( ddy > 0 )
          {
            y_min.emplace_back( yy );
            x_min_pos.emplace_back( X0 + rr );
            i_min_pos.emplace_back( i );
          }
          else if ( ddy < 0 )
          {
            y_max.emplace_back( yy );
            x_max_pos.emplace_back( X0 + rr );
            i_max_pos.emplace_back( i );
          }
        }
        if ( i + 1 >= m_npts ) continue;
        if ( abs( DP1 ) > ( m_X[i + 1] - m_X[i - 1] ) * epsi ) continue;
        real_type const & X2  = m_X[i + 1];
        real_type const & P2  = m_Y[i + 1];
        real_type const & DP2 = m_Yp[i + 1];
        real_type         A1, B1, C1, D1;
        Hermite3_to_poly( X2 - X1, P1, P2, DP1, DP2, A1, B1, C1, D1 );
        real_type const DD = 2 * A * H + B;
        if ( DD >= 0 && B1 >= 0 )
        {
          y_min.emplace_back( P1 );
          x_min_pos.emplace_back( X1 );
          i_min_pos.emplace_back( i );
        }
        else if ( DD <= 0 && B1 <= 0 )
        {
          y_max.emplace_back( P1 );
          x_max_pos.emplace_back( X1 );
          i_max_pos.emplace_back( i );
        }
      }
      if ( m_Yp[m_npts - 1] <= 0 )
      {
        y_min.emplace_back( m_Y[m_npts - 1] );
        x_min_pos.emplace_back( m_X[m_npts - 1] );
        i_min_pos.emplace_back( 0 );
      }
      if ( m_Yp[m_npts - 1] >= 0 )
      {
        y_max.emplace_back( m_Y[m_npts - 1] );
        x_max_pos.emplace_back( m_X[m_npts - 1] );
        i_max_pos.emplace_back( m_npts - 1 );
      }
    }

    // --------------------------- VIRTUALS -----------------------------------

    //!
    //! \name Evaluation
    //!
    ///@{

    real_type eval( real_type const x ) const override
    {
      std::pair<integer, real_type> res( 0, x );
      m_search.find( res );
      return this->id_eval( res.first, res.second );
    }

    real_type D( real_type const x ) const override
    {
      std::pair<integer, real_type> res( 0, x );
      m_search.find( res );
      return this->id_D( res.first, res.second );
    }

    real_type DD( real_type const x ) const override
    {
      std::pair<integer, real_type> res( 0, x );
      m_search.find( res );
      return this->id_DD( res.first, res.second );
    }

    real_type DDD( real_type const x ) const override
    {
      std::pair<integer, real_type> res( 0, x );
      m_search.find( res );
      return this->id_DDD( res.first, res.second );
    }

    real_type DDDD( real_type const ) const override { return 0; }
    real_type DDDDD( real_type const ) const override { return 0; }

    void D( real_type const x, real_type dd[2] ) const override
    {
      std::pair<integer, real_type> res( 0, x );
      m_search.find( res );
      integer const   ni = res.first;
      real_type const X  = res.second;
      real_type       dx = X - m_X[ni];
      real_type       DX = m_X[ni + 1] - m_X[ni];
      real_type       base[4], base_D[4];
      Splines::Hermite3( dx, DX, base );
      Splines::Hermite3_D( dx, DX, base_D );
      dd[0] = base[0] * m_Y[ni] + base[1] * m_Y[ni + 1] + base[2] * m_Yp[ni] + base[3] * m_Yp[ni + 1];
      dd[1] = base_D[0] * m_Y[ni] + base_D[1] * m_Y[ni + 1] + base_D[2] * m_Yp[ni] + base_D[3] * m_Yp[ni + 1];
    }

    void DD( real_type const x, real_type dd[3] ) const override
    {
      std::pair<integer, real_type> res( 0, x );
      m_search.find( res );
      integer const   ni = res.first;
      real_type const X  = res.second;
      real_type       base[4], base_D[4], base_DD[4];
      real_type       dx = X - m_X[ni];
      real_type       DX = m_X[ni + 1] - m_X[ni];
      Splines::Hermite3( dx, DX, base );
      Splines::Hermite3_D( dx, DX, base_D );
      Splines::Hermite3_DD( dx, DX, base_DD );
      dd[0] = base[0] * m_Y[ni] + base[1] * m_Y[ni + 1] + base[2] * m_Yp[ni] + base[3] * m_Yp[ni + 1];
      dd[1] = base_D[0] * m_Y[ni] + base_D[1] * m_Y[ni + 1] + base_D[2] * m_Yp[ni] + base_D[3] * m_Yp[ni + 1];
      dd[2] = base_DD[0] * m_Y[ni] + base_DD[1] * m_Y[ni + 1] + base_DD[2] * m_Yp[ni] + base_DD[3] * m_Yp[ni + 1];
    }

    ///@}

    //!
    //! \name Evaluation when segment is known
    ///@{
    real_type id_eval( integer const ni, real_type const x ) const override
    {
      if ( m_curve_can_extend && m_curve_extended_constant )
      {
        // estendo solo quando esco effettivamente
        if ( x < m_X[0] ) return m_Y[0];
        if ( x > m_X[m_npts - 1] ) return m_Y[m_npts - 1];
      }
      real_type base[4];
      Splines::Hermite3( x - m_X[ni], m_X[ni + 1] - m_X[ni], base );
      return base[0] * m_Y[ni] + base[1] * m_Y[ni + 1] + base[2] * m_Yp[ni] + base[3] * m_Yp[ni + 1];
    }

    real_type id_D( integer const ni, real_type const x ) const override
    {
      if ( m_curve_can_extend && m_curve_extended_constant )
      {
        // estendo solo quando esco effettivamente
        if ( x < m_X[0] || x > m_X[m_npts - 1] ) return 0;
      }
      real_type base_D[4];
      Splines::Hermite3_D( x - m_X[ni], m_X[ni + 1] - m_X[ni], base_D );
      return base_D[0] * m_Y[ni] + base_D[1] * m_Y[ni + 1] + base_D[2] * m_Yp[ni] + base_D[3] * m_Yp[ni + 1];
    }

    real_type id_DD( integer const ni, real_type const x ) const override
    {
      if ( m_curve_can_extend && m_curve_extended_constant )
      {
        // estendo solo quando esco effettivamente
        if ( x < m_X[0] || x > m_X[m_npts - 1] ) return 0;
      }
      real_type base_DD[4];
      Splines::Hermite3_DD( x - m_X[ni], m_X[ni + 1] - m_X[ni], base_DD );
      return base_DD[0] * m_Y[ni] + base_DD[1] * m_Y[ni + 1] + base_DD[2] * m_Yp[ni] + base_DD[3] * m_Yp[ni + 1];
    }

    real_type id_DDD( integer const ni, real_type const x ) const override
    {
      if ( m_curve_can_extend && m_curve_extended_constant )
      {
        // estendo solo quando esco effettivamente
        if ( x < m_X[0] || x > m_X[m_npts - 1] ) return 0;
      }
      real_type base_DDD[4];
      Splines::Hermite3_DDD( x - m_X[ni], m_X[ni + 1] - m_X[ni], base_DDD );
      return base_DDD[0] * m_Y[ni] + base_DDD[1] * m_Y[ni + 1] + base_DDD[2] * m_Yp[ni] + base_DDD[3] * m_Yp[ni + 1];
    }

    real_type id_DDDD( integer const, real_type const ) const override { return 0; }
    real_type id_DDDDD( integer const, real_type const ) const override { return 0; }
    ///@}

#ifdef AUTODIFF_SUPPORT
    autodiff::dual1st eval( autodiff::dual1st const & x ) const override
    {
      using autodiff::dual1st;
      using autodiff::detail::val;
      real_type dd[2];
      D( val( x ), dd );
      dual1st res{ dd[0] };
      res.grad = dd[1] * x.grad;
      return res;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    autodiff::dual2nd eval( autodiff::dual2nd const & x ) const override
    {
      using autodiff::dual2nd;
      using autodiff::detail::val;
      real_type dd[3], xg{ val( x.grad ) };
      DD( val( x ), dd );
      dual2nd res{ dd[0] };
      res.grad      = dd[1] * xg;
      res.grad.grad = dd[1] * x.grad.grad + dd[2] * ( xg * xg );
      return res;
    }

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

    void write_to_stream( ostream_type & s ) const override
    {
      integer const nseg{ m_npts > 0 ? m_npts - 1 : 0 };
      for ( integer i = 0; i < nseg; ++i )
        fmt::print(
          s,
          "segment N.{:4} X:[{:.5},{:.5}] Y:[{:.5},{:.5}] Yp:[{:.5},{:.5}] slope: {}\n",
          i,
          m_X[i],
          m_X[i + 1],
          m_Y[i],
          m_Y[i + 1],
          m_Yp[i],
          m_Yp[i + 1],
          ( m_Y[i + 1] - m_Y[i] ) / ( m_X[i + 1] - m_X[i] ) );
    }

    // --------------------------- VIRTUALS -----------------------------------

    //!
    //! \name Build
    //!
    ///@{

    //!
    //! Build a spline.
    //!
    //! \param[in] x     vector of x-coordinates
    //! \param[in] incx  access elements as x[0], x[incx], x[2*incx],...
    //! \param[in] y     vector of y-coordinates
    //! \param[in] incy  access elements as y[0], y[incy], y[2*incy],...
    //! \param[in] yp    vector of y-defivative
    //! \param[in] incyp access elements as yp[0], yp[incyp], yp[2*incyy],...
    //! \param[in] n     total number of points
    //!
    void build(
      real_type const x[],
      integer const   incx,
      real_type const y[],
      integer const   incy,
      real_type const yp[],
      integer const   incyp,
      integer         n )
    {
      this->reserve( n );
      for ( integer i = 0; i < n; ++i )
      {
        m_X[i]  = x[i * incx];
        m_Y[i]  = y[i * incy];
        m_Yp[i] = yp[i * incyp];
      }
      m_npts = n;
      m_search.must_reset();
    }

    //!
    //! Build a spline.
    //!
    //! \param[in] x  vector of x-coordinates
    //! \param[in] y  vector of y-coordinates
    //! \param[in] yp vector of y'-coordinates
    //! \param[in] n  total number of points
    //!
    inline void build( real_type const x[], real_type const y[], real_type const yp[], integer const n )
    {
      this->build( x, 1, y, 1, yp, 1, n );
    }

    //!
    //! Build a spline.
    //!
    //! \param[in] x  vector of x-coordinates
    //! \param[in] y  vector of y-coordinates
    //! \param[in] yp vector of y'-coordinates
    //!
    void build( vector<real_type> const & x, vector<real_type> const & y, vector<real_type> const & yp )
    {
      integer N{ static_cast<integer>( x.size() ) };
      if ( N > static_cast<integer>( y.size() ) ) N = static_cast<integer>( y.size() );
      if ( N > static_cast<integer>( yp.size() ) ) N = static_cast<integer>( yp.size() );
      this->build( x.data(), 1, y.data(), 1, yp.data(), 1, N );
    }

    void reserve( integer npts ) override
    {
      if ( m_external_alloc && npts <= m_npts_reserved )
      {
        // nothing to do!, already allocated
      }
      else
      {
        m_npts_reserved = npts;
        m_mem_cubic.reallocate( 3 * npts );
        m_X              = m_mem_cubic( npts );
        m_Y              = m_mem_cubic( npts );
        m_Yp             = m_mem_cubic( npts );
        m_external_alloc = false;
      }
      m_npts = 0;
    }

    ///@}

    void clear() override
    {
      if ( !m_external_alloc ) m_mem_cubic.free();
      m_npts = m_npts_reserved = 0;
      m_external_alloc         = false;
      m_X = m_Y = m_Yp = nullptr;
    }

    integer  // order
    coeffs( real_type cfs[], real_type nodes[], bool transpose = false ) const override
    {
      UTILS_ASSERT( m_npts >= 2, "CubicSplineBase::coeffs, npts={} must be >= 2\n", m_npts );

      integer const n = m_npts - 1;
      for ( integer i = 0; i < n; ++i )
      {
        real_type const H{ m_X[i + 1] - m_X[i] };
        real_type       a, b, c, d;

        Hermite3_to_poly( H, m_Y[i], m_Y[i + 1], m_Yp[i], m_Yp[i + 1], a, b, c, d );

        if ( transpose )
        {
          cfs[4 * i + 3] = a;
          cfs[4 * i + 2] = b;
          cfs[4 * i + 1] = c;
          cfs[4 * i + 0] = d;
        }
        else
        {
          cfs[i + 3 * n] = a;
          cfs[i + 2 * n] = b;
          cfs[i + 1 * n] = c;
          cfs[i + 0 * n] = d;
        }
      }
      if ( m_npts > 0 ) std::memcpy( nodes, m_X, m_npts * sizeof( real_type ) );
      return 4;
    }

    integer order() const override { return 4; }

    bool is_monotone() const { return check_cubic_spline_monotonicity( m_X, m_Y, m_Yp, m_npts ); }

#ifdef SPLINES_BACK_COMPATIBILITY
    void      copySpline( CubicSplineBase const & S ) { this->copy_spline( S ); }
    integer   numPoints() const { return m_npts; }
    real_type xNode( integer i ) const { return m_X[i]; }
    real_type yNode( integer i ) const { return m_Y[i]; }
    real_type ypNode( integer i ) const { return this->yp_node( i ); }
    real_type xBegin() const { return m_X[0]; }
    real_type yBegin() const { return m_Y[0]; }
    real_type xEnd() const { return m_X[m_npts - 1]; }
    real_type yEnd() const { return m_Y[m_npts - 1]; }
    real_type xMin() const { return m_X[0]; }
    real_type xMax() const { return m_X[m_npts - 1]; }
    real_type yMin() const { return y_min(); }
    real_type yMax() const { return y_max(); }
#endif
  };

}  // namespace Splines
