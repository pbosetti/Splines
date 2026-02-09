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
 |   _     _                       ____        _ _
 |  | |   (_)_ __   ___  __ _ _ __/ ___| _ __ | (_)_ __   ___
 |  | |   | | '_ \ / _ \/ _` | '__\___ \| '_ \| | | '_ \ / _ \
 |  | |___| | | | |  __/ (_| | |   ___) | |_) | | | | | |  __/
 |  |_____|_|_| |_|\___|\__,_|_|  |____/| .__/|_|_|_| |_|\___|
 |                                      |_|
\*/

namespace Splines
{

  //! Linear spline class
  class LinearSpline : public Spline
  {
    Malloc_real m_mem_linear;
    bool        m_external_alloc = false;

  public:
    using Spline::build;

    //!
    //! Build an empty spline of `LinearSpline` type
    //!
    //! \param name the name of the spline
    //!
    explicit LinearSpline( string_view name = "LinearSpline" )
      : Spline( name ), m_mem_linear( fmt::format( "LinearSpline[{}]", name ) )
    {
      m_curve_extended_constant = true;  // by default linear spline extend constant
    }

    //!
    //! Spline destructor.
    //!
    ~LinearSpline() override {}

    //! Use externally allocated memory for `npts` points
    void reserve_external( integer const n, real_type *& p_x, real_type *& p_y )
    {
      if ( !m_external_alloc ) m_mem_linear.free();
      m_npts           = 0;
      m_npts_reserved  = n;
      m_external_alloc = true;
      m_X              = p_x;
      m_Y              = p_y;
    }

    // --------------------------- VIRTUALS -----------------------------------

    real_type eval( real_type const x ) const override
    {
      std::pair<integer, real_type> res( 0, x );
      m_search.find( res );
      return this->id_eval( res.first, res.second );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    real_type D( real_type const x ) const override
    {
      if ( m_curve_can_extend && m_curve_extended_constant )
      {
        // metto a 0 solo se esco dalla spline
        if ( x < m_X[0] || x > m_X[m_npts - 1] ) return 0;
      }
      std::pair<integer, real_type> res( 0, x );
      m_search.find( res );
      return this->id_D( res.first, res.second );
    }

    real_type DD( real_type const ) const override { return 0; }
    real_type DDD( real_type const ) const override { return 0; }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void D( real_type const x, real_type dd[2] ) const override
    {
      std::pair<integer, real_type> res( 0, x );
      m_search.find( res );
      integer const   ni = res.first;
      real_type const X  = res.second;
      real_type       DX = m_X[ni + 1] - m_X[ni];
      real_type       s  = ( X - m_X[ni] ) / DX;
      dd[0]              = ( 1 - s ) * m_Y[ni] + s * m_Y[ni + 1];
      dd[1]              = ( m_Y[ni + 1] - m_Y[ni] ) / DX;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void DD( real_type const x, real_type dd[3] ) const override
    {
      std::pair<integer, real_type> res( 0, x );
      m_search.find( res );
      integer const   ni = res.first;
      real_type const X  = res.second;
      real_type       DX = m_X[ni + 1] - m_X[ni];
      real_type       s  = ( X - m_X[ni] ) / DX;
      dd[0]              = ( 1 - s ) * m_Y[ni] + s * m_Y[ni + 1];
      dd[1]              = ( m_Y[ni + 1] - m_Y[ni] ) / DX;
      dd[2]              = 0;
    }

    real_type id_eval( integer const ni, real_type const x ) const override
    {
      if ( m_curve_can_extend && m_curve_extended_constant )
      {
        // estendo solo quando esco effettivamente
        if ( x < m_X[0] ) return m_Y[0];
        if ( x > m_X[m_npts - 1] ) return m_Y[m_npts - 1];
      }
      real_type const s = ( x - m_X[ni] ) / ( m_X[ni + 1] - m_X[ni] );
      return ( 1 - s ) * m_Y[ni] + s * m_Y[ni + 1];
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    real_type id_D( integer const i, real_type const x ) const override
    {
      if ( m_curve_can_extend && m_curve_extended_constant )
      {
        // estendo solo quando esco effettivamente
        if ( x < m_X[0] || x > m_X[m_npts - 1] ) return 0;
      }
      return ( m_Y[i + 1] - m_Y[i] ) / ( m_X[i + 1] - m_X[i] );
    }

    real_type id_DD( integer const, real_type const ) const override { return 0; }
    real_type id_DDD( integer const, real_type const ) const override { return 0; }

    void write_to_stream( ostream_type & s ) const override
    {
      integer const nseg = m_npts > 0 ? m_npts - 1 : 0;
      for ( integer i = 0; i < nseg; ++i )
        fmt::print(
          s,
          "segment N.{:4} X:[{:.5},{:.5}] Y:[{:.5},{:.5}] slope: {:.5}\n",
          i,
          m_X[i],
          m_X[i + 1],
          m_Y[i],
          m_Y[i + 1],
          ( m_Y[i + 1] - m_Y[i] ) / ( m_X[i + 1] - m_X[i] ) );
    }

    SplineType1D type() const override { return SplineType1D::LINEAR; }

    // --------------------------- VIRTUALS -----------------------------------

#ifdef AUTODIFF_SUPPORT
    autodiff::dual1st eval( autodiff::dual1st const & x ) const override
    {
      if ( m_curve_can_extend && m_curve_extended_constant )
      {
        // estendo solo quando esco effettivamente
        if ( x.val < m_X[0] )
        {
          autodiff::dual1st res;
          res.val  = m_Y[0];
          res.grad = 0;
          return res;
        }
        if ( x.val > m_X[m_npts - 1] )
        {
          autodiff::dual1st res;
          res.val  = m_Y[m_npts - 1];
          res.grad = 0;
          return res;
        }
      }

      std::pair<integer, real_type> res( 0, x.val );
      m_search.find( res );
      integer ni  = res.first;
      integer ni1 = ni + 1;

      autodiff::dual1st xx;
      xx.val              = res.second;
      xx.grad             = x.grad;
      autodiff::dual1st s = ( xx - m_X[ni] ) / ( m_X[ni1] - m_X[ni] );

      return ( 1 - s ) * m_Y[ni] + s * m_Y[ni1];
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    autodiff::dual2nd eval( autodiff::dual2nd const & x ) const override
    {
      if ( m_curve_can_extend && m_curve_extended_constant )
      {
        // estendo solo quando esco effettivamente
        if ( x.val.val < m_X[0] )
        {
          autodiff::dual2nd res;
          res.val.val   = m_Y[0];
          res.val.grad  = 0;
          res.grad.val  = 0;
          res.grad.grad = 0;
          return res;
        }
        if ( x.val.val > m_X[m_npts - 1] )
        {
          autodiff::dual2nd res;
          res.val.val   = m_Y[m_npts - 1];
          res.val.grad  = 0;
          res.grad.val  = 0;
          res.grad.grad = 0;
          return res;
        }
      }

      std::pair<integer, real_type> res( 0, x.val.val );
      m_search.find( res );
      integer ni  = res.first;
      integer ni1 = ni + 1;

      autodiff::dual2nd xx;
      xx.val.val   = res.second;
      xx.val.grad  = x.val.grad;
      xx.grad.val  = x.val.grad;
      xx.grad.grad = 0;

      // Calcolo di s = (xx - m_X[ni]) / (m_X[ni1] - m_X[ni])
      // Poiché stiamo sottraendo una costante, le derivate di xx rimangono le stesse
      autodiff::dual2nd s = ( xx - m_X[ni] ) / ( m_X[ni1] - m_X[ni] );

      return ( 1 - s ) * m_Y[ni] + s * m_Y[ni1];
    }

#endif  // AUTODIFF_SUPPORT

    void reserve( integer const npts ) override
    {
      if ( m_external_alloc && npts <= m_npts_reserved )
      {
        // nothing to do!, already allocated
      }
      else
      {
        m_mem_linear.reallocate( 2 * npts );
        m_npts_reserved  = npts;
        m_external_alloc = false;
        m_X              = m_mem_linear( npts );
        m_Y              = m_mem_linear( npts );
      }
      m_npts = 0;
    }

    void build() override { m_search.must_reset(); }

    void clear() override
    {
      if ( !m_external_alloc ) m_mem_linear.free();
      m_npts = m_npts_reserved = 0;
      m_external_alloc         = false;
      m_X = m_Y = nullptr;
    }

    integer  // order
    coeffs( real_type cfs[], real_type nodes[], bool const transpose ) const override
    {
      UTILS_ASSERT( m_npts >= 2, "LinearSpline::coeffs, npts={} must be >= 2\n", m_npts );

      integer const n = m_npts - 1;
      for ( integer i = 0; i < n; ++i )
      {
        real_type const a = m_Y[i];
        real_type const b = ( m_Y[i + 1] - m_Y[i] ) / ( m_X[i + 1] - m_X[i] );
        if ( transpose )
        {
          cfs[2 * i + 1] = a;
          cfs[2 * i + 0] = b;
        }
        else
        {
          cfs[n + i] = a;
          cfs[i]     = b;
        }
      }
      if ( m_npts > 0 ) std::memcpy( nodes, m_X, m_npts * sizeof( *nodes ) );
      return 2;
    }

    integer order() const override { return 2; }

    void setup( GenericContainer const & gc ) override
    {
      /*
      // gc["xdata"]
      // gc["ydata"]
      //
      */
      string const where = fmt::format( "LinearSpline[{}]::setup( gc ):", m_name );

      std::set<std::string> keywords;
      for ( auto const & pair : gc.get_map( where ) ) { keywords.insert( pair.first ); }
      keywords.erase( "spline_type" );

      GenericContainer const & gc_x = gc( "xdata", where );
      keywords.erase( "xdata" );
      GenericContainer const & gc_y = gc( "ydata", where );
      keywords.erase( "ydata" );

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

      this->build( x, y );
    }

    void y_min_max(
      integer &   i_min_pos,
      real_type & x_min_pos,
      real_type & y_min,
      integer &   i_max_pos,
      real_type & x_max_pos,
      real_type & y_max ) const override
    {
      UTILS_ASSERT( m_npts > 0, "LinearSpline[{}]::y_min_max() empty spline!", m_name );
      // find max min along the nodes
      i_min_pos = i_max_pos = 0;
      x_min_pos = x_max_pos = m_X[0];
      y_min = y_max = m_Y[0];
      for ( integer i = 1; i < m_npts; ++i )
      {
        real_type const & P1 = m_Y[i];
        if ( P1 > y_max )
        {
          y_max     = P1;
          x_max_pos = m_X[i];
          i_max_pos = i;
        }
        else if ( P1 < y_min )
        {
          y_min     = P1;
          x_min_pos = m_X[i];
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
      i_min_pos.clear();
      i_max_pos.clear();
      x_min_pos.clear();
      x_max_pos.clear();
      y_min.clear();
      y_max.clear();
      UTILS_ASSERT( m_npts > 0, "LinearSpline[{}]::y_min_max() empty spline!", m_name );
      // find max min along the nodes
      for ( integer i = 1; i < m_npts - 1; ++i )
      {
        real_type const & P0 = m_Y[i - 1];
        real_type const & P1 = m_Y[i];
        real_type const & P2 = m_Y[i + 1];
        if ( P1 > P0 && P1 > P2 )
        {
          y_max.emplace_back( P1 );
          x_max_pos.emplace_back( m_X[i] );
          i_max_pos.emplace_back( i );
        }
        else if ( P1 < P0 && P1 < P2 )
        {
          y_min.emplace_back( P1 );
          x_min_pos.emplace_back( m_X[i] );
          i_min_pos.emplace_back( i );
        }
      }
    }
    void copy_spline( LinearSpline const & S )
    {
      LinearSpline::reserve( S.m_npts );
      m_npts = S.m_npts;
      if ( m_npts > 0 )
      {
        size_t const size_bytes = m_npts * sizeof( *m_X );
        std::memcpy( m_X, S.m_X, size_bytes );
        std::memcpy( m_Y, S.m_Y, size_bytes );
      }
      copy_flags( S );
    }
  };

}  // namespace Splines

// EOF: SplineLinear.hxx
