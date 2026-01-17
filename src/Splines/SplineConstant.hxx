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
 |    ____                _              _       ____        _ _
 |   / ___|___  _ __  ___| |_ __ _ _ __ | |_ ___/ ___| _ __ | (_)_ __   ___
 |  | |   / _ \| '_ \/ __| __/ _` | '_ \| __/ __\___ \| '_ \| | | '_ \ / _ \
 |  | |__| (_) | | | \__ \ || (_| | | | | |_\__ \___) | |_) | | | | | |  __/
 |   \____\___/|_| |_|___/\__\__,_|_| |_|\__|___/____/| .__/|_|_|_| |_|\___|
 |                                                    |_|
\*/

namespace Splines
{

  using std::copy_n;

  //! Picewise constants spline class
  class ConstantSpline : public Spline
  {
    Malloc_real m_mem_constant;
    bool        m_external_alloc{ false };

  public:
    using Spline::build;

    //!
    //! Build an empty spline of `ConstantSpline` type
    //!
    //! \param name the name of the spline
    //!
    explicit ConstantSpline( string_view name = "ConstantSpline" )
      : Spline( name ), m_mem_constant( fmt::format( "ConstantSpline[{}]", name ) )
    {
    }

    //!
    //! Spline destructor.
    //!
    ~ConstantSpline() override {}

    //! Use externally allocated memory for `npts` points
    void reserve_external( integer const n, real_type *& p_x, real_type *& p_y )
    {
      if ( !m_external_alloc ) m_mem_constant.free();
      m_npts           = 0;
      m_npts_reserved  = n;
      m_external_alloc = true;
      m_X              = p_x;
      m_Y              = p_y;
    }

    // --------------------------- VIRTUALS -----------------------------------
    //!
    //! \name Build
    //!
    ///@{

    //!
    //! Build the spline with the data stored
    //!
    void build() override { m_search.must_reset(); }

    //!
    //! Build the spline with the data passed as arguments
    //!
    //! \param x    \f$ x \f$ coordinates of the points
    //! \param incx access elements as `x[0]`, `x[incx]`, `x[2*incx]`,...
    //! \param y    \f$ y \f$ coordinates of the points
    //! \param incy access elements as `y[0]`, `y[incx]`, `y[2*incx]`,...
    //! \param n    the number of the points
    //!
    void build( real_type const x[], integer incx, real_type const y[], integer incy, integer const n ) override
    {
      reserve( n );
      for ( integer i = 0; i < n; ++i ) m_X[i] = x[i * incx];
      for ( integer i = 0; i + 1 < n; ++i ) m_Y[i] = y[i * incy];
      m_npts = n;
      build();
    }

    ///@}

    //!
    //! \name Evaluate
    //!
    ///@{
    real_type eval( real_type const x ) const override
    {
      std::pair<integer, real_type> res( 0, x );
      m_search.find( res );
      return m_Y[res.first];
    }

    real_type D( real_type const ) const override { return 0; }
    real_type DD( real_type const ) const override { return 0; }
    real_type DDD( real_type const ) const override { return 0; }

    void D( real_type const x, real_type dd[2] ) const override
    {
      dd[0] = eval( x );
      dd[1] = 0;
    }

    void DD( real_type const x, real_type dd[3] ) const override
    {
      dd[0] = eval( x );
      dd[1] = 0;
      dd[2] = 0;
    }

    ///@}

    //!
    //! \name Evaluation when segment is known
    //!
    ///@{
    real_type id_eval( integer const ni, [[maybe_unused]] real_type const x ) const override { return m_Y[ni]; }
    real_type id_D( [[maybe_unused]] integer const ni, [[maybe_unused]] real_type const x ) const override { return 0; }
    real_type id_DD( [[maybe_unused]] integer const ni, [[maybe_unused]] real_type const x ) const override
    {
      return 0;
    }
    real_type id_DDD( [[maybe_unused]] integer const ni, [[maybe_unused]] real_type const x ) const override
    {
      return 0;
    }
    ///@}

#ifdef AUTODIFF_SUPPORT
    autodiff::dual1st eval( autodiff::dual1st const & x ) const override
    {
      using autodiff::dual1st;
      using autodiff::detail::val;
      dual1st res{ eval( val( x ) ) };
      res.grad = 0;
      return res;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    autodiff::dual2nd eval( autodiff::dual2nd const & x ) const override
    {
      using autodiff::dual2nd;
      using autodiff::detail::val;
      dual2nd res{ eval( val( x ) ) };
      res.grad      = 0;
      res.grad.grad = 0;
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
        fmt::print( s, "segment N. {:4} X:[{:.5},{:.5}] Y:{:.5}\n", i, m_X[i], m_X[i + 1], m_Y[i] );
    }

    SplineType1D type() const override { return SplineType1D::CONSTANT; }

    // --------------------------- VIRTUALS -----------------------------------

    void reserve( integer const npts ) override
    {
      if ( m_external_alloc && npts <= m_npts_reserved )
      {
        // nothing to do!, already allocated
      }
      else
      {
        m_mem_constant.reallocate( 2 * npts );
        m_npts_reserved  = npts;
        m_external_alloc = false;
        m_X              = m_mem_constant( npts );
        m_Y              = m_mem_constant( npts );
      }
      m_npts = 0;
    }

    void clear() override
    {
      if ( !m_external_alloc ) m_mem_constant.free();
      m_npts = m_npts_reserved = 0;
      m_external_alloc         = false;
      m_X = m_Y = nullptr;
    }

    integer  // order
    coeffs( real_type cfs[], real_type nodes[], [[maybe_unused]] bool transpose = false ) const override
    {
      integer const nseg{ m_npts > 0 ? m_npts - 1 : 0 };
      copy_n( m_X, m_npts, nodes );
      copy_n( m_Y, nseg, cfs );
      return 1;
    }

    integer order() const override { return 1; }

    void setup( GenericContainer const & gc ) override
    {
      /*
      // gc["xdata"]
      // gc["ydata"]
      //
      */
      string const where{ fmt::format( "ConstantSpline[{}]::setup( gc ):", m_name ) };

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
      UTILS_ASSERT( m_npts > 0, "ConstantSpline[{}]::y_min_max() empty spline!", m_name );
      // find max min alongh the nodes
      i_min_pos = i_max_pos = 0;
      x_min_pos = x_max_pos = m_X[0];
      y_min = y_max = m_Y[0];
      for ( integer i = 1; i < m_npts - 1; ++i )
      {
        real_type const & P1{ m_Y[i] };
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
      UTILS_ASSERT( m_npts > 0, "ConstantSpline[{}]::y_min_max() empty spline!", m_name );
      // find max min along the nodes
      for ( integer i = 1; i < m_npts - 1; ++i )
      {
        real_type const & P0{ m_Y[i - 1] };
        real_type const & P1{ m_Y[i] };
        real_type const & P2{ m_Y[i + 1] };
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

    void copy_spline( ConstantSpline const & S )
    {
      ConstantSpline::reserve( S.m_npts );
      m_npts = S.m_npts;
      copy_n( S.m_X, m_npts, m_X );
      copy_n( S.m_Y, m_npts - 1, m_Y );
      copy_flags( S );
    }
  };
}  // namespace Splines

// EOF: SplineConstant.hxx
