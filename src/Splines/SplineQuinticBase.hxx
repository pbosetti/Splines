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
 |    ___        _       _   _      ____        _ _            ____
 |   / _ \ _   _(_)_ __ | |_(_) ___/ ___| _ __ | (_)_ __   ___| __ )  __ _ ___  ___
 |  | | | | | | | | '_ \| __| |/ __\___ \| '_ \| | | '_ \ / _ \  _ \ / _` / __|/ _ \
 |  | |_| | |_| | | | | | |_| | (__ ___) | |_) | | | | | |  __/ |_) | (_| \__ \  __/
 |   \__\_\\__,_|_|_| |_|\__|_|\___|____/| .__/|_|_|_| |_|\___|____/ \__,_|___/\___|
 |                                       |_|
\*/

#ifndef SPLINE_QUINTIC_BASE_HH
#define SPLINE_QUINTIC_BASE_HH

namespace Splines
{

  //!
  //! Quintic spline base class
  //!
  class QuinticSplineBase : public Spline
  {
  protected:
    Malloc_real     m_base_quintic;
    real_type *     m_Yp             = nullptr;
    real_type *     m_Ypp            = nullptr;
    bool            m_external_alloc = false;
    Spline_sub_type m_sub_type;

  public:
    //!
    //! \name Constructors
    //!
    ///@{

    using Spline::build;

    //!
    //! Build an empty spline of `QuinticSplineBase` type
    //!
    //! \param name the name of the spline
    //!
    explicit QuinticSplineBase( Spline_sub_type type = Spline_sub_type::PCHIP, string_view name = "QuinticSplineBase" )
      : Spline( name ), m_base_quintic( fmt::format( "QuinticSplineBase[{}]", name ) ), m_sub_type( type )
    {
    }

    //!
    //! Spline destructor.
    //!
    ~QuinticSplineBase() override {}

    ///@}

    //!
    //! Build a copy of spline `S`
    //!
    void copy_spline( QuinticSplineBase const & S )
    {
      QuinticSplineBase::reserve( S.m_npts );
      m_npts = S.m_npts;
      if ( m_npts > 0 )
      {
        size_t const chunk_size = m_npts * sizeof( *m_X );
        std::memcpy( m_X, S.m_X, chunk_size );
        std::memcpy( m_Y, S.m_Y, chunk_size );
        std::memcpy( m_Yp, S.m_Yp, chunk_size );
        std::memcpy( m_Ypp, S.m_Ypp, chunk_size );
      }
      copy_flags( S );
    }

    //!
    //! \name Info
    //!
    ///@{

    //!
    //! Return the pointer of values of yp-nodes.
    //!
    real_type const * yp_nodes() const { return m_Yp; }

    //!
    //! Return the pointer of values of ypp-nodes.
    //!
    real_type const * ypp_nodes() const { return m_Ypp; }

    //!
    //! Return the i-th node of the spline (y' component).
    //!
    real_type yp_node( integer const i ) const { return m_Yp[i]; }

    //!
    //! Return the i-th node of the spline (y'' component).
    //!
    real_type ypp_node( integer const i ) const { return m_Ypp[i]; }

    void y_min_max(
      integer &   i_min_pos,
      real_type & x_min_pos,
      real_type & y_min,
      integer &   i_max_pos,
      real_type & x_max_pos,
      real_type & y_max ) const override
    {
      UTILS_ASSERT( m_npts > 0, "QuinticSplineBase[{}]::y_min_max() empty spline!", m_name );
      // find max min along the nodes
      i_min_pos = i_max_pos = 0;
      x_min_pos = x_max_pos = m_X[0];
      y_min = y_max = m_Y[0];
      PolynomialRoots::Quartic q;
      for ( integer i = 1; i < m_npts; ++i )
      {
        real_type const & X0   = m_X[i - 1];
        real_type const & X1   = m_X[i];
        real_type const & P0   = m_Y[i - 1];
        real_type const & P1   = m_Y[i];
        real_type const & DP0  = m_Yp[i - 1];
        real_type const & DP1  = m_Yp[i];
        real_type const & DDP0 = m_Ypp[i - 1];
        real_type const & DDP1 = m_Ypp[i];
        real_type const   H    = X1 - X0;
        real_type         A, B, C, D, E, F;
        Hermite5_to_poly( H, P0, P1, DP0, DP1, DDP0, DDP1, A, B, C, D, E, F );
        q.setup( 5 * A, 4 * B, 3 * C, 2 * D, E );
        real_type     r[4];
        integer const nr = q.getRootsInOpenRange( 0, H, r );
        for ( integer j = 0; j < nr; ++j )
        {
          real_type const rr = r[j];
          real_type       yy = ( ( ( ( ( A * rr ) + B ) * rr + C ) * rr + D ) * rr + E ) * rr + F;
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
      std::vector<integer> &   i_min_pos,
      std::vector<real_type> & x_min_pos,
      std::vector<real_type> & y_min,
      std::vector<integer> &   i_max_pos,
      std::vector<real_type> & x_max_pos,
      std::vector<real_type> & y_max ) const override
    {
      constexpr real_type epsi = 1e-8;
      i_min_pos.clear();
      i_max_pos.clear();
      x_min_pos.clear();
      x_max_pos.clear();
      y_min.clear();
      y_max.clear();
      UTILS_ASSERT( m_npts > 0, "QuinticSplineBase[{}]::y_min_max() empty spline!", m_name );
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
      PolynomialRoots::Quartic q;
      for ( integer i = 1; i < m_npts; ++i )
      {
        real_type const & X0   = m_X[i - 1];
        real_type const & X1   = m_X[i];
        real_type const & P0   = m_Y[i - 1];
        real_type const & P1   = m_Y[i];
        real_type const & DP0  = m_Yp[i - 1];
        real_type const & DP1  = m_Yp[i];
        real_type const & DDP0 = m_Ypp[i - 1];
        real_type const & DDP1 = m_Ypp[i];
        real_type const   H    = X1 - X0;
        real_type         A, B, C, D, E, F;
        Hermite5_to_poly( H, P0, P1, DP0, DP1, DDP0, DDP1, A, B, C, D, E, F );
        q.setup( 5 * A, 4 * B, 3 * C, 2 * D, E );
        real_type     r[4];
        integer const nr = q.get_roots_in_open_range( 0, H, r );
        for ( integer j = 0; j < nr; ++j )
        {
          real_type const rr  = r[j];
          real_type       yy  = ( ( ( ( ( A * rr ) + B ) * rr + C ) * rr + D ) * rr + E ) * rr + F;
          real_type       ddy = ( ( ( 20 * A * rr ) + 12 * B ) * rr + 6 * C ) * rr + 2 * D;
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
        if ( std::abs( DP1 ) > ( m_X[i + 1] - m_X[i - 1] ) * epsi ) continue;
        real_type const & X2   = m_X[i + 1];
        real_type const & P2   = m_Y[i + 1];
        real_type const & DP2  = m_Yp[i + 1];
        real_type const & DDP2 = m_Ypp[i + 1];
        real_type         A1, B1, C1, D1, E1, F1;
        Hermite5_to_poly( X2 - X1, P1, P2, DP1, DP2, DDP1, DDP2, A1, B1, C1, D1, E1, F1 );
        real_type DD = ( ( ( 20 * A * H ) + 12 * B ) * H + 6 * C ) * H + 2 * D;
        if ( DD >= 0 && D1 >= 0 )
        {
          y_min.emplace_back( P1 );
          x_min_pos.emplace_back( X1 );
          i_min_pos.emplace_back( i );
        }
        else if ( DD <= 0 && D1 <= 0 )
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

    void write_to_stream( std::ostream & s ) const override
    {
      integer const nseg = m_npts > 0 ? m_npts - 1 : 0;
      for ( integer i = 0; i < nseg; ++i )
        fmt::print(
          s,
          "segment N.{:4} X:[{:.5},{:.5}] Y:[{:.5},{:.5}] Yp:[{:.5},{:.5}] Ypp:[{:.5},{:.5}] slope: {:.5}\n",
          i,
          m_X[i],
          m_X[i + 1],
          m_Y[i],
          m_Y[i + 1],
          m_Yp[i],
          m_Yp[i + 1],
          m_Ypp[i],
          m_Ypp[i + 1],
          ( m_Y[i + 1] - m_Y[i] ) / ( m_X[i + 1] - m_X[i] ) );
    }

    SplineType1D type() const override
    {
      switch ( m_sub_type )
      {
        case Spline_sub_type::CUBIC: return SplineType1D::QUINTIC_CUBIC;
        case Spline_sub_type::AKIMA: return SplineType1D::QUINTIC_AKIMA;
        case Spline_sub_type::VANLEER: return SplineType1D::QUINTIC_VANLEER;
        case Spline_sub_type::PCHIP: return SplineType1D::QUINTIC_PCHIP;
      }
    }

    ///@}

    //!
    //! Change X-range of the spline
    //!
    void set_range( real_type const xmin, real_type const xmax )
    {
      UTILS_ASSERT( m_npts > 0, "QuinticSplineBase[{}]::set_range, empty spline", m_name );
      real_type const L    = m_X[m_npts - 1] - m_X[0];
      real_type const newL = xmax - xmin;
      for ( integer i = 0; i < m_npts; ++i )
      {
        real_type const t = ( m_X[i] - m_X[0] ) / L;
        m_X[i]            = xmin + t * newL;
      }
    }

    //!
    //! Use externally allocated memory for `npts` points
    //!
    void reserve_external( integer const n, real_type *& p_x, real_type *& p_y, real_type *& p_Yp, real_type *& p_Ypp )
    {
      if ( !m_external_alloc ) m_base_quintic.free();
      m_npts           = 0;
      m_npts_reserved  = n;
      m_external_alloc = true;
      m_X              = p_x;
      m_Y              = p_y;
      m_Yp             = p_Yp;
      m_Ypp            = p_Ypp;
    }

    // --------------------------- VIRTUALS -----------------------------------

    //!
    //! \name Evaluation Aliases
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

    real_type DDDD( real_type const x ) const override
    {
      std::pair<integer, real_type> res( 0, x );
      m_search.find( res );
      return this->id_DDDD( res.first, res.second );
    }

    real_type DDDDD( real_type const x ) const override
    {
      std::pair<integer, real_type> res( 0, x );
      m_search.find( res );
      return this->id_DDDDD( res.first, res.second );
    }

    void D( real_type const x, real_type dd[2] ) const override
    {
      std::pair<integer, real_type> res( 0, x );
      m_search.find( res );
      integer const   ni = res.first;
      real_type const X  = res.second;
      real_type       dx = X - m_X[ni];
      real_type       DX = m_X[ni + 1] - m_X[ni];
      real_type       base[6], base_D[6];
      Hermite5( dx, DX, base );
      Hermite5_D( dx, DX, base_D );

      dd[0] = base[0] * m_Y[ni] + base[1] * m_Y[ni + 1] + base[2] * m_Yp[ni] + base[3] * m_Yp[ni + 1] +
              base[4] * m_Ypp[ni] + base[5] * m_Ypp[ni + 1];

      dd[1] = base_D[0] * m_Y[ni] + base_D[1] * m_Y[ni + 1] + base_D[2] * m_Yp[ni] + base_D[3] * m_Yp[ni + 1] +
              base_D[4] * m_Ypp[ni] + base_D[5] * m_Ypp[ni + 1];
    }

    void DD( real_type const x, real_type dd[3] ) const override
    {
      std::pair<integer, real_type> res( 0, x );
      m_search.find( res );
      integer const   ni = res.first;
      real_type const X  = res.second;
      real_type       dx = X - m_X[ni];
      real_type       DX = m_X[ni + 1] - m_X[ni];
      real_type       base[6], base_D[6], base_DD[6];
      Hermite5( dx, DX, base );
      Hermite5_D( dx, DX, base_D );
      Hermite5_DD( dx, DX, base_DD );

      dd[0] = base[0] * m_Y[ni] + base[1] * m_Y[ni + 1] + base[2] * m_Yp[ni] + base[3] * m_Yp[ni + 1] +
              base[4] * m_Ypp[ni] + base[5] * m_Ypp[ni + 1];

      dd[1] = base_D[0] * m_Y[ni] + base_D[1] * m_Y[ni + 1] + base_D[2] * m_Yp[ni] + base_D[3] * m_Yp[ni + 1] +
              base_D[4] * m_Ypp[ni] + base_D[5] * m_Ypp[ni + 1];

      dd[2] = base_DD[0] * m_Y[ni] + base_DD[1] * m_Y[ni + 1] + base_DD[2] * m_Yp[ni] + base_DD[3] * m_Yp[ni + 1] +
              base_DD[4] * m_Ypp[ni] + base_DD[5] * m_Ypp[ni + 1];
    }
    ///@}

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
      xx.val  = res.second;  // Potrebbe essere clampato se x.val è fuori dall'intervallo
      xx.grad = x.grad;

      autodiff::dual1st base[6];
      Hermite5<autodiff::dual1st>( xx - m_X[ni], m_X[ni1] - m_X[ni], base );

      return base[0] * m_Y[ni] + base[1] * m_Y[ni1] + base[2] * m_Yp[ni] + base[3] * m_Yp[ni1] + base[4] * m_Ypp[ni] +
             base[5] * m_Ypp[ni1];
    }

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
      xx.val.val   = res.second;  // Potrebbe essere clampato se x.val.val è fuori dall'intervallo
      xx.val.grad  = x.val.grad;
      xx.grad.val  = x.grad.val;
      xx.grad.grad = x.grad.grad;

      autodiff::dual2nd base[6];
      Hermite5<autodiff::dual2nd>( xx - m_X[ni], m_X[ni1] - m_X[ni], base );

      return base[0] * m_Y[ni] + base[1] * m_Y[ni1] + base[2] * m_Yp[ni] + base[3] * m_Yp[ni1] + base[4] * m_Ypp[ni] +
             base[5] * m_Ypp[ni1];
    }
#endif

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
      real_type base[6];
      integer   ni1 = ni + 1;
      Hermite5( x - m_X[ni], m_X[ni1] - m_X[ni], base );
      return base[0] * m_Y[ni] + base[1] * m_Y[ni1] + base[2] * m_Yp[ni] + base[3] * m_Yp[ni1] + base[4] * m_Ypp[ni] +
             base[5] * m_Ypp[ni1];
    }

    real_type id_D( integer const ni, real_type const x ) const override
    {
      if ( m_curve_can_extend && m_curve_extended_constant )
      {
        // estendo solo quando esco effettivamente
        if ( x < m_X[0] || x > m_X[m_npts - 1] ) return 0;
      }
      real_type base_D[6];
      integer   ni1 = ni + 1;
      Hermite5_D( x - m_X[ni], m_X[ni1] - m_X[ni], base_D );
      return base_D[0] * m_Y[ni] + base_D[1] * m_Y[ni1] + base_D[2] * m_Yp[ni] + base_D[3] * m_Yp[ni1] +
             base_D[4] * m_Ypp[ni] + base_D[5] * m_Ypp[ni1];
    }

    real_type id_DD( integer const ni, real_type const x ) const override
    {
      if ( m_curve_can_extend && m_curve_extended_constant )
      {
        // estendo solo quando esco effettivamente
        if ( x < m_X[0] || x > m_X[m_npts - 1] ) return 0;
      }
      real_type base_DD[6];
      integer   ni1 = ni + 1;
      Hermite5_DD( x - m_X[ni], m_X[ni1] - m_X[ni], base_DD );
      return base_DD[0] * m_Y[ni] + base_DD[1] * m_Y[ni1] + base_DD[2] * m_Yp[ni] + base_DD[3] * m_Yp[ni1] +
             base_DD[4] * m_Ypp[ni] + base_DD[5] * m_Ypp[ni1];
    }

    real_type id_DDD( integer const ni, real_type const x ) const override
    {
      if ( m_curve_can_extend && m_curve_extended_constant )
      {
        if ( x < m_X[0] || x > m_X[m_npts - 1] ) return 0;
      }
      real_type base_DDD[6];
      integer   ni1 = ni + 1;
      Hermite5_DDD( x - m_X[ni], m_X[ni1] - m_X[ni], base_DDD );
      return base_DDD[0] * m_Y[ni] + base_DDD[1] * m_Y[ni1] + base_DDD[2] * m_Yp[ni] + base_DDD[3] * m_Yp[ni1] +
             base_DDD[4] * m_Ypp[ni] + base_DDD[5] * m_Ypp[ni1];
    }

    real_type id_DDDD( integer const ni, real_type const x ) const override
    {
      if ( m_curve_can_extend && m_curve_extended_constant )
      {
        // estendo solo quando esco effettivamente
        if ( x < m_X[0] || x > m_X[m_npts - 1] ) return 0;
      }
      real_type base_DDDD[6];
      integer   ni1 = ni + 1;
      Hermite5_DDDD( x - m_X[ni], m_X[ni1] - m_X[ni], base_DDDD );
      return base_DDDD[0] * m_Y[ni] + base_DDDD[1] * m_Y[ni1] + base_DDDD[2] * m_Yp[ni] + base_DDDD[3] * m_Yp[ni1] +
             base_DDDD[4] * m_Ypp[ni] + base_DDDD[5] * m_Ypp[ni1];
    }

    real_type id_DDDDD( integer const ni, real_type const x ) const override
    {
      if ( m_curve_can_extend && m_curve_extended_constant )
      {
        // estendo solo quando esco effettivamente
        if ( x < m_X[0] || x > m_X[m_npts - 1] ) return 0;
      }
      real_type base_DDDDD[6];
      integer   ni1 = ni + 1;
      Hermite5_DDDDD( x - m_X[ni], m_X[ni1] - m_X[ni], base_DDDDD );
      return base_DDDDD[0] * m_Y[ni] + base_DDDDD[1] * m_Y[ni1] + base_DDDDD[2] * m_Yp[ni] + base_DDDDD[3] * m_Yp[ni1] +
             base_DDDDD[4] * m_Ypp[ni] + base_DDDDD[5] * m_Ypp[ni1];
    }
    ///@}

    void reserve( integer npts ) override
    {
      if ( m_external_alloc && npts <= m_npts_reserved )
      {
        // nothing to do!, already allocated
      }
      else
      {
        m_base_quintic.reallocate( 4 * npts );
        m_npts_reserved  = npts;
        m_external_alloc = false;
        m_X              = m_base_quintic( npts );
        m_Y              = m_base_quintic( npts );
        m_Yp             = m_base_quintic( npts );
        m_Ypp            = m_base_quintic( npts );
      }
      m_npts = 0;
    }

    void clear() override
    {
      if ( !m_external_alloc ) m_base_quintic.free();
      m_npts = m_npts_reserved = 0;
      m_external_alloc         = false;
      m_X = m_Y = m_Yp = m_Ypp = nullptr;
    }

    //!
    //! Get the piecewise polinomials of the spline
    //!
    integer  // order
    coeffs( real_type cfs[], real_type nodes[], bool transpose = false ) const override
    {
      UTILS_ASSERT( m_npts >= 2, "QuinticSplineBase, npts={} must be >= 2\n", m_npts );

      integer const n = m_npts - 1;

      for ( integer i = 0; i < n; ++i )
      {
        real_type const H = m_X[i + 1] - m_X[i];
        real_type       a, b, c, d, e, f;
        Hermite5_to_poly( H, m_Y[i], m_Y[i + 1], m_Yp[i], m_Yp[i + 1], m_Ypp[i], m_Ypp[i + 1], a, b, c, d, e, f );
        if ( transpose )
        {
          cfs[6 * i + 5] = a;
          cfs[6 * i + 4] = b;
          cfs[6 * i + 3] = c;
          cfs[6 * i + 2] = d;
          cfs[6 * i + 1] = e;
          cfs[6 * i + 0] = f;
        }
        else
        {
          cfs[i + 5 * n] = a;
          cfs[i + 4 * n] = b;
          cfs[i + 3 * n] = c;
          cfs[i + 2 * n] = d;
          cfs[i + 1 * n] = e;
          cfs[i + 0 * n] = f;
        }
      }
      if ( m_npts > 0 ) std::memcpy( nodes, m_X, m_npts * sizeof( *nodes ) );
      return 6;
    }

    integer order() const override { return 6; }

    bool is_monotone() const { return check_quintic_spline_monotonicity( m_X, m_Y, m_Yp, m_Ypp, m_npts ); }

#ifdef SPLINES_BACK_COMPATIBILITY
    void      copySpline( QuinticSplineBase const & S ) { this->copy_spline( S ); }
    real_type ypNode( integer i ) const { return this->yp_node( i ); }
    real_type yppNode( integer i ) const { return this->ypp_node( i ); }
    void      setRange( real_type xmin, real_type xmax ) { this->set_range( xmin, xmax ); }
#endif
  };

}  // namespace Splines

#endif

// EOF: SplineQuinticBase.hxx
