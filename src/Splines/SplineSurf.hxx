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

#pragma once

#ifndef SPLINE_SURF_HH
#define SPLINE_SURF_HH

namespace Splines
{

  /*\
   |   ____        _ _            ____              __
   |  / ___| _ __ | (_)_ __   ___/ ___| _   _ _ __ / _|
   |  \___ \| '_ \| | | '_ \ / _ \___ \| | | | '__| |_
   |   ___) | |_) | | | | | |  __/___) | |_| | |  |  _|
   |  |____/| .__/|_|_|_| |_|\___|____/ \__,_|_|  |_|
   |        |_|
  \*/

  //!
  //! Spline Management Class
  //!
  class SplineSurf
  {
    Malloc_real m_mem;

  protected:
    string const m_name;
    bool         m_x_closed     = false;
    bool         m_y_closed     = false;
    bool         m_x_can_extend = true;
    bool         m_y_can_extend = true;

    integer m_nx = 0;
    integer m_ny = 0;

    real_type * m_X     = nullptr;
    real_type * m_Y     = nullptr;
    real_type * m_Z_ptr = nullptr;
    Eigen::Map<MatC> mZ{ nullptr, 0, 0 };

    real_type m_Z_min = 0;
    real_type m_Z_max = 0;

    Utils::SearchInterval<real_type, integer> m_search_x;
    Utils::SearchInterval<real_type, integer> m_search_y;

    real_type & z_node_ref( integer const i, integer const j ) { return mZ.coeffRef(i,j); }

    template <typename MAT>
    void
    load_Z( MAT const & Z, bool transposed ) {
      if ( transposed ) {
        for ( integer ix = 0; ix < m_nx; ++ix )
          for ( integer iy = 0; iy < m_ny; ++iy )
            z_node_ref( ix, iy ) = Z( iy, ix );
      } else {
        for ( integer ix = 0; ix < m_nx; ++ix )
          for ( integer iy = 0; iy < m_ny; ++iy )
            z_node_ref( ix, iy ) = Z( ix, iy );
      }
      m_Z_max = Z.maxCoeff();
      m_Z_min = Z.minCoeff();
    }

    void load_Z( real_type const z[], integer const ldZ, bool fortran_storage, bool transposed )
    {
      //
      //  +--------------+
      //  |              | ny = nr
      //  |              |
      //  +--------------+
      //      nx = nc
      //
      integer const tf{ ( transposed ? 1 : 0 ) + ( fortran_storage ? 2 : 0 ) };
      switch ( tf )
      {
        case 0:  // NO transpose NO fortran
        {
          Eigen::Map<const MatC> ZZ( z, m_nx, ldZ );
          load_Z( ZZ, false );
        }
        break;

        case 1:  // YES transpose NO fortran
        {
          Eigen::Map<const MatC> ZZ( z, m_ny, ldZ );
          load_Z( ZZ, true );
        }
        break;

        case 2:  // NO transpose YES fortran
        {
          Eigen::Map<const Mat> ZZ( z, ldZ, m_ny );
          load_Z( ZZ, false );
        }
        break;

        case 3:  // YES transpose YES fortran
        {
          Eigen::Map<const Mat> ZZ( z, ldZ, m_nx );
          load_Z( ZZ, true );
        }
        break;
      }
    }

    virtual void make_spline() = 0;

    void make_derivative_x( Spline_sub_type sub, real_type const z[], real_type dx[] )
    {
      auto interpolate = [this]( CubicSplineBase * S, real_type const z[], real_type dx[] )
      {
        for ( integer j = 0; j < m_ny; ++j )
        {
          S->build( m_X, 1, z + j, m_ny, m_nx );
          for ( integer i = 0; i < m_nx; ++i ) dx[i * m_ny + j] = S->yp_node( i );
        }
      };
      CubicSpline  cs;
      AkimaSpline  ak;
      BesselSpline be;
      PchipSpline  pc;
      switch ( sub )
      {
        case Spline_sub_type::CUBIC: interpolate( &cs, z, dx ); break;
        case Spline_sub_type::AKIMA: interpolate( &ak, z, dx ); break;
        case Spline_sub_type::BESSEL: interpolate( &be, z, dx ); break;
        case Spline_sub_type::PCHIP: interpolate( &pc, z, dx ); break;
      }
    }

    void make_derivative_y( Spline_sub_type sub, real_type const z[], real_type dy[] )
    {
      auto interpolate = [this]( CubicSplineBase * S, real_type const z[], real_type dy[] )
      {
        for ( integer i = 0; i < m_nx; ++i )
        {
          S->build( m_Y, 1, z + i * m_ny, 1, m_ny );
          for ( integer j = 0; j < m_ny; ++j ) dy[i * m_ny + j] = S->yp_node( j );
        }
      };
      CubicSpline  cs;
      AkimaSpline  ak;
      BesselSpline be;
      PchipSpline  pc;
      switch ( sub )
      {
        case Spline_sub_type::CUBIC: interpolate( &cs, z, dy ); break;
        case Spline_sub_type::AKIMA: interpolate( &ak, z, dy ); break;
        case Spline_sub_type::BESSEL: interpolate( &be, z, dy ); break;
        case Spline_sub_type::PCHIP: interpolate( &pc, z, dy ); break;
      }
    }

    void make_derivative_xy( Spline_sub_type sub, real_type const dx[], real_type const dy[], real_type dxy[] )
    {
      auto interpolate = [this]( CubicSplineBase * S, real_type const dx[], real_type const dy[], real_type dxy[] )
      {
        auto minmod = []( real_type a, real_type b ) -> real_type
        {
          if ( a * b <= 0 ) return 0;
          if ( a > 0 ) return std::min( a, b );
          return std::max( a, b );
        };
        for ( integer j = 0; j < m_ny; ++j )
        {
          S->build( m_X, 1, dy + j, m_ny, m_nx );
          for ( integer i = 0; i < m_nx; ++i ) dxy[i * m_ny + j] = S->yp_node( i );
        }

        for ( integer i = 0; i < m_nx; ++i )
        {
          S->build( m_Y, 1, dx + i * m_ny, 1, m_ny );
          for ( integer j = 0; j < m_ny; ++j )
          {
            integer const ij = i * m_ny + j;
            dxy[ij] = minmod( dxy[ij], S->yp_node( j ) );
          }
        }
      };
      CubicSpline  cs;
      AkimaSpline  ak;
      BesselSpline be;
      PchipSpline  pc;
      switch ( sub )
      {
        case Spline_sub_type::CUBIC: interpolate( &cs, dx, dy, dxy ); break;
        case Spline_sub_type::AKIMA: interpolate( &ak, dx, dy, dxy ); break;
        case Spline_sub_type::BESSEL: interpolate( &be, dx, dy, dxy ); break;
        case Spline_sub_type::PCHIP: interpolate( &pc, dx, dy, dxy ); break;
      }
    }

  public:
    SplineSurf( SplineSurf const & )                   = delete;  // block copy constructor
    SplineSurf const & operator=( SplineSurf const & ) = delete;  // block copy method

    //!
    //! Spline constructor
    //!
    explicit SplineSurf( string_view name = "Spline" ) : m_mem( name.data() ), m_name( name )
    {
      m_search_x.setup( &m_name, &m_nx, &m_X, &m_x_closed, &m_x_can_extend );
      m_search_y.setup( &m_name, &m_ny, &m_Y, &m_y_closed, &m_y_can_extend );
    }

    //!
    //! Spline destructor
    //!
    virtual ~SplineSurf() {}

    //!
    //! \name Open/Close
    //!
    ///@{

    //!
    //! Return `true` if the surface is assumed closed in the `x` direction.
    //!
    bool is_x_closed() const { return m_x_closed; }

    //!
    //! Setup the surface as closed in the `x` direction.
    //!
    void make_x_closed() { m_x_closed = true; }

    //!
    //! Setup the surface as open in the `x` direction.
    //!
    void make_x_opened() { m_x_closed = false; }

    //!
    //! Return `true` if the surface is assumed closed in the `y` direction.
    //!
    bool is_y_closed() const { return m_y_closed; }

    //!
    //! Setup the surface as closed in the `y` direction.
    //!
    void make_y_closed() { m_y_closed = true; }

    //!
    //! Setup the surface as open in the `y` direction.
    //!
    void make_y_opened() { m_y_closed = false; }

    //!
    //! Return `true` if the parameter `x` assumed bounded.
    //! If false the spline is estrapolated for `x` values
    //! outside the range.
    //!
    bool is_x_bounded() const { return m_x_can_extend; }

    //!
    //! Make the spline surface unbounded in the `x` direction.
    //!
    void make_x_unbounded() { m_x_can_extend = true; }

    //!
    //! Make the spline surface bounded in the `x` direction.
    //!
    void make_x_bounded() { m_x_can_extend = false; }

    //!
    //! Return `true` if the parameter `y` assumed bounded.
    //! If false the spline is extrapolated for `y` values
    //! outside the range.
    //!
    bool is_y_bounded() const { return m_y_can_extend; }

    //!
    //! Make the spline surface unbounded in the `y` direction
    //!
    void make_y_unbounded() { m_y_can_extend = true; }

    //!
    //! Make the spline surface bounded in the `x` direction.
    //!
    void make_y_bounded() { m_y_can_extend = false; }

    ///@}

    //!
    //! Cancel the support points, empty the spline.
    //!
    void clear()
    {
      m_mem.free();
      m_nx = 0;
      m_ny = 0;

      m_X = nullptr;
      m_Y = nullptr;
      m_Z_ptr = nullptr;
      new (&mZ) Eigen::Map<MatC>{ nullptr, 0, 0 };

      m_Z_min = 0;
      m_Z_max = 0;
    }

    //!
    //! \name Info
    //!
    ///@{

    //!
    //! \return string with the name of the spline
    //!
    string_view name() const { return m_name; }

    //!
    //! Return the number of support points of the spline along x direction.
    //!
    integer num_point_x() const { return m_nx; }

    //!
    //! Return the number of support points of the spline along y direction.
    //!
    integer num_point_y() const { return m_ny; }

    //!
    //! Return the i-th node of the spline (x component).
    //!
    real_type x_node( integer const i ) const { return m_X[i]; }

    //!
    //! Return the i-th node of the spline (y component).
    //!
    real_type y_node( integer const i ) const { return m_Y[i]; }

    //!
    //! Return the i-th node of the spline (y component).
    //!
    real_type z_node( integer const i, integer const j ) const { return mZ.coeff(i,j); }

    //!
    //! Return x-minumum spline value.
    //!
    real_type x_min() const { return m_X[0]; }

    //!
    //! Return x-maximum spline value.
    //!
    real_type x_max() const { return m_X[m_nx - 1]; }

    //!
    //! Return y-minumum spline value.
    //!
    real_type y_min() const { return m_Y[0]; }

    //!
    //! Return y-maximum spline value.
    //!
    real_type y_max() const { return m_Y[m_ny - 1]; }

    //!
    //! Return z-minumum spline value.
    //!
    real_type z_min() const { return m_Z_min; }

    //!
    //! Return z-maximum spline value.
    //!
    real_type z_max() const { return m_Z_max; }

    ///@}

    //!
    //! \name Build Spline
    //!
    ///@{

    //!
    //! Build surface spline
    //!
    //! \param x               vector of `x`-coordinates
    //! \param incx            access elements as `x[0]`, `x[incx]`, `x[2*incx]`,...
    //! \param y               vector of `y`-coordinates
    //! \param incy            access elements as `y[0]`, `y[incy]`, `y[2*incy]`,...
    //! \param z               matrix of `z`-values. Elements are stored
    //!                        by row Z(i,j) = z[i*ny+j] as C-matrix
    //! \param ldZ             leading dimension of `z`
    //! \param nx              number of points in `x` direction
    //! \param ny              number of points in `y` direction
    //! \param fortran_storage if true elements are stored by column
    //!                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
    //! \param transposed      if true matrix Z is stored transposed
    //!
    void build(
      real_type const x[],
      integer const   incx,
      real_type const y[],
      integer const   incy,
      real_type const z[],
      integer const   ldZ,
      integer const   nx,
      integer const   ny,
      bool            fortran_storage = false,
      bool            transposed      = false )
    {
      m_nx = nx;
      m_ny = ny;
      m_mem.reallocate( ( nx + 1 ) * ( ny + 1 ) );
      m_X = m_mem( nx );
      m_Y = m_mem( ny );

      m_Z_ptr = m_mem( nx * ny );
      new (&mZ) Eigen::Map<MatC>( m_Z_ptr, nx, ny );

      for ( integer i = 0; i < nx; ++i ) m_X[i] = x[i * incx];
      for ( integer j = 0; j < ny; ++j ) m_Y[j] = y[j * incy];
      load_Z( z, ldZ, fortran_storage, transposed );
      make_spline();
    }

    //!
    //! Build surface spline
    //!
    //! \param x               vector of x-coordinates, nx = x.size()
    //! \param y               vector of y-coordinates, ny = y.size()
    //! \param Z               matrix of z-values. Elements are stored
    //!                        by row Z(i,j) = z[i*ny+j] as C-matrix
    //! \param transposed      if true matrix Z is stored transposed
    //!
    void build(
      Eigen::Ref<const Vec> x,
      Eigen::Ref<const Vec> y,
      Eigen::Ref<const Mat> Z,
      bool const transposed )
    {
      load_Z( Z, transposed );
      make_spline();
    }

    //!
    //! Build surface spline
    //!
    //! \param x               vector of x-coordinates, nx = x.size()
    //! \param y               vector of y-coordinates, ny = y.size()
    //! \param z               matrix of z-values. Elements are stored
    //!                        by row Z(i,j) = z[i*ny+j] as C-matrix
    //! \param fortran_storage if true elements are stored by column
    //!                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
    //! \param transposed      if true matrix Z is stored transposed
    //!
    void build(
      vector<real_type> const & x,
      vector<real_type> const & y,
      vector<real_type> const & z,
      bool                      fortran_storage = false,
      bool                      transposed      = false )
    {
      this->build(
        x.data(),
        1,
        y.data(),
        1,
        z.data(),
        integer( fortran_storage ? y.size() : x.size() ),
        integer( x.size() ),
        integer( y.size() ),
        fortran_storage,
        transposed );
    }

    void build(
      real_type const z[],
      integer         ldZ,
      integer         nx,
      integer         ny,
      bool            fortran_storage = false,
      bool            transposed      = false )
    {
      m_nx = nx;
      m_ny = ny;
      m_mem.reallocate( ( nx + 1 ) * ( ny + 1 ) );
      m_X = m_mem( nx );
      m_Y = m_mem( ny );

      m_Z_ptr = m_mem( nx * ny );
      new (&mZ) Eigen::Map<MatC>( m_Z_ptr, nx, ny );

      for ( integer i = 0; i < nx; ++i ) m_X[i] = static_cast<real_type>( i );
      for ( integer j = 0; j < ny; ++j ) m_Y[j] = static_cast<real_type>( j );
      load_Z( z, ldZ, fortran_storage, transposed );
      make_spline();
    }

    //!
    //! Build surface spline
    //!
    //! \param z               matrix of z-values. Elements are stored
    //!                        by row Z(i,j) = z[i*ny+j] as C-matrix.
    //!                        ldZ leading dimension of the matrix is ny for C-storage
    //!                        and nx for Fortran storage.
    //! \param nx              x-dimension
    //! \param ny              y-dimension
    //! \param fortran_storage if true elements are stored by column
    //!                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
    //! \param transposed      if true matrix Z is stored transposed
    //!
    void build(
      vector<real_type> const & z,
      integer const             nx,
      integer const             ny,
      bool                      fortran_storage = false,
      bool                      transposed      = false )
    {
      this->build( z.data(), nx, ny, fortran_storage ? nx : ny, fortran_storage, transposed );
    }

    //!
    //! Build spline using data in `gc`
    //!
    void setup( GenericContainer const & gc )
    {
      /*
      // gc["xdata"]
      // gc["ydata"]
      // gc["zdata"]
      //
      */
      string const where{ fmt::format( "SplineSurf[{}]::setup( gc ):", m_name ) };

      std::set<std::string> keywords;
      for ( auto const & pair : gc.get_map( where ) ) { keywords.insert( pair.first ); }

      GenericContainer const & gc_x{ gc( "xdata", where ) };
      keywords.erase( "xdata" );
      GenericContainer const & gc_y{ gc( "ydata", where ) };
      keywords.erase( "ydata" );
      GenericContainer const & gc_z{ gc( "zdata", where ) };
      keywords.erase( "zdata" );

      m_nx = static_cast<integer>( gc_x.get_num_elements() );
      m_ny = static_cast<integer>( gc_y.get_num_elements() );
      m_mem.reallocate( ( m_nx + 1 ) * ( m_ny + 1 ) );
      m_X = m_mem( m_nx );
      m_Y = m_mem( m_ny );

      m_Z_ptr = m_mem( m_nx * m_ny );
      new (&mZ) Eigen::Map<MatC>( m_Z_ptr, m_nx, m_ny );

      for ( integer i = 0; i < m_nx; ++i ) m_X[i] = gc_x.get_number_at( i );
      for ( integer j = 0; j < m_ny; ++j ) m_Y[j] = gc_y.get_number_at( j );

      bool fortran_storage{ gc.get_map_bool( "fortran_storage", where ) };
      keywords.erase( "fortran_storage" );
      bool transposed{ gc.get_map_bool( "transposed", where ) };
      keywords.erase( "transposed" );

      /*
      //     +------+
      //  j ny      | (xi,yj)
      //     +  nx  +
      //         i
      */

      // cosa mi aspetto in lettura
      integer NR, NC;
      if ( transposed )
      {
        NC = m_ny;
        NR = m_nx;
      }
      else
      {
        NC = m_nx;
        NR = m_ny;
      }
      integer const LD{ fortran_storage ? NR : NC };

      if (
        GC_type::MAT_REAL == gc_z.get_type() || GC_type::MAT_INTEGER == gc_z.get_type() ||
        GC_type::MAT_LONG == gc_z.get_type() )
      {
        UTILS_ASSERT(
          static_cast<unsigned>( NR ) == gc_z.num_rows() && static_cast<unsigned>( NC ) == gc_z.num_cols(),
          "{}, field `zdata` is a matrix expected to be of size {} x {}, found: {} x {}\n",
          where,
          NR,
          NC,
          gc_z.num_rows(),
          gc_z.num_cols() );

        if ( GC_type::MAT_REAL == gc_z.get_type() )
        {
          load_Z( gc_z.get_mat_real().data(), LD, fortran_storage, transposed );
        }
        else
        {
          GenericContainer::mat_real_type z_tmp;
          gc_z.copyto_mat_real( z_tmp );
          load_Z( z_tmp.data(), LD, fortran_storage, transposed );
        }
      }
      else if (
        GC_type::VEC_REAL == gc_z.get_type() || GC_type::VEC_INTEGER == gc_z.get_type() ||
        GC_type::VEC_LONG == gc_z.get_type() )
      {
        integer nz{ static_cast<integer>( gc_z.get_num_elements() ) };
        integer nxy{ m_nx * m_ny };
        UTILS_ASSERT(
          nz == nxy,
          "{}, field `zdata` expected to be of size {} = {}x{}, found: `{}`\n",
          where,
          nxy,
          m_nx,
          m_ny,
          nz );

        for ( integer k = 0; k < nxy; ++k )
        {
          integer         i, j;
          real_type const v{ gc_z.get_number_at( k ) };
          if ( fortran_storage )
          {
            i = k % NR;
            j = k / NR;
          }
          else
          {
            j = k % NC;
            i = k / NC;
          }
          if ( transposed )
            z_node_ref( j, i ) = v;
          else
            z_node_ref( i, j ) = v;
        }
      }
      else if ( GC_type::VECTOR == gc_z.get_type() )
      {
        vector_type const & data{ gc_z.get_vector() };
        vec_real_type       tmp;

        if ( transposed )
        {
          UTILS_ASSERT(
            static_cast<size_t>( NC ) == data.size(),
            "{}, field `zdata` (vector of vector transposed) expected of size {} found of size {}\n",
            where,
            NC,
            data.size() );
          for ( integer j = 0; j < NC; ++j )
          {
            GenericContainer const & col{ data[j] };
            string const             msg1{ fmt::format( "{}, reading row {}\n", where, j ) };
            col.copyto_vec_real( tmp, msg1 );
            UTILS_ASSERT(
              static_cast<size_t>( NR ) == tmp.size(),
              "{}, col {}-th of size {}, expected {}\n",
              where,
              j,
              tmp.size(),
              NR );
            for ( integer i = 0; i < NC; ++i ) z_node_ref( i, j ) = tmp[i];
          }
        }
        else
        {
          UTILS_ASSERT(
            static_cast<size_t>( NR ) == data.size(),
            "{}, field `zdata` (vector of vector) expected of size {} found of size {}\n",
            where,
            NR,
            data.size() );
          for ( integer i = 0; i < NR; ++i )
          {
            GenericContainer const & row{ data[i] };
            string const             msg1{ fmt::format( "{}, reading row {}\n", where, i ) };
            row.copyto_vec_real( tmp, msg1 );
            UTILS_ASSERT(
              static_cast<size_t>( NC ) == tmp.size(),
              "{}, row {}-th of size {}, expected {}\n",
              where,
              i,
              tmp.size(),
              NC );
            for ( integer j = 0; j < NC; ++j ) z_node_ref( i, j ) = tmp[j];
          }
        }
      }
      else
      {
        UTILS_ERROR(
          "{}, field `zdata` expected to be of type"
          " `mat_real_type` or  `vec_real_type` or `vector_type` found: `{}`\n",
          where,
          gc_z.get_type_name() );
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

      make_spline();
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

    //!
    //! Build a spline using data in `GenericContainer`
    //!
    void build( GenericContainer const & gc ) { setup( gc ); }

    //!
    //! Build a spline using data in a file `file_name`
    //!
    void build( string const & file_name ) { setup( file_name ); }

    ///@}

    //!
    //! \name Evaluate
    //!
    ///@{

    //!
    //! Evaluate spline value at point \f$ (x,y) \f$.
    //!
    virtual real_type eval( real_type const x, real_type const y ) const = 0;

#ifdef AUTODIFF_SUPPORT
    // Metodi base per dual1st e dual2nd
    virtual autodiff::dual1st eval( autodiff::dual1st const & x, autodiff::dual1st const & y ) const
    {
      using autodiff::dual1st;
      using autodiff::detail::val;

      real_type dd[3];
      D( val( x ), val( y ), dd );

      dual1st res = dd[0];
      res.grad    = dd[1] * x.grad + dd[2] * y.grad;
      return res;
    }

    virtual autodiff::dual2nd eval( autodiff::dual2nd const & x, autodiff::dual2nd const & y ) const
    {
      using autodiff::derivative;
      using autodiff::dual2nd;

      real_type dd[6];
      real_type dx  = val( x.grad );
      real_type dy  = val( y.grad );
      real_type ddx = x.grad.grad;
      real_type ddy = y.grad.grad;
      DD( val( x ), val( y ), dd );

      dual2nd res   = dd[0];
      res.grad      = dd[1] * dx + dd[2] * dy;
      res.grad.grad = dx * dx * dd[3] + 2 * dx * dy * dd[4] + dy * dy * dd[5] + ddx * dd[1] + ddy * dd[2];
      return res;
    }

    // Template per due parametri (x, y) - per SplineSurf
    // Promuove automaticamente a double, dual1st o dual2nd in base ai tipi di input
    template <typename T1, typename T2> auto eval( T1 const & x, T2 const & y ) const
    {
      if constexpr ( std::is_arithmetic<T1>::value && std::is_arithmetic<T2>::value )
      {
        // Entrambi numerici: ritorna real_type
        return eval( static_cast<real_type>( x ), static_cast<real_type>( y ) );
      }
      else
      {
        // Almeno uno è duale: determina il tipo duale appropriato
        constexpr int order1    = std::is_arithmetic<T1>::value ? 0 : autodiff::detail::DualOrder<T1>::value;
        constexpr int order2    = std::is_arithmetic<T2>::value ? 0 : autodiff::detail::DualOrder<T2>::value;
        constexpr int max_order = ( order1 > order2 ) ? order1 : order2;

        if constexpr ( max_order == 1 )
        {
          // Promuovi a dual1st
          autodiff::dual1st X = std::is_arithmetic<T1>::value ? autodiff::dual1st{ static_cast<real_type>( x ) }
                                                              : autodiff::dual1st{ x };
          autodiff::dual1st Y = std::is_arithmetic<T2>::value ? autodiff::dual1st{ static_cast<real_type>( y ) }
                                                              : autodiff::dual1st{ y };
          return eval( X, Y );
        }
        else
        {
          // Promuovi a dual2nd
          autodiff::dual2nd X = std::is_arithmetic<T1>::value ? autodiff::dual2nd{ static_cast<real_type>( x ) }
                                                              : autodiff::dual2nd{ x };
          autodiff::dual2nd Y = std::is_arithmetic<T2>::value ? autodiff::dual2nd{ static_cast<real_type>( y ) }
                                                              : autodiff::dual2nd{ y };
          return eval( X, Y );
        }
      }
    }

    // Operator() per due parametri
    template <typename T1, typename T2> auto operator()( T1 const & x, T2 const & y ) const -> decltype( eval( x, y ) )
    {
      return this->eval( x, y );
    }
#endif

    //!
    //! Value and first derivatives at point \f$ (x,y) \f$:
    //!
    //! - d[0] value of the spline \f$ S(x,y) \f$
    //! - d[1] derivative respect to \f$ x \f$ of the spline: \f$ S_x(x,y) \f$
    //! - d[2] derivative respect to \f$ y \f$ of the spline: \f$ S_y(x,y) \f$
    //!
    virtual void D( real_type const x, real_type const y, real_type d[3] ) const = 0;

    //!
    //! First derivatives respect to \f$ x \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_x(x,y) \f$.
    //!
    virtual real_type Dx( real_type const x, real_type const y ) const = 0;

    //!
    //! First derivatives respect to \f$ y \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_y(x,y) \f$.
    //!
    virtual real_type Dy( real_type const x, real_type const y ) const = 0;

    //!
    //! Value, first and second derivatives at point \f$ (x,y) \f$:
    //!
    //! - dd[0] value of the spline \f$ S(x,y) \f$
    //! - dd[1] derivative respect to \f$ x \f$ of the spline: \f$ S_x(x,y) \f$
    //! - dd[2] derivative respect to \f$ y \f$ of the spline: \f$ S_y(x,y) \f$
    //! - dd[3] second derivative respect to \f$ x \f$ of the spline: \f$ S_{xx}(x,y) \f$
    //! - dd[4] mixed second derivative: \f$ S_{xy}(x,y) \f$
    //! - dd[5] second derivative respect to \f$ y \f$ of the spline: \f$ S_{yy}(x,y) \f$
    //!
    virtual void DD( real_type const x, real_type const y, real_type dd[6] ) const = 0;

    //!
    //! Second derivatives respect to \f$ x \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_{xx}(x,y) \f$.
    //!
    virtual real_type Dxx( real_type const x, real_type const y ) const = 0;

    //!
    //! Mixed second derivatives: \f$ S_{xy}(x,y) \f$.
    //!
    virtual real_type Dxy( real_type const x, real_type const y ) const = 0;

    //!
    //! Second derivatives respect to \f$ y \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_{yy}(x,y) \f$.
    //!
    virtual real_type Dyy( real_type const x, real_type const y ) const = 0;

    //!
    //! Evaluate spline value at point \f$ (x,y) \f$.
    //!
    real_type operator()( real_type const x, real_type const y ) const { return this->eval( x, y ); }

    //!
    //! Alias for `Dx(x,y)`
    //!
    real_type eval_D_1( real_type const x, real_type const y ) const { return this->Dx( x, y ); }

    //!
    //! Alias for `Dy(x,y)`
    //!
    real_type eval_D_2( real_type const x, real_type const y ) const { return this->Dy( x, y ); }

    //!
    //! Alias for `Dxx(x,y)`
    //!
    real_type eval_D_1_1( real_type const x, real_type const y ) const { return this->Dxx( x, y ); }

    //!
    //! Alias for `Dxy(x,y)`
    //!
    real_type eval_D_1_2( real_type const x, real_type const y ) const { return this->Dxy( x, y ); }

    //!
    //! Alias for `Dyy(x,y)`
    //!
    real_type eval_D_2_2( real_type const x, real_type const y ) const { return this->Dyy( x, y ); }

    ///@}

    //!
    //! Print spline coefficients.
    //!
    virtual void write_to_stream( ostream_type & s ) const = 0;

    //!
    //! Return spline type as a string pointer.
    //!
    virtual char const * type_name() const = 0;

    //!
    //! String information of the kind and order of the spline
    //!
    virtual string info() const { return fmt::format( "Bivariate spline [{}] of type = {}", name(), type_name() ); }

    //!
    //! Print information of the kind and order of the spline
    //!
    void info( ostream_type & stream ) const { stream << this->info() << '\n'; }

    //!
    //! Print stored data x, y, and matrix z.
    //!
    void dump_data( ostream_type & s ) const
    {
      s << "X = [ " << m_X[0];
      for ( integer i = 1; i < m_nx; ++i ) s << ", " << m_X[i];
      s << " ]\nY = [ " << m_Y[0];
      for ( integer j = 1; j < m_ny; ++j ) s << ", " << m_Y[j];
      s << " ]\nZ = [\n";
      for ( integer j = 0; j < m_ny; ++j )
      {
        s << "  [ " << z_node( 0, j );
        for ( integer i = 1; i < m_nx; ++i ) s << ", " << z_node( i, j );
        s << " ]\n";
      }
      s << "\n];\n";
    }

#ifdef SPLINES_BACK_COMPATIBILITY
    integer   numPointX() const { return m_nx; }
    integer   numPointY() const { return m_ny; }
    real_type xNode( integer i ) const { return m_X[i]; }
    real_type yNode( integer i ) const { return m_Y[i]; }
    real_type zNode( integer i, integer j ) const { return z_node( i, j ); }
    real_type xMin() const { return this->x_min(); }
    real_type xMax() const { return this->x_max(); }
    real_type yMin() const { return this->y_min(); }
    real_type yMax() const { return this->y_max(); }
    real_type zMin() const { return m_Z_min; }
    real_type zMax() const { return m_Z_max; }
    void      writeToStream( ostream_type & s ) const { write_to_stream( s ); }
#endif
  };

}  // namespace Splines

#endif

// EOF: SplineSurf.hxx
