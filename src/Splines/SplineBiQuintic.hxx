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
 |   ____  _  ___        _       _   _      ____        _ _            ____
 |  | __ )(_)/ _ \ _   _(_)_ __ | |_(_) ___/ ___| _ __ | (_)_ __   ___| __ )  __ _ ___  ___
 |  |  _ \| | | | | | | | | '_ \| __| |/ __\___ \| '_ \| | | '_ \ / _ \  _ \ / _` / __|/ _ \
 |  | |_) | | |_| | |_| | | | | | |_| | (__ ___) | |_) | | | | | |  __/ |_) | (_| \__ \  __/
 |  |____/|_|\__\_\\__,_|_|_| |_|\__|_|\___|____/| .__/|_|_|_| |_|\___|____/ \__,_|___/\___|
 |                                               |_|
\*/

#pragma once

#ifndef SPLINE_BIQUINTIC_SPLINE_BASE_HXX
#define SPLINE_BIQUINTIC_SPLINE_BASE_HXX

namespace Splines
{

  //! Bi-quintic spline base class
  class BiQuinticSplineBase : public SplineSurf
  {
  protected:
    Malloc_real m_mem_biquintic{ "BiQuinticSplineBase" };

    real_type * m_DX = nullptr;
    real_type * m_DY = nullptr;

    real_type * m_DXX = nullptr;
    real_type * m_DYY = nullptr;
    real_type * m_DXY = nullptr;

    real_type * m_DXYY  = nullptr;
    real_type * m_DXXY  = nullptr;
    real_type * m_DXXYY = nullptr;

    using SplineSurf::m_nx;
    using SplineSurf::m_ny;

    using SplineSurf::m_X;
    using SplineSurf::m_Y;
    using SplineSurf::m_Z;

    real_type & Dx_node_ref( integer const i, integer const j ) { return m_DX[ipos_C( i, j )]; }
    real_type & Dy_node_ref( integer const i, integer const j ) { return m_DY[ipos_C( i, j )]; }

    real_type & Dxx_node_ref( integer const i, integer const j ) { return m_DXX[ipos_C( i, j )]; }
    real_type & Dyy_node_ref( integer const i, integer const j ) { return m_DYY[ipos_C( i, j )]; }
    real_type & Dxy_node_ref( integer const i, integer const j ) { return m_DXY[ipos_C( i, j )]; }

    real_type & Dxyy_node_ref( integer const i, integer const j ) { return m_DXYY[ipos_C( i, j )]; }
    real_type & Dxxy_node_ref( integer const i, integer const j ) { return m_DXXY[ipos_C( i, j )]; }
    real_type & Dxxyy_node_ref( integer const i, integer const j ) { return m_DXXYY[ipos_C( i, j )]; }

    void load( integer const i, integer const j, real_type bili5[6][6] ) const;

  public:
    using SplineSurf::eval;

    //! spline constructor
    explicit BiQuinticSplineBase( string_view name = "Spline" ) : SplineSurf( name ) {}

    ~BiQuinticSplineBase() override { m_mem_biquintic.free(); }

    //!
    //! \name Estimated derivatives at interpolation nodes
    //!
    ///@{

    //!
    //! Estimated `x` derivatives at node `(i,j)`
    //!
    real_type Dx_node( integer const i, integer const j ) const { return m_DX[ipos_C( i, j )]; }

    //!
    //! Estimated `y` derivatives at node `(i,j)`
    //!
    real_type Dy_node( integer const i, integer const j ) const { return m_DY[ipos_C( i, j )]; }

    //!
    //! Estimated second `x` second derivatives at node `(i,j)`
    //!
    real_type Dxx_node( integer const i, integer const j ) const { return m_DXX[ipos_C( i, j )]; }

    //!
    //! Estimated second`y` derivatives at node `(i,j)`
    //!
    real_type Dyy_node( integer const i, integer const j ) const { return m_DYY[ipos_C( i, j )]; }

    //!
    //! Estimated mixed `xy` derivatives at node `(i,j)`
    //!
    real_type Dxy_node( integer const i, integer const j ) const { return m_DXY[ipos_C( i, j )]; }

    //!
    //! Estimated `xxy` second derivatives at node `(i,j)`
    //!
    real_type Dxxy_node( integer const i, integer const j ) const { return m_DXXY[ipos_C( i, j )]; }

    //!
    //! Estimated `xyy` derivatives at node `(i,j)`
    //!
    real_type Dxyy_node( integer const i, integer const j ) const { return m_DXYY[ipos_C( i, j )]; }

    //!
    //! Estimated mixed `xxyy` derivatives at node `(i,j)`
    //!
    real_type Dxxyy_node( integer const i, integer const j ) const { return m_DXXYY[ipos_C( i, j )]; }

    ///@}

    //!
    //! \name Evaluate
    //!
    ///@{
    //!
    //! Evaluate spline at point \f$ (x,y) \f$
    //!
    real_type eval( real_type x, real_type y ) const override;

    //!
    //! Evaluate spline with derivative at point \f$ (x,y) \f$
    //!
    //! - `d[0]` the value of the spline
    //! - `d[1]` the value of the spline `x` derivative
    //! - `d[2]` the value of the spline `y` derivative
    //!
    void D( real_type const x, real_type const y, real_type d[3] ) const override;
    //!
    //! Evaluate spline `x`  derivative at point \f$ (x,y) \f$
    //!
    real_type Dx( real_type const x, real_type const y ) const override;
    //!
    //! Evaluate spline `y`  derivative at point \f$ (x,y) \f$
    //!
    real_type Dy( real_type const x, real_type const y ) const override;

    //!
    //! Evaluate spline with derivative at point \f$ (x,y) \f$
    //!
    //! - `d[0]` the value of the spline
    //! - `d[1]` the value of the spline `x` derivative
    //! - `d[2]` the value of the spline `y` derivative
    //! - `d[3]` the value of the spline `x` second derivative
    //! - `d[4]` the value of the spline `y` second derivative
    //! - `d[5]` the value of the spline `xy` mixed derivative
    //!
    void DD( real_type const x, real_type const y, real_type dd[6] ) const override;
    //!
    //! Evaluate spline `x` second derivative at point \f$ (x,y) \f$
    //!
    real_type Dxx( real_type const x, real_type const y ) const override;
    //!
    //! Evaluate spline `xy` mixed derivative at point \f$ (x,y) \f$
    //!
    real_type Dxy( real_type const x, real_type const y ) const override;
    //!
    //! Evaluate spline `y` second derivative at point \f$ (x,y) \f$
    //!
    real_type Dyy( real_type const x, real_type const y ) const override;

    ///@}

#ifdef AUTODIFF_SUPPORT
    //!
    //! \name Autodiff
    //!
    ///@{
    autodiff::dual1st eval( autodiff::dual1st const & x, autodiff::dual1st const & y ) const override
    {
      using autodiff::dual1st;
      using autodiff::detail::val;

      real_type dd[3];
      D( val( x ), val( y ), dd );

      dual1st res{ dd[0] };
      res.grad = dd[1] * x.grad + dd[2] * y.grad;

      return res;
    }

    autodiff::dual2nd eval( autodiff::dual2nd const & x, autodiff::dual2nd const & y ) const override
    {
      using autodiff::derivative;
      using autodiff::dual2nd;

      real_type dd[6], dx{ val( x.grad ) }, dy{ val( y.grad ) }, ddx{ x.grad.grad }, ddy{ y.grad.grad };
      DD( val( x ), val( y ), dd );

      dual2nd res{ dd[0] };
      res.grad      = dd[1] * dx + dd[2] * dy;
      res.grad.grad = dd[3] * dx * dx + 2 * dx * dy * dd[4] + dy * dy * dd[5] + ddx * dd[1] + ddy * dd[2];
      return res;
    }

    template <typename T1, typename T2>
    autodiff::HigherOrderDual<autodiff::detail::DualOrder<T1, T2>::value, real_type> eval( T1 const & x, T2 const & y )
      const
    {
      autodiff::HigherOrderDual<autodiff::detail::DualOrder<T1, T2>::value, real_type> X{ x }, Y{ y };
      return eval( X, Y );
    }
///@}
#endif
  };

  /*\
   |   ____  _  ___        _       _   _      ____        _ _
   |  | __ )(_)/ _ \ _   _(_)_ __ | |_(_) ___/ ___| _ __ | (_)_ __   ___
   |  |  _ \| | | | | | | | | '_ \| __| |/ __\___ \| '_ \| | | '_ \ / _ \
   |  | |_) | | |_| | |_| | | | | | |_| | (__ ___) | |_) | | | | | |  __/
   |  |____/|_|\__\_\\__,_|_|_| |_|\__|_|\___|____/| .__/|_|_|_| |_|\___|
   |                                               |_|
  \*/
  //! cubic spline base class
  class BiQuinticSpline : public BiQuinticSplineBase
  {
    void make_spline() override
    {
      integer const dim{ m_nx * m_ny };
      m_mem_biquintic.reallocate( 8 * dim );
      m_DX    = m_mem_biquintic( dim );
      m_DY    = m_mem_biquintic( dim );
      m_DXY   = m_mem_biquintic( dim );
      m_DXX   = m_mem_biquintic( dim );
      m_DYY   = m_mem_biquintic( dim );
      m_DXYY  = m_mem_biquintic( dim );
      m_DXXY  = m_mem_biquintic( dim );
      m_DXXYY = m_mem_biquintic( dim );

      make_derivative_x( m_Z, m_DX );
      make_derivative_y( m_Z, m_DY );
      make_derivative_xy( m_DX, m_DY, m_DXY );

      make_derivative_x( m_DX, m_DXX );
      make_derivative_y( m_DY, m_DYY );

      make_derivative_y( m_DXX, m_DXXY );
      make_derivative_x( m_DYY, m_DXYY );

      make_derivative_xy( m_DXXY, m_DXYY, m_DXXYY );

      m_search_x.must_reset();
      m_search_y.must_reset();
    }

  public:
    //!
    //! Build an empty spline of `BiQuinticSpline` type
    //!
    //! \param name the name of the spline
    //!
    explicit BiQuinticSpline( string_view name = "BiQuinticSpline" ) : BiQuinticSplineBase( name ) {}

    //!
    //! Spline destructor.
    //!
    ~BiQuinticSpline() override {}

    void write_to_stream( ostream_type & s ) const override
    {
      fmt::print( s, "Nx = {} Ny = {}\n", m_nx, m_ny );
      for ( integer i = 1; i < m_nx; ++i )
      {
        real_type dx{ m_X[i] - m_X[i - 1] };
        for ( integer j = 1; j < m_ny; ++j )
        {
          integer const   i00{ ipos_C( i - 1, j - 1 ) };
          integer const   i10{ ipos_C( i, j - 1 ) };
          integer const   i01{ ipos_C( i - 1, j ) };
          integer const   i11{ ipos_C( i, j ) };
          real_type const dy{ m_Y[j] - m_Y[j - 1] };
          fmt::print(
            s,
            "patch ({},{})\n"
            "  DX    = {:<12.4}  DY    = {:<12.4}\n"
            "  Z00   = {:<12.4}  Z10   = {:<12.4}\n"
            "  Z01   = {:<12.4}  Z11   = {:<12.4}\n"
            "  Dx00  = {:<12.4}  Dx10  = {:<12.4}\n"
            "  Dx01  = {:<12.4}  Dx11  = {:<12.4}\n"
            "  Dy00  = {:<12.4}  Dy10  = {:<12.4}\n"
            "  Dy01  = {:<12.4}  Dy11  = {:<12.4}\n"
            "  Dxy00 = {:<12.4}  Dxy10 = {:<12.4}\n"
            "  Dxy01 = {:<12.4}  Dxy11 = {:<12.4}\n",
            i,
            j,
            dx,
            dy,
            m_Z[i00],
            m_Z[i10],
            m_Z[i01],
            m_Z[i11],
            m_DX[i00],
            m_DX[i10],
            m_DX[i01],
            m_DX[i11],
            m_DY[i00],
            m_DY[i10],
            m_DY[i01],
            m_DY[i11],
            m_DXY[i00],
            m_DXY[i10],
            m_DXY[i01],
            m_DXY[i11] );
        }
      }
    }

    char const * type_name() const override { return "BiQuintic"; }
  };

}  // namespace Splines

#endif

//
// EOF: SplineBiQuintic.hxx
//
