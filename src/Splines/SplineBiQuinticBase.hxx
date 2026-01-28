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
    MatC mDX;
    MatC mDY;

    MatC mDXX;
    MatC mDYY;
    MatC mDXY;

    MatC mDXYY;
    MatC mDXXY;
    MatC mDXXYY;

    Spline_sub_type m_sub_type;

    using SplineSurf::m_nx;
    using SplineSurf::m_ny;

    using SplineSurf::mX;
    using SplineSurf::mY;
    using SplineSurf::mZ;

    real_type & Dx_node_ref( integer const i, integer const j ) { return mDX.coeffRef( i, j ); }
    real_type & Dy_node_ref( integer const i, integer const j ) { return mDY.coeffRef( i, j ); }

    real_type & Dxx_node_ref( integer const i, integer const j ) { return mDXX.coeffRef( i, j ); }
    real_type & Dyy_node_ref( integer const i, integer const j ) { return mDYY.coeffRef( i, j ); }
    real_type & Dxy_node_ref( integer const i, integer const j ) { return mDXY.coeffRef( i, j ); }

    real_type & Dxyy_node_ref( integer const i, integer const j ) { return mDXYY.coeffRef( i, j ); }
    real_type & Dxxy_node_ref( integer const i, integer const j ) { return mDXXY.coeffRef( i, j ); }

    real_type & Dxxyy_node_ref( integer const i, integer const j ) { return mDXXYY.coeffRef( i, j ); }

    void load( integer const i, integer const j, real_type bili5[6][6] ) const
    {
      Eigen::Map<Eigen::Matrix<real_type, 6, 6, Eigen::RowMajor>> bili5_eigen( &bili5[0][0] );
      // Usa Eigen per estrarre i blocchi direttamente
      // Nota: richiede che mZ, mDX, ecc. siano matrici Eigen dense
      bili5_eigen.block<2, 2>( 0, 0 ) = mZ.block( i, j, 2, 2 );
      bili5_eigen.block<2, 2>( 2, 0 ) = mDX.block( i, j, 2, 2 );
      bili5_eigen.block<2, 2>( 4, 0 ) = mDXX.block( i, j, 2, 2 );
      bili5_eigen.block<2, 2>( 0, 2 ) = mDY.block( i, j, 2, 2 );
      bili5_eigen.block<2, 2>( 0, 4 ) = mDYY.block( i, j, 2, 2 );
      bili5_eigen.block<2, 2>( 2, 2 ) = mDXY.block( i, j, 2, 2 );
      bili5_eigen.block<2, 2>( 4, 2 ) = mDXXY.block( i, j, 2, 2 );
      bili5_eigen.block<2, 2>( 2, 4 ) = mDXYY.block( i, j, 2, 2 );
      bili5_eigen.block<2, 2>( 4, 4 ) = mDXXYY.block( i, j, 2, 2 );
    }

  public:
    using SplineSurf::eval;

    //! spline constructor
    explicit BiQuinticSplineBase( Spline_sub_type sub_type, string_view name = "Spline" )
      : SplineSurf( name ), m_sub_type( sub_type )
    {
    }

    ~BiQuinticSplineBase() override = default;

    //!
    //! \name Estimated derivatives at interpolation nodes
    //!
    ///@{

    //!
    //! Estimated `x` derivatives at node `(i,j)`
    //!
    real_type Dx_node( integer const i, integer const j ) const { return mDX.coeff( i, j ); }

    //!
    //! Estimated `y` derivatives at node `(i,j)`
    //!
    real_type Dy_node( integer const i, integer const j ) const { return mDY.coeff( i, j ); }

    //!
    //! Estimated second `x` second derivatives at node `(i,j)`
    //!
    real_type Dxx_node( integer const i, integer const j ) const { return mDXX.coeff( i, j ); }

    //!
    //! Estimated second`y` derivatives at node `(i,j)`
    //!
    real_type Dyy_node( integer const i, integer const j ) const { return mDYY.coeff( i, j ); }

    //!
    //! Estimated mixed `xy` derivatives at node `(i,j)`
    //!
    real_type Dxy_node( integer const i, integer const j ) const { return mDXY.coeff( i, j ); }

    //!
    //! Estimated `xxy` second derivatives at node `(i,j)`
    //!
    real_type Dxxy_node( integer const i, integer const j ) const { return mDXXY.coeff( i, j ); }

    //!
    //! Estimated `xyy` derivatives at node `(i,j)`
    //!
    real_type Dxyy_node( integer const i, integer const j ) const { return mDXYY.coeff( i, j ); }

    //!
    //! Estimated mixed `xxyy` derivatives at node `(i,j)`
    //!
    real_type Dxxyy_node( integer const i, integer const j ) const { return mDXXYY.coeff( i, j ); }

    ///@}

    //!
    //! \name Evaluate
    //!
    ///@{
    //!
    //! Evaluate spline at point \f$ (x,y) \f$
    //!
    real_type eval( real_type x, real_type y ) const override
    {
      real_type bili5[6][6], u[6], v[6];

      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i{ X.first };
      integer const j{ Y.first };

      real_type const dx{ X.second - mX.coeff( i ) };
      real_type const dy{ Y.second - mY.coeff( j ) };
      real_type const DX{ mX.coeff( i + 1 ) - mX.coeff( i ) };
      real_type const DY{ mY.coeff( j + 1 ) - mY.coeff( j ) };

      Hermite5( dx, DX, u );
      Hermite5( dy, DY, v );

      load( i, j, bili5 );

      return bilinear5( u, bili5, v );
    }


    //!
    //! Evaluate spline with derivative at point \f$ (x,y) \f$
    //!
    //! - `d[0]` the value of the spline
    //! - `d[1]` the value of the spline `x` derivative
    //! - `d[2]` the value of the spline `y` derivative
    //!
    void D( real_type const x, real_type const y, real_type d[3] ) const override
    {
      real_type bili5[6][6], u[6], u_D[6], v[6], v_D[6];

      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i{ X.first };
      integer const j{ Y.first };

      real_type const dx{ X.second - mX.coeff( i ) };
      real_type const dy{ Y.second - mY.coeff( j ) };
      real_type const DX{ mX.coeff( i + 1 ) - mX.coeff( i ) };
      real_type const DY{ mY.coeff( j + 1 ) - mY.coeff( j ) };

      Hermite5( dx, DX, u );
      Hermite5_D( dx, DX, u_D );
      Hermite5( dy, DY, v );
      Hermite5_D( dy, DY, v_D );

      load( i, j, bili5 );

      d[0] = bilinear5( u, bili5, v );
      d[1] = bilinear5( u_D, bili5, v );
      d[2] = bilinear5( u, bili5, v_D );
    }

    //!
    //! Evaluate spline `x`  derivative at point \f$ (x,y) \f$
    //!
    real_type Dx( real_type const x, real_type const y ) const override
    {
      real_type bili5[6][6], u_D[6], v[6];

      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i{ X.first };
      integer const j{ Y.first };

      real_type const dx{ X.second - mX.coeff( i ) };
      real_type const dy{ Y.second - mY.coeff( j ) };
      real_type const DX{ mX.coeff( i + 1 ) - mX.coeff( i ) };
      real_type const DY{ mY.coeff( j + 1 ) - mY.coeff( j ) };

      Hermite5_D( dx, DX, u_D );
      Hermite5( dy, DY, v );

      load( i, j, bili5 );

      return bilinear5( u_D, bili5, v );
    }

    //!
    //! Evaluate spline `y`  derivative at point \f$ (x,y) \f$
    //!
    real_type Dy( real_type const x, real_type const y ) const override
    {
      real_type bili5[6][6], u[6], v_D[6];

      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i{ X.first };
      integer const j{ Y.first };

      real_type const dx{ X.second - mX.coeff( i ) };
      real_type const dy{ Y.second - mY.coeff( j ) };
      real_type const DX{ mX.coeff( i + 1 ) - mX.coeff( i ) };
      real_type const DY{ mY.coeff( j + 1 ) - mY.coeff( j ) };

      Hermite5( dx, DX, u );
      Hermite5_D( dy, DY, v_D );

      load( i, j, bili5 );

      return bilinear5( u, bili5, v_D );
    }

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
    void DD( real_type const x, real_type const y, real_type dd[6] ) const override
    {
      real_type bili5[6][6], u[6], u_D[6], u_DD[6], v[6], v_D[6], v_DD[6];

      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i{ X.first };
      integer const j{ Y.first };

      real_type const dx{ X.second - mX.coeff( i ) };
      real_type const dy{ Y.second - mY.coeff( j ) };
      real_type const DX{ mX.coeff( i + 1 ) - mX.coeff( i ) };
      real_type const DY{ mY.coeff( j + 1 ) - mY.coeff( j ) };

      Hermite5( dx, DX, u );
      Hermite5_D( dx, DX, u_D );
      Hermite5_DD( dx, DX, u_DD );
      Hermite5( dy, DY, v );
      Hermite5_D( dy, DY, v_D );
      Hermite5_DD( dy, DY, v_DD );

      load( i, j, bili5 );

      dd[0] = bilinear5( u, bili5, v );
      dd[1] = bilinear5( u_D, bili5, v );
      dd[2] = bilinear5( u, bili5, v_D );
      dd[3] = bilinear5( u_DD, bili5, v );
      dd[4] = bilinear5( u_D, bili5, v_D );
      dd[5] = bilinear5( u, bili5, v_DD );
    }

    //!
    //! Evaluate spline `x` second derivative at point \f$ (x,y) \f$
    //!
    real_type Dxx( real_type const x, real_type const y ) const override
    {
      real_type bili5[6][6], u_DD[6], v[6];

      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i{ X.first };
      integer const j{ Y.first };

      real_type const dx{ X.second - mX.coeff( i ) };
      real_type const dy{ Y.second - mY.coeff( j ) };
      real_type const DX{ mX.coeff( i + 1 ) - mX.coeff( i ) };
      real_type const DY{ mY.coeff( j + 1 ) - mY.coeff( j ) };

      Hermite5_DD( dx, DX, u_DD );
      Hermite5( dy, DY, v );

      load( i, j, bili5 );

      return bilinear5( u_DD, bili5, v );
    }

    //!
    //! Evaluate spline `xy` mixed derivative at point \f$ (x,y) \f$
    //!
    real_type Dxy( real_type const x, real_type const y ) const override
    {
      real_type bili5[6][6], u_D[6], v_D[6];

      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i{ X.first };
      integer const j{ Y.first };

      real_type const dx{ X.second - mX.coeff( i ) };
      real_type const dy{ Y.second - mY.coeff( j ) };
      real_type const DX{ mX.coeff( i + 1 ) - mX.coeff( i ) };
      real_type const DY{ mY.coeff( j + 1 ) - mY.coeff( j ) };

      Hermite5_D( dx, DX, u_D );
      Hermite5_D( dy, DY, v_D );

      load( i, j, bili5 );

      return bilinear5( u_D, bili5, v_D );
    }

    //!
    //! Evaluate spline `y` second derivative at point \f$ (x,y) \f$
    //!
    real_type Dyy( real_type const x, real_type const y ) const override
    {
      real_type bili5[6][6], u[6], v_DD[6];

      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i{ X.first };
      integer const j{ Y.first };

      real_type const dx{ X.second - mX.coeff( i ) };
      real_type const dy{ Y.second - mY.coeff( j ) };
      real_type const DX{ mX.coeff( i + 1 ) - mX.coeff( i ) };
      real_type const DY{ mY.coeff( j + 1 ) - mY.coeff( j ) };

      Hermite5( dx, DX, u );
      Hermite5_DD( dy, DY, v_DD );

      load( i, j, bili5 );

      return bilinear5( u, bili5, v_DD );
    }

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

      real_type dd[6];
      real_type dx  = val( x.grad );
      real_type dy  = val( y.grad );
      real_type ddx = x.grad.grad;
      real_type ddy = y.grad.grad;

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

}  // namespace Splines

#endif

//
// EOF: SplineBiQuintic.hxx
//
