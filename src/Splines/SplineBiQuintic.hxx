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

#ifndef SPLINE_BIQUINTIC_SPLINE_HXX
#define SPLINE_BIQUINTIC_SPLINE_HXX

namespace Splines
{

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
      m_DX_ptr    = m_mem_biquintic( dim );
      m_DY_ptr    = m_mem_biquintic( dim );
      m_DXY_ptr   = m_mem_biquintic( dim );
      m_DXX_ptr   = m_mem_biquintic( dim );
      m_DYY_ptr   = m_mem_biquintic( dim );
      m_DXYY_ptr  = m_mem_biquintic( dim );
      m_DXXY_ptr  = m_mem_biquintic( dim );
      m_DXXYY_ptr = m_mem_biquintic( dim );

      new (&mDX) Eigen::Map<MatC>( m_DX_ptr, m_nx, m_ny );
      new (&mDY) Eigen::Map<MatC>( m_DY_ptr, m_nx, m_ny );
      new (&mDXY) Eigen::Map<MatC>( m_DXY_ptr, m_nx, m_ny );
      new (&mDXX) Eigen::Map<MatC>( m_DXX_ptr, m_nx, m_ny );
      new (&mDYY) Eigen::Map<MatC>( m_DYY_ptr, m_nx, m_ny );
      new (&mDXYY) Eigen::Map<MatC>( m_DXYY_ptr, m_nx, m_ny );
      new (&mDXXY) Eigen::Map<MatC>( m_DXXY_ptr, m_nx, m_ny );
      new (&mDXXYY) Eigen::Map<MatC>( m_DXXYY_ptr, m_nx, m_ny );

      make_derivative_x( m_sub_type, m_Z_ptr, m_DX_ptr );
      make_derivative_y( m_sub_type, m_Z_ptr, m_DY_ptr );
      make_derivative_xy( m_sub_type, m_DX_ptr, m_DY_ptr, m_DXY_ptr );

      make_derivative_x( m_sub_type, m_DX_ptr, m_DXX_ptr );
      make_derivative_y( m_sub_type, m_DY_ptr, m_DYY_ptr );

      make_derivative_y( m_sub_type, m_DXX_ptr, m_DXXY_ptr );
      make_derivative_x( m_sub_type, m_DYY_ptr, m_DXYY_ptr );

      make_derivative_xy( m_sub_type, m_DXXY_ptr, m_DXYY_ptr, m_DXXYY_ptr );

      m_search_x.must_reset();
      m_search_y.must_reset();
    }

  public:
    //!
    //! Build an empty spline of `BiQuinticSpline` type
    //!
    //! \param name the name of the spline
    //!
    explicit BiQuinticSpline( Spline_sub_type sub_type, string_view name = "BiQuinticSpline" )
      : BiQuinticSplineBase( sub_type, name )
    {
    }

    //!
    //! Spline destructor.
    //!
    ~BiQuinticSpline() override {}

    void write_to_stream( ostream_type & s ) const override
    {
      fmt::print( s, "Nx = {} Ny = {}\n", m_nx, m_ny );
      for ( integer i = 1; i < m_nx; ++i )
      {
        real_type dx = m_X[i] - m_X[i - 1];
        for ( integer j = 1; j < m_ny; ++j )
        {
          real_type const dy = m_Y[j] - m_Y[j - 1];
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
            mZ.coeff( i - 1, j - 1 ),
            mZ.coeff( i, j - 1 ),
            mZ.coeff( i, j + 1 ),
            mZ.coeff( i + 1, j + 1 ),
            mDX.coeff( i - 1, j - 1 ),
            mDX.coeff( i, j - 1 ),
            mDX.coeff( i, j + 1 ),
            mDX.coeff( i + 1, j + 1 ),
            mDY.coeff( i - 1, j - 1 ),
            mDY.coeff( i, j - 1 ),
            mDY.coeff( i, j + 1 ),
            mDY.coeff( i + 1, j + 1 ),
            mDXY.coeff( i - 1, j - 1 ),
            mDXY.coeff( i, j - 1 ),
            mDXY.coeff( i, j + 1 ),
            mDXY.coeff( i + 1, j + 1 ) );
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
