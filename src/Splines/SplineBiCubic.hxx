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

#ifndef SPLINE_BICUBIC_HXX
#define SPLINE_BICUBIC_HXX

namespace Splines
{

  /*\
   |   ____  _  ____      _     _      ____        _ _
   |  | __ )(_)/ ___|   _| |__ (_) ___/ ___| _ __ | (_)_ __   ___
   |  |  _ \| | |  | | | | '_ \| |/ __\___ \| '_ \| | | '_ \ / _ \
   |  | |_) | | |__| |_| | |_) | | (__ ___) | |_) | | | | | |  __/
   |  |____/|_|\____\__,_|_.__/|_|\___|____/| .__/|_|_|_| |_|\___|
   |                                        |_|
  \*/
  //!
  //! Cubic spline base class
  //!
  class BiCubicSpline : public BiCubicSplineBase
  {
    void make_spline() override
    {
      integer const nn{ m_nx * m_ny };
      m_mem_bicubic.reallocate( 3 * nn );
      m_DX  = m_mem_bicubic( nn );
      m_DY  = m_mem_bicubic( nn );
      m_DXY = m_mem_bicubic( nn );

      make_derivative_x( m_sub_type, m_Z, m_DX );
      make_derivative_y( m_sub_type, m_Z, m_DY );
      make_derivative_xy( m_sub_type, m_DX, m_DY, m_DXY );

      m_search_x.must_reset();
      m_search_y.must_reset();
    }

    using BiCubicSplineBase::m_DX;
    using BiCubicSplineBase::m_DXY;
    using BiCubicSplineBase::m_DY;
    using BiCubicSplineBase::m_mem_bicubic;

  public:
    using BiCubicSplineBase::eval;

    //!
    //! Build an empty spline of `BiCubicSpline` type
    //!
    //! \param name the name of the spline
    //!
    explicit BiCubicSpline( Spline_sub_type sub_type, string_view name = "BiCubicSpline" )
      : BiCubicSplineBase( sub_type, name )
    {
    }

    //!
    //! Spline destructor.
    //!
    ~BiCubicSpline() override {}

    void write_to_stream( ostream_type & s ) const override
    {
      fmt::print( "Nx = {} Ny = {}\n", m_nx, m_ny );
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

    char const * type_name() const override { return "BiCubic"; }
  };

}  // namespace Splines

#endif

//
// EOF: SplineBiCubic.hxx
//
