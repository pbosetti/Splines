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
