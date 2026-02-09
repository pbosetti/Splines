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
      mDX.resize( m_nx, m_ny );
      mDY.resize( m_nx, m_ny );
      mDXY.resize( m_nx, m_ny );
      mDXX.resize( m_nx, m_ny );
      mDYY.resize( m_nx, m_ny );
      mDXYY.resize( m_nx, m_ny );
      mDXXY.resize( m_nx, m_ny );
      mDXXYY.resize( m_nx, m_ny );

      CubicSpline   cs;
      AkimaSpline   ak;
      VanLeerSpline vl;
      PchipSpline   pc;

      CubicSplineBase * S;
      switch ( m_sub_type )
      {
        case Spline_sub_type::CUBIC: S = &cs; break;
        case Spline_sub_type::AKIMA: S = &ak; break;
        case Spline_sub_type::VANLEER: S = &vl; break;
        case Spline_sub_type::PCHIP: S = &pc; break;
      }

      make_derivative_x( S, mZ, mDX );
      make_derivative_y( S, mZ, mDY );
      make_derivative_xy( S, mDX, mDY, mDXY );

      make_derivative_x( S, mDX, mDXX );
      make_derivative_y( S, mDY, mDYY );

      make_derivative_y( S, mDXX, mDXXY );
      make_derivative_x( S, mDYY, mDXYY );

      make_derivative_xy( S, mDXXY, mDXYY, mDXXYY );

      m_search_x.must_reset();
      m_search_y.must_reset();
    }

  public:
    //!
    //! Build an empty spline of `BiQuinticSpline` type
    //!
    //! \param type spline type
    //! \param name the name of the spline
    //!
    explicit BiQuinticSpline( Spline_sub_type sub_type, string_view name = "BiQuinticSpline" ) : BiQuinticSplineBase( sub_type, name ) { }

    //!
    //! Build an empty spline of `BiQuinticSpline` type
    //!
    //! \param name the name of the spline
    //!
    explicit BiQuinticSpline( string_view name ) : BiQuinticSplineBase( Spline_sub_type::PCHIP, name ) { }

    //!
    //! Spline destructor.
    //!
    ~BiQuinticSpline() override {}

    void write_to_stream( ostream_type & s ) const override
    {
      fmt::print( s, "Nx = {} Ny = {}\n", m_nx, m_ny );
      for ( integer i = 1; i < m_nx; ++i )
      {
        integer   i1 = i - 1;
        real_type dx = mX.coeff( i ) - mX.coeff( i1 );
        for ( integer j = 1; j < m_ny; ++j )
        {
          integer         j1 = j - 1;
          real_type const dy = mY.coeff( j ) - mY.coeff( j1 );
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
            mZ.coeff( i1, j1 ),
            mZ.coeff( i, j1 ),
            mZ.coeff( i1, j ),
            mZ.coeff( i, j ),
            mDX.coeff( i1, j1 ),
            mDX.coeff( i, j1 ),
            mDX.coeff( i1, j ),
            mDX.coeff( i, j ),
            mDY.coeff( i1, j1 ),
            mDY.coeff( i, j1 ),
            mDY.coeff( i1, j ),
            mDY.coeff( i, j ),
            mDXY.coeff( i1, j1 ),
            mDXY.coeff( i, j1 ),
            mDXY.coeff( i1, j ),
            mDXY.coeff( i, j ) );
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
