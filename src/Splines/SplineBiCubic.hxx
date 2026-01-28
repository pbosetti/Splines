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
      m_DX_ptr  = m_mem_bicubic( nn );
      m_DY_ptr  = m_mem_bicubic( nn );
      m_DXY_ptr = m_mem_bicubic( nn );
      
      new (&mDX)  Eigen::Map<MatC>( m_DX_ptr, m_nx, m_ny );
      new (&mDY)  Eigen::Map<MatC>( m_DY_ptr, m_nx, m_ny );
      new (&mDXY) Eigen::Map<MatC>( m_DXY_ptr, m_nx, m_ny );

      make_derivative_x( m_sub_type, mZ.data(), m_DX_ptr );
      make_derivative_y( m_sub_type, mZ.data(), m_DY_ptr );
      make_derivative_xy( m_sub_type, m_DX_ptr, m_DY_ptr, m_DXY_ptr );

      m_search_x.must_reset();
      m_search_y.must_reset();
    }

    using BiCubicSplineBase::mDX;
    using BiCubicSplineBase::mDXY;
    using BiCubicSplineBase::mDY;
    using BiCubicSplineBase::m_DX_ptr;
    using BiCubicSplineBase::m_DXY_ptr;
    using BiCubicSplineBase::m_DY_ptr;
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
      // Stampa intestazione
      fmt::print( s, "Nx = {} Ny = {}\n", m_nx, m_ny );

      // Map dei vettori delle coordinate
      Eigen::Map<Vec> X( m_X, m_nx );
      Eigen::Map<Vec> Y( m_Y, m_ny );

      for ( integer i = 1; i < m_nx; ++i )
      {
        // Accesso diretto al vettore mappato
        real_type const dx = X[i] - X[i - 1];

        for ( integer j = 1; j < m_ny; ++j )
        {
          real_type const dy = Y[j] - Y[j - 1];

          // Indici base per la patch corrente
          integer const r0 = i - 1;
          integer const r1 = i;
          integer const c0 = j - 1;
          integer const c1 = j;

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
            // Z values
            mZ( r0, c0 ),
            mZ( r1, c0 ),
            mZ( r0, c1 ),
            mZ( r1, c1 ),
            // DX values
            mDX( r0, c0 ),
            mDX( r1, c0 ),
            mDX( r0, c1 ),
            mDX( r1, c1 ),
            // DY values
            mDY( r0, c0 ),
            mDY( r1, c0 ),
            mDY( r0, c1 ),
            mDY( r1, c1 ),
            // DXY values
            mDXY( r0, c0 ),
            mDXY( r1, c0 ),
            mDXY( r0, c1 ),
            mDXY( r1, c1 ) );
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
