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
      // Stampa intestazione
      fmt::print( s, "Nx = {} Ny = {}\n", m_nx, m_ny );

      using MatrixView = Eigen::Map<const Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>;
      using VectorView = Eigen::Map<const Eigen::VectorX<real_type>>;

      // Map dei vettori delle coordinate
      VectorView X( m_X, m_nx );
      VectorView Y( m_Y, m_ny );

      // Map delle matrici dei coefficienti
      // Assumiamo che m_Z e gli altri abbiano dimensione m_nx * m_ny
      MatrixView Z( m_Z, m_nx, m_ny );
      MatrixView DX( m_DX, m_nx, m_ny );
      MatrixView DY( m_DY, m_nx, m_ny );
      MatrixView DXY( m_DXY, m_nx, m_ny );

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
            Z( r0, c0 ),
            Z( r1, c0 ),
            Z( r0, c1 ),
            Z( r1, c1 ),
            // DX values
            DX( r0, c0 ),
            DX( r1, c0 ),
            DX( r0, c1 ),
            DX( r1, c1 ),
            // DY values
            DY( r0, c0 ),
            DY( r1, c0 ),
            DY( r0, c1 ),
            DY( r1, c1 ),
            // DXY values
            DXY( r0, c0 ),
            DXY( r1, c0 ),
            DXY( r0, c1 ),
            DXY( r1, c1 ) );
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
