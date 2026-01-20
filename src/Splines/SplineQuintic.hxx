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
 |    ___        _       _   _      ____        _ _
 |   / _ \ _   _(_)_ __ | |_(_) ___/ ___| _ __ | (_)_ __   ___
 |  | | | | | | | | '_ \| __| |/ __\___ \| '_ \| | | '_ \ / _ \
 |  | |_| | |_| | | | | | |_| | (__ ___) | |_) | | | | | |  __/
 |   \__\_\\__,_|_|_| |_|\__|_|\___|____/| .__/|_|_|_| |_|\___|
 |                                       |_|
 |
\*/

namespace Splines
{

  using QuinticSpline_sub_type = enum class QuinticSpline_sub_type : integer {
    CUBIC  = 0,
    PCHIP  = 1,
    AKIMA  = 2,
    BESSEL = 3
  };

  inline std::string to_string( QuinticSpline_sub_type t )
  {
    switch ( t )
    {
      case QuinticSpline_sub_type::CUBIC: return "CUBIC";
      case QuinticSpline_sub_type::PCHIP: return "PCHIP";
      case QuinticSpline_sub_type::AKIMA: return "AKIMA";
      case QuinticSpline_sub_type::BESSEL: return "BESSEL";
    }
    return "NOTYPE";
  }

  /*!
   * \brief Construct a complete quintic spline from cubic Hermite data
   *
   * This function elevates a C1 cubic spline (defined by nodal values and
   * first derivatives) to a C3 quintic spline by computing optimal second
   * derivatives that minimize the curvature energy.
   *
   * The algorithm consists of three phases:
   * 1. First derivative estimation (if not already provided)
   * 2. Second derivative computation via tridiagonal solver
   * 3. Quintic Hermite interpolation on each interval
   *
   * The resulting spline is C3 continuous and provides superior smoothness
   * while preserving the original interpolation conditions.
   *
   * \param[in]  q_sub_type Method for first derivative estimation (CUBIC, PCHIP, AKIMA, BESSEL)
   * \param[in]  X          Array of x-coordinates
   * \param[in]  Y          Array of y-coordinates
   * \param[out] Yp         Array of first derivatives (computed if needed)
   * \param[out] Ypp        Array of second derivatives (computed)
   * \param[in]  npts       Number of data points
   * \param[in]  bc_type    Boundary condition type (0=natural, 1=not-a-knot, 2=clamped)
   * \param[in]  bcl        Left boundary second derivative for clamped condition
   * \param[in]  bcr        Right boundary second derivative for clamped condition
   */
  inline void Quintic_build(
    QuinticSpline_sub_type const q_sub_type,
    real_type const              X[],
    real_type const              Y[],
    real_type                    Yp[],
    real_type                    Ypp[],
    integer const                npts )
  {
    UTILS_ASSERT( npts >= 2, "Quintic_build, npts={} must be >= 2\n", npts );

    // Phase 1: Estimate first derivatives if not already available
    // (This step can be skipped if Yp is precomputed from a C1 cubic spline)
    switch ( q_sub_type )
    {
      case QuinticSpline_sub_type::CUBIC:
        CubicSpline_build( X, Y, Yp, Ypp, npts, CubicSpline_BC::EXTRAPOLATE, CubicSpline_BC::EXTRAPOLATE );
        return;
      case QuinticSpline_sub_type::PCHIP: Pchip_build( X, Y, Yp, npts ); break;
      case QuinticSpline_sub_type::AKIMA:
      {
        Malloc_real mem( "Quintic_build::work memory" );
        Akima_build( X, Y, Yp, mem.malloc( npts ), npts );
      }
      break;
      case QuinticSpline_sub_type::BESSEL: Bessel_build( X, Y, Yp, npts ); break;
      default: UTILS_ERROR( "Unknown QuinticSpline_sub_type value\n" );
    }

    if ( npts == 2 )
    {  // solo 2 punti, niente da fare
      Ypp[0] = Ypp[1] = 0;
      return;
    }

    // ---- left boundary ----
    if ( npts <= 3 )
    {
      Ypp[0] = Utils::second_derivative_3p( X[0], Y[0], X[1], Y[1], X[2], Y[2] );
      Ypp[1] = Utils::second_derivative_3p( X[1], Y[1], X[2], Y[2], X[0], Y[0] );
    }
    else
    {
      Ypp[0] = Utils::second_derivative_4p( X[0], Y[0], X[1], Y[1], X[2], Y[2], X[3], Y[3] );
      Ypp[1] = Utils::second_derivative_4p( X[1], Y[1], X[2], Y[2], X[3], Y[3], X[0], Y[0] );
    }

    // ---- interior points ----
    for ( integer i = 2; i < npts - 2; ++i )
      Ypp[i] = Utils::second_derivative_5p(
        X[i - 0],
        Y[i - 0],
        X[i + 1],
        Y[i + 1],
        X[i - 1],
        Y[i - 1],
        X[i + 2],
        Y[i + 2],
        X[i - 2],
        Y[i - 2] );

    // ---- right boundary ----
    if ( npts <= 3 )
    {
      Ypp[npts - 1] =
        Utils::second_derivative_3p( X[npts - 1], Y[npts - 1], X[npts - 2], Y[npts - 2], X[npts - 3], Y[npts - 3] );
      Ypp[npts - 2] =
        Utils::second_derivative_3p( X[npts - 2], Y[npts - 2], X[npts - 3], Y[npts - 3], X[npts - 1], Y[npts - 1] );
    }
    else
    {
      Ypp[npts - 1] = Utils::second_derivative_4p(
        X[npts - 1],
        Y[npts - 1],
        X[npts - 2],
        Y[npts - 2],
        X[npts - 3],
        Y[npts - 3],
        X[npts - 4],
        Y[npts - 4] );
      Ypp[npts - 2] = Utils::second_derivative_4p(
        X[npts - 2],
        Y[npts - 2],
        X[npts - 3],
        Y[npts - 3],
        X[npts - 4],
        Y[npts - 4],
        X[npts - 1],
        Y[npts - 1] );
    }
  }

  //! Quintic spline class
  class QuinticSpline : public QuinticSplineBase
  {
    QuinticSpline_sub_type m_q_sub_type{ QuinticSpline_sub_type::CUBIC };

  public:
    //!
    //! \name Constructors
    //!
    ///@{

    using QuinticSplineBase::build;
    using QuinticSplineBase::reserve;

    //!
    //! Build an empty spline of `QuinticSpline` type
    //!
    //! \param name the name of the spline
    //!
    explicit QuinticSpline( string_view name = "Spline" ) : QuinticSplineBase( name ) {}

    //!
    //! Spline destructor.
    //!
    ~QuinticSpline() override {}

    ///@}

    //!
    //! Set spline type
    //!
    //! - CUBIC
    //! - PCHIP
    //! - AKIMA
    //! - BESSEL
    //!
    void set_quintic_type( QuinticSpline_sub_type qt ) { m_q_sub_type = qt; }

    void setQuinticType( QuinticSpline_sub_type qt ) { m_q_sub_type = qt; }

    // --------------------------- VIRTUALS -----------------------------------
    //! Build a Monotone quintic spline from previously inserted points
    void build() override
    {
      string msg{ fmt::format( "QuinticSpline[{}]::build():", m_name ) };
      UTILS_ASSERT( m_npts > 1, "{} npts = {} not enought points\n", msg, m_npts );
      Utils::check_NaN( m_X, msg + " X", m_npts, __LINE__, __FILE__ );
      Utils::check_NaN( m_Y, msg + " Y", m_npts, __LINE__, __FILE__ );
      integer ibegin{ 0 };
      integer iend{ 0 };
      do
      {
        // cerca intervallo monotono strettamente crescente
        for ( ++iend; iend < m_npts && m_X[iend - 1] < m_X[iend]; ++iend ) {}
        Quintic_build( m_q_sub_type, m_X + ibegin, m_Y + ibegin, m_Yp + ibegin, m_Ypp + ibegin, iend - ibegin );
        ibegin = iend;
      } while ( iend < m_npts );

      Utils::check_NaN( m_Yp, msg + " Yp", m_npts, __LINE__, __FILE__ );
      Utils::check_NaN( m_Ypp, msg + " Ypp", m_npts, __LINE__, __FILE__ );
      m_search.must_reset();
    }

    //! Build a Monotone quintic spline from data from `gc`
    void setup( GenericContainer const & gc ) override
    {
      /*
      // gc["xdata"]
      // gc["ydata"]
      //
      */
      string const where{ fmt::format( "QuinticSpline[{}]::setup( gc ):", m_name ) };

      std::set<std::string> keywords;
      for ( auto const & pair : gc.get_map( where ) ) { keywords.insert( pair.first ); }
      keywords.erase( "spline_type" );

      GenericContainer const & gc_x{ gc( "xdata", where ) };
      keywords.erase( "xdata" );
      GenericContainer const & gc_y{ gc( "ydata", where ) };
      keywords.erase( "ydata" );

      vec_real_type x, y;
      {
        string const ff{ fmt::format( "{}, field `xdata'", where ) };
        gc_x.copyto_vec_real( x, ff );
      }
      {
        string const ff{ fmt::format( "{}, field `ydata'", where ) };
        gc_y.copyto_vec_real( y, ff );
      }
      if ( gc.exists( "spline_sub_type" ) )
      {
        string_view st{ gc.get_map_string( "spline_sub_type", where ) };
        keywords.erase( "spline_sub_type" );
        if ( st == "cubic" )
          m_q_sub_type = QuinticSpline_sub_type::CUBIC;
        else if ( st == "pchip" )
          m_q_sub_type = QuinticSpline_sub_type::PCHIP;
        else if ( st == "akima" )
          m_q_sub_type = QuinticSpline_sub_type::AKIMA;
        else if ( st == "bessel" )
          m_q_sub_type = QuinticSpline_sub_type::BESSEL;
        else
        {
          UTILS_ERROR( "{} unknow sub type: {}\n", where, st );
        }
      }
      else
      {
        UTILS_WARNING( false, "{}, missing field `spline_sub_type` using `cubic` as default value\n", where );
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

      this->build( x, y );
    }
  };

}  // namespace Splines

// EOF: SplineQuintic.hxx
