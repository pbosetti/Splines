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
 |   ____      _     _      ____        _ _
 |  |  _ \ ___| |__ (_)_ __/ ___| _ __ | (_)_ __   ___
 |  | |_) / __| '_ \| | '_ \___ \| '_ \| | | '_ \ / _ \
 |  |  __/ (__| | | | | |_) |__) | |_) | | | | | |  __/
 |  |_|   \___|_| |_|_| .__/____/| .__/|_|_|_| |_|\___|
 |                    |_|        |_|
\*/

#ifndef SPLINE_PCHIP_HH
#define SPLINE_PCHIP_HH

namespace Splines
{
  /*!
   * \brief Computes derivatives for a Monotone Piecewise Cubic Interpolation (PCHIP).
   *
   * This function constructs the first derivatives \f$ d_k \f$ at the data points such that
   * the resulting piecewise cubic Hermite spline preserves the monotonicity of the given data.
   * Unlike a natural cubic spline, PCHIP avoids "overshoots" and oscillations locally.
   *
   * \details
   * The algorithm calculates the derivatives \f$ Y'_k \f$ based on the secant slopes between points.
   *
   * <b>Interior Points:</b>
   * For interior points \f$ k \f$, the derivative is computed using a weighted harmonic mean
   * of the adjacent secant slopes \f$ \delta_{k-1} \f$ and \f$ \delta_k \f$.
   * - If \f$ \delta_{k-1} \f$ and \f$ \delta_k \f$ have opposite signs (or one is zero),
   * the derivative \f$ d_k \f$ is set to 0 (local extremum).
   * - If they have the same sign, the derivative is determined by the Fritsch-Butland formula
   * (modified by Brodlie) to ensure strict monotonicity.
   *
   * The formula for the weighted harmonic mean is:
   * \f[
   * d_k = \frac{\delta_{k-1} \cdot \delta_k}{\alpha \delta_k + \beta \delta_{k-1}}
   * \f]
   * where weights \f$ \alpha, \beta \f$ depend on the interval lengths \f$ h_{k-1}, h_k \f$.
   *
   * <b>Boundary Points:</b>
   * Derivatives at the endpoints (start and end) are computed using a non-centered three-point
   * formula that is subsequently adjusted (clamped) to enforce monotonicity constraints if necessary.
   *
   * \param[in]  X    Array of x-coordinates (must be strictly increasing).
   * \param[in]  Y    Array of y-coordinates.
   * \param[out] Yp   Output array for the computed derivatives (must be allocated with size npts).
   * \param[in]  npts Number of points (dimension of arrays). Must be >= 2.
   *
   * \note
   * This implementation uses a sliding window approach to minimize memory access overhead
   * and redundant calculations.
   *
   * \references
   * - F.N. Fritsch, R.E. Carlson: <i>Monotone Piecewise Cubic Interpolation</i>,
   * SIAM J. Numer. Anal. Vol 17, No. 2, April 1980.
   * - F.N. Fritsch and J. Butland: <i>A method for constructing local monotone piecewise cubic interpolants</i>,
   * SIAM Journal on Scientific and Statistical Computing 5, 2 (June 1984), pp. 300-304.
   */
  inline void Pchip_build( real_type const X[], real_type const Y[], real_type Yp[], integer const npts )
  {
    UTILS_ASSERT( npts >= 2, "Pchip_build, npts={} must be >= 2\n", npts );

    integer const n = npts - 1;  // Index of the last point

    // Helper lambda: returns 1 if signs match, -1 if opposite, 0 if either is zero.
    auto sign_test = []( real_type const a, real_type const b ) -> int
    {
      if ( a * b > 0 ) return 1;
      return ( a == 0 || b == 0 ) ? 0 : -1;
    };

    // Initialize first interval variables
    real_type h_cur   = X[1] - X[0];
    real_type del_cur = ( Y[1] - Y[0] ) / h_cur;

    // Special case: only 2 points (Linear Interpolation)
    if ( n == 1 )
    {
      Yp[0] = Yp[1] = del_cur;
      return;
    }

    /* -----------------------------------------------------------
     * 1. Left Boundary (Start Point)
     * Use non-centered 3-point formula, then adjust for shape.
     * ----------------------------------------------------------- */
    real_type h_next   = X[2] - X[1];
    real_type del_next = ( Y[2] - Y[1] ) / h_next;

    real_type h_sum = h_cur + h_next;
    // Weights for the non-centered formula
    real_type w1 = ( h_cur + h_sum ) / h_sum;
    real_type w2 = -h_cur / h_sum;

    Yp[0] = w1 * del_cur + w2 * del_next;

    // Check monotonicity constraints for the start point
    if ( sign_test( Yp[0], del_cur ) <= 0 ) { Yp[0] = 0; }
    else if ( sign_test( del_cur, del_next ) < 0 )
    {
      // If convexity changes, clamp the derivative to prevent overshoot
      real_type const dmax = 3 * del_cur;
      if ( std::abs( Yp[0] ) > std::abs( dmax ) ) Yp[0] = dmax;
    }

    /* -----------------------------------------------------------
     * 2. Interior Points Loop
     * Uses a sliding window: prev -> cur -> next
     * ----------------------------------------------------------- */
    real_type h_prev   = h_cur;
    real_type del_prev = del_cur;

    // Shift current to next for the start of the loop
    h_cur   = h_next;
    del_cur = del_next;

    for ( integer i = 1; i < n; ++i )
    {
      // Prefetch 'next' values if not at the second-to-last point
      if ( i < n - 1 )
      {
        h_next   = X[i + 1] - X[i];
        del_next = ( Y[i + 1] - Y[i] ) / h_next;
      }

      // Check the relationship between adjacent slopes
      int const sign_prod = sign_test( del_prev, del_cur );

      if ( sign_prod > 0 )
      {
        // Case: Strict Monotonicity (Same signs)
        // Use Brodlie modification of Fritsch-Butland formula (Weighted Harmonic Mean)
        real_type const sum_h = h_prev + h_cur;

        // Calculate weights based on interval lengths
        real_type const wa = ( h_prev + 2 * h_cur ) / ( 3 * sum_h );
        real_type const wb = ( 2 * h_prev + h_cur ) / ( 3 * sum_h );

        // Harmonic mean avoids division by zero implicitly
        Yp[i] = ( del_prev * del_cur ) / ( wa * del_cur + wb * del_prev );
      }
      else
      {
        // Case: Local Extremum (Opposite signs) or Plateau (Zero slope)
        // Derivative must be zero to preserve monotonicity
        Yp[i] = 0;
      }

      // Slide the window for the next iteration
      if ( i < n - 1 )
      {
        h_prev   = h_cur;
        del_prev = del_cur;
        h_cur    = h_next;
        del_cur  = del_next;
      }
    }

    /* -----------------------------------------------------------
     * 3. Right Boundary (End Point)
     * Re-evaluation of last intervals required for the 3-point formula.
     * Note: h_prev and del_prev currently hold values for interval [n-1, n].
     * We need interval [n-2, n-1] as well.
     * ----------------------------------------------------------- */

    // Retriev last interval (already in h_prev, but explicit for clarity)
    real_type const h_last   = X[n] - X[n - 1];
    real_type const del_last = ( Y[n] - Y[n - 1] ) / h_last;

    // Retrieve second-to-last interval
    real_type const h_prev_last   = X[n - 1] - X[n - 2];
    real_type const del_prev_last = ( Y[n - 1] - Y[n - 2] ) / h_prev_last;

    h_sum = h_last + h_prev_last;
    w1    = -h_last / h_sum;
    w2    = ( h_last + h_sum ) / h_sum;

    Yp[n] = w1 * del_prev_last + w2 * del_last;

    // Check monotonicity constraints for the end point
    if ( sign_test( Yp[n], del_last ) <= 0 ) { Yp[n] = 0; }
    else if ( sign_test( del_prev_last, del_last ) < 0 )
    {
      real_type const dmax = 3 * del_last;
      if ( std::abs( Yp[n] ) > std::abs( dmax ) ) Yp[n] = dmax;
    }
  }

  //! Pchip (Piecewise Cubic Hermite Interpolating Polynomial) spline class
  class PchipSpline : public CubicSplineBase
  {
  public:
    using CubicSplineBase::build;
    using CubicSplineBase::reserve;

    //!
    //! Build an empty spline of `PchipSpline` type
    //!
    //! \param name the name of the spline
    //!
    explicit PchipSpline( string_view name = "PchipSpline" ) : CubicSplineBase( name ) {}

    //!
    //! Spline destructor.
    //!
    ~PchipSpline() override {}

    //! Return spline type (as number)
    SplineType1D type() const override { return SplineType1D::PCHIP; }

    // --------------------------- VIRTUALS -----------------------------------

    //! Build a Monotone spline from previously inserted points
    void build() override
    {
      string msg{ fmt::format( "PchipSpline[{}]::build():", m_name ) };
      UTILS_ASSERT( m_npts > 1, "{} npts = {} not enought points\n", msg, m_npts );

      Utils::check_NaN( m_X, msg + " X", m_npts, __LINE__, __FILE__ );
      Utils::check_NaN( m_Y, msg + " Y", m_npts, __LINE__, __FILE__ );

      integer ibegin{ 0 };
      integer iend{ 0 };

      do
      {
        // cerca intervallo monotono strettamente crescente
        for ( ++iend; iend < m_npts && m_X[iend - 1] < m_X[iend]; ++iend ) {}
        Pchip_build( m_X + ibegin, m_Y + ibegin, m_Yp + ibegin, iend - ibegin );
        ibegin = iend;
      } while ( iend < m_npts );

      Utils::check_NaN( m_Yp, msg + " Yp", m_npts, __LINE__, __FILE__ );
      m_search.must_reset();
    }

    void setup( GenericContainer const & gc ) override
    {
      /*
      // gc["xdata"]
      // gc["ydata"]
      //
      */
      string const where{ fmt::format( "PchipSpline[{}]::setup( gc ):", m_name ) };

      std::set<std::string> keywords;
      for ( auto const & pair : gc.get_map( where ) ) { keywords.insert( pair.first ); }
      keywords.erase( "spline_type" );

      GenericContainer const & gc_x{ gc( "xdata", where ) };
      keywords.erase( "xdata" );
      GenericContainer const & gc_y{ gc( "ydata", where ) };
      keywords.erase( "ydata" );

      std::vector<real_type> x, y;
      {
        string const ff{ fmt::format( "{}, field `xdata'", where ) };
        gc_x.copyto_vec_real( x, ff );
      }
      {
        string const ff{ fmt::format( "{}, field `ydata'", where ) };
        gc_y.copyto_vec_real( y, ff );
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

#endif

// EOF: SplinePchip.hxx
