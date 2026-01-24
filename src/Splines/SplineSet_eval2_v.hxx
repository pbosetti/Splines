
/**
 * \brief Evaluate all splines using a specific spline as independent variable (C array version)
 *
 * \param[in] spl Index of independent spline
 * \param[in] zeta Target value of independent spline
 * \param[out] x Computed x-value such that spline[spl](x) = zeta
 * \param[out] vals Output array to store spline values
 * \param[in] incy Stride between consecutive values (default: 1)
 */
void eval2( integer const spl, real_type const zeta, real_type & x, real_type vals[], integer const incy = 1 ) const
{
  intersect( spl, zeta, x );
  size_t ii{ 0 };
  for ( integer i = 0; i < m_nspl; ++i, ii += incy ) vals[ii] = m_splines[i]->eval( x );
}

/**
 * \brief Evaluate first derivatives using a specific spline as independent variable (C array version)
 *
 * Computes derivatives with respect to the independent spline value, not x.
 *
 * \param[in] spl Index of independent spline
 * \param[in] zeta Target value of independent spline
 * \param[out] x Computed x-value
 * \param[out] vals Output array for derivatives
 * \param[in] incy Stride between consecutive values (default: 1)
 *
 * \par Mathematical Formula:
 * \f[
 * \frac{dy_i}{d\zeta} = \frac{dy_i/dx}{ds/dx}
 * \f]
 * where s is the independent spline.
 */
void eval2_D( integer const spl, real_type const zeta, real_type & x, real_type vals[], integer const incy = 1 ) const
{
  Spline const *  S{ intersect( spl, zeta, x ) };
  real_type const ds{ S->D( x ) };
  size_t          ii{ 0 };
  for ( integer i = 0; i < m_nspl; ++i, ii += incy ) vals[ii] = m_splines[i]->D( x ) / ds;
}

/**
 * \brief Evaluate second derivatives using a specific spline as independent variable (C array version)
 *
 * \param[in] spl Index of independent spline
 * \param[in] zeta Target value of independent spline
 * \param[out] x Computed x-value
 * \param[out] vals Output array for second derivatives
 * \param[in] incy Stride between consecutive values (default: 1)
 *
 * \par Mathematical Formula:
 * \f[
 * \frac{d^2y_i}{d\zeta^2} = \frac{d^2y_i/dx^2}{(ds/dx)^2} - \frac{dy_i/dx \cdot d^2s/dx^2}{(ds/dx)^3}
 * \f]
 */
void eval2_DD( integer const spl, real_type const zeta, real_type & x, real_type vals[], integer const incy = 1 ) const
{
  Spline const *  S{ intersect( spl, zeta, x ) };
  real_type const dt{ 1 / S->D( x ) };
  real_type const dt2{ dt * dt };
  real_type const ddt{ -S->DD( x ) * ( dt * dt2 ) };
  size_t          ii{ 0 };
  for ( integer i = 0; i < m_nspl; ++i, ii += incy )
  {
    auto & Si{ m_splines[i] };
    vals[ii] = Si->DD( x ) * dt2 + Si->D( x ) * ddt;
  }
}

/**
 * \brief Evaluate third derivatives using a specific spline as independent variable (C array version)
 *
 * \param[in] spl Index of independent spline
 * \param[in] zeta Target value of independent spline
 * \param[out] x Computed x-value
 * \param[out] vals Output array for third derivatives
 * \param[in] incy Stride between consecutive values (default: 1)
 *
 * \par Mathematical Formula:
 * Uses chain rule for third derivatives with respect to parameterized variable.
 */
void eval2_DDD( integer const spl, real_type const zeta, real_type & x, real_type vals[], integer const incy = 1 ) const
{
  Spline const *  S{ intersect( spl, zeta, x ) };
  real_type const dt{ 1 / S->D( x ) };
  real_type const dt3{ dt * dt * dt };
  real_type const ddt{ -S->DD( x ) * dt3 };
  real_type const dddt{ 3 * ( ddt * ddt ) / dt - S->DDD( x ) * ( dt * dt3 ) };
  size_t          ii{ 0 };
  for ( integer i = 0; i < m_nspl; ++i, ii += incy )
  {
    auto & Si{ m_splines[i] };
    vals[ii] = Si->DDD( x ) * dt3 + 3 * Si->DD( x ) * dt * ddt + Si->D( x ) * dddt;
  }
}
