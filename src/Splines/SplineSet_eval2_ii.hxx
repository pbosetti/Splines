/**
 * \brief Evaluate a specific spline using another spline as independent variable (index version)
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline
 * \param[out] x Computed x-value
 * \param[in] spl Index of spline to evaluate
 * \return Value of spline `spl` at computed x
 */
real_type eval2( real_type const zeta, integer const indep, real_type & x, integer const spl ) const
{
  intersect( indep, zeta, x );
  return m_splines[spl]->eval( x );
}

/**
 * \brief Evaluate first derivative using another spline as independent variable (index version)
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline
 * \param[out] x Computed x-value
 * \param[in] spl Index of spline to evaluate
 * \return First derivative with respect to independent spline value
 */
real_type eval2_D( real_type const zeta, integer const indep, real_type & x, integer const spl ) const
{
  Spline const * S{ intersect( indep, zeta, x ) };
  return m_splines[spl]->D( x ) / S->D( x );
}

/**
 * \brief Evaluate second derivative using another spline as independent variable (index version)
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline
 * \param[out] x Computed x-value
 * \param[in] spl Index of spline to evaluate
 * \return Second derivative with respect to independent spline value
 */
real_type eval2_DD( real_type const zeta, integer const indep, real_type & x, integer const spl ) const
{
  Spline const *  S{ intersect( indep, zeta, x ) };
  real_type const dt{ 1 / S->D( x ) };
  real_type const dt2{ dt * dt };
  real_type const ddt{ -S->DD( x ) * ( dt * dt2 ) };
  auto &          SPL{ m_splines[spl] };
  return SPL->DD( x ) * dt2 + SPL->D( x ) * ddt;
}

/**
 * \brief Evaluate third derivative using another spline as independent variable (index version)
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline
 * \param[out] x Computed x-value
 * \param[in] spl Index of spline to evaluate
 * \return Third derivative with respect to independent spline value
 */
real_type eval2_DDD( real_type const zeta, integer const indep, real_type & x, integer const spl ) const
{
  Spline const *  S{ intersect( indep, zeta, x ) };
  real_type const dt{ 1 / S->D( x ) };
  real_type const dt3{ dt * dt * dt };
  real_type const ddt{ -S->DD( x ) * dt3 };
  real_type const dddt{ 3 * ( ddt * ddt ) / dt - S->DDD( x ) * ( dt * dt3 ) };

  auto & SPL{ m_splines[spl] };
  return SPL->DDD( x ) * dt3 + 3 * SPL->DD( x ) * dt * ddt + SPL->D( x ) * dddt;
}

/**
 * \brief Evaluate a specific spline using another spline as independent variable (index version)
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline
 * \param[in] spl Index of spline to evaluate
 * \return Value of spline `spl` at computed x
 */
real_type eval2( real_type const zeta, integer const indep, integer const spl ) const
{
  real_type x;
  return eval2( zeta, indep, x, spl );
}

/**
 * \brief Evaluate first derivative using another spline as independent variable (index version)
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline
 * \param[in] spl Index of spline to evaluate
 * \return First derivative with respect to independent spline value
 */
real_type eval2_D( real_type const zeta, integer const indep, integer const spl ) const
{
  real_type x;
  return eval2_D( zeta, indep, x, spl );
}

/**
 * \brief Evaluate second derivative using another spline as independent variable (index version)
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline
 * \param[in] spl Index of spline to evaluate
 * \return Second derivative with respect to independent spline value
 */
real_type eval2_DD( real_type const zeta, integer const indep, integer const spl ) const
{
  real_type x;
  return eval2_DD( zeta, indep, x, spl );
}

/**
 * \brief Evaluate third derivative using another spline as independent variable (index version)
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline
 * \param[in] spl Index of spline to evaluate
 * \return Third derivative with respect to independent spline value
 */
real_type eval2_DDD( real_type const zeta, integer const indep, integer const spl ) const
{
  real_type x;
  return eval2_DDD( zeta, indep, x, spl );
}

#ifdef AUTODIFF_SUPPORT

/**
 * \brief Evaluate with first-order automatic differentiation (index version)
 *
 * \param[in] zeta Dual number target value
 * \param[in] indep Index of independent spline
 * \param[out] x Computed x-value (real type)
 * \param[in] spl Index of spline to evaluate
 * \return Dual number containing value and derivative
 */
autodiff::dual1st eval2( autodiff::dual1st const & zeta, integer const indep, real_type & x, integer const spl ) const
{
  using autodiff::derivative;
  using autodiff::dual1st;
  real_type zv{ val( zeta ) };
  dual1st   res{ eval2( zv, indep, x, spl ) };
  res.grad = eval_D( x, spl ) * zeta.grad;  // x gia calcolato
  return res;
}

/**
 * \brief Evaluate with first-order automatic differentiation (index version)
 *
 * \param[in] zeta Dual number target value
 * \param[in] indep Index of independent spline
 * \param[in] spl Index of spline to evaluate
 * \return Dual number containing value and derivative
 */
autodiff::dual1st eval2( autodiff::dual1st const & zeta, integer const indep, integer const spl ) const
{
  real_type x;
  return eval2( zeta, indep, x, spl );
}

/**
 * \brief Evaluate with second-order automatic differentiation (index version)
 *
 * \param[in] zeta Dual number target value
 * \param[in] indep Index of independent spline
 * \param[out] x Computed x-value (real type)
 * \param[in] spl Index of spline to evaluate
 * \return Dual number containing value, first and second derivatives
 */
autodiff::dual2nd eval2( autodiff::dual2nd const & zeta, integer const indep, real_type & x, integer const spl ) const
{
  using autodiff::derivative;
  using autodiff::dual2nd;

  real_type zv{ val( zeta ) };
  real_type zg{ val( zeta.grad ) };
  real_type dfx{ eval2_D( zv, indep, x, spl ) };
  real_type dxx{ eval_DD( x, spl ) };
  dual2nd   res{ eval( x, spl ) };

  res.grad      = dfx * zg;
  res.grad.grad = dfx * zeta.grad.grad + dxx * ( zg * zg );
  return res;
}

/**
 * \brief Evaluate with second-order automatic differentiation (index version)
 *
 * \param[in] zeta Dual number target value
 * \param[in] indep Index of independent spline
 * \param[in] spl Index of spline to evaluate
 * \return Dual number containing value, first and second derivatives
 */
autodiff::dual2nd eval2( autodiff::dual2nd const & zeta, integer const indep, integer const spl ) const
{
  real_type x;
  return eval2( zeta, indep, x, spl );
}

/**
 * \brief Generic template for automatic differentiation (index version)
 *
 * \tparam T Input type
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline
 * \param[out] x Computed x-value (real type)
 * \param[in] spl Index of spline to evaluate
 * \return Result matching input type
 */
template <typename T> auto eval2( T const & zeta, integer const indep, real_type & x, integer const spl ) const
{
  if constexpr ( std::is_arithmetic<T>::value )
  {
    // Se T è un tipo numerico (int, float, double, etc.), promuovi a real_type
    return eval2( static_cast<real_type>( zeta ), indep, x, spl );
  }
  else
  {
    // Altrimenti deduce automaticamente il tipo duale appropriato
    return eval2( autodiff::detail::to_dual( zeta ), indep, x, spl );
  }
}

/**
 * \brief Generic template for automatic differentiation (index version)
 *
 * \tparam T Input type
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline
 * \param[in] spl Index of spline to evaluate
 * \return Result matching input type
 */
template <typename T> auto eval2( T const & zeta, integer const indep, integer const spl ) const
{
  real_type x;
  return eval2( zeta, indep, x, spl );
}

#endif
