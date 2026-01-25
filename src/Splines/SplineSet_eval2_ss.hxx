
/**
 * \brief Evaluate a specific spline using another spline as independent variable
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value
 * \param[in] name Name of spline to evaluate
 * \return Value of spline `name` at computed x
 */
real_type eval2( real_type const zeta, string_view const indep, real_type & x, string_view const name ) const
{
  return this->eval2( zeta, this->get_position( indep ), x, this->get_position( name ) );
}

/**
 * \brief Evaluate first derivative using another spline as independent variable
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value
 * \param[in] name Name of spline to evaluate
 * \return First derivative of spline `name` with respect to independent spline value
 */
real_type eval2_D( real_type const zeta, string_view const indep, real_type & x, string_view const name ) const
{
  return this->eval2_D( zeta, this->get_position( indep ), x, this->get_position( name ) );
}

/**
 * \brief Evaluate second derivative using another spline as independent variable
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value
 * \param[in] name Name of spline to evaluate
 * \return Second derivative of spline `name` with respect to independent spline value
 */
real_type eval2_DD( real_type const zeta, string_view const indep, real_type & x, string_view const name ) const
{
  return this->eval2_DD( zeta, this->get_position( indep ), x, this->get_position( name ) );
}

/**
 * \brief Evaluate third derivative using another spline as independent variable
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value
 * \param[in] name Name of spline to evaluate
 * \return Third derivative of spline `name` with respect to independent spline value
 */
real_type eval2_DDD( real_type const zeta, string_view const indep, real_type & x, string_view const name ) const
{
  return this->eval2_DDD( zeta, this->get_position( indep ), x, this->get_position( name ) );
}

/**
 * \brief Evaluate a specific spline using another spline as independent variable
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[in] name Name of spline to evaluate
 * \return Value of spline `name` at computed x
 */
real_type eval2( real_type const zeta, string_view const indep, string_view const name ) const
{
  real_type x;
  return this->eval2( zeta, indep, x, name );
}

/**
 * \brief Evaluate first derivative using another spline as independent variable
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[in] name Name of spline to evaluate
 * \return First derivative of spline `name` with respect to independent spline value
 */
real_type eval2_D( real_type const zeta, string_view const indep, string_view const name ) const
{
  real_type x;
  return this->eval2_D( zeta, indep, x, name );
}

/**
 * \brief Evaluate second derivative using another spline as independent variable
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[in] name Name of spline to evaluate
 * \return Second derivative of spline `name` with respect to independent spline value
 */
real_type eval2_DD( real_type const zeta, string_view const indep, string_view const name ) const
{
  real_type x;
  return this->eval2_DD( zeta, indep, x, name );
}

/**
 * \brief Evaluate third derivative using another spline as independent variable
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[in] name Name of spline to evaluate
 * \return Third derivative of spline `name` with respect to independent spline value
 */
real_type eval2_DDD( real_type const zeta, string_view const indep, string_view const name ) const
{
  real_type x;
  return this->eval2_DDD( zeta, indep, x, name );
}

#ifdef AUTODIFF_SUPPORT

/**
 * \brief Evaluate with first-order automatic differentiation using independent spline
 *
 * \param[in] zeta Dual number target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value (real type)
 * \param[in] name Name of spline to evaluate
 * \return Dual number containing value and derivative
 */
autodiff::dual1st eval2(
  autodiff::dual1st const & zeta,
  string_view const         indep,
  real_type &               x,
  string_view const         name ) const
{
  using autodiff::derivative;
  using autodiff::dual1st;
  real_type zv{ val( zeta ) };
  dual1st   res{ eval2( zv, indep, x, name ) };
  res.grad = eval_D( x, name ) * zeta.grad;  // x gia calcolato
  return res;
}

/**
 * \brief Evaluate with first-order automatic differentiation using independent spline
 *
 * \param[in] zeta Dual number target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[in] name Name of spline to evaluate
 * \return Dual number containing value and derivative
 */
autodiff::dual1st eval2(
  autodiff::dual1st const & zeta,
  string_view const         indep,
  string_view const         name ) const
{
  real_type x;
  return eval2( zeta, indep, x, name );
}

/**
 * \brief Evaluate with second-order automatic differentiation using independent spline
 *
 * \param[in] zeta Dual number target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value (real type)
 * \param[in] name Name of spline to evaluate
 * \return Dual number containing value, first and second derivatives
 */
autodiff::dual2nd eval2(
  autodiff::dual2nd const & zeta,
  string_view const         indep,
  real_type &               x,
  string_view const         name ) const
{
  using autodiff::derivative;
  using autodiff::dual2nd;

  real_type zv{ val( zeta ) };
  real_type zg{ val( zeta.grad ) };
  real_type dfx{ eval2_D( zv, indep, x, name ) };
  real_type dxx{ eval_DD( x, name ) };  // x gia calcolato
  dual2nd   res{ eval( x, name ) };

  res.grad      = dfx * zg;
  res.grad.grad = dfx * zeta.grad.grad + dxx * ( zg * zg );
  return res;
}

/**
 * \brief Evaluate with second-order automatic differentiation using independent spline
 *
 * \param[in] zeta Dual number target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[in] name Name of spline to evaluate
 * \return Dual number containing value, first and second derivatives
 */
autodiff::dual2nd eval2(
  autodiff::dual2nd const & zeta,
  string_view const         indep,
  string_view const         name ) const
{
  real_type x;
  return eval2( zeta, indep, x, name );
}

/**
 * \brief Generic template for automatic differentiation using independent spline
 *
 * \tparam T Input type
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value (real type)
 * \param[in] name Name of spline to evaluate
 * \return Result matching input type
 */
template <typename T> auto eval2( T const & zeta, string_view const indep, real_type & x, string_view const name ) const
{
  if constexpr ( std::is_arithmetic<T>::value )
  {
    // Se T è un tipo numerico (int, float, double, etc.), promuovi a real_type
    return eval2( static_cast<real_type>( zeta ), indep, x, name );
  }
  else
  {
    // Altrimenti deduce automaticamente il tipo duale appropriato
    return eval2( autodiff::detail::to_dual( zeta ), indep, x, name );
  }
}

/**
 * \brief Generic template for automatic differentiation using independent spline
 *
 * \tparam T Input type
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[in] name Name of spline to evaluate
 * \return Result matching input type
 */
template <typename T> auto eval2( T const & zeta, string_view const indep, string_view const name ) const
{
  real_type x;
  return eval2<T>( zeta, indep, x, name );
}

#endif
