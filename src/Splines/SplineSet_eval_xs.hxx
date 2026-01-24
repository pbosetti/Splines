/**
 * \brief Evaluate a spline by name at given x
 *
 * \param[in] x Parameter value
 * \param[in] name Name of the spline
 * \return Value of spline `name` at parameter `x`
 *
 * \throw UTILS_ASSERT if spline name not found
 */
real_type eval( real_type const x, string_view name ) const
{
  Spline const * S{ this->get_spline( name ) };
  return S->eval( x );
}

/**
 * \brief Evaluate first derivative of a spline by name
 *
 * \param[in] x Parameter value
 * \param[in] name Name of the spline
 * \return First derivative of spline `name` at `x`
 */
real_type eval_D( real_type const x, string_view name ) const
{
  Spline const * S{ this->get_spline( name ) };
  return S->D( x );
}

/**
 * \brief Evaluate second derivative of a spline by name
 *
 * \param[in] x Parameter value
 * \param[in] name Name of the spline
 * \return Second derivative of spline `name` at `x`
 */
real_type eval_DD( real_type const x, string_view name ) const
{
  Spline const * S{ this->get_spline( name ) };
  return S->DD( x );
}

/**
 * \brief Evaluate third derivative of a spline by name
 *
 * \param[in] x Parameter value
 * \param[in] name Name of the spline
 * \return Third derivative of spline `name` at `x`
 */
real_type eval_DDD( real_type const x, string_view name ) const
{
  Spline const * S{ this->get_spline( name ) };
  return S->DDD( x );
}

/**
 * \brief Evaluate fourth derivative of a spline by name
 *
 * \param[in] x Parameter value
 * \param[in] name Name of the spline
 * \return Fourth derivative of spline `name` at `x`
 */
real_type eval_DDDD( real_type const x, string_view name ) const
{
  Spline const * S{ this->get_spline( name ) };
  return S->DDDD( x );
}

/**
 * \brief Evaluate fifth derivative of a spline by name
 *
 * \param[in] x Parameter value
 * \param[in] name Name of the spline
 * \return Fifth derivative of spline `name` at `x`
 */
real_type eval_DDDDD( real_type const x, string_view name ) const
{
  Spline const * S{ this->get_spline( name ) };
  return S->DDDDD( x );
}

#ifdef AUTODIFF_SUPPORT
/**
 * \brief Evaluate spline by name with first-order automatic differentiation
 *
 * \param[in] x Dual number parameter
 * \param[in] name Name of the spline
 * \return Dual number containing value and gradient
 */
autodiff::dual1st eval( autodiff::dual1st const & x, string_view name ) const
{
  using autodiff::derivative;
  using autodiff::dual1st;
  real_type xv{ val( x ) };
  dual1st   res{ eval( xv, name ) };
  res.grad = eval_D( xv, name ) * x.grad;
  return res;
}

/**
 * \brief Evaluate spline by name with second-order automatic differentiation
 *
 * \param[in] x Dual number parameter (second order)
 * \param[in] name Name of the spline
 * \return Dual number containing value, first and second derivatives
 */
autodiff::dual2nd eval( autodiff::dual2nd const & x, string_view name ) const
{
  using autodiff::derivative;
  using autodiff::dual2nd;

  real_type xv{ val( x ) };
  real_type xg{ val( x.grad ) };
  real_type dfx{ eval_D( xv, name ) };
  real_type dxx{ eval_DD( xv, name ) };
  dual2nd   res{ eval( xv, name ) };

  res.grad      = dfx * xg;
  res.grad.grad = dfx * x.grad.grad + dxx * ( xg * xg );
  return res;
}

/**
 * \brief Generic template for automatic differentiation by name
 *
 * \tparam T Input type
 * \param[in] x Parameter value
 * \param[in] name Name of the spline
 * \return Result matching input type
 */
template <typename T> auto eval( T const & x, string_view name ) const
{
  if constexpr ( std::is_arithmetic<T>::value )
  {
    // Se T è un tipo numerico (int, float, double, etc.), promuovi a real_type
    return eval( static_cast<real_type>( x ), name );
  }
  else
  {
    // Altrimenti deduce automaticamente il tipo duale appropriato
    return eval( autodiff::detail::to_dual( x ), name );
  }
}

/**
 * \brief Generic function call operator for autodiff by name
 *
 * \tparam T Parameter type
 * \param[in] x Parameter value
 * \param[in] name Name of the spline
 * \return Result matching input type
 */
template <typename T> auto operator()( T const & x, string_view name ) const -> decltype( eval( x, name ) )
{
  return eval( x, name );
}

#endif
