
/**
 * \brief Evaluate all splines using independent spline and store in GenericContainer
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline
 * \param[out] x Computed x-value
 * \param[out] gc GenericContainer to store results
 */
void eval2( real_type const zeta, integer const indep, real_type & x, GenericContainer & gc ) const
{
  map_type & vals = gc.set_map();
  intersect( indep, zeta, x );
  for ( auto const & [fst, snd] : m_header_to_position ) vals[fst] = m_splines[snd]->eval( x );
}

/**
 * \brief Evaluate first derivatives using independent spline
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline
 * \param[out] x Computed x-value
 * \param[out] gc GenericContainer with derivative values
 */
void eval2_D( real_type const zeta, integer const indep, real_type & x, GenericContainer & gc ) const
{
  map_type & vals = gc.set_map();
  intersect( indep, zeta, x );
  for ( auto const & [fst, snd] : m_header_to_position ) vals[fst] = m_splines[snd]->eval_D( x );
}

/**
 * \brief Evaluate second derivatives using independent spline
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline
 * \param[out] x Computed x-value
 * \param[out] gc GenericContainer with second derivative values
 */
void eval2_DD( real_type const zeta, integer const indep, real_type & x, GenericContainer & gc ) const
{
  map_type & vals = gc.set_map();
  intersect( indep, zeta, x );
  for ( auto const & [fst, snd] : m_header_to_position ) vals[fst] = m_splines[snd]->eval_DD( x );
}

/**
 * \brief Evaluate third derivatives using independent spline
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline
 * \param[out] x Computed x-value
 * \param[out] gc GenericContainer with third derivative values
 */
void eval2_DDD( real_type const zeta, integer const indep, real_type & x, GenericContainer & gc ) const
{
  map_type & vals = gc.set_map();
  intersect( indep, zeta, x );
  for ( auto const & [fst, snd] : m_header_to_position ) vals[fst] = m_splines[snd]->eval_DDD( x );
}

/**
 * \brief Evaluate all splines at multiple zeta values using independent spline
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Index of independent spline
 * \param[out] gc GenericContainer with results
 */
void eval2( vec_real_type const & zetas, integer const indep, GenericContainer & gc ) const
{
  integer const npts = static_cast<integer>( zetas.size() );
  map_type &    vals = gc.set_map();

  // pre-allocation
  for ( auto const & [fst, snd] : m_header_to_position ) vals[fst].set_vec_real( npts );

  for ( integer i = 0; i < npts; ++i )
  {
    real_type x;
    intersect( indep, zetas[i], x );
    for ( auto const & [fst, snd] : m_header_to_position )
    {
      vec_real_type & v{ vals[fst].get_vec_real() };
      v[i] = m_splines[snd]->eval( x );
    }
  }
}

/**
 * \brief Evaluate first derivatives at multiple zeta values using independent spline
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Index of independent spline
 * \param[out] gc GenericContainer with vectors of derivatives
 */
void eval2_D( vec_real_type const & zetas, integer const indep, GenericContainer & gc ) const
{
  integer const npts = static_cast<integer>( zetas.size() );
  map_type &    vals = gc.set_map();

  // pre-allocation
  for ( auto const & [fst, snd] : m_header_to_position ) vals[fst].set_vec_real( npts );

  for ( integer i = 0; i < npts; ++i )
  {
    real_type x;
    intersect( indep, zetas[i], x );
    for ( auto const & [fst, snd] : m_header_to_position )
    {
      vec_real_type & v = vals[fst].get_vec_real();
      v[i]              = m_splines[snd]->eval_D( x );
    }
  }
}

/**
 * \brief Evaluate second derivatives at multiple zeta values using independent spline
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Index of independent spline
 * \param[out] gc GenericContainer with vectors of second derivatives
 */
void eval2_DD( vec_real_type const & zetas, integer const indep, GenericContainer & gc ) const
{
  integer const npts = static_cast<integer>( zetas.size() );
  map_type &    vals = gc.set_map();

  // pre-allocation
  for ( auto const & [fst, snd] : m_header_to_position ) vals[fst].set_vec_real( npts );

  for ( integer i = 0; i < npts; ++i )
  {
    real_type x;
    intersect( indep, zetas[i], x );
    for ( auto const & [fst, snd] : m_header_to_position )
    {
      vec_real_type & v = vals[fst].get_vec_real();
      v[i]              = m_splines[snd]->eval_DD( x );
    }
  }
}


/**
 * \brief Evaluate third derivatives at multiple zeta values using independent spline
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Index of independent spline
 * \param[out] gc GenericContainer with vectors of third derivatives
 */
void eval2_DDD( vec_real_type const & zetas, integer const indep, GenericContainer & gc ) const
{
  integer const npts = static_cast<integer>( zetas.size() );
  map_type &    vals = gc.set_map();

  // pre-allocation
  for ( auto const & [fst, snd] : m_header_to_position ) vals[fst].set_vec_real( npts );

  for ( integer i = 0; i < npts; ++i )
  {
    real_type x;
    intersect( indep, zetas[i], x );
    for ( auto const & [fst, snd] : m_header_to_position )
    {
      vec_real_type & v = vals[fst].get_vec_real();
      v[i]              = m_splines[snd]->eval_DDD( x );
    }
  }
}


/**
 * \brief Evaluate specific splines using independent spline and store in GenericContainer
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline
 * \param[out] x Computed x-value
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer to store results
 */
void eval2(
  real_type const         zeta,
  integer const           indep,
  real_type &             x,
  vec_string_type const & columns,
  GenericContainer &      gc ) const
{
  map_type & vals = gc.set_map();
  intersect( indep, zeta, x );
  for ( auto const & S : columns )
  {
    Spline const * p_spl = get_spline( S );
    vals[S]              = p_spl->eval( x );
  }
}

/**
 * \brief Evaluate first derivatives of specific splines using independent spline
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline
 * \param[out] x Computed x-value
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with derivative values
 */
void eval2_D(
  real_type const         zeta,
  integer const           indep,
  real_type &             x,
  vec_string_type const & columns,
  GenericContainer &      gc ) const
{
  map_type & vals = gc.set_map();
  intersect( indep, zeta, x );
  for ( auto const & S : columns )
  {
    Spline const * p_spl = get_spline( S );
    vals[S]              = p_spl->eval_D( x );
  }
}

/**
 * \brief Evaluate second derivatives of specific splines using independent spline
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline
 * \param[out] x Computed x-value
 * \param[in] columns Names of splines to evaluate
 * \param[out] vals GenericContainer with second derivative values
 */
void eval2_DD(
  real_type const         zeta,
  integer const           indep,
  real_type &             x,
  vec_string_type const & columns,
  GenericContainer &      gc ) const
{
  map_type & vals = gc.set_map();
  intersect( indep, zeta, x );
  for ( auto const & S : columns )
  {
    Spline const * p_spl = get_spline( S );
    vals[S]              = p_spl->eval_DD( x );
  }
}

//!
//! Evaluate all the splines at `zeta` and fill a map of values
//! in a GenericContainer with keys in `columns` and `indep`
//! as independent spline
//!
void eval2_DDD(
  real_type const         zeta,
  integer const           indep,
  real_type &             x,
  vec_string_type const & columns,
  GenericContainer &      gc ) const
{
  map_type & vals = gc.set_map();
  intersect( indep, zeta, x );
  for ( auto const & S : columns )
  {
    Spline const * p_spl = get_spline( S );
    vals[S]              = p_spl->eval_DDD( x );
  }
}


/**
 * \brief Evaluate specific splines at multiple zeta values using independent spline
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Index of independent spline
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with results
 */
void eval2( vec_real_type const & zetas, integer const indep, vec_string_type const & columns, GenericContainer & gc )
  const
{
  integer const npts = static_cast<integer>( zetas.size() );
  map_type &    vals = gc.set_map();

  // pre-allocation
  for ( auto const & S : columns ) vals[S].set_vec_real( npts );

  for ( integer i = 0; i < npts; ++i )
  {
    real_type x;
    intersect( indep, zetas[i], x );
    for ( auto const & S : columns )
    {
      vec_real_type & v{ vals[S].get_vec_real() };
      Spline const *  p_spl{ get_spline( S ) };
      v[i] = p_spl->eval( x );
    }
  }
}

/**
 * \brief Evaluate first derivatives of specific splines at multiple zeta values
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Index of independent spline
 * \param[in] columns Names of splines to evaluate
 * \param[out] vals GenericContainer with vectors of derivatives
 */
void eval2_D( vec_real_type const & zetas, integer const indep, vec_string_type const & columns, GenericContainer & gc )
  const
{
  integer const npts = static_cast<integer>( zetas.size() );
  map_type &    vals = gc.set_map();

  // pre-allocation
  for ( auto const & S : columns ) vals[S].set_vec_real( npts );

  for ( integer i = 0; i < npts; ++i )
  {
    real_type x;
    intersect( indep, zetas[i], x );
    for ( auto const & S : columns )
    {
      vec_real_type & v     = vals[S].get_vec_real();
      Spline const *  p_spl = get_spline( S );
      v[i]                  = p_spl->eval_D( x );
    }
  }
}

/**
 * \brief Evaluate second derivatives of specific splines at multiple zeta values using independent spline
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Index of independent spline
 * \param[in] columns Names of splines to evaluate
 * \param[out] vals GenericContainer with vectors of second derivatives
 */
void eval2_DD(
  vec_real_type const &   zetas,
  integer const           indep,
  vec_string_type const & columns,
  GenericContainer &      gc ) const
{
  integer const npts = static_cast<integer>( zetas.size() );
  map_type &    vals = gc.set_map();

  // pre-allocation
  for ( auto const & S : columns ) vals[S].set_vec_real( npts );

  for ( integer i = 0; i < npts; ++i )
  {
    real_type x;
    intersect( indep, zetas[i], x );
    for ( auto const & S : columns )
    {
      vec_real_type & v     = vals[S].get_vec_real();
      Spline const *  p_spl = get_spline( S );
      v[i]                  = p_spl->eval_DD( x );
    }
  }
}

/**
 * \brief Evaluate third derivatives of specific splines at multiple zeta values using independent spline
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Index of independent spline
 * \param[in] columns Names of splines to evaluate
 * \param[out] vals GenericContainer with vectors of third derivatives
 */
void eval2_DDD(
  vec_real_type const &   zetas,
  integer const           indep,
  vec_string_type const & columns,
  GenericContainer &      gc ) const
{
  integer const npts = static_cast<integer>( zetas.size() );
  map_type &    vals = gc.set_map();

  // pre-allocation
  for ( auto const & S : columns ) vals[S].set_vec_real( npts );

  for ( integer i = 0; i < npts; ++i )
  {
    real_type x;
    intersect( indep, zetas[i], x );
    for ( auto const & S : columns )
    {
      vec_real_type & v     = vals[S].get_vec_real();
      Spline const *  p_spl = get_spline( S );
      v[i]                  = p_spl->eval_DDD( x );
    }
  }
}

/**
 * \brief Evaluate all splines using independent spline by name
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value
 * \param[out] vals GenericContainer to store results
 */
void eval2( real_type const zeta, string_view indep, real_type & x, GenericContainer & vals ) const
{
  this->eval2( zeta, this->get_position( indep ), x, vals );
}

/**
 * \brief Evaluate first derivatives using independent spline by name
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value
 * \param[out] vals GenericContainer with derivative values
 */
void eval2_D( real_type const zeta, string_view indep, real_type & x, GenericContainer & vals ) const
{
  this->eval2_D( zeta, this->get_position( indep ), x, vals );
}

/**
 * \brief Evaluate second derivatives using independent spline by name
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value
 * \param[out] vals GenericContainer with second derivative values
 */
void eval2_DD( real_type const zeta, string_view indep, real_type & x, GenericContainer & vals ) const
{
  this->eval2_DD( zeta, this->get_position( indep ), x, vals );
}

/**
 * \brief Evaluate third derivatives using independent spline by name
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value
 * \param[out] vals GenericContainer with third derivative values
 */
void eval2_DDD( real_type const zeta, string_view indep, real_type & x, GenericContainer & vals ) const
{
  this->eval2_DDD( zeta, this->get_position( indep ), x, vals );
}


/**
 * \brief Evaluate all splines at multiple zeta values using independent spline by name
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Name of independent spline
 * \param[out] vals GenericContainer with results
 */
void eval2( vec_real_type const & zetas, string_view indep, GenericContainer & vals ) const
{
  this->eval2( zetas, this->get_position( indep ), vals );
}

/**
 * \brief Evaluate first derivatives at multiple zeta values using independent spline by name
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Name of independent spline
 * \param[out] vals GenericContainer with vectors of derivatives
 */
void eval2_D( vec_real_type const & zetas, string_view indep, GenericContainer & vals ) const
{
  this->eval2_D( zetas, this->get_position( indep ), vals );
}

/**
 * \brief Evaluate second derivatives at multiple zeta values using independent spline by name
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Name of independent spline
 * \param[out] vals GenericContainer with vectors of second derivatives
 */
void eval2_DD( vec_real_type const & zetas, string_view indep, GenericContainer & vals ) const
{
  this->eval2_DD( zetas, this->get_position( indep ), vals );
}

/**
 * \brief Evaluate third derivatives at multiple zeta values using independent spline by name
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Name of independent spline
 * \param[out] vals GenericContainer with vectors of third derivatives
 */
void eval2_DDD( vec_real_type const & zetas, string_view indep, GenericContainer & vals ) const
{
  this->eval2_DDD( zetas, this->get_position( indep ), vals );
}


/**
 * \brief Evaluate specific splines using independent spline by name
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value
 * \param[in] columns Names of splines to evaluate
 * \param[out] vals GenericContainer to store results
 */
void eval2(
  real_type const         zeta,
  string_view             indep,
  real_type &             x,
  vec_string_type const & columns,
  GenericContainer &      vals ) const
{
  this->eval2( zeta, this->get_position( indep ), x, columns, vals );
}

/**
 * \brief Evaluate first derivatives of specific splines using independent spline by name
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value
 * \param[in] columns Names of splines to evaluate
 * \param[out] vals GenericContainer with derivative values
 */
void eval2_D(
  real_type const         zeta,
  string_view             indep,
  real_type &             x,
  vec_string_type const & columns,
  GenericContainer &      vals ) const
{
  this->eval2_D( zeta, this->get_position( indep ), x, columns, vals );
}

/**
 * \brief Evaluate second derivatives of specific splines using independent spline by name
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value
 * \param[in] columns Names of splines to evaluate
 * \param[out] vals GenericContainer with second derivative values
 */
void eval2_DD(
  real_type const         zeta,
  string_view             indep,
  real_type &             x,
  vec_string_type const & columns,
  GenericContainer &      vals ) const
{
  this->eval2_DD( zeta, this->get_position( indep ), x, columns, vals );
}

/**
 * \brief Evaluate third derivatives of specific splines using independent spline by name
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value
 * \param[in] columns Names of splines to evaluate
 * \param[out] vals GenericContainer with third derivative values
 */
void eval2_DDD(
  real_type const         zeta,
  string_view             indep,
  real_type &             x,
  vec_string_type const & columns,
  GenericContainer &      vals ) const
{
  this->eval2_DDD( zeta, this->get_position( indep ), x, columns, vals );
}


/**
 * \brief Evaluate specific splines at multiple zeta values using independent spline by name
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Name of independent spline
 * \param[in] columns Names of splines to evaluate
 * \param[out] vals GenericContainer with results
 */
void eval2( vec_real_type const & zetas, string_view indep, vec_string_type const & columns, GenericContainer & vals )
  const
{
  this->eval2( zetas, this->get_position( indep ), columns, vals );
}

/**
 * \brief Evaluate first derivatives of specific splines at multiple zeta values using independent spline by name
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Name of independent spline
 * \param[in] columns Names of splines to evaluate
 * \param[out] vals GenericContainer with vectors of derivatives
 */
void eval2_D( vec_real_type const & zetas, string_view indep, vec_string_type const & columns, GenericContainer & vals )
  const
{
  this->eval2_D( zetas, this->get_position( indep ), columns, vals );
}

/**
 * \brief Evaluate second derivatives of specific splines at multiple zeta values using independent spline by name
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Name of independent spline
 * \param[in] columns Names of splines to evaluate
 * \param[out] vals GenericContainer with vectors of second derivatives
 */
void eval2_DD(
  vec_real_type const &   zetas,
  string_view             indep,
  vec_string_type const & columns,
  GenericContainer &      vals ) const
{
  this->eval2_DD( zetas, this->get_position( indep ), columns, vals );
}

/**
 * \brief Evaluate third derivatives of specific splines at multiple zeta values using independent spline by name
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Name of independent spline
 * \param[in] columns Names of splines to evaluate
 * \param[out] vals GenericContainer with vectors of second derivatives
 */
void eval2_DDD(
  vec_real_type const &   zetas,
  string_view             indep,
  vec_string_type const & columns,
  GenericContainer &      vals ) const
{
  this->eval2_DDD( zetas, this->get_position( indep ), columns, vals );
}
