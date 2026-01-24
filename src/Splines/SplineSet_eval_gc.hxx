
/**
 * \brief Evaluate all splines at x and store results in GenericContainer
 *
 * Results are stored as a map where keys are spline names and values
 * are the evaluated spline values.
 *
 * \param[in] x Parameter value
 * \param[out] gc GenericContainer to store results
 */
void eval( real_type const x, GenericContainer & gc ) const
{
  map_type & vals = gc.set_map();
  for ( auto const & [fst, snd] : m_header_to_position ) vals[fst] = m_splines[snd]->eval( x );
}

/**
 * \brief Evaluate first derivatives of all splines and store in GenericContainer
 *
 * \param[in] x Parameter value
 * \param[out] gc GenericContainer with derivative values
 */
void eval_D( real_type const x, GenericContainer & gc ) const
{
  map_type & vals = gc.set_map();
  for ( auto const & [fst, snd] : m_header_to_position ) vals[fst] = m_splines[snd]->eval_D( x );
}

/**
 * \brief Evaluate second derivatives of all splines and store in GenericContainer
 *
 * \param[in] x Parameter value
 * \param[out] gc GenericContainer with second derivative values
 */
void eval_DD( real_type const x, GenericContainer & gc ) const
{
  map_type & vals = gc.set_map();
  for ( auto const & [fst, snd] : m_header_to_position ) vals[fst] = m_splines[snd]->eval_DD( x );
}

/**
 * \brief Evaluate third derivatives of all splines and store in GenericContainer
 *
 * \param[in] x Parameter value
 * \param[out] gc GenericContainer with third derivative values
 */
void eval_DDD( real_type const x, GenericContainer & gc ) const
{
  map_type & vals{ gc.set_map() };
  for ( auto const & [fst, snd] : m_header_to_position ) vals[fst] = m_splines[snd]->eval_DDD( x );
}

/**
 * \brief Evaluate all splines at multiple x-values and store in GenericContainer
 *
 * For each spline, creates a vector of values corresponding to each x in `vec`.
 *
 * \param[in] vec Vector of x-values
 * \param[out] gc GenericContainer with spline names as keys and vectors as values
 */
void eval( vec_real_type const & vec, GenericContainer & gc ) const
{
  integer const npts = static_cast<integer>( vec.size() );
  map_type &    vals = gc.set_map();
  for ( auto const & [fst, snd] : m_header_to_position )
  {
    vec_real_type & v     = vals[fst].set_vec_real( snd );
    Spline const *  p_spl = m_splines[snd].get();
    for ( integer i = 0; i < npts; ++i ) v[i] = p_spl->eval( vec[i] );
  }
}

/**
 * \brief Evaluate first derivatives of all splines at multiple x-values
 *
 * \param[in] vec Vector of x-values
 * \param[out] gc GenericContainer with vectors of derivatives
 */
void eval_D( vec_real_type const & vec, GenericContainer & gc ) const
{
  integer const npts = static_cast<integer>( vec.size() );
  map_type &    vals = gc.set_map();
  for ( auto const & [fst, snd] : m_header_to_position )
  {
    vec_real_type & v     = vals[fst].set_vec_real( npts );
    Spline const *  p_spl = m_splines[snd].get();
    for ( integer i = 0; i < npts; ++i ) v[i] = p_spl->eval_D( vec[i] );
  }
}

/**
 * \brief Evaluate second derivatives of all splines at multiple x-values
 *
 * \param[in] vec Vector of x-values
 * \param[out] gc GenericContainer with vectors of second derivatives
 */
void eval_DD( vec_real_type const & vec, GenericContainer & gc ) const
{
  integer const npts = static_cast<integer>( vec.size() );
  map_type &    vals = gc.set_map();
  for ( auto const & [fst, snd] : m_header_to_position )
  {
    vec_real_type & v     = vals[fst].set_vec_real( npts );
    auto &          p_spl = m_splines[snd];
    for ( integer i = 0; i < npts; ++i ) v[i] = p_spl->eval_DD( vec[i] );
  }
}

/**
 * \brief Evaluate third derivatives of all splines at multiple x-values
 *
 * \param[in] vec Vector of x-values
 * \param[out] gc GenericContainer with vectors of third derivatives
 */
void eval_DDD( vec_real_type const & vec, GenericContainer & gc ) const
{
  integer const npts = static_cast<integer>( vec.size() );
  map_type &    vals = gc.set_map();
  for ( auto const & [fst, snd] : m_header_to_position )
  {
    vec_real_type & v     = vals[fst].set_vec_real( npts );
    Spline const *  p_spl = m_splines[snd].get();
    for ( integer i = 0; i < npts; ++i ) v[i] = p_spl->eval_DDD( vec[i] );
  }
}

/**
 * \brief Evaluate specific splines at x and store in GenericContainer
 *
 * \param[in] x Parameter value
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with requested spline values
 */
void eval( real_type const x, vec_string_type const & columns, GenericContainer & gc ) const
{
  map_type & vals = gc.set_map();
  for ( auto const & S : columns )
  {
    Spline const * p_spl = get_spline( S );
    vals[S]              = p_spl->eval( x );
  }
}

/**
 * \brief Evaluate first derivatives of specific splines
 *
 * \param[in] x Parameter value
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with derivative values
 */
void eval_D( real_type const x, vec_string_type const & columns, GenericContainer & gc ) const
{
  map_type & vals = gc.set_map();
  for ( auto const & S : columns )
  {
    Spline const * p_spl = get_spline( S );
    vals[S]              = p_spl->eval_D( x );
  }
}

/**
 * \brief Evaluate second derivatives of specific splines
 *
 * \param[in] x Parameter value
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with second derivative values
 */
void eval_DD( real_type const x, vec_string_type const & columns, GenericContainer & gc ) const
{
  map_type & vals = gc.set_map();
  for ( auto const & S : columns )
  {
    Spline const * p_spl = get_spline( S );
    vals[S]              = p_spl->eval_DD( x );
  }
}

/**
 * \brief Evaluate third derivatives of specific splines
 *
 * \param[in] x Parameter value
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with third derivative values
 */
void eval_DDD( real_type const x, vec_string_type const & columns, GenericContainer & gc ) const
{
  map_type & vals = gc.set_map();
  for ( auto const & S : columns )
  {
    Spline const * p_spl = get_spline( S );
    vals[S]              = p_spl->eval_DDD( x );
  }
}

/**
 * \brief Evaluate specific splines at multiple x-values and store in GenericContainer
 *
 * \param[in] vec Vector of x-values
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with vectors of values for requested splines
 */
void eval( vec_real_type const & vec, vec_string_type const & columns, GenericContainer & gc ) const
{
  integer const npts = static_cast<integer>( vec.size() );
  map_type &    vals = gc.set_map();
  for ( auto const & S : columns )
  {
    Spline const *  p_spl = get_spline( S );
    vec_real_type & v     = vals[S].set_vec_real( npts );
    for ( integer i = 0; i < npts; ++i ) v[i] = p_spl->eval( vec[i] );
  }
}

/**
 * \brief Evaluate first derivatives of specific splines at multiple x-values
 *
 * \param[in] vec Vector of x-values
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with vectors of derivatives
 */
void eval_D( vec_real_type const & vec, vec_string_type const & columns, GenericContainer & gc ) const
{
  integer const npts = static_cast<integer>( vec.size() );
  map_type &    vals = gc.set_map();
  for ( auto const & S : columns )
  {
    Spline const *  p_spl = get_spline( S );
    vec_real_type & v     = vals[S].set_vec_real( npts );
    for ( integer i = 0; i < npts; ++i ) v[i] = p_spl->eval_D( vec[i] );
  }
}
/**
 * \brief Evaluate second derivatives of specific splines at multiple x-values
 *
 * \param[in] vec Vector of x-values
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with vectors of second derivatives
 */
void eval_DD( vec_real_type const & vec, vec_string_type const & columns, GenericContainer & gc ) const
{
  integer const npts = static_cast<integer>( vec.size() );
  map_type &    vals = gc.set_map();
  for ( auto const & S : columns )
  {
    Spline const *  p_spl = get_spline( S );
    vec_real_type & v     = vals[S].set_vec_real( npts );
    for ( integer i = 0; i < npts; ++i ) v[i] = p_spl->eval_DD( vec[i] );
  }
}

/**
 * \brief Evaluate third derivatives of specific splines at multiple x-values
 *
 * \param[in] vec Vector of x-values
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with vectors of third derivatives
 */
void eval_DDD( vec_real_type const & vec, vec_string_type const & columns, GenericContainer & gc ) const
{
  integer const npts = static_cast<integer>( vec.size() );
  map_type &    vals = gc.set_map();
  for ( auto const & S : columns )
  {
    Spline const *  p_spl = get_spline( S );
    vec_real_type & v     = vals[S].set_vec_real( npts );
    for ( integer i = 0; i < npts; ++i ) v[i] = p_spl->eval_DDD( vec[i] );
  }
}
