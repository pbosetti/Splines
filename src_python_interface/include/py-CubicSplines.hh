/**
 * PYTHON Wrapper for Splines
 *
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich, Enrico Bertolazzi
 */

#ifndef PY_CUBIC_SPLINES_HH
#define PY_CUBIC_SPLINES_HH

#include <string>

#include <Splines.hh>
#include <pybind11/pybind11.h>

namespace pySpline
{
  using Splines::integer;
  using Splines::ostream_type;
  using Splines::real_type;

  using Splines::CubicSplineBase;

  using pybind11::module;

  class PythonicCubicSplineBase : public CubicSplineBase
  {
  public:
    PythonicCubicSplineBase( std::string const & name = "CubicSplineBase" ) : CubicSplineBase( name ) {}

    real_type eval( real_type x ) const override { PYBIND11_OVERLOAD( real_type, CubicSplineBase, eval, x ); }

    real_type D( real_type x ) const override { PYBIND11_OVERLOAD( real_type, CubicSplineBase, D, x ); }

    real_type DD( real_type x ) const override { PYBIND11_OVERLOAD( real_type, CubicSplineBase, DD, x ); }

    real_type DDD( real_type x ) const override { PYBIND11_OVERLOAD( real_type, CubicSplineBase, DDD, x ); }

    integer coeffs( real_type cfs[], real_type nodes[], bool transpose = false ) const override
    { PYBIND11_OVERLOAD( integer, CubicSplineBase, coeffs, cfs, nodes, transpose ); }

    integer order() const override { PYBIND11_OVERLOAD( integer, CubicSplineBase, order ); }

    void write_to_stream( ostream_type & s ) const override
    { PYBIND11_OVERLOAD( void, CubicSplineBase, write_to_stream, s ); }

    Splines::SplineType1D type() const override
    { PYBIND11_OVERLOAD_PURE( Splines::SplineType1D, CubicSplineBase, type ); }

    void reserve( integer npts ) override { PYBIND11_OVERLOAD( void, CubicSplineBase, reserve, npts ); }

    void build() override { PYBIND11_OVERLOAD_PURE( void, CubicSplineBase, build ); }

    void D( real_type x, real_type dd[2] ) const override { PYBIND11_OVERLOAD( void, CubicSplineBase, D, x, dd ); }

    void DD( real_type x, real_type dd[3] ) const override { PYBIND11_OVERLOAD( void, CubicSplineBase, DD, x, dd ); }

    real_type id_eval( integer ni, real_type x ) const override
    { PYBIND11_OVERLOAD( real_type, CubicSplineBase, id_eval, ni, x ); }

    real_type id_D( integer ni, real_type x ) const override
    { PYBIND11_OVERLOAD( real_type, CubicSplineBase, id_D, ni, x ); }

    real_type id_DD( integer ni, real_type x ) const override
    { PYBIND11_OVERLOAD( real_type, CubicSplineBase, id_DD, ni, x ); }

    real_type id_DDD( integer ni, real_type x ) const override
    { PYBIND11_OVERLOAD( real_type, CubicSplineBase, id_DDD, ni, x ); }
  };

  void python_register_cubic_splines_base_class( module & m );
  void python_register_cubic_splines_class( module & m );

}  // namespace pySpline

#endif /* PY_CUBIC_SPLINES_HH */
