/**
 * PYTHON Wrapper for Splines
 *
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich, Enrico Bertolazzi
 */

#ifndef PY_SPLINES_HH
#define PY_SPLINES_HH

#include <string>
#include <vector>

#include <Splines.hh>
#include <pybind11/pybind11.h>

namespace pySpline
{
  using Splines::integer;
  using Splines::ostream_type;
  using Splines::real_type;

  using Splines::integer;
  using Splines::real_type;
  using Splines::Spline;

  using pybind11::module;

  using GC_namespace::GenericContainer;

  class PythonicSpline : public Spline
  {
  public:
    PythonicSpline( std::string const & name = "Spline" ) : Spline( name ) {}

    void reserve( integer npts ) override { PYBIND11_OVERLOAD_PURE( void, Spline, reserve, npts ); }

    void build() override { PYBIND11_OVERLOAD_PURE( void, Spline, build ); }

    void setup( GenericContainer const & gc ) override { PYBIND11_OVERLOAD_PURE( void, Spline, setup, gc ); }

    void clear() override { PYBIND11_OVERLOAD_PURE( void, Spline, clear ); }

    Splines::SplineType1D type() const override { PYBIND11_OVERLOAD_PURE( Splines::SplineType1D, Spline, type ); }

    real_type eval( real_type x ) const override { PYBIND11_OVERLOAD_PURE( real_type, Spline, eval, x ); }

    real_type D( real_type x ) const override { PYBIND11_OVERLOAD_PURE( real_type, Spline, D, x ); }

    real_type DD( real_type x ) const override { PYBIND11_OVERLOAD_PURE( real_type, Spline, DD, x ); }

    real_type DDD( real_type x ) const override { PYBIND11_OVERLOAD_PURE( real_type, Spline, DDD, x ); }

    real_type DDDD( real_type x ) const override { PYBIND11_OVERLOAD_PURE( real_type, Spline, DDDD, x ); }

    real_type DDDDD( real_type x ) const override { PYBIND11_OVERLOAD_PURE( real_type, Spline, DDDDD, x ); }

    void D( real_type x, real_type dd[2] ) const override { PYBIND11_OVERLOAD_PURE( void, Spline, D, x, dd ); }

    void DD( real_type x, real_type dd[3] ) const override { PYBIND11_OVERLOAD_PURE( void, Spline, DD, x, dd ); }

    real_type id_eval( integer ni, real_type x ) const override { PYBIND11_OVERLOAD_PURE( real_type, Spline, id_eval, ni, x ); }

    real_type id_D( integer ni, real_type x ) const override { PYBIND11_OVERLOAD_PURE( real_type, Spline, id_D, ni, x ); }

    real_type id_DD( integer ni, real_type x ) const override { PYBIND11_OVERLOAD_PURE( real_type, Spline, id_DD, ni, x ); }

    real_type id_DDD( integer ni, real_type x ) const override { PYBIND11_OVERLOAD_PURE( real_type, Spline, id_DDD, ni, x ); }

    integer coeffs( real_type cfs[], real_type nodes[], bool transpose = false ) const override
    { PYBIND11_OVERLOAD_PURE( integer, Spline, coeffs, cfs, nodes, transpose ); }

    integer order() const override { PYBIND11_OVERLOAD_PURE( integer, Spline, order ); }

    void write_to_stream( ostream_type & s ) const override
    { PYBIND11_OVERLOAD_PURE( void, Spline, write_to_stream, s ); }
  };

  void python_register_splines_class( module & m );
  void python_register_hermite_functions( module & m );
}  // namespace pySpline

#endif /* PY_SPLINES_HH */
