/**
 * PYTHON Wrapper for Splines
 *
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich, Enrico Bertolazzi
 */

#include <string>

#include <Splines.hh>
#include "py-ConstantSplines.hh"

namespace pySpline
{
  namespace py = pybind11;

  using Splines::ConstantSpline;
  using Splines::Spline;

  void python_register_constant_splines_class( py::module & m )
  {
    py::class_<ConstantSpline, Spline>( m, "ConstantSpline" )
      .def( py::init<std::string const &>(), py::arg( "name" ) = "ConstantSpline" );
  }

}  // namespace pySpline