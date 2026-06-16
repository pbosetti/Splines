/**
 * PYTHON Wrapper for Splines
 *
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich, Enrico Bertolazzi
 */

#include <string>

#include <Splines.hh>
#include "py-VanLeerSplines.hh"

namespace pySpline
{
  namespace py = pybind11;

  using Splines::Spline;
  using Splines::VanLeerSpline;

  void python_register_vanleer_splines_class( py::module & m )
  {
    py::class_<VanLeerSpline, Spline>( m, "VanLeerSpline" )
      .def( py::init<std::string const &>(), py::arg( "name" ) = "VanLeerSpline" );
  }

}  // namespace pySpline
