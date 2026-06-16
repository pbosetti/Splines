/**
 * PYTHON Wrapper for Splines
 *
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich, Enrico Bertolazzi
 */

#include "py-Splines.hh"
#include "py-CubicSplines.hh"
#include "py-ConstantSplines.hh"
#include "py-AkimaSplines.hh"
#include "py-VanLeerSplines.hh"
#include "py-HermiteSplines.hh"
#include "py-LinearSplines.hh"
#include "py-PchipSplines.hh"
#include "py-QuinticSplines.hh"

#include <pybind11/stl.h>
#include <tuple>

namespace pySpline
{
  using Splines::integer;
  using Splines::real_type;

  using Splines::CubicSpline_BC;
  using Splines::Spline;
  using Splines::Spline_sub_type;
  using Splines::SplineType1D;
  using Splines::SplineType2D;

  namespace py = pybind11;

  void python_register_splines_class( py::module & m )
  {
    py::enum_<SplineType1D>( m, "SplineType1D" )
      .value( "CONSTANT_TYPE", SplineType1D::CONSTANT )
      .value( "LINEAR_TYPE", SplineType1D::LINEAR )
      .value( "CUBIC_TYPE", SplineType1D::CUBIC )
      .value( "AKIMA_TYPE", SplineType1D::AKIMA )
      .value( "VANLEER_TYPE", SplineType1D::VANLEER )
      .value( "PCHIP_TYPE", SplineType1D::PCHIP )
      .value( "QUINTIC_CUBIC_TYPE", SplineType1D::QUINTIC_CUBIC )
      .value( "QUINTIC_AKIMA_TYPE", SplineType1D::QUINTIC_AKIMA )
      .value( "QUINTIC_VANLEER_TYPE", SplineType1D::QUINTIC_VANLEER )
      .value( "QUINTIC_PCHIP_TYPE", SplineType1D::QUINTIC_PCHIP )
      .value( "HERMITE_TYPE", SplineType1D::HERMITE )
      .value( "SPLINE_SET_TYPE", SplineType1D::SPLINE_SET )
      .value( "SPLINE_VEC_TYPE", SplineType1D::SPLINE_VEC )
      .export_values();

    py::enum_<SplineType2D>( m, "SplineType2D" )
      .value( "BILINEAR_TYPE", SplineType2D::BILINEAR )
      .value( "BICUBIC_CUBIC_TYPE", SplineType2D::BICUBIC_CUBIC )
      .value( "BICUBIC_AKIMA_TYPE", SplineType2D::BICUBIC_AKIMA )
      .value( "BICUBIC_VANLEER_TYPE", SplineType2D::BICUBIC_VANLEER )
      .value( "BICUBIC_PCHIP_TYPE", SplineType2D::BICUBIC_PCHIP )
      .value( "BIQUINTIC_CUBIC_TYPE", SplineType2D::BIQUINTIC_CUBIC )
      .value( "BIQUINTIC_AKIMA_TYPE", SplineType2D::BIQUINTIC_AKIMA )
      .value( "BIQUINTIC_VANLEER_TYPE", SplineType2D::BIQUINTIC_VANLEER )
      .value( "BIQUINTIC_PCHIP_TYPE", SplineType2D::BIQUINTIC_PCHIP )
      .export_values();

    py::enum_<CubicSpline_BC>( m, "CubicSpline_BC" )
      .value( "EXTRAPOLATE", CubicSpline_BC::EXTRAPOLATE )
      .value( "NATURAL", CubicSpline_BC::NATURAL )
      .value( "PARABOLIC_RUNOUT", CubicSpline_BC::PARABOLIC_RUNOUT )
      .value( "NOT_A_KNOT", CubicSpline_BC::NOT_A_KNOT )
      .export_values();

    py::enum_<Spline_sub_type>( m, "Spline_sub_type" )
      .value( "CUBIC", Spline_sub_type::CUBIC )
      .value( "PCHIP", Spline_sub_type::PCHIP )
      .value( "AKIMA", Spline_sub_type::AKIMA )
      .value( "VANLEER", Spline_sub_type::VANLEER )
      .export_values();

    py::class_<Spline, PythonicSpline>( m, "Spline" )
      .def( "name", &Spline::name )
      .def( "is_closed", &Spline::is_closed )
      .def( "make_closed", &Spline::make_closed )
      .def( "make_opened", &Spline::make_opened )
      .def( "is_bounded", &Spline::is_bounded )
      .def( "make_unbounded", &Spline::make_unbounded )
      .def( "make_bounded", &Spline::make_bounded )
      .def( "num_points", &Spline::num_points )
      .def( "x_node", &Spline::x_node )
      .def( "y_node", &Spline::y_node )
      .def( "x_begin", &Spline::x_begin )
      .def( "y_begin", &Spline::y_begin )
      .def( "x_end", &Spline::x_end )
      .def( "y_end", &Spline::y_end )
      .def( "reserve", &Spline::reserve )
      .def( "push_back", &Spline::push_back )
      .def( "drop_back", &Spline::drop_back )
      .def( "build", []( Spline & self ) { self.build(); } )
      .def(
        "build",
        []( Spline & self, std::vector<real_type> const & x, std::vector<real_type> const & y ) { self.build( x, y ); } )
      .def( "clear", &Spline::clear )
      .def( "x_min", &Spline::x_min )
      .def( "x_max", &Spline::x_max )
      .def( "y_min", &Spline::y_min )
      .def( "y_max", &Spline::y_max )
      .def( "type", &Spline::type )
      .def( "set_origin", &Spline::set_origin )
      .def( "set_range", &Spline::set_range )
      .def( "__call__", []( Spline const & self, real_type x ) { return self( x ); } )
      .def( "D", []( Spline const & self, real_type x ) { return self.D( x ); } )
      .def( "DD", []( Spline const & self, real_type x ) { return self.DD( x ); } )
      .def( "DDD", &Spline::DDD )
      .def( "DDDD", &Spline::DDDD )
      .def( "DDDDD", &Spline::DDDDD )
      .def( "eval", []( Spline const & self, real_type x ) { return self.eval( x ); } )
      .def( "eval_D", &Spline::eval_D )
      .def( "eval_DD", &Spline::eval_DD )
      .def( "eval_DDD", &Spline::eval_DDD )
      .def( "eval_DDDD", &Spline::eval_DDDD )
      .def( "eval_DDDDD", &Spline::eval_DDDDD )
      .def(
        "coeffs",
        []( Spline & self, bool transpose ) -> std::tuple<std::vector<real_type>, std::vector<real_type> >
        {
          integer                n_nodes = self.num_points() - 1;
          integer                n_order = self.order();
          std::vector<real_type> cfs( n_order * n_nodes );
          std::vector<real_type> nodes( n_nodes );
          self.coeffs( cfs.data(), nodes.data(), transpose );
          return std::make_tuple( cfs, nodes );
        },
        py::arg( "transpose" ) = false )
      .def( "order", &Spline::order )
      .def( "type_name", &Spline::type_name )
      .def(
        "__str__",
        []( const Spline & self ) -> std::string
        {
          std::ostringstream str;
          self.info( str );
          return str.str();
        } )
      .def(
        "write_to_stream",
        []( const Spline & self ) -> std::string
        {
          std::ostringstream str;
          self.write_to_stream( str );
          return str.str();
        } )
      .def(
        "dump",
        []( const Spline & self, integer nintervals, std::string header ) -> std::string
        {
          std::ostringstream str;
          self.dump( str, nintervals, header.c_str() );
          return str.str();
        },
        py::arg( "nintervals" ),
        py::arg( "header" ) = "x\ty" );
  }

  void python_register_hermite_functions( module & m )
  {
    m.def(
       "Hermite3",
       []( real_type const x, real_type const H ) -> std::vector<real_type>
       {
         std::vector<real_type> result( 4 );
         Splines::Hermite3( x, H, result.data() );
         return result;
       } )
      .def(
        "Hermite3_D",
        []( real_type const x, real_type const H ) -> std::vector<real_type>
        {
          std::vector<real_type> result( 4 );
          Splines::Hermite3_D( x, H, result.data() );
          return result;
        } )
      .def(
        "Hermite3_DD",
        []( real_type const x, real_type const H ) -> std::vector<real_type>
        {
          std::vector<real_type> result( 4 );
          Splines::Hermite3_DD( x, H, result.data() );
          return result;
        } )
      .def(
        "Hermite3_DDD",
        []( real_type const x, real_type const H ) -> std::vector<real_type>
        {
          std::vector<real_type> result( 4 );
          Splines::Hermite3_DDD( x, H, result.data() );
          return result;
        } )
      .def(
        "Hermite5",
        []( real_type const x, real_type const H ) -> std::vector<real_type>
        {
          std::vector<real_type> result( 6 );
          Splines::Hermite5( x, H, result.data() );
          return result;
        } )
      .def(
        "Hermite5_D",
        []( real_type const x, real_type const H ) -> std::vector<real_type>
        {
          std::vector<real_type> result( 6 );
          Splines::Hermite5_D( x, H, result.data() );
          return result;
        } )
      .def(
        "Hermite5_DD",
        []( real_type const x, real_type const H ) -> std::vector<real_type>
        {
          std::vector<real_type> result( 6 );
          Splines::Hermite5_DD( x, H, result.data() );
          return result;
        } )
      .def(
        "Hermite5_DDD",
        []( real_type const x, real_type const H ) -> std::vector<real_type>
        {
          std::vector<real_type> result( 6 );
          Splines::Hermite5_DDD( x, H, result.data() );
          return result;
        } )
      .def(
        "Hermite5_DDDD",
        []( real_type const x, real_type const H ) -> std::vector<real_type>
        {
          std::vector<real_type> result( 6 );
          Splines::Hermite5_DDDD( x, H, result.data() );
          return result;
        } )
      .def(
        "Hermite5_DDDDD",
        []( real_type const x, real_type const H ) -> std::vector<real_type>
        {
          std::vector<real_type> result( 6 );
          Splines::Hermite5_DDDDD( x, H, result.data() );
          return result;
        } );
  }

}  // namespace pySpline

PYBIND11_MODULE( Splines, m )
{
  pySpline::python_register_splines_class( m );
  pySpline::python_register_cubic_splines_base_class( m );
  pySpline::python_register_cubic_splines_class( m );
  pySpline::python_register_constant_splines_class( m );
  pySpline::python_register_akima_splines_class( m );
  pySpline::python_register_vanleer_splines_class( m );
  pySpline::python_register_hermite_splines_class( m );
  pySpline::python_register_linear_splines_class( m );
  pySpline::python_register_pchip_splines_class( m );
  pySpline::python_register_quintic_splines_base_class( m );
  pySpline::python_register_quintic_splines_class( m );
  pySpline::python_register_hermite_functions( m );
}
