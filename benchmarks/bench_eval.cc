/*\
 |  Splines micro-benchmark — scalar evaluation hot paths.
 |
 |  Baseline harness for the "faster" work in MODERNIZATION.md (Tier 3):
 |  it times the current per-point virtual-dispatch path (operator(), D, DD)
 |  over a large array of query points. When the batch-eval API (3a) and
 |  reciprocal hoisting (3b) land, re-run this and compare — the scalar loop
 |  here is exactly the baseline those changes must beat / not regress.
 |
 |  Usage: bench_eval [N_knots] [M_queries] [reps]
 |         defaults: 1000 knots, 2,000,000 queries, best of 5
\*/

#include "Splines/Splines.hh"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <random>
#include <vector>

using Splines::integer;
using Splines::real_type;

namespace
{
  using clock_type = std::chrono::steady_clock;

  // Run `fn` (evaluates all queries, returns a checksum) `reps` times and
  // return the best (min) ns-per-op. `sink` accumulates the checksum so the
  // optimizer cannot eliminate the work.
  template <typename F>
  double best_ns_per_op( F const & fn, std::size_t n_ops, int reps, double & sink )
  {
    double best = std::numeric_limits<double>::infinity();
    for ( int r = 0; r < reps; ++r )
    {
      auto const   t0 = clock_type::now();
      sink += fn();
      auto const   t1 = clock_type::now();
      double const ns = std::chrono::duration<double, std::nano>( t1 - t0 ).count();
      best            = std::min( best, ns / static_cast<double>( n_ops ) );
    }
    return best;
  }

  void bench_one(
    char const *                   name,
    Splines::Spline &              s,
    std::vector<real_type> const & xq,
    int                            reps,
    double &                       sink )
  {
    std::size_t const M = xq.size();

    double const ns_eval = best_ns_per_op(
      [&] { double a = 0; for ( real_type const x : xq ) a += s( x ); return a; }, M, reps, sink );
    double const ns_D = best_ns_per_op(
      [&] { double a = 0; for ( real_type const x : xq ) a += s.D( x ); return a; }, M, reps, sink );
    double const ns_DD = best_ns_per_op(
      [&] { double a = 0; for ( real_type const x : xq ) a += s.DD( x ); return a; }, M, reps, sink );

    std::printf(
      "%-10s %11.2f %11.2f %11.2f %13.1f\n",
      name, ns_eval, ns_D, ns_DD, 1000.0 / ns_eval );
  }
}  // namespace

int main( int argc, char ** argv )
{
  integer     N    = 1000;       // number of knots
  std::size_t M    = 2'000'000;  // number of query points
  int         reps = 5;          // best-of-reps

  if ( argc > 1 ) N    = std::atoi( argv[1] );
  if ( argc > 2 ) M    = static_cast<std::size_t>( std::atoll( argv[2] ) );
  if ( argc > 3 ) reps = std::atoi( argv[3] );

  // Knots: strictly increasing x with random spacing, smooth y.
  std::vector<real_type>                     X( static_cast<std::size_t>( N ) );
  std::vector<real_type>                     Y( static_cast<std::size_t>( N ) );
  std::mt19937_64                            rng( 12345 );
  std::uniform_real_distribution<real_type>  gap( 0.2, 1.0 );
  real_type                                  acc = 0;
  for ( integer i = 0; i < N; ++i )
  {
    X[static_cast<std::size_t>( i )] = acc;
    Y[static_cast<std::size_t>( i )] = std::sin( 0.3 * acc ) + 0.1 * acc;
    acc += gap( rng );
  }

  // Query points: uniform over the knot range, generated up-front so the
  // RNG is never inside a timed loop.
  std::vector<real_type>                    xq( M );
  std::uniform_real_distribution<real_type> q( X.front(), X.back() );
  for ( real_type & v : xq ) v = q( rng );

  std::printf(
    "Splines micro-benchmark  (N=%d knots, M=%zu queries, best of %d, ns/op)\n",
    N, M, reps );
  std::printf( "%-10s %11s %11s %11s %13s\n", "type", "eval", "D", "DD", "Meval/s" );

  double sink = 0;

  Splines::LinearSpline li;
  li.build( X, Y );
  bench_one( "Linear", li, xq, reps, sink );

  Splines::CubicSpline cs;
  cs.build( X, Y );
  bench_one( "Cubic", cs, xq, reps, sink );

  Splines::AkimaSpline ak;
  ak.build( X, Y );
  bench_one( "Akima", ak, xq, reps, sink );

  Splines::PchipSpline pc;
  pc.build( X, Y );
  bench_one( "Pchip", pc, xq, reps, sink );

  Splines::QuinticSpline qs( Splines::Spline_sub_type::CUBIC );
  qs.build( X, Y );
  bench_one( "Quintic", qs, xq, reps, sink );

  std::printf( "checksum: %.6g\n", sink );
  return 0;
}
