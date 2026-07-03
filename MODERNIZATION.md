# Splines — modernization plan (`src/` + paired headers)

Scope: the library code under [`src/`](src/) and its public headers under
[`include/Splines/`](include/Splines/). Goal: a **safer** and **faster**
library without changing the numerical results or the public API surface
(additions are fine; silent behavior changes are not).

This document is the agreed scope for a staged effort. Each tier is
independently shippable and gated by the existing 19-test suite
(`ctest`, `test00`..`test18`) on both the Linux VM (`orb -m mads`, `lbuild/`)
and macOS (`build/`).

---

## What is already solid (do **not** churn)

- **Interval search is already well-optimized.** `Utils::SearchInterval::find()`
  ([`Utils_search_intervals2.hh:544`](include/Splines/Splines.hh)) has a
  lock-free fast path (atomic `m_ready`, mutex taken only on the first call /
  after `must_reset`), table bucketing, and a bounded binary search. The
  per-evaluation cost is **not** dominated by locking or by the search — do
  not "optimize" it.
- **The C API is exception-safe at the boundary.** Every entry point routes
  through `c_api_call_int/real/cstr/ptr`, which catch `...`
  ([`SplinesCinterface.cc:56-106`](src/SplinesCinterface.cc)).
- **The pool allocator is RAII-owned and not leaking.** `Utils::Malloc`
  members free themselves; build paths carry broad assertion coverage
  (130 `UTILS_ASSERT`/`UTILS_ERROR` across `src/`).

---

## Tier 1 — Correctness / portability bugs (do first, low risk)

### 1a. Integer `abs()` applied to a `double`
- **Where:** [`src/SplineCubicBase.cc:262`](src/SplineCubicBase.cc#L262)
  ```cpp
  if ( abs( DP1 ) > ( m_X[i+1] - m_X[i-1] ) * epsi ) continue;  // DP1 is real_type&
  ```
- **Problem:** no `using std::abs` / `using namespace std` is in scope in this
  file, and `Splines.hh` does not pull in `std::abs`. Bare `abs` therefore
  binds to C's `::abs(int)`: the derivative is truncated to an integer before
  the magnitude is taken, so any `|DP1| < 1` collapses to `0` and this
  min/max-detection guard silently misfires. It only "works" today because
  libstdc++/libc++ happen to leak `::abs(double)` into the global namespace —
  precisely the class of construct that diverges on MSVC.
- **Fix:** use `std::abs`. Grep the whole tree for the same trap while here.
- **Risk:** minimal; changes a currently-wrong comparison to the intended one.
  Confirm the min/max extrema tests still pass (this code is the
  extrema/`y_min_max` path).

### 1b. `std::ofstream( string_view.data() )`
- **Where (7 sites):** [`include/Splines/Splines.hh:1077`](include/Splines/Splines.hh#L1077),
  [`SplineCubicBase.hxx:181,187,193`](include/Splines/SplineCubicBase.hxx#L181),
  [`SplineQuinticBase.hxx:122,128,134`](include/Splines/SplineQuinticBase.hxx#L122)
  ```cpp
  void dump( string_view fname, ... ) { std::ofstream file( fname.data() ); ... }
  ```
- **Problem:** `string_view::data()` is **not** guaranteed NUL-terminated. A
  true substring view opens a filename that runs past its intended end. Works
  today only because callers pass string literals / `std::string`.
- **Fix:** construct `std::string{ fname }` for the stream (or change the
  overload to take `std::string const&`). Keep the `string_view` public
  signature to avoid an API break; just materialize a terminated string
  internally.
- **Risk:** minimal.

---

## Tier 2 — Safety hardening

### 2a. Owning raw pointers → span-backed views
- **Where:** `m_X`, `m_Y` in [`Splines.hh:680-681`](include/Splines/Splines.hh#L680);
  `m_Yp` in [`SplineCubicBase.hxx:39`](include/Splines/SplineCubicBase.hxx#L39);
  the `m_external_alloc` flag + `reserve_external(...)` bookkeeping; `memcpy`
  copies in [`SplineCubicBase.cc:105`](src/SplineCubicBase.cc#L105).
- **Proposal:** keep the `Utils::Malloc` pool (it backs SplineSet's
  single-block sub-allocation), but stop passing bare `real_type*` around
  internally — expose the block through `std::span<real_type>` views so
  indexing is bounds-checkable in debug builds and the manual pointer juggling
  / `m_external_alloc` branching shrinks.
- **Risk:** moderate (touches every derived class's `reserve`); purely internal.

### 2b. Re-enable narrowing diagnostics
- **Context:** `integer = int` ([`Splines.hh:68`](include/Splines/Splines.hh#L68))
  is used for all indexing. The CMake rewrite dropped the old
  `-Wconversion -Wsign-conversion` flags (now only `-Wall -Wextra`).
- **Proposal:** add a **strict** CI build variant re-enabling
  `-Wconversion -Wsign-conversion` (GCC/Clang) to surface `int`/`size_t`
  boundary narrowing, then fix the hits. Decide per-site whether to keep `int`
  or move to `std::ptrdiff_t`.
- **Risk:** low; diagnostic-driven.

### 2c. `[[nodiscard]]` on pure observers
- **Context:** currently **0** occurrences in the library.
- **Proposal:** annotate `eval`, `D`, `DD`, `DDD`, `name`, `x_min/x_max`,
  `is_closed`, … so dropped results are caught at compile time.
- **Risk:** minimal; may reveal (and want to fix) existing misuse.

### 2d. C-API global mutable state
- **Where:** `spline_stored` map + `head`
  ([`SplinesCinterface.cc:132-133`](src/SplinesCinterface.cc#L132)).
- **Problem:** process-global, non-reentrant; concurrent C-API calls race.
- **Proposal:** at minimum document "not thread-safe"; better, guard the map
  with a mutex, or add an explicit context-handle variant of the API.
- **Risk:** low (documentation) to moderate (new API surface).

### 2e. `std::span` overloads on buffer APIs
- **Proposal:** offer `std::span<const real_type>` overloads alongside the
  existing `(ptr, n)` build/eval signatures (C++20 is already required).
- **Risk:** minimal; additive.

---

## Tier 3 — Performance (measurement-driven)

> There is **no** benchmark harness today (`test09` uses `TicToc` but is not
> one). **Step 0 of this tier is to add a micro-benchmark** so every change
> below is validated against a baseline on the VM and macOS.

### 3a. Batch evaluation API — DONE
- **Shipped:** `eval`/`D`/`DD( std::span<const real_type>, std::span<real_type> )`
  on the base `Spline` ([`Splines.hh`](include/Splines/Splines.hh)), one
  implementation, no per-leaf boilerplate. Bit-identical to the scalar path
  (asserted in [`benchmarks/bench_eval.cc`](benchmarks/bench_eval.cc) across
  exact knots, out-of-domain points, sorted input, and a closed/periodic
  spline).
- **Finding:** amortizing the *dispatch* alone was ~1.03x — measurement showed
  the outer virtual call is **not** the cost; the per-point interval `find()`
  is. The win comes from skipping redundant searches: for **sorted/monotone**
  input the batch reuses the current segment (strict half-open guard, falls
  back to `find()` at boundaries / out-of-domain / wrapping), giving **3-4x**.
  An up-front `is_sorted` gate keeps **random** input at ~1.0x (no regression,
  no footgun).
- **Left for later:** the remaining per-point cost is `find()` + the virtual
  `id_eval`. Devirtualizing `id_eval` (leaf overrides, leaning on the Tier-3c
  `final`) and/or a SIMD inner loop are the next levers if profiles justify.

### 3b. Hoist reciprocals / precompute inverse interval widths
- **Where:** `id_eval` and the Hermite basis recompute
  `1/(m_X[ni+1]-m_X[ni])`, `1/H`, `1/H²` on every call
  ([`Splines.hh:210-249`](include/Splines/Splines.hh#L210),
  [`SplineLinear.cc:84`](src/SplineLinear.cc#L84)).
- **Proposal:** precompute inverse interval widths at build time; pass the
  reciprocal into the Hermite helpers.
- **Risk:** low–moderate; must confirm bit-for-bit-close numerics (division
  vs. multiply-by-reciprocal changes rounding — verify test tolerances hold).

### 3c. Devirtualization
- **Proposal:** ensure `final` on leaf spline classes (24 already have it) so
  the compiler can devirtualize when the static type is known; measure whether
  it moves the needle for scalar loops (3a may make this moot).
- **Risk:** minimal.

---

## Tier 4 — Idiom cleanup (low value; fold into the tiers above)

- `typedef` → `using` ([`Splines.hh:67-68`](include/Splines/Splines.hh#L67)).
- `noexcept` on non-throwing observers (only 4 in the library today).
- Drop `using namespace std;` from a translation unit
  ([`SplinesCinterface.cc:52`](src/SplinesCinterface.cc#L52)).
- `constexpr` on the Hermite basis helpers where feasible.

---

## Suggested sequencing

1. **Tier 1** (bug fixes) — ship gated by the test suite; unambiguous.
2. **Tier 3 step 0** (benchmark harness) — establishes the baseline.
3. **Tier 2a/2c/2e** (span views, `[[nodiscard]]`, span overloads) — the bulk
   of the "safer" payoff.
4. **Tier 3a/3b** (batch eval, reciprocal hoisting) — the "faster" payoff,
   each accepted only if the benchmark shows a real gain and numerics hold.
5. **Tier 2b/2d + Tier 4** — cleanup and CI hardening.

## Verification (every tier)
- Build + `ctest` (19/19) on the Linux VM (`orb -m mads`, `lbuild/`) **and**
  macOS (`build/`).
- For Tier 3, the new micro-benchmark must show no regression (3b) or a
  measured improvement (3a) on both platforms.
- Numerical-equivalence check for any hot-path change (3a/3b): compare against
  the current scalar path within existing test tolerances.
