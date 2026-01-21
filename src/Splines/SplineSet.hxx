/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2016                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#pragma once

#ifndef SPLINESET_HXX
#define SPLINESET_HXX

#include "Splines.hh"
#include "Utils_fmt.hh"
#include "Utils_AlgoHNewton.hh"

#include <limits>
#include <cmath>
#include <set>
#include <algorithm>

namespace Splines
{

  using std::abs;
  using std::copy_n;
  using std::sqrt;

  /*\
   |   ____        _ _            ____       _
   |  / ___| _ __ | (_)_ __   ___/ ___|  ___| |_
   |  \___ \| '_ \| | | '_ \ / _ \___ \ / _ \ __|
   |   ___) | |_) | | | | | |  __/___) |  __/ |_
   |  |____/| .__/|_|_|_| |_|\___|____/ \___|\__|
   |        |_|
  \*/

  /**
   * \class SplineSet
   * \brief A collection of splines sharing the same independent variable (x-nodes)
   *
   * The SplineSet class manages multiple 1D splines that share the same
   * independent variable values (x-nodes). This is particularly useful for
   * representing multiple dependent variables that are all functions of the
   * same independent variable.
   *
   * \par Key Features:
   * - Manages multiple splines of different types (linear, cubic, quintic, etc.)
   * - All splines share the same x-nodes
   * - Efficient evaluation of all splines at once
   * - Support for using any spline as the independent variable
   * - Automatic monotonicity checking
   * - Integration with GenericContainer for I/O
   * - Automatic differentiation support
   *
   * \par Architecture:
   * The class stores splines as unique pointers in a vector and maintains
   * additional data structures for fast lookup by name and efficient evaluation.
   *
   * \par Memory Management:
   * Uses custom allocators for efficient memory handling and cache locality.
   * All spline data is stored in contiguous memory blocks.
   *
   * \par Use Cases:
   * - Multiple physical quantities as functions of time
   * - Curve families with shared parameterization
   * - Multi-dimensional function approximation
   * - Inverse function evaluation (using spline as independent variable)
   *
   * \note All splines in the set must have the same number of points and
   *       share the same x-values.
   *
   * \see Spline, ConstantSpline, LinearSpline, CubicSpline, AkimaSpline,
   *      BesselSpline, PchipSpline, HermiteSpline, QuinticSpline
   */
  class SplineSet
  {
    SplineSet( SplineSet const & )                   = delete;
    SplineSet const & operator=( SplineSet const & ) = delete;

    /**
     * \class BinarySearch
     * \brief Internal helper class for efficient name-to-index lookup
     *
     * This class implements a binary search data structure for mapping
     * spline names to their indices in the spline set. It maintains
     * a sorted vector of (name, index) pairs for O(log n) lookups.
     *
     * \par Design:
     * - Uses std::vector for storage with pre-allocated capacity
     * - Maintains sorted order by spline name for binary search
     * - Supports insertion and lookup operations
     *
     * \par Complexity:
     * - Lookup: O(log n)
     * - Insertion: O(n) due to shifting elements
     *
     * \note This is an internal class and not intended for direct use.
     */
    class BinarySearch
    {
    public:
      /// Data type stored in the search structure: (name, index) pair
      typedef std::pair<std::string, integer> DATA_TYPE;

    private:
      /// Sorted vector of (name, index) pairs
      mutable std::vector<DATA_TYPE> data;

    public:
      /**
       * \brief Default constructor
       *
       * Initializes an empty search structure with pre-allocated capacity.
       */
      BinarySearch()
      {
        data.clear();
        data.reserve( 256 );
      }

      /**
       * \brief Destructor
       *
       * Cleans up the internal data structure.
       */
      ~BinarySearch() { data.clear(); }

      /**
       * \brief Clear all entries from the search structure
       *
       * Removes all stored (name, index) pairs and resets capacity.
       */
      void clear()
      {
        data.clear();
        data.reserve( 256 );
      }

      /**
       * \brief Get the number of elements in the search structure
       *
       * \return Number of stored (name, index) pairs
       */
      integer n_elem() const { return integer( data.size() ); }

      /**
       * \brief Get element at specified position
       *
       * \param[in] i Index of element to retrieve
       * \return Constant reference to the (name, index) pair at position i
       *
       * \throw May throw std::out_of_range if i is out of bounds
       */
      DATA_TYPE const & get_elem( integer i ) const { return data[i]; }

      /**
       * \brief Search for a spline by name using binary search
       *
       * Finds the index of a spline given its name. The search is
       * case-sensitive and performs binary search on the sorted name list.
       *
       * \param[in] id Name of the spline to find
       * \return Index of the spline if found, -1 otherwise
       *
       * \par Algorithm:
       * 1. Perform binary search on sorted name list
       * 2. Compare strings lexicographically
       * 3. Return associated index if found, -1 otherwise
       *
       * \par Complexity:
       * Time: O(log n) where n is number of splines
       * Space: O(1)
       */
      integer search( string_view id ) const
      {
        size_t U{ data.size() };
        size_t L{ 0 };
        while ( U - L > 1 )
        {
          size_t const pos{ ( L + U ) >> 1 };  // se L=U+1 --> (L+U)/2 ==> L
          string_view  id_pos{ data[pos].first };
          if ( id_pos < id )
            L = pos;
          else
            U = pos;
        }
        if ( data[L].first == id ) return data[L].second;
        if ( U < data.size() )
          if ( data[U].first == id ) return data[U].second;
        return -1;  // non trovato
      }

      /**
       * \brief Insert a new (name, index) pair into the search structure
       *
       * Inserts a new entry while maintaining the sorted order by name.
       * If the name already exists, it will be updated with the new index.
       *
       * \param[in] id Name of the spline
       * \param[in] position Index associated with the name
       *
       * \par Algorithm:
       * 1. Append new element at the end
       * 2. Bubble up while previous element has greater name
       * 3. Overwrite existing entry if name already exists
       *
       * \par Complexity:
       * Time: O(n) worst-case due to element shifting
       * Space: O(1) additional
       */
      void insert( string_view const id, integer const position )
      {
        size_t pos{ data.size() };
        data.emplace_back( id, position );
        while ( pos > 0 )
        {
          size_t const pos1{ pos - 1 };
          data[pos].first  = data[pos1].first;
          data[pos].second = data[pos1].second;
          if ( data[pos1].first < id ) break;
          pos = pos1;
        }
        data[pos] = DATA_TYPE( id, position );
      }
    };

  protected:
    /// Name of the spline set for identification and error messages
    string const m_name;

    /// Custom memory allocator for real-type values
    Utils::Malloc<real_type> m_mem;

    /// Custom memory allocator for real-type pointers
    Utils::Malloc<real_type *> m_mem_p;

    /// Custom memory allocator for integer values
    Utils::Malloc<int> m_mem_int;

    /// Vector of unique pointers to individual splines
    vector<std::unique_ptr<Spline>> m_splines;

    /// Number of points in each spline (shared across all splines)
    integer m_npts{ 0 };

    /// Number of splines in the set
    integer m_nspl{ 0 };

    /// Array of x-nodes (independent variable values), size = m_npts
    real_type * m_X{ nullptr };

    /// Array of pointers to y-values for each spline, size = m_nspl × m_npts
    real_type ** m_Y{ nullptr };

    /// Array of pointers to first derivatives for each spline (if available)
    real_type ** m_Yp{ nullptr };

    /// Array of pointers to second derivatives for each spline (if available)
    real_type ** m_Ypp{ nullptr };

    /// Minimum y-value for each spline, size = m_nspl
    real_type * m_Ymin{ nullptr };

    /// Maximum y-value for each spline, size = m_nspl
    real_type * m_Ymax{ nullptr };

    /// Monotonicity information for each spline, size = m_nspl
    /// - -2: non-monotone data
    /// - -1: not monotone
    /// - 0: monotone (non-strict)
    /// - 1: strictly monotone
    int * m_is_monotone{ nullptr };

    /// Map from spline name to its index for O(1) lookup
    std::map<string, integer> m_header_to_position;

  private:
    /**
     * \brief Find x-value where monotone spline intersects given y-value
     *
     * For a monotone spline, finds the x-value such that spline(x) = zeta.
     * Uses a combination of interval search and Halley's method for root finding.
     *
     * \param[in] spl Index of the spline to use as independent variable
     * \param[in] zeta Target y-value to intersect
     * \param[out] x Computed x-value such that spline[spl](x) = zeta
     * \return Pointer to the spline used as independent variable
     *
     * \throw UTILS_ASSERT if spline index is invalid
     * \throw UTILS_ASSERT if spline is not monotone
     * \throw UTILS_ASSERT if zeta is outside spline range
     *
     * \par Algorithm:
     * 1. Verify spline is monotone
     * 2. Find interval containing zeta using binary search
     * 3. Apply Halley's method (AlgoHNewton) to solve spline(x) - zeta = 0
     *
     * \note Only works with monotone splines (m_is_monotone[spl] > 0)
     *
     * \see Utils::AlgoHNewton
     */
    Spline const * intersect( integer const spl, real_type const zeta, real_type & x ) const
    {
      string msg = fmt::format( "SplineSet[{}]::intersect(...):", m_name );
      UTILS_ASSERT( spl >= 0 && spl < m_nspl, "{}\nSpline n.{} is not in SplineSet", msg, spl );
      UTILS_ASSERT(
        m_is_monotone[spl] > 0,
        "{}\nSpline n.{} is not monotone and can't be used as independent",
        msg,
        spl );
      Spline * S = m_splines[spl].get();
      // cerco intervallo intersezione
      real_type const * Y{ m_Y[spl] };
      UTILS_ASSERT(
        zeta >= Y[0] && zeta <= Y[m_npts - 1],
        "{} evaluation at zeta = {} is out of range: [{},{}]\n",
        msg,
        zeta,
        Y[0],
        Y[m_npts - 1] );

      integer interval{ static_cast<integer>( std::lower_bound( Y, Y + m_npts, zeta ) - Y ) };
      if ( interval > 0 ) --interval;
      if ( Utils::is_zero( Y[interval] - Y[interval + 1] ) ) ++interval;  // degenerate interval for duplicated nodes
      if ( interval >= m_npts - 1 ) interval = m_npts - 2;

      Utils::AlgoHNewton<real_type> HN;
      class HNfun : public Utils::AlgoHNewton_base_fun<real_type>
      {
        Spline *  m_S;
        real_type m_zeta;

      public:
        explicit HNfun( Spline * S, real_type zeta ) : m_S( S ), m_zeta( zeta ) {}
        real_type eval( real_type x ) const { return m_S->eval( x ) - m_zeta; }
        real_type D( real_type x ) const { return m_S->D( x ); }
      };

      HNfun fun( S, zeta );

      // compute intersection
      x = HN.eval( m_X[interval], m_X[interval + 1], &fun );
      return S;
    }

  public:
    //!
    //! \name Constructors
    //!
    ///@{

    /**
     * \brief Construct an empty SplineSet
     *
     * Creates an uninitialized spline set with no splines. Use build()
     * or setup() methods to populate the spline set with data.
     *
     * \param[in] name Identifier for the spline set (default: "SplineSet")
     *
     * \par Memory:
     * Initializes custom allocators but allocates no memory until build() is called.
     */
    explicit SplineSet( string_view name = "SplineSet" )
      : m_name( name )
      , m_mem( fmt::format( "SplineSet[{}]::m_mem", name ) )
      , m_mem_p( fmt::format( "SplineSet[{}]::m_mem_p", name ) )
      , m_mem_int( fmt::format( "SplineSet[{}]::m_mem_int", name ) )
    {
    }

    /**
     * \brief Destructor
     *
     * Releases all allocated memory including spline objects and internal arrays.
     * Automatically handles cleanup of all resources.
     */
    virtual ~SplineSet()
    {
      m_mem.free();
      m_mem_p.free();
      m_mem_int.free();
    }
    ///@}

    //!
    //! \name Info
    //!
    ///@{

    /**
     * \brief Get the name of the spline set
     *
     * \return String view containing the spline set identifier
     */
    string_view name() const { return m_name; }

    /**
     * \brief Get the name of a specific spline in the set
     *
     * \param[in] i Index of the spline (0 ≤ i < m_nspl)
     * \return Name of the i-th spline
     *
     * \throw UTILS_ASSERT if i is out of bounds
     */
    string_view header( integer const i ) const { return m_splines[i]->name(); }

    /**
     * \brief Get all spline names as a vector of strings
     *
     * \param[out] names Vector to be filled with spline names
     *
     * \par Post-conditions:
     * - names.size() == m_nspl
     * - names[i] == header(i) for all i
     */
    void get_headers( std::vector<std::string> & names ) const
    {
      names.clear();
      names.reserve( m_nspl );
      for ( integer i = 0; i < m_nspl; ++i ) names.emplace_back( m_splines[i]->name() );
    }

    /**
     * \brief Get a formatted string listing all spline names
     *
     * \return String in format: "[ 'name1' 'name2' ... 'nameN' ]"
     *
     * \par Example:
     * \code
     * "[ 'velocity' 'acceleration' 'position' ]"
     * \endcode
     */
    string name_list() const
    {
      string tmp{ "[ " };
      for ( integer i = 0; i < m_nspl; ++i ) tmp += fmt::format( "'{}' ", m_splines[i]->name() );
      tmp += "]";
      return tmp;
    }

    /**
     * \brief Get monotonicity information for a spline
     *
     * \param[in] i Index of the spline
     * \return Monotonicity code:
     *         - +1 = strictly monotone
     *         - 0 = weakly monotone (constant segments allowed)
     *         - -1 = non-monotone
     *
     * \note For cubic splines, this is determined by checking if
     *       the derivative doesn't change sign.
     */
    int is_monotone( integer const i ) const { return m_is_monotone[i]; }

    /**
     * \brief Get number of points in each spline
     *
     * \return Number of x-nodes (same for all splines in the set)
     */
    integer num_points() const { return m_npts; }

    /**
     * \brief Get number of splines in the set
     *
     * \return Number of splines managed by this SplineSet
     */
    integer num_splines() const { return m_nspl; }

    /**
     * \brief Find index of spline by name
     *
     * \param[in] hdr Name of the spline to find
     * \return Index of the spline (0 ≤ index < m_nspl)
     *
     * \throw UTILS_ASSERT if spline name not found
     */
    integer get_position( string_view hdr ) const
    {
      integer const pos{ m_header_to_position.at( hdr.data() ) };
      UTILS_ASSERT(
        pos >= 0 && pos < m_nspl,
        "SplineSet[{}]::get_position(\"{}\") not found!\n"
        "available keys: {}\n",
        m_name,
        hdr,
        name_list() );
      return pos;
    }

    /**
     * \brief Get pointer to x-nodes array
     *
     * \return Constant pointer to array of x-values (size = m_npts)
     */
    real_type const * x_nodes() const { return m_X; }

    /**
     * \brief Get pointer to y-nodes array for a specific spline
     *
     * \param[in] i Index of the spline
     * \return Constant pointer to array of y-values for spline i (size = m_npts)
     *
     * \throw UTILS_ASSERT if i is out of bounds
     */
    real_type const * y_nodes( integer const i ) const
    {
      UTILS_ASSERT(
        i >= 0 && i < m_nspl,
        "SplineSet[{}]::y_nodes({}) argument out of range [0,{}]\n",
        m_name,
        i,
        m_nspl - 1 );
      return m_Y[i];
    }

    /**
     * \brief Get a specific x-node value
     *
     * \param[in] npt Index of the x-node (0 ≤ npt < m_npts)
     * \return X-value at the specified node
     */
    real_type x_node( integer const npt ) const { return m_X[npt]; }

    /**
     * \brief Get a specific y-node value
     *
     * \param[in] npt Index of the node point
     * \param[in] spl Index of the spline
     * \return Y-value of spline `spl` at node `npt`
     */
    real_type y_node( integer const npt, integer const spl ) const { return m_Y[spl][npt]; }

    /**
     * \brief Get minimum x-value in the domain
     *
     * \return Smallest x-node value (first element of m_X)
     */
    real_type x_min() const { return m_X[0]; }

    /**
     * \brief Get maximum x-value in the domain
     *
     * \return Largest x-node value (last element of m_X)
     */
    real_type x_max() const { return m_X[m_npts - 1]; }

    /**
     * \brief Get minimum y-value of a specific spline
     *
     * \param[in] spl Index of the spline
     * \return Minimum y-value attained by spline `spl` over its domain
     */
    real_type y_min( integer const spl ) const { return m_Ymin[spl]; }

    /**
     * \brief Get maximum y-value of a specific spline
     *
     * \param[in] spl Index of the spline
     * \return Maximum y-value attained by spline `spl` over its domain
     */
    real_type y_max( integer const spl ) const { return m_Ymax[spl]; }

    /**
     * \brief Get minimum y-value of a spline by name
     *
     * \param[in] spl Name of the spline
     * \return Minimum y-value of the specified spline
     *
     * \throw UTILS_ASSERT if spline name not found
     */
    real_type y_min( string_view spl ) const
    {
      integer idx{ this->get_position( spl ) };
      return m_Ymin[idx];
    }

    /**
     * \brief Get maximum y-value of a spline by name
     *
     * \param[in] spl Name of the spline
     * \return Maximum y-value of the specified spline
     *
     * \throw UTILS_ASSERT if spline name not found
     */
    real_type y_max( string_view spl ) const
    {
      integer idx{ this->get_position( spl ) };
      return m_Ymax[idx];
    }

    ///@}

    //! \name Access splines
    ///@{

    /**
     * \brief Get pointer to a spline by index
     *
     * \param[in] i Index of the spline (0 ≤ i < m_nspl)
     * \return Pointer to the i-th spline (non-const)
     *
     * \throw UTILS_ASSERT if i is out of bounds
     *
     * \note The returned pointer remains valid until the SplineSet is modified.
     */
    Spline * get_spline( integer const i ) const
    {
      UTILS_ASSERT(
        i >= 0 && i < m_nspl,
        "SplineSet[{}]::get_spline({}) argument out of range [0,{}]\n",
        m_name,
        i,
        m_nspl - 1 );
      return m_splines[i].get();
    }

    /**
     * \brief Get pointer to a spline by name
     *
     * \param[in] hdr Name of the spline
     * \return Pointer to the requested spline (non-const)
     *
     * \throw UTILS_ASSERT if spline name not found
     */
    Spline * get_spline( string_view hdr ) const
    {
      integer idx{ this->get_position( hdr ) };
      return m_splines[idx].get();
    }

    ///@}

    //! \name Evaluate
    ///@{

    /**
     * \brief Evaluate a spline by index at given x
     *
     * \param[in] x Parameter value (in the domain of x-nodes)
     * \param[in] spl Index of the spline to evaluate
     * \return Value of spline `spl` at parameter `x`
     *
     * \par Complexity:
     * - Time: O(log npts) for interval search plus spline evaluation
     * - Space: O(1)
     */
    real_type operator()( real_type const x, integer const spl ) const
    {
      Spline const * S{ this->get_spline( spl ) };
      return S->eval( x );
    }

    /**
     * \brief Evaluate a spline by index at given x
     *
     * \copydetails operator()()
     * This is an alias for operator()().
     */
    real_type eval( real_type const x, integer const spl ) const
    {
      Spline const * S{ this->get_spline( spl ) };
      return S->eval( x );
    }

    /**
     * \brief Evaluate first derivative of a spline by index
     *
     * \param[in] x Parameter value
     * \param[in] spl Index of the spline
     * \return First derivative of spline `spl` at `x`
     */
    real_type D( real_type const x, integer const spl ) const
    {
      Spline const * S{ this->get_spline( spl ) };
      return S->D( x );
    }

    /**
     * \brief Evaluate first derivative of a spline by index
     *
     * \copydetails D()
     * This is an alias for D().
     */
    real_type eval_D( real_type const x, integer const spl ) const
    {
      Spline const * S{ this->get_spline( spl ) };
      return S->D( x );
    }

    /**
     * \brief Evaluate second derivative of a spline by index
     *
     * \param[in] x Parameter value
     * \param[in] spl Index of the spline
     * \return Second derivative of spline `spl` at `x`
     */
    real_type DD( real_type const x, integer const spl ) const
    {
      Spline const * S{ this->get_spline( spl ) };
      return S->DD( x );
    }

    /**
     * \brief Evaluate second derivative of a spline by index
     *
     * \copydetails DD()
     * This is an alias for DD().
     */
    real_type eval_DD( real_type const x, integer const spl ) const
    {
      Spline const * S{ this->get_spline( spl ) };
      return S->DD( x );
    }

    /**
     * \brief Evaluate third derivative of a spline by index
     *
     * \param[in] x Parameter value
     * \param[in] spl Index of the spline
     * \return Third derivative of spline `spl` at `x`
     */
    real_type DDD( real_type const x, integer const spl ) const
    {
      Spline const * S{ this->get_spline( spl ) };
      return S->DDD( x );
    }

    /**
     * \brief Evaluate third derivative of a spline by index
     *
     * \copydetails DDD()
     * This is an alias for DDD().
     */
    real_type eval_DDD( real_type const x, integer const spl ) const
    {
      Spline const * S{ this->get_spline( spl ) };
      return S->DDD( x );
    }

    /**
     * \brief Evaluate fourth derivative of a spline by index
     *
     * \param[in] x Parameter value
     * \param[in] spl Index of the spline
     * \return Fourth derivative of spline `spl` at `x`
     *
     * \note For cubic splines, returns zero beyond third derivative.
     */
    real_type DDDD( real_type const x, integer const spl ) const
    {
      Spline const * S{ this->get_spline( spl ) };
      return S->DDDD( x );
    }

    /**
     * \brief Evaluate fourth derivative of a spline by index
     *
     * \copydetails DDDD()
     * This is an alias for DDDD().
     */
    real_type eval_DDDD( real_type const x, integer const spl ) const
    {
      Spline const * S{ this->get_spline( spl ) };
      return S->DDDD( x );
    }

    /**
     * \brief Evaluate fifth derivative of a spline by index
     *
     * \param[in] x Parameter value
     * \param[in] spl Index of the spline
     * \return Fifth derivative of spline `spl` at `x`
     *
     * \note For cubic splines, returns zero beyond third derivative.
     */
    real_type DDDDD( real_type const x, integer const spl ) const
    {
      Spline const * S{ this->get_spline( spl ) };
      return S->DDDDD( x );
    }

    /**
     * \brief Evaluate fifth derivative of a spline by index
     *
     * \copydetails DDDDD()
     * This is an alias for DDDDD().
     */
    real_type eval_DDDDD( real_type const x, integer const spl ) const
    {
      Spline const * S{ this->get_spline( spl ) };
      return S->DDDDD( x );
    }

#ifdef AUTODIFF_SUPPORT
    //!
    //! \name Autodiff
    //!
    ///@{
    /**
     * \brief Evaluate spline with first-order automatic differentiation
     *
     * \param[in] x Dual number parameter (contains value and gradient)
     * \param[in] spl Index of the spline
     * \return Dual number containing function value and derivative
     *
     * \par Usage:
     * \code{.cpp}
     * autodiff::dual1st t = 1.5;
     * autodiff::dual1st y = spline_set.eval(t, 0);
     * // y.val = value at t=1.5
     * // y.grad = derivative dY/dt at t=1.5
     * \endcode
     */
    autodiff::dual1st eval( autodiff::dual1st const & x, integer const spl ) const
    {
      Spline const * S{ this->get_spline( spl ) };
      return S->eval( x );
    }

    /**
     * \brief Evaluate spline with second-order automatic differentiation
     *
     * \param[in] x Dual number parameter (second order)
     * \param[in] spl Index of the spline
     * \return Dual number containing function value, first and second derivatives
     */
    autodiff::dual2nd eval( autodiff::dual2nd const & x, integer const spl ) const
    {
      Spline const * S{ this->get_spline( spl ) };
      return S->eval( x );
    }

    /**
     * \brief Generic template for automatic differentiation
     *
     * Automatically selects the appropriate evaluation method based on input type.
     *
     * \tparam T Input type (arithmetic or autodiff dual type)
     * \param[in] x Parameter value
     * \param[in] spl Index of the spline
     * \return Result matching input type
     *
     * \par Type Detection:
     * - Arithmetic types (int, float, double): promotes to real_type
     * - Autodiff dual types: uses appropriate autodiff evaluation
     */
    template <typename T> auto eval( T const & x, integer const spl ) const
    {
      if constexpr ( std::is_arithmetic<T>::value )
      {
        // Se T è un tipo numerico (int, float, double, etc.), promuovi a real_type
        return eval( static_cast<real_type>( x ), spl );
      }
      else
      {
        // Altrimenti deduce automaticamente il tipo duale appropriato
        return eval( autodiff::detail::to_dual( x ), spl );
      }
    }

    /**
     * \brief Generic function call operator for autodiff
     *
     * \tparam T Parameter type
     * \param[in] x Parameter value
     * \param[in] spl Index of the spline
     * \return Result matching input type
     */
    template <typename T> auto operator()( T const & x, integer const spl ) const -> decltype( eval( x, spl ) )
    {
      return eval( x, spl );
    }
///@}
#endif

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

    ///@}

#ifdef AUTODIFF_SUPPORT
    //!
    //! \name Autodiff
    //!
    ///@{
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
///@}
#endif

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    //!
    //! \name Evaluate to std vector
    //!
    ///@{

    /**
     * \brief Evaluate all splines at given x and store in std::vector
     *
     * \param[in] x Parameter value
     * \param[out] vals Vector to store results (automatically resized to m_nspl)
     *
     * \par Post-conditions:
     * - vals.size() == m_nspl
     * - vals[i] == value of i-th spline at x
     */
    void eval( real_type const x, vector<real_type> & vals ) const
    {
      vals.resize( m_nspl );
      for ( integer i = 0; i < m_nspl; ++i ) vals[i] = m_splines[i]->eval( x );
    }

    /**
     * \brief Evaluate first derivatives of all splines and store in std::vector
     *
     * \param[in] x Parameter value
     * \param[out] vals Vector to store derivatives (automatically resized)
     */
    void eval_D( real_type const x, vector<real_type> & vals ) const
    {
      vals.resize( m_nspl );
      for ( integer i = 0; i < m_nspl; ++i ) vals[i] = m_splines[i]->D( x );
    }

    /**
     * \brief Evaluate second derivatives of all splines and store in std::vector
     *
     * \param[in] x Parameter value
     * \param[out] vals Vector to store second derivatives (automatically resized)
     */
    void eval_DD( real_type const x, vector<real_type> & vals ) const
    {
      vals.resize( m_nspl );
      for ( integer i = 0; i < m_nspl; ++i ) vals[i] = m_splines[i]->DD( x );
    }

    /**
     * \brief Evaluate third derivatives of all splines and store in std::vector
     *
     * \param[in] x Parameter value
     * \param[out] vals Vector to store third derivatives (automatically resized)
     */
    void eval_DDD( real_type const x, vector<real_type> & vals ) const
    {
      vals.resize( m_nspl );
      for ( integer i = 0; i < m_nspl; ++i ) vals[i] = m_splines[i]->DDD( x );
    }

    /**
     * \brief Evaluate fourth derivatives of all splines and store in std::vector
     *
     * \param[in] x Parameter value
     * \param[out] vals Vector to store fourth derivatives (automatically resized)
     */
    void eval_DDDD( real_type const x, vector<real_type> & vals ) const
    {
      vals.resize( m_nspl );
      for ( integer i = 0; i < m_nspl; ++i ) vals[i] = m_splines[i]->DDDD( x );
    }

    /**
     * \brief Evaluate fifth derivatives of all splines and store in std::vector
     *
     * \param[in] x Parameter value
     * \param[out] vals Vector to store fifth derivatives (automatically resized)
     */
    void eval_DDDDD( real_type const x, vector<real_type> & vals ) const
    {
      vals.resize( m_nspl );
      for ( integer i = 0; i < m_nspl; ++i ) vals[i] = m_splines[i]->DDDDD( x );
    }

    ///@}

    //! \name Evaluate to vector
    ///@{

    /**
     * \brief Evaluate all splines at given x and store in C array with stride
     *
     * \param[in] x Parameter value
     * \param[out] vals Output array (must have space for at least m_nspl*incy elements)
     * \param[in] incy Stride between consecutive values in output array (default: 1)
     *
     * \par Memory Layout:
     * - incy = 1: vals[0..m_nspl-1] = spline values
     * - incy = k: vals[0, k, 2k, ...] = spline values
     */
    void eval( real_type const x, real_type vals[], integer const incy = 1 ) const
    {
      integer ii{ 0 };
      for ( integer i = 0; i < m_nspl; ++i, ii += incy ) vals[ii] = m_splines[i]->eval( x );
    }

    /**
     * \brief Evaluate first derivatives of all splines to C array with stride
     *
     * \param[in] x Parameter value
     * \param[out] vals Output array for derivatives
     * \param[in] incy Stride between consecutive values (default: 1)
     */
    void eval_D( real_type const x, real_type vals[], integer const incy = 1 ) const
    {
      size_t ii{ 0 };
      for ( integer i = 0; i < m_nspl; ++i, ii += incy ) vals[ii] = m_splines[i]->D( x );
    }

    /**
     * \brief Evaluate second derivatives of all splines to C array with stride
     *
     * \param[in] x Parameter value
     * \param[out] vals Output array for second derivatives
     * \param[in] incy Stride between consecutive values (default: 1)
     */
    void eval_DD( real_type const x, real_type vals[], integer const incy = 1 ) const
    {
      size_t ii{ 0 };
      for ( integer i = 0; i < m_nspl; ++i, ii += incy ) vals[ii] = m_splines[i]->DD( x );
    }

    /**
     * \brief Evaluate third derivatives of all splines to C array with stride
     *
     * \param[in] x Parameter value
     * \param[out] vals Output array for third derivatives
     * \param[in] incy Stride between consecutive values (default: 1)
     */
    void eval_DDD( real_type const x, real_type vals[], integer const incy = 1 ) const
    {
      size_t ii{ 0 };
      for ( integer i = 0; i < m_nspl; ++i, ii += incy ) vals[ii] = m_splines[i]->DDD( x );
    }

    /**
     * \brief Evaluate fourth derivatives of all splines to C array with stride
     *
     * \param[in] x Parameter value
     * \param[out] vals Output array for fourth derivatives
     * \param[in] incy Stride between consecutive values (default: 1)
     */
    void eval_DDDD( real_type const x, real_type vals[], integer const incy = 1 ) const
    {
      size_t ii{ 0 };
      for ( integer i = 0; i < m_nspl; ++i, ii += incy ) vals[ii] = m_splines[i]->DDDD( x );
    }

    /**
     * \brief Evaluate fifth derivatives of all splines to C array with stride
     *
     * \param[in] x Parameter value
     * \param[out] vals Output array for fifth derivatives
     * \param[in] incy Stride between consecutive values (default: 1)
     */
    void eval_DDDDD( real_type const x, real_type vals[], integer const incy = 1 ) const
    {
      size_t ii{ 0 };
      for ( integer i = 0; i < m_nspl; ++i, ii += incy ) vals[ii] = m_splines[i]->DDDDD( x );
    }

    ///@}

    //! \name Evaluate using another spline as independent
    ///@{

    /**
     * \brief Evaluate all splines using a specific spline as independent variable
     *
     * Instead of evaluating at a given x, evaluates at the x where spline `spl`
     * equals `zeta`. This is useful for inverse evaluation or cross-section analysis.
     *
     * \param[in] spl Index of spline to use as independent variable (must be monotone)
     * \param[in] zeta Target value of the independent spline
     * \param[out] x Computed x-value such that spline[spl](x) = zeta
     * \param[out] vals Vector to store all spline values at computed x
     *
     * \throw UTILS_ASSERT if spline `spl` is not monotone
     * \throw UTILS_ASSERT if zeta is outside range of spline `spl`
     */
    void eval2( integer const spl, real_type const zeta, real_type & x, vector<real_type> & vals ) const
    {
      vals.resize( m_nspl );
      this->eval2( spl, zeta, x, vals.data(), 1 );
    }

    /**
     * \brief Evaluate first derivatives using a specific spline as independent variable
     *
     * \param[in] spl Index of independent spline
     * \param[in] zeta Target value of independent spline
     * \param[out] x Computed x-value
     * \param[out] vals Vector to store first derivatives at computed x
     */
    void eval2_D( integer const spl, real_type const zeta, real_type & x, vector<real_type> & vals ) const
    {
      vals.resize( m_nspl );
      this->eval2_D( spl, zeta, x, vals.data(), 1 );
    }

    /**
     * \brief Evaluate second derivatives using a specific spline as independent variable
     *
     * \param[in] spl Index of independent spline
     * \param[in] zeta Target value of independent spline
     * \param[out] x Computed x-value
     * \param[out] vals Vector to store second derivatives at computed x
     */
    void eval2_DD( integer const spl, real_type const zeta, real_type & x, vector<real_type> & vals ) const
    {
      vals.resize( m_nspl );
      this->eval2_DD( spl, zeta, x, vals.data(), 1 );
    }

    /**
     * \brief Evaluate third derivatives using a specific spline as independent variable
     *
     * \param[in] spl Index of independent spline
     * \param[in] zeta Target value of independent spline
     * \param[out] x Computed x-value
     * \param[out] vals Vector to store third derivatives at computed x
     */
    void eval2_DDD( integer const spl, real_type const zeta, real_type & x, vector<real_type> & vals ) const
    {
      vals.resize( m_nspl );
      this->eval2_DDD( spl, zeta, x, vals.data(), 1 );
    }

    ///@}

    //! \name Evaluate using another spline as independent
    ///@{

    /**
     * \brief Evaluate all splines using a specific spline as independent variable (C array version)
     *
     * \param[in] spl Index of independent spline
     * \param[in] zeta Target value of independent spline
     * \param[out] x Computed x-value such that spline[spl](x) = zeta
     * \param[out] vals Output array to store spline values
     * \param[in] incy Stride between consecutive values (default: 1)
     */
    void eval2( integer const spl, real_type const zeta, real_type & x, real_type vals[], integer const incy = 1 ) const
    {
      intersect( spl, zeta, x );
      size_t ii{ 0 };
      for ( integer i = 0; i < m_nspl; ++i, ii += incy ) vals[ii] = m_splines[i]->eval( x );
    }

    /**
     * \brief Evaluate first derivatives using a specific spline as independent variable (C array version)
     *
     * Computes derivatives with respect to the independent spline value, not x.
     *
     * \param[in] spl Index of independent spline
     * \param[in] zeta Target value of independent spline
     * \param[out] x Computed x-value
     * \param[out] vals Output array for derivatives
     * \param[in] incy Stride between consecutive values (default: 1)
     *
     * \par Mathematical Formula:
     * \f[
     * \frac{dy_i}{d\zeta} = \frac{dy_i/dx}{ds/dx}
     * \f]
     * where s is the independent spline.
     */
    void eval2_D( integer const spl, real_type const zeta, real_type & x, real_type vals[], integer const incy = 1 )
      const
    {
      Spline const *  S{ intersect( spl, zeta, x ) };
      real_type const ds{ S->D( x ) };
      size_t          ii{ 0 };
      for ( integer i = 0; i < m_nspl; ++i, ii += incy ) vals[ii] = m_splines[i]->D( x ) / ds;
    }

    /**
     * \brief Evaluate second derivatives using a specific spline as independent variable (C array version)
     *
     * \param[in] spl Index of independent spline
     * \param[in] zeta Target value of independent spline
     * \param[out] x Computed x-value
     * \param[out] vals Output array for second derivatives
     * \param[in] incy Stride between consecutive values (default: 1)
     *
     * \par Mathematical Formula:
     * \f[
     * \frac{d^2y_i}{d\zeta^2} = \frac{d^2y_i/dx^2}{(ds/dx)^2} - \frac{dy_i/dx \cdot d^2s/dx^2}{(ds/dx)^3}
     * \f]
     */
    void eval2_DD( integer const spl, real_type const zeta, real_type & x, real_type vals[], integer const incy = 1 )
      const
    {
      Spline const *  S{ intersect( spl, zeta, x ) };
      real_type const dt{ 1 / S->D( x ) };
      real_type const dt2{ dt * dt };
      real_type const ddt{ -S->DD( x ) * ( dt * dt2 ) };
      size_t          ii{ 0 };
      for ( integer i = 0; i < m_nspl; ++i, ii += incy )
      {
        auto & Si{ m_splines[i] };
        vals[ii] = Si->DD( x ) * dt2 + Si->D( x ) * ddt;
      }
    }

    /**
     * \brief Evaluate third derivatives using a specific spline as independent variable (C array version)
     *
     * \param[in] spl Index of independent spline
     * \param[in] zeta Target value of independent spline
     * \param[out] x Computed x-value
     * \param[out] vals Output array for third derivatives
     * \param[in] incy Stride between consecutive values (default: 1)
     *
     * \par Mathematical Formula:
     * Uses chain rule for third derivatives with respect to parameterized variable.
     */
    void eval2_DDD( integer const spl, real_type const zeta, real_type & x, real_type vals[], integer const incy = 1 )
      const
    {
      Spline const *  S{ intersect( spl, zeta, x ) };
      real_type const dt{ 1 / S->D( x ) };
      real_type const dt3{ dt * dt * dt };
      real_type const ddt{ -S->DD( x ) * dt3 };
      real_type const dddt{ 3 * ( ddt * ddt ) / dt - S->DDD( x ) * ( dt * dt3 ) };
      size_t          ii{ 0 };
      for ( integer i = 0; i < m_nspl; ++i, ii += incy )
      {
        auto & Si{ m_splines[i] };
        vals[ii] = Si->DDD( x ) * dt3 + 3 * Si->DD( x ) * dt * ddt + Si->D( x ) * dddt;
      }
    }

    ///@}

    // \name Evaluate using another spline as independent
    ///@{

    /**
     * \brief Evaluate a specific spline using another spline as independent variable
     *
     * \param[in] zeta Target value of independent spline
     * \param[in] indep Name of independent spline
     * \param[out] x Computed x-value
     * \param[in] name Name of spline to evaluate
     * \return Value of spline `name` at computed x
     */
    real_type eval2( real_type const zeta, string_view const indep, real_type & x, string_view const name ) const
    {
      return this->eval2( zeta, this->get_position( indep ), x, this->get_position( name ) );
    }

    /**
     * \brief Evaluate first derivative using another spline as independent variable
     *
     * \param[in] zeta Target value of independent spline
     * \param[in] indep Name of independent spline
     * \param[out] x Computed x-value
     * \param[in] name Name of spline to evaluate
     * \return First derivative of spline `name` with respect to independent spline value
     */
    real_type eval2_D( real_type const zeta, string_view const indep, real_type & x, string_view const name ) const
    {
      return this->eval2_D( zeta, this->get_position( indep ), x, this->get_position( name ) );
    }

    /**
     * \brief Evaluate second derivative using another spline as independent variable
     *
     * \param[in] zeta Target value of independent spline
     * \param[in] indep Name of independent spline
     * \param[out] x Computed x-value
     * \param[in] name Name of spline to evaluate
     * \return Second derivative of spline `name` with respect to independent spline value
     */
    real_type eval2_DD( real_type const zeta, string_view const indep, real_type & x, string_view const name ) const
    {
      return this->eval2_DD( zeta, this->get_position( indep ), x, this->get_position( name ) );
    }

    /**
     * \brief Evaluate third derivative using another spline as independent variable
     *
     * \param[in] zeta Target value of independent spline
     * \param[in] indep Name of independent spline
     * \param[out] x Computed x-value
     * \param[in] name Name of spline to evaluate
     * \return Third derivative of spline `name` with respect to independent spline value
     */
    real_type eval2_DDD( real_type const zeta, string_view const indep, real_type & x, string_view const name ) const
    {
      return this->eval2_DDD( zeta, this->get_position( indep ), x, this->get_position( name ) );
    }

    ///@}

#ifdef AUTODIFF_SUPPORT
    //!
    //! \name Autodiff
    //!
    ///@{
    /**
     * \brief Evaluate with first-order automatic differentiation using independent spline
     *
     * \param[in] zeta Dual number target value of independent spline
     * \param[in] indep Name of independent spline
     * \param[out] x Computed x-value (real type)
     * \param[in] name Name of spline to evaluate
     * \return Dual number containing value and derivative
     */
    autodiff::dual1st eval2(
      autodiff::dual1st const & zeta,
      string_view const         indep,
      real_type &               x,
      string_view const         name ) const
    {
      using autodiff::derivative;
      using autodiff::dual1st;
      real_type zv{ val( zeta ) };
      dual1st   res{ eval2( zv, indep, x, name ) };
      res.grad = eval_D( x, name ) * zeta.grad;  // x gia calcolato
      return res;
    }

    /**
     * \brief Evaluate with second-order automatic differentiation using independent spline
     *
     * \param[in] zeta Dual number target value of independent spline
     * \param[in] indep Name of independent spline
     * \param[out] x Computed x-value (real type)
     * \param[in] name Name of spline to evaluate
     * \return Dual number containing value, first and second derivatives
     */
    autodiff::dual2nd eval2(
      autodiff::dual2nd const & zeta,
      string_view const         indep,
      real_type &               x,
      string_view const         name ) const
    {
      using autodiff::derivative;
      using autodiff::dual2nd;

      real_type zv{ val( zeta ) };
      real_type zg{ val( zeta.grad ) };
      real_type dfx{ eval2_D( zv, indep, x, name ) };
      real_type dxx{ eval_DD( x, name ) };  // x gia calcolato
      dual2nd   res{ eval( x, name ) };

      res.grad      = dfx * zg;
      res.grad.grad = dfx * zeta.grad.grad + dxx * ( zg * zg );
      return res;
    }

    /**
     * \brief Generic template for automatic differentiation using independent spline
     *
     * \tparam T Input type
     * \param[in] zeta Target value of independent spline
     * \param[in] indep Name of independent spline
     * \param[out] x Computed x-value (real type)
     * \param[in] name Name of spline to evaluate
     * \return Result matching input type
     */
    template <typename T>
    auto eval2( T const & zeta, string_view const indep, real_type & x, string_view const name ) const
    {
      if constexpr ( std::is_arithmetic<T>::value )
      {
        // Se T è un tipo numerico (int, float, double, etc.), promuovi a real_type
        return eval2( static_cast<real_type>( zeta ), indep, x, name );
      }
      else
      {
        // Altrimenti deduce automaticamente il tipo duale appropriato
        return eval2( autodiff::detail::to_dual( zeta ), indep, x, name );
      }
    }
///@}
#endif

    //! \name Evaluate using another spline as independent
    ///@{

    /**
     * \brief Evaluate a specific spline using another spline as independent variable (index version)
     *
     * \param[in] zeta Target value of independent spline
     * \param[in] indep Index of independent spline
     * \param[out] x Computed x-value
     * \param[in] spl Index of spline to evaluate
     * \return Value of spline `spl` at computed x
     */
    real_type eval2( real_type const zeta, integer const indep, real_type & x, integer const spl ) const
    {
      intersect( indep, zeta, x );
      return m_splines[spl]->eval( x );
    }

    /**
     * \brief Evaluate first derivative using another spline as independent variable (index version)
     *
     * \param[in] zeta Target value of independent spline
     * \param[in] indep Index of independent spline
     * \param[out] x Computed x-value
     * \param[in] spl Index of spline to evaluate
     * \return First derivative with respect to independent spline value
     */
    real_type eval2_D( real_type const zeta, integer const indep, real_type & x, integer const spl ) const
    {
      Spline const * S{ intersect( indep, zeta, x ) };
      return m_splines[spl]->D( x ) / S->D( x );
    }

    /**
     * \brief Evaluate second derivative using another spline as independent variable (index version)
     *
     * \param[in] zeta Target value of independent spline
     * \param[in] indep Index of independent spline
     * \param[out] x Computed x-value
     * \param[in] spl Index of spline to evaluate
     * \return Second derivative with respect to independent spline value
     */
    real_type eval2_DD( real_type const zeta, integer const indep, real_type & x, integer const spl ) const
    {
      Spline const *  S{ intersect( indep, zeta, x ) };
      real_type const dt{ 1 / S->D( x ) };
      real_type const dt2{ dt * dt };
      real_type const ddt{ -S->DD( x ) * ( dt * dt2 ) };
      auto &          SPL{ m_splines[spl] };
      return SPL->DD( x ) * dt2 + SPL->D( x ) * ddt;
    }

    /**
     * \brief Evaluate third derivative using another spline as independent variable (index version)
     *
     * \param[in] zeta Target value of independent spline
     * \param[in] indep Index of independent spline
     * \param[out] x Computed x-value
     * \param[in] spl Index of spline to evaluate
     * \return Third derivative with respect to independent spline value
     */
    real_type eval2_DDD( real_type const zeta, integer const indep, real_type & x, integer const spl ) const
    {
      Spline const *  S{ intersect( indep, zeta, x ) };
      real_type const dt{ 1 / S->D( x ) };
      real_type const dt3{ dt * dt * dt };
      real_type const ddt{ -S->DD( x ) * dt3 };
      real_type const dddt{ 3 * ( ddt * ddt ) / dt - S->DDD( x ) * ( dt * dt3 ) };

      auto & SPL{ m_splines[spl] };
      return SPL->DDD( x ) * dt3 + 3 * SPL->DD( x ) * dt * ddt + SPL->D( x ) * dddt;
    }

    ///@}

#ifdef AUTODIFF_SUPPORT
    //!
    //! \name Autodiff
    //!
    ///@{
    /**
     * \brief Evaluate with first-order automatic differentiation (index version)
     *
     * \param[in] zeta Dual number target value
     * \param[in] indep Index of independent spline
     * \param[out] x Computed x-value (real type)
     * \param[in] spl Index of spline to evaluate
     * \return Dual number containing value and derivative
     */
    autodiff::dual1st eval2( autodiff::dual1st const & zeta, integer const indep, real_type & x, integer const spl )
      const
    {
      using autodiff::derivative;
      using autodiff::dual1st;
      real_type zv{ val( zeta ) };
      dual1st   res{ eval2( zv, indep, x, spl ) };
      res.grad = eval_D( x, spl ) * zeta.grad;  // x gia calcolato
      return res;
    }

    /**
     * \brief Evaluate with second-order automatic differentiation (index version)
     *
     * \param[in] zeta Dual number target value
     * \param[in] indep Index of independent spline
     * \param[out] x Computed x-value (real type)
     * \param[in] spl Index of spline to evaluate
     * \return Dual number containing value, first and second derivatives
     */
    autodiff::dual2nd eval2( autodiff::dual2nd const & zeta, integer const indep, real_type & x, integer const spl )
      const
    {
      using autodiff::derivative;
      using autodiff::dual2nd;

      real_type zv{ val( zeta ) };
      real_type zg{ val( zeta.grad ) };
      real_type dfx{ eval2_D( zv, indep, x, spl ) };
      real_type dxx{ eval_DD( x, spl ) };
      dual2nd   res{ eval( x, spl ) };

      res.grad      = dfx * zg;
      res.grad.grad = dfx * zeta.grad.grad + dxx * ( zg * zg );
      return res;
    }

    /**
     * \brief Generic template for automatic differentiation (index version)
     *
     * \tparam T Input type
     * \param[in] zeta Target value of independent spline
     * \param[in] indep Index of independent spline
     * \param[out] x Computed x-value (real type)
     * \param[in] spl Index of spline to evaluate
     * \return Result matching input type
     */
    template <typename T> auto eval2( T const & zeta, integer const indep, real_type & x, integer const spl ) const
    {
      if constexpr ( std::is_arithmetic<T>::value )
      {
        // Se T è un tipo numerico (int, float, double, etc.), promuovi a real_type
        return eval2( static_cast<real_type>( zeta ), indep, x, spl );
      }
      else
      {
        // Altrimenti deduce automaticamente il tipo duale appropriato
        return eval2( autodiff::detail::to_dual( zeta ), indep, x, spl );
      }
    }
///@}
#endif

    //! \name Evaluate onto a vector
    ///@{

    /**
     * \brief Evaluate all splines at x and store results in GenericContainer
     *
     * Results are stored as a map where keys are spline names and values
     * are the evaluated spline values.
     *
     * \param[in] x Parameter value
     * \param[out] vals GenericContainer to store results
     */
    void eval( real_type const x, GenericContainer & vals ) const;

    /**
     * \brief Evaluate all splines at multiple x-values and store in GenericContainer
     *
     * For each spline, creates a vector of values corresponding to each x in `vec`.
     *
     * \param[in] vec Vector of x-values
     * \param[out] vals GenericContainer with spline names as keys and vectors as values
     */
    void eval( vec_real_type const & vec, GenericContainer & vals ) const;

    /**
     * \brief Evaluate specific splines at x and store in GenericContainer
     *
     * \param[in] x Parameter value
     * \param[in] columns Names of splines to evaluate
     * \param[out] vals GenericContainer with requested spline values
     */
    void eval( real_type const x, vec_string_type const & columns, GenericContainer & vals ) const;

    /**
     * \brief Evaluate specific splines at multiple x-values and store in GenericContainer
     *
     * \param[in] vec Vector of x-values
     * \param[in] columns Names of splines to evaluate
     * \param[out] vals GenericContainer with vectors of values for requested splines
     */
    void eval( vec_real_type const & vec, vec_string_type const & columns, GenericContainer & vals ) const;

    ///@}

    //! \name Evaluate to GenericContainer using another spline as independent
    ///@{

    /**
     * \brief Evaluate all splines using independent spline and store in GenericContainer
     *
     * \param[in] zeta Target value of independent spline
     * \param[in] indep Index of independent spline
     * \param[out] x Computed x-value
     * \param[out] vals GenericContainer to store results
     */
    void eval2( real_type const zeta, integer const indep, real_type & x, GenericContainer & vals ) const;

    /**
     * \brief Evaluate all splines at multiple zeta values using independent spline
     *
     * \param[in] zetas Vector of target values for independent spline
     * \param[in] indep Index of independent spline
     * \param[out] vals GenericContainer with results
     */
    void eval2( vec_real_type const & zetas, integer const indep, GenericContainer & vals ) const;

    /**
     * \brief Evaluate specific splines using independent spline and store in GenericContainer
     *
     * \param[in] zeta Target value of independent spline
     * \param[in] indep Index of independent spline
     * \param[out] x Computed x-value
     * \param[in] columns Names of splines to evaluate
     * \param[out] vals GenericContainer to store results
     */
    void eval2(
      real_type const         zeta,
      integer const           indep,
      real_type &             x,
      vec_string_type const & columns,
      GenericContainer &      vals ) const;

    /**
     * \brief Evaluate specific splines at multiple zeta values using independent spline
     *
     * \param[in] zetas Vector of target values for independent spline
     * \param[in] indep Index of independent spline
     * \param[in] columns Names of splines to evaluate
     * \param[out] vals GenericContainer with results
     */
    void eval2(
      vec_real_type const &   zetas,
      integer const           indep,
      vec_string_type const & columns,
      GenericContainer &      vals ) const;

    ///@}

    //! \name Evaluate to GenericContainer using another spline as independent
    ///@{

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
     * \brief Evaluate specific splines at multiple zeta values using independent spline by name
     *
     * \param[in] zetas Vector of target values for independent spline
     * \param[in] indep Name of independent spline
     * \param[in] columns Names of splines to evaluate
     * \param[out] vals GenericContainer with results
     */
    void eval2(
      vec_real_type const &   zetas,
      string_view             indep,
      vec_string_type const & columns,
      GenericContainer &      vals ) const
    {
      this->eval2( zetas, this->get_position( indep ), columns, vals );
    }
    ///@}

    //! \name Evaluate derivative into a Generic Container
    ///@{

    /**
     * \brief Evaluate first derivatives of all splines and store in GenericContainer
     *
     * \param[in] x Parameter value
     * \param[out] vals GenericContainer with derivative values
     */
    void eval_D( real_type const x, GenericContainer & vals ) const;

    /**
     * \brief Evaluate first derivatives of all splines at multiple x-values
     *
     * \param[in] vec Vector of x-values
     * \param[out] vals GenericContainer with vectors of derivatives
     */
    void eval_D( vec_real_type const & vec, GenericContainer & vals ) const;

    /**
     * \brief Evaluate first derivatives of specific splines
     *
     * \param[in] x Parameter value
     * \param[in] columns Names of splines to evaluate
     * \param[out] vals GenericContainer with derivative values
     */
    void eval_D( real_type const x, vec_string_type const & columns, GenericContainer & vals ) const;

    /**
     * \brief Evaluate first derivatives of specific splines at multiple x-values
     *
     * \param[in] vec Vector of x-values
     * \param[in] columns Names of splines to evaluate
     * \param[out] vals GenericContainer with vectors of derivatives
     */
    void eval_D( vec_real_type const & vec, vec_string_type const & columns, GenericContainer & vals ) const;

    ///@}

    //! \name Evaluate derivative to GenericContainer using another spline as independent
    ///@{

    /**
     * \brief Evaluate first derivatives using independent spline
     *
     * \param[in] zeta Target value of independent spline
     * \param[in] indep Index of independent spline
     * \param[out] x Computed x-value
     * \param[out] vals GenericContainer with derivative values
     */
    void eval2_D( real_type const zeta, integer const indep, real_type & x, GenericContainer & vals ) const;

    /**
     * \brief Evaluate first derivatives at multiple zeta values using independent spline
     *
     * \param[in] zetas Vector of target values for independent spline
     * \param[in] indep Index of independent spline
     * \param[out] vals GenericContainer with vectors of derivatives
     */
    void eval2_D( vec_real_type const & zetas, integer const indep, GenericContainer & vals ) const;

    /**
     * \brief Evaluate first derivatives of specific splines using independent spline
     *
     * \param[in] zeta Target value of independent spline
     * \param[in] indep Index of independent spline
     * \param[out] x Computed x-value
     * \param[in] columns Names of splines to evaluate
     * \param[out] vals GenericContainer with derivative values
     */
    void eval2_D(
      real_type const         zeta,
      integer const           indep,
      real_type &             x,
      vec_string_type const & columns,
      GenericContainer &      vals ) const;

    /**
     * \brief Evaluate first derivatives of specific splines at multiple zeta values
     *
     * \param[in] zetas Vector of target values for independent spline
     * \param[in] indep Index of independent spline
     * \param[in] columns Names of splines to evaluate
     * \param[out] vals GenericContainer with vectors of derivatives
     */
    void eval2_D(
      vec_real_type const &   zetas,
      integer const           indep,
      vec_string_type const & columns,
      GenericContainer &      vals ) const;

    ///@}

    //! \name Evaluate derivative to GenericContainer using another spline as independent
    ///@{

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
     * \brief Evaluate first derivatives of specific splines at multiple zeta values using independent spline by name
     *
     * \param[in] zetas Vector of target values for independent spline
     * \param[in] indep Name of independent spline
     * \param[in] columns Names of splines to evaluate
     * \param[out] vals GenericContainer with vectors of derivatives
     */
    void eval2_D(
      vec_real_type const &   zetas,
      string_view             indep,
      vec_string_type const & columns,
      GenericContainer &      vals ) const
    {
      this->eval2_D( zetas, this->get_position( indep ), columns, vals );
    }

    ///@}

    //! \name Evaluate second derivative to GenericContainer using another spline as independent
    ///@{

    /**
     * \brief Evaluate second derivatives of all splines and store in GenericContainer
     *
     * \param[in] x Parameter value
     * \param[out] vals GenericContainer with second derivative values
     */
    void eval_DD( real_type const x, GenericContainer & vals ) const;

    /**
     * \brief Evaluate second derivatives of all splines at multiple x-values
     *
     * \param[in] vec Vector of x-values
     * \param[out] vals GenericContainer with vectors of second derivatives
     */
    void eval_DD( vec_real_type const & vec, GenericContainer & vals ) const;

    /**
     * \brief Evaluate second derivatives of specific splines
     *
     * \param[in] x Parameter value
     * \param[in] columns Names of splines to evaluate
     * \param[out] vals GenericContainer with second derivative values
     */
    void eval_DD( real_type const x, vec_string_type const & columns, GenericContainer & vals ) const;

    /**
     * \brief Evaluate second derivatives of specific splines at multiple x-values
     *
     * \param[in] vec Vector of x-values
     * \param[in] columns Names of splines to evaluate
     * \param[out] vals GenericContainer with vectors of second derivatives
     */
    void eval_DD( vec_real_type const & vec, vec_string_type const & columns, GenericContainer & vals ) const;

    /**
     * \brief Evaluate second derivatives using independent spline
     *
     * \param[in] zeta Target value of independent spline
     * \param[in] indep Index of independent spline
     * \param[out] x Computed x-value
     * \param[out] vals GenericContainer with second derivative values
     */
    void eval2_DD( real_type const zeta, integer const indep, real_type & x, GenericContainer & vals ) const;

    /**
     * \brief Evaluate second derivatives at multiple zeta values using independent spline
     *
     * \param[in] zetas Vector of target values for independent spline
     * \param[in] indep Index of independent spline
     * \param[out] vals GenericContainer with vectors of second derivatives
     */
    void eval2_DD( vec_real_type const & zetas, integer const indep, GenericContainer & vals ) const;

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
      GenericContainer &      vals ) const;

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
      GenericContainer &      vals ) const;

    ///@}

    //! \name Evaluate second derivative to GenericContainer using another spline as independent
    ///@{

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
    ///@}

    //! \name Evaluate third derivative to GenericContainer
    ///@{

    /**
     * \brief Evaluate third derivatives of all splines and store in GenericContainer
     *
     * \param[in] x Parameter value
     * \param[out] vals GenericContainer with third derivative values
     */
    void eval_DDD( real_type const x, GenericContainer & vals ) const;

    /**
     * \brief Evaluate third derivatives of all splines at multiple x-values
     *
     * \param[in] vec Vector of x-values
     * \param[out] vals GenericContainer with vectors of third derivatives
     */
    void eval_DDD( vec_real_type const & vec, GenericContainer & vals ) const;

    /**
     * \brief Evaluate third derivatives of specific splines
     *
     * \param[in] x Parameter value
     * \param[in] columns Names of splines to evaluate
     * \param[out] vals GenericContainer with third derivative values
     */
    void eval_DDD( real_type const x, vec_string_type const & columns, GenericContainer & vals ) const;

    /**
     * \brief Evaluate third derivatives of specific splines at multiple x-values
     *
     * \param[in] vec Vector of x-values
     * \param[in] columns Names of splines to evaluate
     * \param[out] vals GenericContainer with vectors of third derivatives
     */
    void eval_DDD( vec_real_type const & vec, vec_string_type const & columns, GenericContainer & vals ) const;
    ///@}

    //! \name Evaluate third derivative to GenericContainer using another spline as independent
    ///@{

    /**
     * \brief Evaluate third derivatives using independent spline
     *
     * \param[in] zeta Target value of independent spline
     * \param[in] indep Index of independent spline
     * \param[out] x Computed x-value
     * \param[out] vals GenericContainer with third derivative values
     */
    void eval2_DDD( real_type const zeta, integer const indep, real_type & x, GenericContainer & vals ) const;

    /**
     * \brief Evaluate third derivatives at multiple zeta values using independent spline
     *
     * \param[in] zetas Vector of target values for independent spline
     * \param[in] indep Index of independent spline
     * \param[out] vals GenericContainer with vectors of third derivatives
     */
    void eval2_DDD( vec_real_type const & zetas, integer const indep, GenericContainer & vals ) const;

    /**
     * \brief Evaluate third derivatives of specific splines using independent spline
     *
     * \param[in] zeta Target value of independent spline
     * \param[in] indep Index of independent spline
     * \param[out] x Computed x-value
     * \param[in] columns Names of splines to evaluate
     * \param[out] vals GenericContainer with third derivative values
     */
    void eval2_DDD(
      real_type const         zeta,
      integer const           indep,
      real_type &             x,
      vec_string_type const & columns,
      GenericContainer &      vals ) const;

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
      GenericContainer &      vals ) const;
    ///@}

    //! \name Evaluate third derivative to GenericContainer using another spline as independent
    ///@{

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
     * \brief Evaluate third derivatives of specific splines at multiple zeta values using independent spline by name
     *
     * \param[in] zetas Vector of target values for independent spline
     * \param[in] indep Name of independent spline
     * \param[in] columns Names of splines to evaluate
     * \param[out] vals GenericContainer with vectors of third derivatives
     */
    void eval2_DDD(
      vec_real_type const &   zetas,
      string_view             indep,
      vec_string_type const & columns,
      GenericContainer &      vals ) const
    {
      this->eval2_DDD( zetas, this->get_position( indep ), columns, vals );
    }

    ///@}

    //! \name Build Spline
    ///@{

    ///////////////////////////////////////////////////////////////////////////
    /**
     * \brief Build the spline set from raw data arrays
     *
     * Constructs multiple splines sharing the same x-nodes. Each spline can be
     * of a different type (linear, cubic, etc.) and may have optional derivative data.
     *
     * \param[in] nspl    Number of splines
     * \param[in] npts    Number of points per spline
     * \param[in] headers Array of spline names (size = nspl)
     * \param[in] stype   Array of spline types (size = nspl)
     * \param[in] X       Array of x-node values (size = npts)
     * \param[in] Y       Array of pointers to y-value arrays (size = nspl)
     * \param[in] Yp      Optional array of pointers to derivative arrays (size = nspl, default: nullptr)
     *
     * \par Memory Management:
     * Allocates internal storage and copies input data. Input arrays can be freed
     * after calling build().
     *
     * \par Supported Spline Types:
     * - CONSTANT: Piecewise constant
     * - LINEAR: Linear interpolation
     * - CUBIC: Cubic spline (natural)
     * - AKIMA: Akima spline
     * - BESSEL: Bessel spline
     * - PCHIP: Piecewise Cubic Hermite Interpolating Polynomial
     * - HERMITE: Hermite spline (requires Yp)
     * - QUINTIC: Quintic spline
     *
     * \throw UTILS_ASSERT if nspl ≤ 0 or npts ≤ 1
     * \throw UTILS_ERROR for unsupported spline types
     *
     * \note For HERMITE splines, Yp must be provided.
     * \note Automatically checks monotonicity for applicable spline types.
     */
    void build(
      integer const           nspl,
      integer const           npts,
      char const * const      headers[],
      SplineType1D const      stype[],
      real_type const         X[],
      real_type const * const Y[],
      real_type const * const Yp[] = nullptr )
    {
      string msg{ fmt::format( "SplineSet[{}]::build(...):", m_name ) };
      UTILS_ASSERT( nspl > 0, "{} expected positive nspl = {}\n", msg, nspl );
      UTILS_ASSERT( npts > 1, "{} expected npts = {} greater than 1", msg, npts );
      m_nspl = nspl;
      m_npts = npts;
      // allocate memory
      m_splines.resize( m_nspl );
      m_is_monotone = m_mem_int.realloc( m_nspl );

      m_header_to_position.clear();

      integer mem{ npts };
      for ( integer spl{ 0 }; spl < nspl; ++spl )
      {
        switch ( stype[spl] )
        {
          case SplineType1D::QUINTIC_CUBIC:
          case SplineType1D::QUINTIC_AKIMA:
          case SplineType1D::QUINTIC_BESSEL:
          case SplineType1D::QUINTIC_PCHIP:
            mem += npts;  // Y, Yp, Ypp
            [[fallthrough]];
          case SplineType1D::CUBIC:
          case SplineType1D::AKIMA:
          case SplineType1D::BESSEL:
          case SplineType1D::PCHIP:
          case SplineType1D::HERMITE:
            mem += npts;  // Y, Yp
            [[fallthrough]];
          case SplineType1D::CONSTANT:
          case SplineType1D::LINEAR: mem += npts; break;
          case SplineType1D::SPLINE_SET:
          case SplineType1D::SPLINE_VEC:
            // default:
            UTILS_ERROR(
              "{} At spline n.{} named {} cannot be done for type = {}\n",
              msg,
              spl,
              headers[spl],
              to_string( stype[spl] ) );
        }
      }

      m_mem.reallocate( mem + 2 * nspl );
      m_mem_p.reallocate( 3 * nspl );

      m_Y    = m_mem_p( m_nspl );
      m_Yp   = m_mem_p( m_nspl );
      m_Ypp  = m_mem_p( m_nspl );
      m_X    = m_mem( m_npts );
      m_Ymin = m_mem( m_nspl );
      m_Ymax = m_mem( m_nspl );

      copy_n( X, npts, m_X );
      for ( integer spl{ 0 }; spl < nspl; ++spl )
      {
        real_type *& pY{ m_Y[spl] };
        real_type *& pYp{ m_Yp[spl] };
        real_type *& pYpp{ m_Ypp[spl] };
        pY = m_mem( m_npts );
        copy_n( Y[spl], npts, pY );
        if ( stype[spl] == SplineType1D::CONSTANT )
        {
          m_Ymin[spl] = *std::min_element( pY, pY + npts - 1 );
          m_Ymax[spl] = *std::max_element( pY, pY + npts - 1 );
        }
        else
        {
          m_Ymin[spl] = *std::min_element( pY, pY + npts );
          m_Ymax[spl] = *std::max_element( pY, pY + npts );
        }
        pYpp = pYp = nullptr;
        switch ( stype[spl] )
        {
          case SplineType1D::QUINTIC_CUBIC:
          case SplineType1D::QUINTIC_AKIMA:
          case SplineType1D::QUINTIC_BESSEL:
          case SplineType1D::QUINTIC_PCHIP: pYpp = m_mem( m_npts ); [[fallthrough]];
          case SplineType1D::CUBIC:
          case SplineType1D::AKIMA:
          case SplineType1D::BESSEL:
          case SplineType1D::PCHIP:
          case SplineType1D::HERMITE:
            pYp = m_mem( m_npts );
            if ( stype[spl] == SplineType1D::HERMITE )
            {
              UTILS_ASSERT(
                Yp != nullptr && Yp[spl] != nullptr,
                "{} At spline n.{} named {}\n"
                "expect to find derivative values",
                msg,
                spl,
                headers[spl] );
              copy_n( Yp[spl], npts, pYp );
            }
            [[fallthrough]];
          case SplineType1D::CONSTANT:
          case SplineType1D::LINEAR:
          case SplineType1D::SPLINE_SET:
          case SplineType1D::SPLINE_VEC:
            // default:
            break;
        }
        string                    h{ headers[spl] };
        std::unique_ptr<Spline> & s{ m_splines[spl] };

        m_is_monotone[spl] = -1;
        switch ( stype[spl] )
        {
          case SplineType1D::CONSTANT:
          {
            auto S = std::make_unique<ConstantSpline>( h );
            S->reserve_external( m_npts, m_X, pY );
            S->m_npts = m_npts;
            S->build();
            s = std::move( S );
          }
          break;

          case SplineType1D::LINEAR:
          {
            auto S = std::make_unique<LinearSpline>( h );
            S->reserve_external( m_npts, m_X, pY );
            S->m_npts = m_npts;
            S->build();
            // check monotonicity of data
            integer flag{ 1 };
            for ( integer j = 1; j < m_npts; ++j )
            {
              if ( pY[j - 1] > pY[j] )
              {
                flag = -1;
                break;
              }  // non monotone data
              if ( Utils::is_zero( pY[j - 1] - pY[j] ) && m_X[j - 1] < m_X[j] ) flag = 0;  // non strict monotone
            }
            m_is_monotone[spl] = flag;
            s                  = std::move( S );
          }
          break;

          case SplineType1D::CUBIC:
          {
            auto S = std::make_unique<CubicSpline>( h );
            S->reserve_external( m_npts, m_X, pY, pYp );
            S->m_npts = m_npts;
            S->build();
            m_is_monotone[spl] = S->is_monotone();
            s                  = std::move( S );
          }
          break;

          case SplineType1D::AKIMA:
          {
            auto S = std::make_unique<AkimaSpline>( h );
            S->reserve_external( m_npts, m_X, pY, pYp );
            S->m_npts = m_npts;
            S->build();
            m_is_monotone[spl] = S->is_monotone();
            s                  = std::move( S );
          }
          break;

          case SplineType1D::BESSEL:
          {
            auto S = std::make_unique<BesselSpline>( h );
            S->reserve_external( m_npts, m_X, pY, pYp );
            S->m_npts = m_npts;
            S->build();
            m_is_monotone[spl] = S->is_monotone();
            s                  = std::move( S );
          }
          break;

          case SplineType1D::PCHIP:
          {
            auto S = std::make_unique<PchipSpline>( h );
            S->reserve_external( m_npts, m_X, pY, pYp );
            S->m_npts = m_npts;
            S->build();
            m_is_monotone[spl] = S->is_monotone();
            s                  = std::move( S );
          }
          break;

          case SplineType1D::HERMITE:
          {
            auto S = std::make_unique<HermiteSpline>( h );
            S->reserve_external( m_npts, m_X, pY, pYp );
            S->m_npts = m_npts;
            S->build();
            m_is_monotone[spl] = S->is_monotone();
            s                  = std::move( S );
          }
          break;

          case SplineType1D::QUINTIC_CUBIC:
          {
            auto S = std::make_unique<QuinticSpline>( Spline_sub_type::CUBIC, h );
            S->reserve_external( m_npts, m_X, pY, pYp, pYpp );
            S->m_npts = m_npts;
            S->build();
            m_is_monotone[spl] = S->is_monotone();
            s                  = std::move( S );
          }
          break;

          case SplineType1D::QUINTIC_AKIMA:
          {
            auto S = std::make_unique<QuinticSpline>( Spline_sub_type::AKIMA, h );
            S->reserve_external( m_npts, m_X, pY, pYp, pYpp );
            S->m_npts = m_npts;
            S->build();
            m_is_monotone[spl] = S->is_monotone();
            s                  = std::move( S );
          }
          break;

          case SplineType1D::QUINTIC_BESSEL:
          {
            auto S = std::make_unique<QuinticSpline>( Spline_sub_type::BESSEL, h );
            S->reserve_external( m_npts, m_X, pY, pYp, pYpp );
            S->m_npts = m_npts;
            S->build();
            m_is_monotone[spl] = S->is_monotone();
            s                  = std::move( S );
          }
          break;

          case SplineType1D::QUINTIC_PCHIP:
          {
            auto S = std::make_unique<QuinticSpline>( Spline_sub_type::PCHIP, h );
            S->reserve_external( m_npts, m_X, pY, pYp, pYpp );
            S->m_npts = m_npts;
            S->build();
            m_is_monotone[spl] = S->is_monotone();
            s                  = std::move( S );
          }
          break;

          case SplineType1D::SPLINE_SET:
          case SplineType1D::SPLINE_VEC:
            // default:
            UTILS_ERROR(
              "{} At spline n.{} named {}\n"
              "{} not allowed as spline type\n"
              "in SplineSet::build for {}-th spline\n",
              msg,
              spl,
              headers[spl],
              to_string( stype[spl] ),
              spl );
        }
        m_header_to_position.insert( { s->name().data(), static_cast<integer>( spl ) } );
      }

      m_mem.must_be_empty( "SplineSet::build, baseValue" );
      m_mem_p.must_be_empty( "SplineSet::build, basePointer" );
    }

    /**
     * \brief Create a deep copy of this spline set into another
     *
     * Copies all spline data, names, and types to the target SplineSet.
     * The target is completely overwritten and will be identical to this one.
     *
     * \param[out] S Target SplineSet to receive the copy
     *
     * \par Note:
     * Also copies flags (extrapolation, closed curve, etc.) from each spline.
     */
    void deep_copy_to( SplineSet & S ) const
    {
      std::vector<char const *> headers( m_nspl );
      std::vector<SplineType1D> stype( m_nspl );

      for ( integer i = 0; i < m_nspl; ++i )
      {
        headers[i] = m_splines[i]->name().data();
        stype[i]   = m_splines[i]->type();
      }

      S.build( m_nspl, m_npts, headers.data(), stype.data(), m_X, m_Y, m_Yp );

      // propagate flags
      for ( integer i = 0; i < m_nspl; ++i ) S.m_splines[i]->copy_flags( *m_splines[i] );
    }

    /**
     * \brief Build spline set from GenericContainer
     *
     * \param[in] gc GenericContainer containing spline set configuration
     *
     * \par Expected GC Structure:
     * Typically contains matrices or named arrays for x and y values,
     * plus metadata about spline types and names.
     *
     * \see setup() for detailed format specification
     */
    void setup( GenericContainer const & gc );

    /**
     * \brief Build spline set from GenericContainer
     *
     * \copydetails setup()
     * This is an alias for setup().
     */
    void build( GenericContainer const & gc ) { this->setup( gc ); }

    ///@}

    /**
     * \brief Get the type identifier of this spline set
     *
     * \return SplineType1D::SPLINE_SET enum value
     */
    SplineType1D type() const { return SplineType1D::SPLINE_SET; }

    /**
     * \brief Generate detailed informational string about the spline set
     *
     * \return Multi-line string containing:
     *         - Name of spline set
     *         - Number of points and splines
     *         - For each spline: index, monotonicity, and spline-specific info
     *
     * \par Example Output:
     * \code
     * SplineSet[my_set] n.points=10 n.splines=3
     * Spline n.0 is strictly monotone
     * ConstantSpline[name0] n.points = 10
     * Spline n.1 is NOT monotone
     * CubicSpline[name1] n.points = 10
     * ...
     * \endcode
     */
    string info() const
    {
      string res{ fmt::format( "SplineSet[{}] n.points={} n.splines={}", name(), m_npts, m_nspl ) };

      for ( integer i = 0; i < m_nspl; ++i )
      {
        res += fmt::format( "\nSpline n.{} ", i );
        switch ( m_is_monotone[i] )
        {
          case -2: res += " with NON monotone data\n"; break;
          case -1: res += " is NOT monotone\n"; break;
          case 0: res += " is monotone\n"; break;
          case 1: res += " is strictly monotone\n"; break;
          default: UTILS_ERROR( "SplineSet::info classification: {} not in range {-2,-1,0,1}\n", m_is_monotone[i] );
        }
        res += m_splines[i]->info();
      }
      return res;
    }

    /**
     * \brief Print informational string to output stream
     *
     * \param[out] stream Output stream (e.g., std::cout, std::ofstream)
     *
     * \see info() for output format
     */
    void info( ostream_type & stream ) const { stream << this->info() << '\n'; }

    /**
     * \brief Generate a table of spline values for plotting
     *
     * Samples all splines uniformly over the x-domain and writes a
     * tab-separated table to the output stream.
     *
     * \param[out] stream Output stream for the table
     * \param[in] num_points Number of sample points (including endpoints)
     *
     * \par Output Format:
     * ```
     * s    spline0    spline1    ...    splineN
     * 0.0  val00      val01      ...    val0N
     * 0.1  val10      val11      ...    val1N
     * ...  ...        ...        ...    ...
     * 1.0  valM0      valM1      ...    valMN
     * ```
     * where M = num_points-1, N = m_nspl-1
     */
    void dump_table( ostream_type & stream, integer const num_points ) const
    {
      Malloc_real mem( "SplineSet::dump_table" );
      mem.allocate( m_nspl );
      real_type * vals{ mem( m_nspl ) };
      stream << 's';
      for ( integer i = 0; i < m_nspl; ++i ) stream << '\t' << header( i );
      stream << '\n';

      for ( integer j = 0; j < num_points; ++j )
      {
        real_type const s{ x_min() + ( ( x_max() - x_min() ) * j ) / ( num_points - 1 ) };
        this->eval( s, vals, 1 );
        stream << s;
        for ( integer i = 0; i < m_nspl; ++i ) stream << '\t' << vals[i];
        stream << '\n';
      }
    }

#ifdef SPLINES_BACK_COMPATIBILITY
    /**
     * \brief Backward compatibility: get monotonicity information
     * \param[in] i Spline index
     * \return Monotonicity code
     */
    int isMonotone( integer i ) const { return m_is_monotone[i]; }

    /**
     * \brief Backward compatibility: get number of points
     * \return Number of points per spline
     */
    integer numPoints() const { return m_npts; }

    /**
     * \brief Backward compatibility: get number of splines
     * \return Number of splines in set
     */
    integer numSplines() const { return m_nspl; }

    /**
     * \brief Backward compatibility: find spline index by name
     * \param[in] hdr Spline name
     * \return Spline index
     */
    integer getPosition( string_view hdr ) const { return get_position( hdr ); }

    /**
     * \brief Backward compatibility: get x-nodes array
     * \return Pointer to x-values
     */
    real_type const * xNodes() const { return m_X; }

    /**
     * \brief Backward compatibility: get y-nodes array for spline
     * \param[in] i Spline index
     * \return Pointer to y-values for spline i
     */
    real_type const * yNodes( integer i ) const;

    /**
     * \brief Backward compatibility: get specific x-node
     * \param[in] npt Node index
     * \return X-value at node npt
     */
    real_type xNode( integer npt ) const { return this->x_node( npt ); }

    /**
     * \brief Backward compatibility: get specific y-node
     * \param[in] npt Node index
     * \param[in] spl Spline index
     * \return Y-value of spline spl at node npt
     */
    real_type yNode( integer npt, integer spl ) const { return this->y_node( npt, spl ); }

    /**
     * \brief Backward compatibility: get minimum x
     * \return Minimum x-value
     */
    real_type xMin() const { return this->x_min(); }

    /**
     * \brief Backward compatibility: get maximum x
     * \return Maximum x-value
     */
    real_type xMax() const { return this->x_max(); }

    /**
     * \brief Backward compatibility: get minimum y for spline
     * \param[in] spl Spline index
     * \return Minimum y-value of spline spl
     */
    real_type yMin( integer spl ) const { return this->y_min( spl ); }

    /**
     * \brief Backward compatibility: get maximum y for spline
     * \param[in] spl Spline index
     * \return Maximum y-value of spline spl
     */
    real_type yMax( integer spl ) const { return this->y_max( spl ); }

    /**
     * \brief Backward compatibility: get minimum y by name
     * \param[in] spl Spline name
     * \return Minimum y-value of named spline
     */
    real_type yMin( string_view spl ) const { return this->y_min( spl ); }

    /**
     * \brief Backward compatibility: get maximum y by name
     * \param[in] spl Spline name
     * \return Maximum y-value of named spline
     */
    real_type yMax( string_view spl ) const { return this->y_max( spl ); }

    /**
     * \brief Backward compatibility: get spline by index
     * \param[in] i Spline index
     * \return Pointer to spline i
     */
    Spline * getSpline( integer i ) const { return this->get_spline( i ); }

    /**
     * \brief Backward compatibility: get spline by name
     * \param[in] hdr Spline name
     * \return Pointer to named spline
     */
    Spline * getSpline( string_view hdr ) const { return this->get_spline( hdr ); }
#endif
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

  using GC_namespace::GC_type;
  using GC_namespace::map_type;
  using GC_namespace::mat_real_type;
  using GC_namespace::vec_int_type;
  using GC_namespace::vec_real_type;
  using GC_namespace::vec_string_type;
  using GC_namespace::vector_type;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline void SplineSet::setup( GenericContainer const & gc )
  {
    string const where{ fmt::format( "SplineSet[{}]::setup( gc ): ", m_name ) };

    std::set<std::string> keywords;
    for ( auto const & pair : gc.get_map( where ) ) { keywords.insert( pair.first ); }

    /*
    // gc["spline_type"]
    // gc["xdata"]
    // gc["ydata"]
    // gc["headers"] (opzionale)
    */
    vec_string_type       spline_type_vec;
    vec_real_type         X;
    vector<SplineType1D>  stype;
    vec_string_type       headers;
    vector<vec_real_type> Y, Yp;

    GenericContainer const & gc_stype{ gc( "spline_type", where ) };
    keywords.erase( "spline_type" );
    GenericContainer const & gc_xdata{ gc( "xdata", where ) };
    keywords.erase( "xdata" );
    GenericContainer const & gc_ydata{ gc( "ydata", where ) };
    keywords.erase( "ydata" );

    //
    // gc["spline_type"]
    //
    gc_stype.copyto_vec_string( spline_type_vec, ( where + "in reading `spline_type'\n" ) );
    m_nspl = static_cast<integer>( spline_type_vec.size() );
    stype.resize( m_nspl );
    for ( integer spl{ 0 }; spl < m_nspl; ++spl ) stype[spl] = string_to_splineType1D( spline_type_vec[spl] );

    // se non tipo MAP deve esserci headers
    //
    // gc["headers"] (opzionale)
    //
    if ( GC_type::MAP != gc_ydata.get_type() )
    {
      GenericContainer const & gc_headers{ gc( "headers", where ) };
      keywords.erase( "headers" );
      gc_headers.copyto_vec_string( headers, ( where + ", reading `headers'" ) );
      UTILS_ASSERT(
        headers.size() == static_cast<size_t>( m_nspl ),
        "{}, field `headers` expected to be of size {} found of size {}\n",
        where,
        m_nspl,
        headers.size() );
    }

    //
    // gc["xdata"]
    //
    gc_xdata.copyto_vec_real( X, ( where + "reading `xdata'" ) );
    m_npts = static_cast<integer>( X.size() );

    // allocate for _nspl splines
    Y.resize( m_nspl );
    Yp.resize( m_nspl );

    switch ( gc_ydata.get_type() )
    {
      case GC_type::VEC_BOOL:
      case GC_type::VEC_INTEGER:
      case GC_type::VEC_LONG:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      {
        // leggo vettore come matrice di una sola colonna
        vec_real_type data;
        gc_ydata.copyto_vec_real( data, ( where + "reading `ydata'" ) );
        UTILS_ASSERT(
          m_nspl == 1,
          "{}, number of splines [{}]\n"
          "is incompatible with the type of `ydata` a vector (or a matrix with 1 column) in data\n",
          where,
          m_nspl );
        UTILS_ASSERT(
          static_cast<size_t>( m_npts ) == data.size(),
          "{}, number of points [{}]\n"
          "differs from the number of `ydata` rows [{}] in data\n",
          where,
          m_npts,
          data.size() );
        Y[0].reserve( data.size() );
        std::copy( data.begin(), data.end(), std::back_inserter( Y[0] ) );
      }
      break;
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      {
        mat_real_type data;
        gc_ydata.copyto_mat_real( data, ( where + "reading `ydata'" ) );
        UTILS_ASSERT(
          static_cast<size_t>( m_nspl ) == data.num_cols(),
          "{}, number of splines [{}]\n"
          "differs from the number of `ydata` columns [{}] in data\n",
          where,
          m_nspl,
          data.num_cols() );
        UTILS_ASSERT(
          static_cast<size_t>( m_npts ) == data.num_rows(),
          "{}, number of points [{}]\n"
          "differs from the number of `ydata` rows [{}] in data\n",
          where,
          m_npts,
          data.num_rows() );
        for ( integer i = 0; i < m_nspl; ++i ) data.get_column( static_cast<unsigned>( i ), Y[i] );
      }
      break;
      case GC_type::VECTOR:
      {
        vector_type const & data{ gc_ydata.get_vector() };
        UTILS_ASSERT(
          static_cast<size_t>( m_nspl ) == data.size(),
          "{}, field `ydata` expected of size {} found of size {}\n",
          where,
          m_nspl,
          data.size() );
        string msg1{ where + " reading `ydata` columns" };
        for ( integer spl = 0; spl < m_nspl; ++spl )
        {
          GenericContainer const & datai{ data[spl] };
          integer                  nrow{ m_npts };
          if ( stype[spl] == SplineType1D::CONSTANT ) --nrow;  // constant spline uses n-1 points
          datai.copyto_vec_real( Y[spl], msg1 );
          UTILS_ASSERT(
            static_cast<size_t>( nrow ) == Y[spl].size(),
            "{}, column {} of `ydata` of type `{}` expected of size {} found of size {}\n",
            where,
            spl,
            spline_type_vec[spl],
            nrow,
            Y[spl].size() );
        }
      }
      break;
      case GC_type::MAP:
      {
        map_type const & data{ gc_ydata.get_map() };
        UTILS_ASSERT(
          data.size() == static_cast<size_t>( m_nspl ),
          "{}, field `ydata` expected of size {} found of size {}\n",
          where,
          m_nspl,
          data.size() );
        headers.clear();
        headers.reserve( data.size() );
        auto   im{ data.begin() };
        string msg1{ where + " reading `ydata` columns" };
        for ( integer spl{ 0 }; im != data.end(); ++im, ++spl )
        {
          headers.emplace_back( im->first );
          GenericContainer const & datai{ im->second };
          integer                  nrow{ m_npts };
          if ( stype[spl] == SplineType1D::CONSTANT ) --nrow;  // constant spline uses n-1 points
          datai.copyto_vec_real( Y[spl], msg1 );
          UTILS_ASSERT(
            static_cast<size_t>( nrow ) == Y[spl].size(),
            "{}, column `{}` of `ydata` ot type `{}` expected of size {} found of size {}\n",
            where,
            im->first,
            spline_type_vec[spl],
            nrow,
            Y[spl].size() );
        }
      }
      break;
      default:
        UTILS_ERROR(
          "{}, field `data` expected\n"
          "to be of type `vec_[int/long/real]_type`, `mat_[int/long/real]_type`, `vector_type` or `map_type'\n"
          "found: `{}`\n",
          where,
          gc_ydata.get_type_name() );
        break;
    }

    if ( gc.exists( "ypdata" ) )
    {  // yp puo' essere solo tipo map
      GenericContainer const & gc_ypdata{ gc( "ypdata", where ) };
      keywords.erase( "ypdata" );
      UTILS_ASSERT(
        GC_type::MAP == gc_ypdata.get_type(),
        "{}, field `ypdata` expected to be of type `map_type` found: `{}`\n",
        where,
        gc_ypdata.get_type_name() );

      std::map<string, integer>       h_to_pos;
      vec_string_type::const_iterator is{ headers.begin() };
      for ( integer idx{ 0 }; is != headers.end(); ++is ) h_to_pos[*is] = idx++;

      string           msg1{ where + " reading `ypdata` columns" };
      map_type const & data{ gc_ypdata.get_map() };
      for ( auto const & [fst, snd] : data )
      {
        // cerca posizione
        auto is_pos{ h_to_pos.find( fst ) };
        UTILS_ASSERT( is_pos != h_to_pos.end(), "{}, column `{}` of `ypdata` not found\n", where, fst );
        integer spl{ is_pos->second };

        GenericContainer const & datai{ snd };
        integer                  nrow{ m_npts };
        if ( stype[spl] == SplineType1D::CONSTANT ) --nrow;  // constant spline uses n-1 points
        datai.copyto_vec_real( Yp[spl], msg1 );
        UTILS_ASSERT(
          static_cast<size_t>( nrow ) == Y[spl].size(),
          "{}, column `{}` of `ypdata` or type `{}` expected of size {} found of size {}\n",
          where,
          fst,
          spline_type_vec[spl],
          nrow,
          Y[spl].size() );
      }
    }

    Utils::Malloc<void *> mem( where );
    mem.allocate( 3 * m_nspl );

    void ** pp__headers = mem( m_nspl );
    void ** pp__Y       = mem( m_nspl );
    void ** pp__Yp      = mem( m_nspl );

    for ( integer spl{ 0 }; spl < m_nspl; ++spl )
    {
      pp__headers[spl] = const_cast<void *>( reinterpret_cast<void const *>( headers[spl].data() ) );
      pp__Y[spl]       = reinterpret_cast<void *>( Y[spl].data() );
      pp__Yp[spl]      = reinterpret_cast<void *>( Yp[spl].empty() ? nullptr : Yp[spl].data() );
    }

    this->build(
      m_nspl,
      m_npts,
      reinterpret_cast<char const **>( const_cast<void const **>( pp__headers ) ),
      stype.data(),
      X.data(),
      reinterpret_cast<real_type const **>( const_cast<void const **>( pp__Y ) ),
      reinterpret_cast<real_type const **>( const_cast<void const **>( pp__Yp ) ) );

    if ( gc.exists( "boundary" ) )
    {
      GenericContainer const & gc_boundary{ gc( "boundary" ) };
      keywords.erase( "boundary" );
      integer ne{ static_cast<integer>( gc_boundary.get_num_elements() ) };
      UTILS_ASSERT(
        ne == m_nspl,
        "{}, field `boundary` expected a"
        " generic vector of size: {} but is of size: {}\n",
        where,
        ne,
        m_nspl );

      for ( integer ispl{ 0 }; ispl < ne; ++ispl )
      {
        Spline *                 S{ m_splines[ispl].get() };
        GenericContainer const & item{ gc_boundary( ispl, "SplineSet boundary data" ) };

        std::set<std::string> keywords2;
        for ( auto const & pair : item.get_map( where ) ) { keywords2.insert( pair.first ); }


        bool is_closed{ false };
        item.get_if_exists( "closed", is_closed );
        keywords2.erase( "closed" );
        if ( is_closed ) { S->make_closed(); }
        else
        {
          keywords2.erase( "extend" );
          keywords2.erase( "can_extend" );
          S->make_opened();
          bool can_extend{ false };
          if ( !item.get_if_exists( "extend", can_extend ) ) item.get_if_exists( "can_extend", can_extend );
          if ( can_extend )
          {
            S->make_unbounded();
            bool extend_constant{ false };
            keywords2.erase( "extend_constant" );
            item.get_if_exists( "extend_constant", extend_constant );
            if ( extend_constant )
              S->make_extended_constant();
            else
              S->make_extended_not_constant();
          }
          else
          {
            S->make_bounded();
          }
        }
        UTILS_WARNING(
          keywords2.empty(),
          "{}, spline N.{} of {} unused keys\n{}\n",
          where,
          ispl,
          ne,
          [&keywords2]() -> string
          {
            string res;
            for ( auto const & it : keywords2 )
            {
              res += it;
              res += ' ';
            };
            return res;
          }() );
      }
    }

    UTILS_WARNING(
      keywords.empty(),
      "{}: unused keys\n{}\n",
      where,
      [&keywords]() -> string
      {
        string res;
        for ( auto const & it : keywords )
        {
          res += it;
          res += ' ';
        };
        return res;
      }() );
  }

  //!
  //! Evaluate all the splines at `x` and
  //! fill a `map of values` in a GenericContainer
  //!
  inline void SplineSet::eval( real_type const x, GenericContainer & gc ) const
  {
    map_type & vals{ gc.set_map() };
    for ( auto const & [fst, snd] : m_header_to_position ) vals[fst] = m_splines[snd]->eval( x );
  }

  //!
  //! Evaluate all the splines at `x` values contained in vec
  //! and fill a `map of vector` in a GenericContainer
  //!
  inline void SplineSet::eval( vec_real_type const & vec, GenericContainer & gc ) const
  {
    integer const npts{ static_cast<integer>( vec.size() ) };
    map_type &    vals{ gc.set_map() };
    for ( auto const & [fst, snd] : m_header_to_position )
    {
      vec_real_type & v{ vals[fst].set_vec_real( snd ) };
      Spline const *  p_spl{ m_splines[snd].get() };
      for ( integer i = 0; i < npts; ++i ) v[i] = p_spl->eval( vec[i] );
    }
  }

  //!
  //! Evaluate all the splines at `x` and fill a map of values
  //! in a GenericContainer with keys in `columns`
  //!
  inline void SplineSet::eval( real_type const x, vec_string_type const & columns, GenericContainer & gc ) const
  {
    map_type & vals{ gc.set_map() };
    for ( auto const & S : columns )
    {
      Spline const * p_spl{ get_spline( S ) };
      vals[S] = p_spl->eval( x );
    }
  }

  //!
  //! Evaluate all the splines at `x` values contained in vec
  //! and fill a map of vector in a GenericContainer with keys in `columns`
  //!
  inline void SplineSet::eval( vec_real_type const & vec, vec_string_type const & columns, GenericContainer & gc ) const
  {
    integer const npts{ static_cast<integer>( vec.size() ) };
    map_type &    vals{ gc.set_map() };
    for ( auto const & S : columns )
    {
      Spline const *  p_spl{ get_spline( S ) };
      vec_real_type & v{ vals[S].set_vec_real( npts ) };
      for ( integer i = 0; i < npts; ++i ) v[i] = p_spl->eval( vec[i] );
    }
  }

  //!
  //! Evaluate all the splines at `zeta` and fill a `map of values`
  //! in a GenericContainer and `indep` as independent spline
  //!
  inline void SplineSet::eval2( real_type const zeta, integer const indep, real_type & x, GenericContainer & gc ) const
  {
    map_type & vals{ gc.set_map() };
    intersect( indep, zeta, x );
    for ( auto const & [fst, snd] : m_header_to_position ) vals[fst] = m_splines[snd]->eval( x );
  }

  //!
  //! Evaluate all the splines at `zeta` values contained in vec
  //! and fill a `map of vector` in a GenericContainer and `indep`
  //! as independent spline
  //!
  inline void SplineSet::eval2( vec_real_type const & zetas, integer const indep, GenericContainer & gc ) const
  {
    integer const npts{ static_cast<integer>( zetas.size() ) };
    map_type &    vals{ gc.set_map() };

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

  //!
  //! Evaluate all the splines at `zeta` and fill a map of values
  //! in a GenericContainer with keys in `columns` and `indep`
  //! as independent spline
  //!
  inline void SplineSet::eval2(
    real_type const         zeta,
    integer const           indep,
    real_type &             x,
    vec_string_type const & columns,
    GenericContainer &      gc ) const
  {
    map_type & vals{ gc.set_map() };
    intersect( indep, zeta, x );
    for ( auto const & S : columns )
    {
      Spline const * p_spl{ get_spline( S ) };
      vals[S] = p_spl->eval( x );
    }
  }

  //!
  //! Evaluate all the splines at `zeta` values contained
  //! in vec and fill a map of vector in a GenericContainer
  //! with keys in `columns` and `indep` as independent spline
  //!
  inline void SplineSet::eval2(
    vec_real_type const &   zetas,
    integer const           indep,
    vec_string_type const & columns,
    GenericContainer &      gc ) const
  {
    integer const npts{ static_cast<integer>( zetas.size() ) };
    map_type &    vals{ gc.set_map() };

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

  /*
  //    __ _          _         _           _            _   _
  //   / _(_)_ __ ___| |_    __| | ___ _ __(_)_   ____ _| |_(_)_   _____
  //  | |_| | '__/ __| __|  / _` |/ _ \ '__| \ \ / / _` | __| \ \ / / _ \
  //  |  _| | |  \__ \ |_  | (_| |  __/ |  | |\ V / (_| | |_| |\ V /  __/
  //  |_| |_|_|  |___/\__|  \__,_|\___|_|  |_| \_/ \__,_|\__|_| \_/ \___|
  */
  inline void SplineSet::eval_D( real_type const x, GenericContainer & gc ) const
  {
    map_type & vals{ gc.set_map() };
    for ( auto const & [fst, snd] : m_header_to_position ) vals[fst] = m_splines[snd]->eval_D( x );
  }

  //!
  //! Evaluate all the splines at `x` values contained in vec and
  //! fill a `map of vector` in a GenericContainer
  //!
  inline void SplineSet::eval_D( vec_real_type const & vec, GenericContainer & gc ) const
  {
    integer const npts{ static_cast<integer>( vec.size() ) };
    map_type &    vals{ gc.set_map() };
    for ( auto const & [fst, snd] : m_header_to_position )
    {
      vec_real_type & v{ vals[fst].set_vec_real( npts ) };
      Spline const *  p_spl{ m_splines[snd].get() };
      for ( integer i = 0; i < npts; ++i ) v[i] = p_spl->eval_D( vec[i] );
    }
  }

  //!
  //! Evaluate all the splines at `x` and fill a map of values in
  //! a GenericContainer with keys in `columns`
  //!
  inline void SplineSet::eval_D( real_type const x, vec_string_type const & columns, GenericContainer & gc ) const
  {
    map_type & vals{ gc.set_map() };
    for ( auto const & S : columns )
    {
      Spline const * p_spl{ get_spline( S ) };
      vals[S] = p_spl->eval_D( x );
    }
  }

  //!
  //! Evaluate all the splines at `x` values contained in vec
  //! and fill a map of vector in a GenericContainer with keys in `columns`
  //!
  inline void SplineSet::eval_D( vec_real_type const & vec, vec_string_type const & columns, GenericContainer & gc )
    const
  {
    integer const npts{ static_cast<integer>( vec.size() ) };
    map_type &    vals{ gc.set_map() };
    for ( auto const & S : columns )
    {
      Spline const *  p_spl{ get_spline( S ) };
      vec_real_type & v{ vals[S].set_vec_real( npts ) };
      for ( integer i = 0; i < npts; ++i ) v[i] = p_spl->eval_D( vec[i] );
    }
  }

  //!
  //! Evaluate all the splines at `zeta` and fill a `map of values`
  //! in a GenericContainer and `indep` as independent spline
  //!
  inline void SplineSet::eval2_D( real_type const zeta, integer const indep, real_type & x, GenericContainer & gc )
    const
  {
    map_type & vals{ gc.set_map() };
    intersect( indep, zeta, x );
    for ( auto const & [fst, snd] : m_header_to_position ) vals[fst] = m_splines[snd]->eval_D( x );
  }

  //!
  //! Evaluate all the splines at `zeta` values contained in vec
  //! and fill a `map of vector` in a GenericContainer and `indep`
  //! as independent spline
  //!
  inline void SplineSet::eval2_D( vec_real_type const & zetas, integer const indep, GenericContainer & gc ) const
  {
    integer const npts{ static_cast<integer>( zetas.size() ) };
    map_type &    vals{ gc.set_map() };

    // pre-allocation
    for ( auto const & [fst, snd] : m_header_to_position ) vals[fst].set_vec_real( npts );

    for ( integer i = 0; i < npts; ++i )
    {
      real_type x;
      intersect( indep, zetas[i], x );
      for ( auto const & [fst, snd] : m_header_to_position )
      {
        vec_real_type & v{ vals[fst].get_vec_real() };
        v[i] = m_splines[snd]->eval_D( x );
      }
    }
  }

  //!
  //! Evaluate all the splines at `zeta` and fill a map of values
  //! in a GenericContainer with keys in `columns` and `indep`
  //! as independent spline
  //!
  inline void SplineSet::eval2_D(
    real_type const         zeta,
    integer const           indep,
    real_type &             x,
    vec_string_type const & columns,
    GenericContainer &      gc ) const
  {
    map_type & vals{ gc.set_map() };
    intersect( indep, zeta, x );
    for ( auto const & S : columns )
    {
      Spline const * p_spl{ get_spline( S ) };
      vals[S] = p_spl->eval_D( x );
    }
  }

  //!
  //! Evaluate all the splines at `zeta` values contained in vec
  //! and fill a map of vector in a GenericContainer with keys
  //! in `columns` and `indep` as independent spline
  //!
  inline void SplineSet::eval2_D(
    vec_real_type const &   zetas,
    integer const           indep,
    vec_string_type const & columns,
    GenericContainer &      gc ) const
  {
    integer const npts{ static_cast<integer>( zetas.size() ) };
    map_type &    vals{ gc.set_map() };

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
        v[i] = p_spl->eval_D( x );
      }
    }
  }

  /*                                _       _           _            _   _
  //   ___  ___  ___ ___  _ __   __| |   __| | ___ _ __(_)_   ____ _| |_(_)_   _____
  //  / __|/ _ \/ __/ _ \| '_ \ / _` |  / _` |/ _ \ '__| \ \ / / _` | __| \ \ / / _ \
  //  \__ \  __/ (_| (_) | | | | (_| | | (_| |  __/ |  | |\ V / (_| | |_| |\ V /  __/
  //  |___/\___|\___\___/|_| |_|\__,_|  \__,_|\___|_|  |_| \_/ \__,_|\__|_| \_/ \___|
  */
  inline void SplineSet::eval_DD( real_type const x, GenericContainer & gc ) const
  {
    map_type & vals{ gc.set_map() };
    for ( auto const & [fst, snd] : m_header_to_position ) vals[fst] = m_splines[snd]->eval_DD( x );
  }

  //!
  //! Evaluate all the splines at `x` values contained in vec
  //! and fill a `map of vector` in a GenericContainer
  //!
  inline void SplineSet::eval_DD( vec_real_type const & vec, GenericContainer & gc ) const
  {
    integer const npts{ static_cast<integer>( vec.size() ) };
    map_type &    vals{ gc.set_map() };
    for ( auto const & [fst, snd] : m_header_to_position )
    {
      vec_real_type & v = vals[fst].set_vec_real( npts );
      auto &          p_spl{ m_splines[snd] };
      for ( integer i = 0; i < npts; ++i ) v[i] = p_spl->eval_DD( vec[i] );
    }
  }

  //!
  //! Evaluate all the splines at `x` and fill a map of values in
  //! a GenericContainer with keys in `columns`
  //!
  inline void SplineSet::eval_DD( real_type const x, vec_string_type const & columns, GenericContainer & gc ) const
  {
    map_type & vals{ gc.set_map() };
    for ( auto const & S : columns )
    {
      Spline const * p_spl{ get_spline( S ) };
      vals[S] = p_spl->eval_DD( x );
    }
  }

  //!
  //! Evaluate all the splines at `x` values contained in vec and
  //! fill a map of vector in a GenericContainer with keys in `columns`
  //!
  inline void SplineSet::eval_DD( vec_real_type const & vec, vec_string_type const & columns, GenericContainer & gc )
    const
  {
    integer const npts{ static_cast<integer>( vec.size() ) };
    map_type &    vals{ gc.set_map() };
    for ( auto const & S : columns )
    {
      Spline const *  p_spl{ get_spline( S ) };
      vec_real_type & v{ vals[S].set_vec_real( npts ) };
      for ( integer i = 0; i < npts; ++i ) v[i] = p_spl->eval_DD( vec[i] );
    }
  }

  //!
  //! Evaluate all the splines at `zeta` and fill a `map of values`
  //! in a GenericContainer and `indep` as independent spline
  //!
  inline void SplineSet::eval2_DD( real_type const zeta, integer const indep, real_type & x, GenericContainer & gc )
    const
  {
    map_type & vals = gc.set_map();
    intersect( indep, zeta, x );
    for ( auto const & [fst, snd] : m_header_to_position ) vals[fst] = m_splines[snd]->eval_DD( x );
  }

  //!
  //! Evaluate all the splines at `zeta` values contained in vec
  //! and fill a `map of vector` in a GenericContainer and `indep`
  //! as independent spline
  //!
  inline void SplineSet::eval2_DD( vec_real_type const & zetas, integer const indep, GenericContainer & gc ) const
  {
    integer const npts{ static_cast<integer>( zetas.size() ) };
    map_type &    vals{ gc.set_map() };

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

  //!
  //! Evaluate all the splines at `zeta` and fill a map of values
  //! in a GenericContainer with keys in `columns` and `indep`
  //! as independent spline
  //!
  inline void SplineSet::eval2_DD(
    real_type const         zeta,
    integer const           indep,
    real_type &             x,
    vec_string_type const & columns,
    GenericContainer &      gc ) const
  {
    map_type & vals{ gc.set_map() };
    intersect( indep, zeta, x );
    for ( auto const & S : columns )
    {
      Spline const * p_spl{ get_spline( S ) };
      vals[S] = p_spl->eval_DD( x );
    }
  }

  //!
  //! Evaluate all the splines at `zeta` values contained in vec
  //! and fill a map of vector in a GenericContainer with keys
  //! in `columns` and `indep` as independent spline
  //!
  inline void SplineSet::eval2_DD(
    vec_real_type const &   zetas,
    integer const           indep,
    vec_string_type const & columns,
    GenericContainer &      gc ) const
  {
    integer const npts{ static_cast<integer>( zetas.size() ) };
    map_type &    vals{ gc.set_map() };

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
        v[i] = p_spl->eval_DD( x );
      }
    }
  }

  /*
  //   _   _     _         _       _           _            _   _
  //  | |_| |__ (_)_ __ __| |   __| | ___ _ __(_)_   ____ _| |_(_)_   _____
  //  | __| '_ \| | '__/ _` |  / _` |/ _ \ '__| \ \ / / _` | __| \ \ / / _ \
  //  | |_| | | | | | | (_| | | (_| |  __/ |  | |\ V / (_| | |_| |\ V /  __/
  //   \__|_| |_|_|_|  \__,_|  \__,_|\___|_|  |_| \_/ \__,_|\__|_| \_/ \___|
  */
  inline void SplineSet::eval_DDD( real_type const x, GenericContainer & gc ) const
  {
    map_type & vals{ gc.set_map() };
    for ( auto const & [fst, snd] : m_header_to_position ) vals[fst] = m_splines[snd]->eval_DDD( x );
  }

  //!
  //! Evaluate all the splines at `x` values contained in vec
  //! and fill a `map of vector` in a GenericContainer
  //!
  inline void SplineSet::eval_DDD( vec_real_type const & vec, GenericContainer & gc ) const
  {
    integer const npts{ static_cast<integer>( vec.size() ) };
    map_type &    vals{ gc.set_map() };
    for ( auto const & [fst, snd] : m_header_to_position )
    {
      vec_real_type & v{ vals[fst].set_vec_real( npts ) };
      Spline const *  p_spl{ m_splines[snd].get() };
      for ( integer i = 0; i < npts; ++i ) v[i] = p_spl->eval_DDD( vec[i] );
    }
  }

  //!
  //! Evaluate all the splines at `x` and fill a map of values
  //! in a GenericContainer with keys in `columns`
  //!
  inline void SplineSet::eval_DDD( real_type const x, vec_string_type const & columns, GenericContainer & gc ) const
  {
    map_type & vals{ gc.set_map() };
    for ( auto const & S : columns )
    {
      Spline const * p_spl{ get_spline( S ) };
      vals[S] = p_spl->eval_DDD( x );
    }
  }

  //!
  //! Evaluate all the splines at `x` values contained in vec and
  //! fill a map of vector in a GenericContainer with keys in `columns`
  //!
  inline void SplineSet::eval_DDD( vec_real_type const & vec, vec_string_type const & columns, GenericContainer & gc )
    const
  {
    integer const npts{ static_cast<integer>( vec.size() ) };
    map_type &    vals{ gc.set_map() };
    for ( auto const & S : columns )
    {
      Spline const *  p_spl{ get_spline( S ) };
      vec_real_type & v{ vals[S].set_vec_real( npts ) };
      for ( integer i = 0; i < npts; ++i ) v[i] = p_spl->eval_DDD( vec[i] );
    }
  }

  //!
  //! Evaluate all the splines at `zeta` and fill a `map of values`
  //! in a GenericContainer and `indep` as independent spline
  //!
  inline void SplineSet::eval2_DDD( real_type const zeta, integer const indep, real_type & x, GenericContainer & gc )
    const
  {
    map_type & vals{ gc.set_map() };
    intersect( indep, zeta, x );
    for ( auto const & [fst, snd] : m_header_to_position ) vals[fst] = m_splines[snd]->eval_DDD( x );
  }

  //!
  //! Evaluate all the splines at `zeta` values contained in vec
  //! and fill a `map of vector` in a GenericContainer and `indep`
  //! as independent spline
  //!
  inline void SplineSet::eval2_DDD( vec_real_type const & zetas, integer const indep, GenericContainer & gc ) const
  {
    integer const npts{ static_cast<integer>( zetas.size() ) };
    map_type &    vals{ gc.set_map() };

    // pre-allocation
    for ( auto const & [fst, snd] : m_header_to_position ) vals[fst].set_vec_real( npts );

    for ( integer i = 0; i < npts; ++i )
    {
      real_type x;
      intersect( indep, zetas[i], x );
      for ( auto const & [fst, snd] : m_header_to_position )
      {
        vec_real_type & v{ vals[fst].get_vec_real() };
        v[i] = m_splines[snd]->eval_DDD( x );
      }
    }
  }

  //!
  //! Evaluate all the splines at `zeta` and fill a map of values
  //! in a GenericContainer with keys in `columns` and `indep`
  //! as independent spline
  //!
  inline void SplineSet::eval2_DDD(
    real_type const         zeta,
    integer const           indep,
    real_type &             x,
    vec_string_type const & columns,
    GenericContainer &      gc ) const
  {
    map_type & vals{ gc.set_map() };
    intersect( indep, zeta, x );
    for ( auto const & S : columns )
    {
      Spline const * p_spl{ get_spline( S ) };
      vals[S] = p_spl->eval_DDD( x );
    }
  }

  //!
  //! Evaluate all the splines at `zeta` values contained in vec
  //! and fill a map of vector in a GenericContainer with keys
  //! in `columns` and `indep` as independent spline
  //!
  inline void SplineSet::eval2_DDD(
    vec_real_type const &   zetas,
    integer const           indep,
    vec_string_type const & columns,
    GenericContainer &      gc ) const
  {
    integer const npts{ static_cast<integer>( zetas.size() ) };
    map_type &    vals{ gc.set_map() };

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
        v[i] = p_spl->eval_DDD( x );
      }
    }
  }

}  // namespace Splines

// EOF: SplineSet.hxx
#endif
