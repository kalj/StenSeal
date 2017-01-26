/* -*- c-basic-offset:2; tab-width:2; indent-tabs-mode:nil -*-
 *
 *
 */

#ifndef _STENCIL_H
#define _STENCIL_H

#include <array>
#include <utility>

#include <deal.II/lac/vector.h>

namespace stenseal
{
  // forward declaration
  template <int n>
  struct Sum;

  /**
   * A class representing a single row in an operator, i.e. a weighted
   * combination of a number elements in the input vector. The template
   * parameter `width` is the number of non-zeros in the row. It is initialized
   * using expression templates based on a variable of the `Symbol` type, e.g.
   *
   *     const stenseal::Symbol u;
   *     constexpr stenseal::Stencil<2> interior((-1.0)*u[-1] + 1.0*u[0]);
   *
   * Here, the indexing into the symbol `u` should be understood as relative
   * offsets from the diagonal.
   *
   * For performance reasons, all constituents of the expression, i.e. indices
   * and weights, should be literals or other `constexpr` entities, which can be
   * enforced by declaring the `Stencil` itself to be `constexpr` as above. If
   * this is done, very aggressive inlining optimizations can be performed by
   * the compiler such that no indices or weights must be read from memory
   * during the stencil application.
   */

  template <int width>
  class Stencil
  {
  private:
    const std::array<double,width> weights;
    const std::array<int,width> offsets;
  public:

    /**
     * Constructor. Called with an expression of matching length `width` as
     * input, which for performance reasons should always be a `constexpr`.
     */
    constexpr  Stencil(const Sum<width> s);

    /**
     * Function for applying this stencil to a vector `u` at position
     * `i`. Returns the resulting value.
     */
    inline double apply(const dealii::Vector<double> &u, int i) const;

    /**
     * Function for applying the flipped/reversed stencil.
     */
    inline double apply_fliplr(const dealii::Vector<double> &u, int i) const;
  };

  //=============================================================================
  // Implementations
  //=============================================================================

  template <int width>
  constexpr  Stencil<width>::Stencil(const Sum<width> s)
    : weights(s.weights), offsets(s.offsets)
  {}

  template <int width>
  inline double Stencil<width>::apply(const dealii::Vector<double> &u, int i) const
  {
    double t=0;
    for(int j = 0; j < width; ++j) {
      t+=u[i+offsets[j]]*weights[j];
    }
    return t;
  }

  template <int width>
  inline double Stencil<width>::apply_fliplr(const dealii::Vector<double> &u, int i) const
  {
    double t=0;
    for(int j = width-1; j >= 0; --j) {
      t += u[i-offsets[j]]*weights[j];
    }
    return t;
  }

  //=============================================================================
  // append to std::array
  //=============================================================================

  template <typename T, std::size_t N, std::size_t... I>
  constexpr std::array<T, N + 1> append_aux(const std::array<T, N> a, T t,
                                            const std::index_sequence<I...>)
  {
    return std::array<T, N + 1>{ a[I]..., t };
  }

  template <typename T, std::size_t N>
  constexpr std::array<T, N+1> append(const std::array<T, N> a, const T t)
  {
    return append_aux(a, t, std::make_index_sequence<N>());
  }


  //=============================================================================
  // infrastructure for expression template
  //=============================================================================

  struct Accessor {
    const int offset;
  };

  struct Symbol
  {
    constexpr Accessor operator[](int i) const
    {
      return Accessor{i};
    }
  };

  struct WeightedTerm
  {
    const double weight;
    const int offset;

    constexpr WeightedTerm(double a, Accessor u)
      : weight(a), offset(u.offset)
    {}
  };

  constexpr WeightedTerm operator*(double a, Accessor u)
  {
    return WeightedTerm(a,u);
  }


  template <int n>
  struct Sum
  {
    const std::array<double,n> weights;
    const std::array<int,n> offsets;

    constexpr Sum(WeightedTerm x, WeightedTerm y)
      : weights({x.weight,y.weight}), offsets({x.offset,y.offset})
    {}

    constexpr Sum(Sum<n-1> x, WeightedTerm y)
      : weights(append(x.weights,y.weight)), offsets(append(x.offsets,y.offset))
    {}
  };

  constexpr Sum<2> operator+(WeightedTerm x, WeightedTerm y)
  {
    return Sum<2>(x,y);
  }

  template <int n>
  constexpr Sum<n+1> operator+(Sum<n> x, WeightedTerm y)
  {
    return Sum<n+1>(x,y);
  }

}

#endif /* _STENCIL_H */
