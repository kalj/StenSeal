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
     * Constructor. Called explicitly with weights and offsets of matching
     * length `width`, which for performance reasons should always be
     * `constexpr`.
     */
    constexpr  Stencil(const std::array<double,width> w, const std::array<int,width> o);

    /**
     * Constructor. Only weights are specified; contiguous centered indexing is
     * assumed. For even width, the center is assumed to be the left one of the
     * two center points. For performance reasons the input should always be
     * `constexpr`.
     */
    constexpr  Stencil(const std::array<double,width> w);

    /**
     * Constructor. Called with an expression of matching length `width` as
     * input, which for performance reasons should always be `constexpr`.
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
    inline double apply_flip(const dealii::Vector<double> &u, int i) const;

    constexpr inline double apply(const std::array<double,width> &u) const;
    constexpr inline double apply_flip(const std::array<double,width> &u) const;

    std::array<double,width> get_weights() const;

    std::array<int,width> get_offsets() const;

  };

  //=============================================================================
  // Implementations
  //=============================================================================

  namespace internal
  {
    template <int m>
    struct gen_offsets;
  }


  template <int width>
  constexpr  Stencil<width>::Stencil(const std::array<double,width> w, const std::array<int,width> o)
    : weights(w), offsets(o)
  {}

  template <int width>
  constexpr  Stencil<width>::Stencil(const std::array<double,width> w)
    : weights(w), offsets(internal::gen_offsets<width>::value)
  {}

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
  inline double Stencil<width>::apply_flip(const dealii::Vector<double> &u, int i) const
  {
    double t=0;
    for(int j = width-1; j >= 0; --j) {
      t += u[i-offsets[j]]*weights[j];
    }
    return t;
  }

  // helper for compile-time computing of dot product of two std::arrays
  template<bool flip,typename T, std::size_t N>
  constexpr T compile_time_dot_product(std::array<T, N> const &A, std::array<T, N> const &B, int const i = 0) {
    return (i < N)? A[flip?N-1-i:i]*B[i] + compile_time_dot_product<flip>(A, B, i + 1) : T(0);
  }

  template <int width>
  constexpr inline double Stencil<width>::apply(const std::array<double,width> &u) const
  {
    return compile_time_dot_product<false>(weights,u);
  }

  template <int width>
  constexpr inline double Stencil<width>::apply_flip(const std::array<double,width> &u) const
  {
    return compile_time_dot_product<true>(weights,u);
  }

  template <int width>
  std::array<double,width> Stencil<width>::get_weights() const
  {
    return weights;
  }

  template <int width>
  std::array<int,width> Stencil<width>::get_offsets() const
  {
    return offsets;
  }

  //=============================================================================
  // append/prepend to std::array
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

  template <typename T, std::size_t N, std::size_t... I>
  constexpr std::array<T, N + 1> prepend_aux(T t,const std::array<T, N> a,
                                             const std::index_sequence<I...>)
  {
    return std::array<T, N + 1>{t, a[I]... };
  }

  template <typename T, std::size_t N>
  constexpr std::array<T, N+1> prepend(const T t, const std::array<T, N> a)
  {
    return prepend_aux(t,a, std::make_index_sequence<N>());
  }

  namespace internal
  {
    // generate contiguous list of offsets centered around 0
    template <>
    struct gen_offsets<1>
    {
      static constexpr std::array<int,1> value{0};
    };

    template <>
    struct gen_offsets<2>
    {
      static constexpr std::array<int,2> value{0,1};
    };

    template <int m>
    struct gen_offsets
    {
      static constexpr std::array<int,m> value = prepend((m%2==0)-m/2, append(gen_offsets<m-2>::value, m/2));
    };

    template <int m>
    constexpr std::array<int,m> gen_offsets<m>::value;

    constexpr std::array<int,1> gen_offsets<1>::value;
    constexpr std::array<int,2> gen_offsets<2>::value;
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
