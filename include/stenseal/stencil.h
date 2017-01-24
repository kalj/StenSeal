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

  //=============================================================================
  // append to std::array
  //=============================================================================

  template <typename T, std::size_t N>
  constexpr std::array<T, N+1> append(const std::array<T, N> a, const T t);

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

  template <int width>
  struct Stencil
  {
    const std::array<double,width> weights;
    const std::array<int,width> offsets;

    constexpr  Stencil(const Sum<width> s)
      : weights(s.weights), offsets(s.offsets)
    {}

    inline double apply(const dealii::Vector<double> &u, int i) const
    {
      double t=0;
      for(int j = 0; j < width; ++j) {
        t+=u[i+offsets[j]]*weights[j];
      }
      return t;
    }

    inline double apply_fliplr(const dealii::Vector<double> &u, int i) const
    {
      double t=0;
      for(int j = width-1; j >= 0; --j) {
        t += u[i-offsets[j]]*weights[j];
      }
      return t;
    }

  };

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


}

#endif /* _STENCIL_H */
