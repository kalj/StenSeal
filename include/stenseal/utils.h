/* -*- c-basic-offset:2; tab-width:2; indent-tabs-mode:nil -*-
 *
 */

#ifndef _UTILS_H
#define _UTILS_H

#include <array>
#include <utility>

namespace stenseal
{
  namespace internal
  {
    template <std::size_t n, typename T>
    constexpr std::array<T,n> repeat_value(const T val);

    template <int m>
    struct gen_offsets;
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

    template <std::size_t n, typename T, std::size_t... I>
    constexpr std::array<T,n> repeat_value_impl(const T val, const std::index_sequence<I...>)
    {
      return { (I,val)...};
    }

    template <std::size_t n, typename T>
    constexpr std::array<T,n> repeat_value(const T val)
    {
      return repeat_value_impl<n>(val,std::make_index_sequence<n>());
    }

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



}

#endif /* _UTILS_H */
