/* -*- c-basic-offset:2; tab-width:2; indent-tabs-mode:nil -*-
 *
 *
 */

#ifndef _STENCIL_TENSOR_H
#define _STENCIL_TENSOR_H

#include "stenseal/stencil.h"

namespace stenseal
{
  template <int n, typename Inner>
  class StencilTensor;

  template <int n, int m>
  using StencilTensor2D = StencilTensor<n,Stencil<m>>;

  template <int n, int m, int l>
  using StencilTensor3D = StencilTensor<n,StencilTensor<m,Stencil<l>>>;

  template <int n, int m, int l, int k>
  using StencilTensor4D = StencilTensor<n,StencilTensor<m,StencilTensor<l,Stencil<k>>>>;

  template <int n, typename Inner>
  class StencilTensor
  {
  private:
    const std::array<Inner,n> stencils;

  public:
    const static int inner_dim = Inner::inner_dim;
    typedef StencilTensor<n,typename Inner::result_type> result_type;

    template <typename... Ss>
    constexpr  StencilTensor(const Ss... s);

    /**
     * Access a specific row of the block.
     */
    inline constexpr const Inner& operator[](int row) const;

    inline constexpr result_type apply_inner(const std::array<double,inner_dim> &src) const;

    inline constexpr result_type apply_inner_flip(const std::array<double,inner_dim> &src) const;
  };

  template <int n, int m>
  class StencilTensor<n,Stencil<m>>
  {
  private:
    const std::array<Stencil<m>,n> stencils;

  public:
    typedef Stencil<n> result_type;
    const static int inner_dim = m;

    template <typename... Ss>
    constexpr StencilTensor(const Ss... s);

    /**
     * Access a specific row of the block.
     */
    constexpr inline const Stencil<m>& operator[](int row) const;

    constexpr inline const result_type apply_inner(const std::array<double,m> &src) const;

    constexpr inline const result_type apply_inner_flip(const std::array<double,m> &src) const;

  };


  //---------------------------------------------------------------------------
  // Implementations
  //---------------------------------------------------------------------------


  // helper for compile-time evaluation of an array of stencils on a single input
  template <bool flip,std::size_t n, typename S, typename T, std::size_t... I>
  constexpr auto for_each_apply_impl( const std::array<S,n> &stencils,
                                      const T input, const std::index_sequence<I...>)
    -> std::array<decltype(stencils[0].apply(input)),n>
  {
    return { flip?stencils[n-1-I].apply_flip(input):stencils[I].apply(input)... };
  }

  template <std::size_t n, typename S, typename T>
  constexpr auto for_each_apply(const std::array<S,n> &stencils, const T input)
    -> std::array<decltype(stencils[0].apply(input)),n>
  {
    return for_each_apply_impl<false,n>(stencils, input, std::make_index_sequence<n>());
  }

  template <std::size_t n, typename S, typename T>
  constexpr auto for_each_apply_flip(const std::array<S,n> &stencils, const T input)
    -> std::array<decltype(stencils[0].apply(input)),n>
  {
    return for_each_apply_impl<true,n>(stencils, input, std::make_index_sequence<n>());
  }

  // general version
  template <int n,  typename Inner>
  template <typename... Ss>
  constexpr StencilTensor<n,Inner>::StencilTensor(const Ss... s)
    : stencils{s...}
  {
    static_assert(sizeof...(Ss) == n, "wrong number of arguments");
  }

  template <int n, typename Inner>
  inline constexpr const Inner& StencilTensor<n,Inner>::operator[](int row) const
  {
    return stencils[row];
  }

  template <int n, typename Inner>
  inline constexpr typename StencilTensor<n,Inner>::result_type
  StencilTensor<n,Inner>::apply_inner(const std::array<double,inner_dim> &src) const
  {
    return result_type{for_each_apply(stencils,src)};
  }

  template <int n, typename Inner>
  inline constexpr typename StencilTensor<n,Inner>::result_type
  StencilTensor<n,Inner>::apply_inner_flip(const std::array<double,inner_dim> &src) const
  {
    return result_type{for_each_apply_flip(stencils,src)};
  }


  // base case explicit specialization
  template <int n, int m>
  template <typename... Ss>
  constexpr StencilTensor<n,Stencil<m>>::StencilTensor(const Ss... s)
    : stencils{s...}
  {
    static_assert(sizeof...(Ss) == n, "wrong number of arguments");
  }

  template <int n, int m>
  inline constexpr const Stencil<m>& StencilTensor<n,Stencil<m>>::operator[](int row) const
  {
    return stencils[row];
  }

  template <int n, int m>
  inline constexpr const Stencil<n> StencilTensor<n,Stencil<m>>::apply_inner(const std::array<double,m> &src) const
  {
    return Stencil<n>{for_each_apply(stencils,src)};
  }

  template <int n, int m>
  inline constexpr const Stencil<n> StencilTensor<n,Stencil<m>>::apply_inner_flip(const std::array<double,m> &src) const
  {
    return Stencil<n>{for_each_apply_flip(stencils,src)};
  }

}

#endif /* _STENCIL_TENSOR_H */
