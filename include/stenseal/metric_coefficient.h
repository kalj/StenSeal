/* -*- c-basic-offset:2; tab-width:2; indent-tabs-mode:nil -*-
 *
 */

#ifndef _METRIC_COEFFICIENT_H
#define _METRIC_COEFFICIENT_H

#include <array>

#include "stenseal/geometry.h"

namespace stenseal
{
  namespace internal
  {
    template <std::size_t n, typename T>
    constexpr std::array<T,n> repeat_value(const T val);
  }

  template <int dim, typename Geometry>
  struct MetricCoefficient;

  template <int dim>
  struct MetricCoefficient<dim,CartesianGeometry<dim>>
  {
    template <typename dummy>
    constexpr MetricCoefficient(const CartesianGeometry<dim>&, const dummy&) {}

    constexpr double get(int) const
    {
      return 1.0;
    }

    template <int n>
      constexpr std::array<double, n> get_centered_array(int) const
    {
      return internal::repeat_value<n>(1.0);
    }

    template <int n>
      constexpr std::array<double, n> get_left_boundary_array() const
    {
      return internal::repeat_value<n>(1.0);
    }

    template <int n>
      constexpr std::array<double, n> get_right_boundary_array() const
    {
      return internal::repeat_value<n>(1.0);
    }
  };



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
  }

}

#endif /* _METRIC_COEFFICIENT_H */
