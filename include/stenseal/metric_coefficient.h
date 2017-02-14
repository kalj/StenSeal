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
      constexpr MetricCoefficient(const CartesianGeometry<dim>&, const dummy&);

    constexpr double get(int) const;

    template <int n>
      constexpr std::array<double, n> get_centered_array(int) const;

    template <int n>
      constexpr std::array<double, n> get_left_boundary_array() const;

    template <int n>
      constexpr std::array<double, n> get_right_boundary_array() const;
  };


  template <int dim>
  struct MetricCoefficient<dim,GeneralGeometry<dim>>
  {
  private:
    dealii::Vector<double> c;
    unsigned int n_points;

  public:
    template <typename DmOp>
    MetricCoefficient(const GeneralGeometry<dim> &g, const DmOp &op);

    double get(int i) const;

    template <int n>
      std::array<double, n> get_centered_array(int i) const;

    template <int n>
      std::array<double, n> get_left_boundary_array() const;

    template <int n>
      std::array<double, n> get_right_boundary_array() const;

  };


  //---------------------------------------------------------------------------
  // implementations
  //---------------------------------------------------------------------------

  // for cartesian geometry
  template <int dim>
  template <typename dummy>
  constexpr MetricCoefficient<dim,CartesianGeometry<dim>>
    ::MetricCoefficient(const CartesianGeometry<dim>&,
                        const dummy&)
  {}

  template <int dim>
  constexpr double MetricCoefficient<dim,CartesianGeometry<dim>>::get(int) const
  {
    return 1.0;
  }

  template <int dim>
  template <int n>
  constexpr std::array<double, n> MetricCoefficient<dim,CartesianGeometry<dim>>
    ::get_centered_array(int) const
  {
    return internal::repeat_value<n>(1.0);
  }

  template <int dim>
  template <int n>
  constexpr std::array<double, n> MetricCoefficient<dim,CartesianGeometry<dim>>
    ::get_left_boundary_array() const
  {
    return internal::repeat_value<n>(1.0);
  }

  template <int dim>
  template <int n>
  constexpr std::array<double, n> MetricCoefficient<dim,CartesianGeometry<dim>>
    ::get_right_boundary_array() const
  {
    return internal::repeat_value<n>(1.0);
  }


  // for general geometry
  template <int dim>
  template <typename DmOp>
  MetricCoefficient<dim,GeneralGeometry<dim>>::MetricCoefficient(const GeneralGeometry<dim> &g, const DmOp &op)
  {
    const std::vector<dealii::Point<dim>> &pts = g.get_node_points();
    n_points = pts.size();

    if(dim==1) {
      dealii::Vector<double> xvals(n_points);

      for(int i=0; i<n_points; ++i)
        xvals[i] = pts[i](0);


      c.reinit(n_points);
      op.apply(c,xvals,n_points);
      c /= g.get_mapped_h(0);

    }
    else {
      AssertThrow(false,dealii::ExcNotImplemented());
    }
  }


  template <int dim>
  double MetricCoefficient<dim,GeneralGeometry<dim>>
    ::get(int i) const
  {
    return c[i];
  }

  template <int dim>
  template <int n>
  std::array<double, n> MetricCoefficient<dim,GeneralGeometry<dim>>
    ::get_centered_array(int i) const
  {
    std::array<double, n> res;
    const std::array<int, n> &offsets = internal::gen_offsets<n>::value;
    for(int j = 0; j < n; ++j)
      res[j] = c[i+offsets[j]];
    return res;
  }

  template <int dim>
  template <int n>
  std::array<double, n> MetricCoefficient<dim,GeneralGeometry<dim>>
    ::get_left_boundary_array() const
  {
    std::array<double, n> res;
    for(int j = 0; j < n; ++j)
      res[j] = c[j];
    return res;
  }

  template <int dim>
  template <int n>
  std::array<double, n> MetricCoefficient<dim,GeneralGeometry<dim>>
    ::get_right_boundary_array() const
  {
    std::array<double, n> res;
    for(int j = 0; j < n; ++j)
      res[j] = c[n_points-n+j];
    return res;
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
  }

}

#endif /* _METRIC_COEFFICIENT_H */
