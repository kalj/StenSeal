/* -*- c-basic-offset:2; tab-width:2; indent-tabs-mode:nil -*-
 *
 *
 */

#ifndef _GEOMETRY_H
#define _GEOMETRY_H

#include <deal.II/base/point.h>

namespace stenseal
{

  namespace internal
  {
    template <std::size_t n, typename T>
    constexpr std::array<T,n> repeat_value(const T val);
  }

  template <int dim>
  class CartesianGeometry {
  public:
    struct MetricCoefficient {
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
  private:
    double h[dim];
    unsigned int n_nodes[dim];
    unsigned int n_nodes_total;
    dealii::Point<dim> lower_left;
    MetricCoefficient coefficient;

  public:
    // FIXME: add default parameters which are [0,0,0,0....] and [1,1,1,1,...]
    CartesianGeometry(unsigned int n_nodes[dim],
                      const dealii::Point<dim> lower_left,
                      const dealii::Point<dim> upper_right);

    inline unsigned get_n_nodes(int d) const
    {
      return n_nodes[d];
    }

    inline unsigned get_n_nodes_total() const
    {
      return n_nodes_total;
    }

    inline double get_h(int d) const
    {
      return h[d];
    }


    inline constexpr const MetricCoefficient& get_metric_coefficient() const
    {
      return coefficient;
    }

    inline double get_lower_left(int d) const
    {
      return lower_left(d);
    }

  };

  /*
  struct TransfiniteInterpolationGeometry {
    double h;
    inline double get_c(int i) const
    {
      return h;
    }
  };

  struct GeneralGeometry {
    double *c;
    inline double get_c(int i) const
    {
      return c[i];
    }
    };
  */

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

#endif /* _GEOMETRY_H */
