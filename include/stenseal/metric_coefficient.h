/* -*- c-basic-offset:2; tab-width:2; indent-tabs-mode:nil -*-
 *
 */

#ifndef _METRIC_COEFFICIENT_H
#define _METRIC_COEFFICIENT_H

#include <array>

#include "stenseal/utils.h"
#include "stenseal/geometry.h"

namespace stenseal
{
  // metric coefficients

  // forward declaration
  template <int dim, typename Geometry>
  struct MetricCoefficient;


  //---------------------------------------------------------------------------
  // coefficient classes
  //---------------------------------------------------------------------------

  /**
   * A coefficient equal to 1.
   */

  class OneCoefficient
  {
  public:
    constexpr OneCoefficient();

    constexpr double get(int) const;

    template <int n>
      constexpr std::array<double, n> get_centered_array(int) const;

    template <int n>
      constexpr std::array<double, n> get_left_boundary_array() const;

    template <int n>
      constexpr std::array<double, n> get_right_boundary_array() const;
  };

  /**
   * A general coefficient
   */

  class VectorCoefficient
  {
  private:
    dealii::Vector<double> vec;

  public:
    VectorCoefficient();

    void init(const dealii::Vector<double> &v);

    double get(int i) const;

    template <int n>
      std::array<double, n> get_centered_array(int i) const;

    template <int n>
      std::array<double, n> get_left_boundary_array() const;

    template <int n>
      std::array<double, n> get_right_boundary_array() const;
  };


  /**
   * MetricCoefficient for cartesian geometry
   */

  template <int dim>
  struct MetricCoefficient<dim,CartesianGeometry<dim>>
  {
  private:
    OneCoefficient coeff;

  public:
    template <typename dummy>
      constexpr MetricCoefficient(const CartesianGeometry<dim>&, const dummy&);

    constexpr const OneCoefficient& inverse_jacobian() const;

  };


  /**
   * MetricCoefficient for general geometry
   */

  template <int dim>
  struct MetricCoefficient<dim,GeneralGeometry<dim>>
  {
  private:
    VectorCoefficient coeff;
    unsigned int n_points;

  public:
    template <typename DmOp>
    MetricCoefficient(const GeneralGeometry<dim> &g, const DmOp &op);

    const VectorCoefficient& inverse_jacobian() const;

  };


  //---------------------------------------------------------------------------
  // implementations
  //---------------------------------------------------------------------------

  // unit coefficient
  constexpr OneCoefficient::OneCoefficient()
  {}

  constexpr double OneCoefficient::get(int) const
  {
    return 1.0;
  }

  template <int n>
  constexpr std::array<double, n> OneCoefficient::get_centered_array(int) const
  {
    return internal::repeat_value<n>(1.0);
  }

  template <int n>
  constexpr std::array<double, n> OneCoefficient::get_left_boundary_array() const
  {
    return internal::repeat_value<n>(1.0);
  }

  template <int n>
  constexpr std::array<double, n> OneCoefficient::get_right_boundary_array() const
  {
    return internal::repeat_value<n>(1.0);
  }


  // general coefficient
  VectorCoefficient::VectorCoefficient()
  {}

  void VectorCoefficient::init(const dealii::Vector<double> &v)
  {
    vec = v;
  }

  double VectorCoefficient::get(int i) const
  {
    return vec[i];
  }

  template <int n>
  std::array<double, n> VectorCoefficient::get_centered_array(int i) const
  {
    std::array<double, n> res;
    const std::array<int, n> &offsets = internal::gen_offsets<n>::value;
    for(int j = 0; j < n; ++j)
      res[j] = vec[i+offsets[j]];
    return res;
  }

  template <int n>
  std::array<double, n> VectorCoefficient::get_left_boundary_array() const
  {
    std::array<double, n> res;
    for(int j = 0; j < n; ++j)
      res[j] = vec[j];
    return res;
  }

  template <int n>
  std::array<double, n> VectorCoefficient::get_right_boundary_array() const
  {
    const unsigned int n_points = vec.size();

    std::array<double, n> res;
    for(int j = 0; j < n; ++j)
      res[j] = vec[n_points-n+j];
    return res;
  }



  // for cartesian geometry
  template <int dim>
  template <typename dummy>
  constexpr MetricCoefficient<dim,CartesianGeometry<dim>>
    ::MetricCoefficient(const CartesianGeometry<dim>&,
                        const dummy&)
  {}

  template <int dim>
  constexpr const OneCoefficient& MetricCoefficient<dim,CartesianGeometry<dim>>::inverse_jacobian() const
  {
    return coeff;
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

      dealii::Vector<double> c(n_points);
      op.apply(c,xvals,n_points);

      // divide by step size, and invert
      for(int i=0; i<n_points; ++i) {
        c[i] = g.get_mapped_h(0)/c[i];
      }

      coeff.init(c);
    }
    else {
      AssertThrow(false,dealii::ExcNotImplemented());
    }
  }

  template <int dim>
  const VectorCoefficient& MetricCoefficient<dim,GeneralGeometry<dim>>::inverse_jacobian() const
  {
    return coeff;
  }

}

#endif /* _METRIC_COEFFICIENT_H */
