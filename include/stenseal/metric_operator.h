/* -*- c-basic-offset:2; tab-width:2; indent-tabs-mode:nil -*-
-*
-*/

#ifndef _METRIC_OPERATOR_H
#define _METRIC_OPERATOR_H

#include <vector>

#include <deal.II/lac/vector.h>

#include "stenseal/stencil.h"
#include "stenseal/stencil_tensor.h"

namespace stenseal
{
  template <int width_interior,
            int width_boundary,
            int height_boundary>
  class MetricOperator
  {
  private:
    const StencilTensor2D<width_interior,width_interior> interior;
    const StencilTensor3D<height_boundary,width_boundary,width_boundary> boundary;

  public:
    constexpr MetricOperator(const StencilTensor2D<width_interior,width_interior> i,
                             const StencilTensor3D<height_boundary,width_boundary,width_boundary> b);

    template <typename MetricCoefficient>
    void apply(dealii::Vector<double> &dst,
               const dealii::Vector<double> &src,
               const MetricCoefficient &coeff,
               const unsigned int n) const;
  };

  //---------------------------------------------------------------------------
  // Implementations
  //---------------------------------------------------------------------------

  template <int width_interior,
            int width_boundary,
            int height_boundary>
  constexpr MetricOperator<width_interior,width_boundary,height_boundary>
  ::MetricOperator(const StencilTensor2D<width_interior,width_interior> i,
                   const StencilTensor3D<height_boundary,width_boundary,width_boundary> b)
    : interior(i), boundary(b)
  {}

  template <int width_interior,
            int width_boundary,
            int height_boundary>
  template <typename MetricCoefficient>
  void MetricOperator<width_interior,width_boundary,height_boundary>
  ::apply(dealii::Vector<double> &dst,
          const dealii::Vector<double> &src,
          const MetricCoefficient &coeff,
          const unsigned int n) const
  {

    for(int i = 0; i < height_boundary; ++i) {
      const auto b = boundary[i].apply_inner(coeff.template get_left_boundary_array<width_boundary>());
      dst[i] = b.apply(src,i);
    }

    for(int i = height_boundary; i < n-height_boundary; ++i) {
      const Stencil<width_interior> s = interior.apply_inner(coeff.template get_centered_array<width_interior>(i));
      dst[i] = s.apply(src,i);
    }

    for(int i = n-height_boundary; i < n; ++i) {
      const auto b = boundary[i-n+height_boundary].apply_inner_flip(coeff.template get_right_boundary_array<width_boundary>());
      dst[i] = b.apply_flip(src,i);
    }

  }

}

#endif /* _METRIC_OPERATOR_H */
