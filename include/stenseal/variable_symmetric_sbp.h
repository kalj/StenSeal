/* -*- c-basic-offset:2; tab-width:2; indent-tabs-mode:nil -*-
-*
-*/

#ifndef _VARIABLE_SYMMETRIC_SBP_H
#define _VARIABLE_SYMMETRIC_SBP_H

#include <vector>

#include <deal.II/lac/vector.h>

#include "stenseal/stencil.h"
#include "stenseal/stencil_tensor.h"

namespace stenseal
{
  template <int width_interior,
            int width_boundary,
            int height_boundary>
  class VariableSymmetricSBP
  {
  private:
    const StencilTensor2D<width_interior,width_interior> interior;
    const StencilTensor3D<height_boundary,width_boundary,width_boundary> boundary;

  public:
    const static int width_b   = width_boundary;
    const static int height_b  = height_boundary;
    const static int width_i   = width_interior;

    constexpr VariableSymmetricSBP(const StencilTensor2D<width_interior,width_interior> i,
                             const StencilTensor3D<height_boundary,width_boundary,width_boundary> b);

    template <typename MetricCoefficient>
    void apply(dealii::Vector<double> &dst,
               const dealii::Vector<double> &src,
               const MetricCoefficient &coeff,
               const unsigned int n) const;

    const StencilTensor3D<height_boundary,width_boundary,width_boundary>& get_boundary() const;

    const StencilTensor2D<width_interior,width_interior>& get_interior() const;
  };

  //---------------------------------------------------------------------------
  // Implementations
  //---------------------------------------------------------------------------

  template <int width_interior,
            int width_boundary,
            int height_boundary>
  constexpr VariableSymmetricSBP<width_interior,width_boundary,height_boundary>
  ::VariableSymmetricSBP(const StencilTensor2D<width_interior,width_interior> i,
                   const StencilTensor3D<height_boundary,width_boundary,width_boundary> b)
    : interior(i), boundary(b)
  {}

  template <int width_interior,
            int width_boundary,
            int height_boundary>
  template <typename MetricCoefficient>
  void VariableSymmetricSBP<width_interior,width_boundary,height_boundary>
  ::apply(dealii::Vector<double> &dst,
          const dealii::Vector<double> &src,
          const MetricCoefficient &coeff,
          const unsigned int n) const
  {

    for(int i = 0; i < height_boundary; ++i) {
      const auto b = boundary[i].apply_inner(coeff.template get_left_boundary_array<width_boundary>());
      // FIXME: temporary fix to apply generated stencil in a non-centered fashion
      dst[i] = b.apply(src,(width_boundary-1)/2);
      // dst[i] = b.apply(src,i);
    }

    for(int i = height_boundary; i < n-height_boundary; ++i) {
      const Stencil<width_interior> s = interior.apply_inner(coeff.template get_centered_array<width_interior>(i));
      dst[i] = s.apply(src,i);
    }

    for(int i = n-height_boundary; i < n; ++i) {
      const auto b = boundary[n-1-i].apply_inner_flip(coeff.template get_right_boundary_array<width_boundary>());

      // FIXME: temporary fix to apply generated stencil in a non-centered fashion
      dst[i] = b.apply(src,n-width_boundary+(width_boundary-1)/2);
      // dst[i] = b.apply(src,i);
    }

  }

  template <int width_interior,
            int width_boundary,
            int height_boundary>
  const StencilTensor3D<height_boundary,width_boundary,width_boundary>& VariableSymmetricSBP<width_interior,width_boundary,height_boundary>
  ::get_boundary() const
  {
    return boundary;
  }

  template <int width_interior,
            int width_boundary,
            int height_boundary>
  const StencilTensor2D<width_interior,width_interior>& VariableSymmetricSBP<width_interior,width_boundary,height_boundary>
  ::get_interior() const
  {
    return interior;
  }

  template <int width_interior,
            int width_boundary,
            int height_boundary>
  const int VariableSymmetricSBP<width_interior,width_boundary,height_boundary>::height_b;

  template <int width_interior,
            int width_boundary,
            int height_boundary>
  const int VariableSymmetricSBP<width_interior,width_boundary,height_boundary>::width_b;

  template <int width_interior,
            int width_boundary,
            int height_boundary>
  const int VariableSymmetricSBP<width_interior,width_boundary,height_boundary>::width_i;

}

#endif /* _VARIABLE_SYMMETRIC_SBP_H */
