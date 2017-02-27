/* -*- c-basic-offset:2; tab-width:2; indent-tabs-mode:nil -*-
 *
 */

#ifndef _OPERATOR_H
#define _OPERATOR_H

#include <vector>

#include <deal.II/lac/vector.h>

#include "stenseal/stencil.h"
#include "stenseal/stencil_tensor.h"

//=============================================================================
// this is an operator in 1D
//=============================================================================
namespace stenseal
{
  /**
   * Class representing a 1D SBP operator, width an interior stencil with width
   * `width_interior`, and distinct boundary stencil block wich is the same for
   * both boundaries with `height_boundary` rows of width `width_boundary`.
   *
   */
  template <int width_interior,
            int width_boundary,
            int height_boundary>
  class SymmetricSBP
  {
  private:
    const StencilTensor2D<height_boundary,width_boundary> boundary;
    const Stencil<width_interior> interior;
  public:
    const static int height_b = height_boundary;
    const static int width_b = width_boundary_r;
    const static int width_i = width_interior;

    /**
     * Constructor. Takes the interior `Stencil` `i`, and the `StencilTensor2D`
     * `b` for the boundary .
     */
    constexpr Operator(const Stencil<width_interior> i,
                       const StencilTensor2D<height_boundary,width_boundary> b);

    /**
     * Apply this `Operator` to the vector `src` and write the result into
     * `dst`, both of which have length `n`.
     */
    void apply(dealii::Vector<double> &dst,
               const dealii::Vector<double> &src,
               const unsigned int n) const;

    Stencil<width_interior> get_interior_stencil() const;

    StencilTensor2D<height_boundary_l,width_boundary_l> get_boundary_stencil() const;

  };

  //---------------------------------------------------------------------------
  // Implementations
  //---------------------------------------------------------------------------

  template <int width_interior,
            int width_boundary,
            int height_boundary>
  constexpr Operator<width_interior,width_boundary,height_boundary>
  ::Operator(const Stencil<width_interior> i,
             const StencilTensor2D<height_boundary,width_boundary> b)
    :interior(i), lboundary(l), rboundary(r)
  {}

  template <int width_interior,
            int width_boundary,
            int height_boundary>
  void Operator<width_interior,width_boundary,height_boundary>
  ::apply(dealii::Vector<double> &dst,
          const dealii::Vector<double> &src,
          const unsigned int n) const
  {
    for(int i = 0; i < height_boundary; ++i) {
      dst[i] = boundary[i].apply(src,i);
    }

    for(int i = height_boundary; i < n-height_boundary; ++i) {
      dst[i] = interior.apply(src,i);
    }

    for(int i = n-height_boundary; i < n; ++i) {
      dst[i] = boundary[i-n+height_boundary_r].apply_flip(src,i);
    }
  }

  template <int width_interior,
            int width_boundary_l,
            int height_boundary_l>
  Stencil<width_interior> Operator<width_interior,width_boundary,height_boundary>
  ::get_interior_stencil() const{
    return interior;
  }

  template <int width_interior,
            int width_boundary,
            int height_boundary>
  StencilTensor2D<height_boundary,width_boundary>  Operator<width_interior,width_boundary,height_boundary>
  ::get_left_boundary_stencil() const{
    return boundary;
  }

}

#endif /* _OPERATOR_H */