/* -*- c-basic-offset:2; tab-width:2; indent-tabs-mode:nil -*-
-*
-*/

#ifndef _OPERATOR_H
#define _OPERATOR_H

#include <vector>

#include <deal.II/lac/vector.h>

#include "stenseal/stencil.h"
#include "stenseal/stencil_array.h"

//=============================================================================
// this is an operator in 1D
//=============================================================================
namespace stenseal
{
  /**
   * Class representing a 1D SBP operator, width an interior stencil with width
   * `width_interior`, and distinct boundary stencil blocks on the left and
   * right boundary, both with `height_boundary` rows of width `width_boundary`.
   *
   */
  template <int width_interior,
            int width_boundary_l,
            int height_boundary_l,
            int width_boundary_r,
            int height_boundary_r>
  class Operator
  {
  private:
    const StencilArray<width_boundary_l,height_boundary_l> lboundary;
    const StencilArray<width_boundary_r,height_boundary_r> rboundary;
    const Stencil<width_interior> interior;
  public:
    const static int height_r = width_boundary_l;
    const static int height_l = height_boundary_l;
    const static int width_i = width_interior;

    /**
     * Constructor. Takes the interior `Stencil` `i`, and the `StencilArray`s
     * `l` and `r` for the left and right boundaries respectively.
     */
    constexpr Operator(const Stencil<width_interior> i,
                       const StencilArray<width_boundary_l,height_boundary_l> l,
                       const StencilArray<width_boundary_r,height_boundary_r> r);

    /**
     * Apply this `Operator` to the vector `src` and write the result into
     * `dst`, both of which have length `n`.
     */
    void apply(dealii::Vector<double> &dst,
               const dealii::Vector<double> &src,
               const unsigned int n) const;
    /**
     * Apply the associated "D+" operator of this object, i.e. the operator
     * defined by reflecting and negating the corresponding matrix.
     */
    void apply_dp(dealii::Vector<double> &dst,
      const dealii::Vector<double> &src,
      const unsigned int n) const;
  };


  //---------------------------------------------------------------------------
  // Implementations
  //---------------------------------------------------------------------------

  template <int width_interior,
            int width_boundary_l,
            int height_boundary_l,
            int width_boundary_r,
            int height_boundary_r>
  constexpr Operator<width_interior,width_boundary_l,height_boundary_l,width_boundary_r,height_boundary_r>
  ::Operator(const Stencil<width_interior> i,
             const StencilArray<width_boundary_l,height_boundary_l> l,
             const StencilArray<width_boundary_r,height_boundary_r> r)
      :interior(i), lboundary(l), rboundary(r)
  {}

  template <int width_interior,
            int width_boundary_l,
            int height_boundary_l,
            int width_boundary_r,
            int height_boundary_r>
  void Operator<width_interior,width_boundary_l,height_boundary_l, width_boundary_r, height_boundary_r>
  ::apply(dealii::Vector<double> &dst,
        const dealii::Vector<double> &src,
        const unsigned int n) const
  {
    for(int i = 0; i < height_boundary_l; ++i) {
      dst[i] = lboundary[i].apply(src,i);
    }

    for(int i = height_boundary_l; i < n-height_boundary_r; ++i) {
      dst[i] = interior.apply(src,i);
    }

    for(int i = n-height_boundary_r; i < n; ++i) {
      dst[i] = rboundary[i-n+height_boundary_r].apply(src,i);
    }
  }

  template <int width_interior,
            int width_boundary_l,
            int height_boundary_l,
            int width_boundary_r,
            int height_boundary_r>
  void Operator<width_interior,width_boundary_l,height_boundary_l, width_boundary_r, height_boundary_r>
  ::apply_dp(dealii::Vector<double> &dst,
           const dealii::Vector<double> &src,
           const unsigned int n) const
  {
    for(int i = 0; i < height_boundary_r; ++i) {
      dst[i] = -rboundary[height_boundary_r-1-i].apply_fliplr(src,i);
    }

    for(int i = height_boundary_r; i < n-height_boundary_l; ++i) {
      dst[i] = -interior.apply_fliplr(src,i);
    }

    for(int i = n-height_boundary_l; i < n; ++i) {
      dst[i] = -lboundary[n-1-i].apply_fliplr(src,i);
    }
  }

}

#endif /* _OPERATOR_H */
