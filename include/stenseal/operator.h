/* -*- c-basic-offset:2; tab-width:2; indent-tabs-mode:nil -*-
-*
-*/

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
    const StencilTensor2D<height_boundary_l,width_boundary_l> lboundary;
    const StencilTensor2D<height_boundary_r,width_boundary_r> rboundary;
    const Stencil<width_interior> interior;
  public:
    const static int height_r = height_boundary_r;
    const static int height_l = height_boundary_l;
    const static int width_r = width_boundary_r;
    const static int width_l = width_boundary_l;
    const static int width_i = width_interior;

    /**
     * Constructor. Takes the interior `Stencil` `i`, and the `StencilTensor2D`s
     * `l` and `r` for the left and right boundaries respectively.
     */
    constexpr Operator(const Stencil<width_interior> i,
                       const StencilTensor2D<height_boundary_l,width_boundary_l> l,
                       const StencilTensor2D<height_boundary_r,width_boundary_r> r);

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

    Stencil<width_interior> get_interior_stencil() const;

    StencilTensor2D<height_boundary_l,width_boundary_l> get_left_boundary_stencil() const;

    StencilTensor2D<height_boundary_r,width_boundary_r> get_right_boundary_stencil() const;


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
             const StencilTensor2D<height_boundary_l,width_boundary_l> l,
             const StencilTensor2D<height_boundary_r,width_boundary_r> r)
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
      dst[i] = -rboundary[height_boundary_r-1-i].apply_flip(src,i);
    }

    for(int i = height_boundary_r; i < n-height_boundary_l; ++i) {
      dst[i] = -interior.apply_flip(src,i);
    }

    for(int i = n-height_boundary_l; i < n; ++i) {
      dst[i] = -lboundary[n-1-i].apply_flip(src,i);
    }
  }

 template <int width_interior,
            int width_boundary_l,
            int height_boundary_l,
            int width_boundary_r,
            int height_boundary_r>
  Stencil<width_interior> Operator<width_interior,width_boundary_l,height_boundary_l, width_boundary_r, height_boundary_r>
  ::get_interior_stencil() const{
   return interior;
  }

 template <int width_interior,
            int width_boundary_l,
            int height_boundary_l,
            int width_boundary_r,
            int height_boundary_r>
  StencilTensor2D<height_boundary_l,width_boundary_l>  Operator<width_interior,width_boundary_l,height_boundary_l, width_boundary_r, height_boundary_r>
  ::get_left_boundary_stencil() const{
    return lboundary;
  }

   template <int width_interior,
            int width_boundary_l,
            int height_boundary_l,
            int width_boundary_r,
            int height_boundary_r>

  StencilTensor2D<height_boundary_r,width_boundary_r> Operator<width_interior,width_boundary_l,height_boundary_l, width_boundary_r, height_boundary_r>
  ::get_right_boundary_stencil() const{
    return rboundary;
  }


}

#endif /* _OPERATOR_H */
