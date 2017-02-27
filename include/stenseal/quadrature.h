/* -*- c-basic-offset:2; tab-width:2; indent-tabs-mode:nil -*-
 *
 */

#ifndef _QUADRATURE_H
#define _QUADRATURE_H

#include <vector>

#include <deal.II/lac/vector.h>


//=============================================================================
// this is a quadrature in 1D
//=============================================================================
namespace stenseal
{
  /**
   * Class representing a  quadrature
   *
   */
  template <int width_boundary>
  class Quadrature
  {
  private:
    const std::array<double , width_boundary> quad;
  public:
    const static int width = width_boundary;


    /**
     * Constructor.
     */
    constexpr Quadrature(const std::array<double,width_boundary> h);

    /**
     * Apply this `Quadrature` to the vector `src` and write the result into
     * `dst`
     */
    void apply(dealii::Vector<double> &dst,
               const dealii::Vector<double> &src,
               const unsigned int n) const;
  };

  //---------------------------------------------------------------------------
  // Implementations
  //---------------------------------------------------------------------------

  template <int width_boundary>
  constexpr Quadrature<width_boundary>
  ::Quadrature(const std::array<double,width_boundary> h)
    :quad(h)
  {}

  template <int width_boundary>
  void Quadrature<width_boundary>
  ::apply(dealii::Vector<double> &dst,
          const dealii::Vector<double> &src,
          const unsigned int n) const
  {
    dst = src;

    for(int i = 0; i < width_boundary; ++i) {
      dst[i] = quad[i]*src(i);
    }

    for(int i = n-width_boundary; i < n; ++i) {
      dst[i] = quad[n - i -1]*src(i);
    }
  }

}

#endif /* _QUADRATURE_H */
