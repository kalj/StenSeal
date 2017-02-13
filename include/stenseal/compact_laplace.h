/* -*- c-basic-offset:2; tab-width:2; indent-tabs-mode:nil -*-
 *
 */

#ifndef _COMPACT_LAPLACE_H
#define _COMPACT_LAPLACE_H

namespace stenseal
{
  /**
   * Class representing an Compact Laplace operator on a `dim` dimensional grid
   * block with geometry defined by the type `Geometry`. The operator is defined
   * by the following stencils:
   *
   * - The compact first derivative Operator D1, including an interior Stencil
   *   and a boundary StencilArray
   * - The compact second derivative D2 defined by metric-coefficient
   *   StencilArrays for the interior and for the boundary (array of StencilArrays).
   */
  template <int dim, typename D2Operator, typename D1Operator, typename Geometry>
  class CompactLaplace
  {
  private:
    const D2Operator D2;
    const D1Operator D1;
    const Geometry geometry;
  public:
    /**
   * Constructor. Takes the 1D compact second-derivative SBP operator `d2`, the
   * 1D compact first-derivative SBP operator d1, and the geometry descriptor
   * `g` of this grid block.
   */
    CompactLaplace(const D2Operator d2, const D1Operator d1, const Geometry geom);

  /**
   * Apply this operator to the vector `src` writing the result into `dst`.
   */
    void apply(dealii::Vector<double> &dst, const dealii::Vector<double> &src) const;
  };



  //---------------------------------------------------------------------------
  // implementations
  //---------------------------------------------------------------------------


  template <int dim, typename D2Operator, typename D1Operator, typename Geometry>
  CompactLaplace<dim,D2Operator,D1Operator,Geometry>
  ::CompactLaplace(const D2Operator d2, const D1Operator d1, const Geometry geom)
    : D2(d2), D1(d1), geometry(geom)
  {}


  template <int dim, typename D2Operator, typename D1Operator, typename Geometry>
  void CompactLaplace<dim,D2Operator,D1Operator,Geometry>
  ::apply(dealii::Vector<double> &dst, const dealii::Vector<double> &src) const
  {
    if(dim == 1) {
      const unsigned int n = geometry.get_n_nodes(0);

      D2.apply(dst,src,geometry.get_metric_coefficient(),n);

      // divide by h^2
      const double h2 = geometry.get_h(0)*geometry.get_h(0);
      dst /= h2;
    }
    else {
      AssertThrow(false,dealii::ExcNotImplemented());
    }
  }
}

#endif /* _COMPACT_LAPLACE_H */
