/* -*- c-basic-offset:2; tab-width:2; indent-tabs-mode:nil -*-
 *
 *
 */

#ifndef _UPWIND_BLOCK_OPERATOR_H
#define _UPWIND_BLOCK_OPERATOR_H

#include <deal.II/lac/vector.h>
#include "stenseal/geometry.h"

namespace stenseal
{

  /**
   * Class representing a Upwind Laplace operator defined by the 1D 1st
   * derivative "D-minus" operator DmT, on a `dim` dimensional grid block with
   * geometry defined by the type `Geometry`.
   */
  template <int dim, typename DmT, typename Geometry>
  class UpwindBlockOperator
  {
  private:
    const DmT Dm;
    const Geometry geometry;
    unsigned int n_nodes_tot;

  public:
  /**
   * Constructor. Takes the 1D D-minus Upwind SBP operator `dm`, and the
   * geometry descriptor `g` of this grid block.
   */
    UpwindBlockOperator(const DmT dm, const Geometry &g);

  /**
   * Apply this operator to the vector `src` writing the result into `dst`.
   */
    void apply(dealii::Vector<double> &dst, const dealii::Vector<double> &src) const;
  };



  //---------------------------------------------------------------------------
  // implementations
  //---------------------------------------------------------------------------

  template <int dim, typename DmT, typename Geometry>
  UpwindBlockOperator<dim,DmT,Geometry>
  ::UpwindBlockOperator(const DmT dm, const Geometry &g)
    : Dm(dm), geometry(g)
  {}


  template <int dim, typename DmT, typename Geometry>
  void UpwindBlockOperator<dim,DmT,Geometry>
  ::apply(dealii::Vector<double> &dst,
          const dealii::Vector<double> &src) const
  {


    /*
      Ideas:
      - blocking
      - use several index variables stepping up/down to avoid complicated shift and index computations
    */

    if(dim==1) {

      const unsigned int n = geometry.get_n_nodes(0);

      // for now, use full temporary array
      // FIXME: this is stupid!! instead use blocking
      dealii::Vector<double> tmp (geometry.get_n_nodes_total());

      //-------------------------------------------------------------------------
      // apply Dm
      //-------------------------------------------------------------------------

      Dm.apply(tmp,src,n);

      // stupid inefficient way of multiplying with c and 1/h^2
      // FIXME: merge with above.
      const auto &coeff = geometry.get_metric_coefficient();
      for(int i = 0; i<n; ++i) {
        tmp[i] *= coeff.get(i);
      }

      //-------------------------------------------------------------------------
      // apply Dp
      //-------------------------------------------------------------------------
      Dm.apply_dp(dst,tmp,n);

      const double h2 = geometry.get_h(0)*geometry.get_h(0);
      dst /= h2;

    }
    else {
      AssertThrow(false,dealii::ExcNotImplemented());
    }

  }

}

#endif /* _UPWIND_BLOCK_OPERATOR_H */
