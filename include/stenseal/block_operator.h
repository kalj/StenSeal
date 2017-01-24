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

  template <int dim, typename DmT, typename Geometry>
  class UpwindBlockOperator
  {
  private:
    const DmT Dm;
    const Geometry geometry;

    unsigned int n_nodes_tot;

  public:
    UpwindBlockOperator(const DmT dm, const Geometry &g);

    void apply(dealii::Vector<double> &dst, const dealii::Vector<double> &src) const;

  };

  template <int dim, typename DmT, typename Geometry>
  UpwindBlockOperator<dim,DmT,Geometry>
  ::UpwindBlockOperator(const DmT dm, const Geometry &g)
    : Dm(dm), geometry(g)
  {
    n_nodes_tot = 1;
    for(int d = 0; d < dim; ++d) {
      n_nodes_tot *= g.n_nodes[d];
    }
  }


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

      const unsigned int n = geometry.n_nodes[0];

      // for now, use full temporary array
      // FIXME: this is stupid!! instead use blocking
      dealii::Vector<double> tmp (n_nodes_tot);

      //-------------------------------------------------------------------------
      // apply Dm
      //-------------------------------------------------------------------------

      Dm.apply(tmp,src,n);

      // stupid inefficient way of multiplying with c and 1/h^2
      // FIXME: merge with above.
      for(int i = 0; i<n; ++i) {
        tmp[i] *= geometry.get_c(i);
      }

      //-------------------------------------------------------------------------
      // apply Dp
      //-------------------------------------------------------------------------
      Dm.apply_dp(dst,tmp,n);


    }
    else {
      AssertThrow(false,dealii::ExcNotImplemented());
    }

  }

}

#endif /* _UPWIND_BLOCK_OPERATOR_H */
