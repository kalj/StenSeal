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

  template <int dim, typename StencilDm, typename StencilDp, typename Geometry>
  class UpwindBlockOperator
  {
  private:
    Geometry geometry;

    unsigned int n_nodes_tot;

  public:
    UpwindBlockOperator(const Geometry &g);

    void apply(dealii::Vector<double> &dst, const dealii::Vector<double> &src) const;

  };

  template <int dim, typename StencilDm, typename StencilDp, typename Geometry>
  UpwindBlockOperator<dim,StencilDm,StencilDp,Geometry>
  ::UpwindBlockOperator(const Geometry &g)
    : geometry(g)
  {
    n_nodes_tot = 1;
    for(int d = 0; d < dim; ++d) {
      n_nodes_tot *= g.n_nodes[d];
    }
  }


  template <int dim, typename StencilDm, typename StencilDp, typename Geometry>
  void UpwindBlockOperator<dim,StencilDm,StencilDp,Geometry>
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
      // FIXME: this is stupid!!
      dealii::Vector<double> tmp (n_nodes_tot);

      //-------------------------------------------------------------------------
      // apply Dm
      //-------------------------------------------------------------------------

      // left boundary
      for(int i = 0;
          i < StencilDm::boundary_height;
          ++i) {

        double t = 0;
        for(int j = 0; j < StencilDm::boundary_width; ++j)
          t += src[0+j]*StencilDm::boundary_stencil[i*StencilDm::boundary_width + j];

        tmp[i] = geometry.get_c(i)*t;
      }

      // interior
      for(int i = StencilDm::boundary_height;
          i < (n - StencilDm::boundary_height) ;
          ++i) {

        double t = 0;
        for(int j = 0;
            j < StencilDm::interior_width;
            ++j)
          t += src[i+StencilDm::interior_start_idx+j]*StencilDm::interior_stencil[j];

        tmp[i] = geometry.get_c(i)*t;

      }

      // right boundary
      for(int i = 0;
          i < StencilDm::boundary_height;
          ++i) {

        double t = 0;
        for(int j = 0; j < StencilDm::boundary_width; ++j)

          t += src[n - StencilDm::boundary_width+j] * (-1) * StencilDm::boundary_stencil[(StencilDm::boundary_height-1 - i)*StencilDm::boundary_width
                                                                                         + (StencilDm::boundary_width-1 - j)];

        tmp[n - StencilDm::boundary_height + i] = geometry.get_c(n - StencilDm::boundary_height + i)*t;
      }


      //-------------------------------------------------------------------------
      // apply Dp
      //-------------------------------------------------------------------------

      // left boundary
      for(int i = 0;
          i < StencilDp::boundary_height;
          ++i) {

        double t = 0;
        for(int j = 0; j < StencilDp::boundary_width; ++j)
          t += tmp[0+j]*StencilDp::boundary_stencil[i*StencilDp::boundary_width + j];

        dst[i] = t;
      }

      // interior
      for(int i = StencilDp::boundary_height;
          i < (n - StencilDp::boundary_height) ;
          ++i) {

        double t = 0;
        for(int j = 0;
            j < StencilDp::interior_width;
            ++j)
          t += tmp[i+StencilDp::interior_start_idx+j]*StencilDp::interior_stencil[j];

        dst[i] = t;

      }

      // right boundary
      for(int i = 0;
          i < StencilDp::boundary_height;
          ++i) {

        double t = 0;
        for(int j = 0; j < StencilDp::boundary_width; ++j)

          t += tmp[n - StencilDp::boundary_width+j] * (-1) * StencilDp::boundary_stencil[(StencilDp::boundary_height-1 - i)*StencilDp::boundary_width
                                                                                         + (StencilDp::boundary_width-1 - j)];

        dst[n - StencilDp::boundary_height + i] = t;
      }



    }
    else {
      AssertThrow(false,dealii::ExcNotImplemented());
    }

  }

}

#endif /* _UPWIND_BLOCK_OPERATOR_H */
