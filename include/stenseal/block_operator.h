/* -*- c-basic-offset:2; tab-width:2; indent-tabs-mode:nil -*-
 *
 * @(#)block_operator.h
 * @author Karl Ljungkvist <k.ljungkvist@gmail.com>
 *
 */

#ifndef _UPWIND_BLOCK_OPERATOR_H
#define _UPWIND_BLOCK_OPERATOR_H

#include <deal.II/lac/vector.h>

namespace stenseal
{

  struct Geometry {
    double get_c(int i) const {
      return 1.0;
    }
  };

  template <int dim, typename StencilDm, typename StencilDp>
  class UpwindBlockOperator
  {
  private:
    Geometry geometry;

    unsigned int n_nodes_1d[dim];
    unsigned int n_nodes_tot;

  public:
    UpwindBlockOperator(int n_nodes[dim]);

    void apply(dealii::Vector<double> &dst, const dealii::Vector<double> &src) const;

  };

  template <int dim, typename StencilDm, typename StencilDp>
  UpwindBlockOperator<dim,StencilDm,StencilDp>::UpwindBlockOperator(int n_nodes[dim])
  {
    n_nodes_tot = 1;
    for(int d = 0; d < dim; ++d) {
      n_nodes_1d[d] = n_nodes[d];
      n_nodes_tot *= n_nodes[d];
    }
  }


  template <int dim, typename StencilDm, typename StencilDp>
  void UpwindBlockOperator<dim,StencilDm,StencilDp>::apply(dealii::Vector<double> &dst,
                                                           const dealii::Vector<double> &src) const
  {


    /*
      Ideas:
      - blocking
      - use several index variables stepping up/down to avoid complicated shift and index computations
    */

    if(dim==1) {

      // for now, use full temporary array
      // FIXME: this is stupid!!
      dealii::Vector<double> tmp (n_nodes_1d[0]);

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
          i < (n_nodes_1d[0]-StencilDm::boundary_height) ;
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

          t += src[n_nodes_1d[0]-StencilDm::boundary_width+j] * StencilDm::boundary_stencil[(StencilDm::boundary_height-1 - i)*StencilDm::boundary_width
                                                                                            + (StencilDm::boundary_width-1 - j)];

        tmp[n_nodes_1d[0]-StencilDm::boundary_height + i] = geometry.get_c(n_nodes_1d[0]-StencilDm::boundary_height + i)*t;
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
          i < (n_nodes_1d[0]-StencilDp::boundary_height) ;
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

          t += tmp[n_nodes_1d[0]-StencilDp::boundary_width+j] * StencilDp::boundary_stencil[(StencilDp::boundary_height-1 - i)*StencilDp::boundary_width
                                                                                            + (StencilDp::boundary_width-1 - j)];

        dst[n_nodes_1d[0]-StencilDp::boundary_height + i] = t;
      }



    }
    else {
      AssertThrow(false,dealii::ExcNotImplemented());
    }

  }

}

#endif /* _UPWIND_BLOCK_OPERATOR_H */
