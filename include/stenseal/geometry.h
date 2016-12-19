/* -*- c-basic-offset:2; tab-width:2; indent-tabs-mode:nil -*-
 *
 *
 */

#ifndef _GEOMETRY_H
#define _GEOMETRY_H

#include <deal.II/base/point.h>

namespace stenseal
{

  template <int dim>
  class CartesianGeometry {
  public:
    double h[dim];
    unsigned int n_nodes[dim];

    // FIXME: add default parameters which are [0,0,0,0....] and [1,1,1,1,...]
    CartesianGeometry(unsigned int n_nodes[dim],
                      const dealii::Point<dim> lower_left,
                      const dealii::Point<dim> upper_right);


    inline double get_c(int i) const
    {
      if(dim==1)
        return h[0];
      else {
        // not implemented yet
        return -1;
      }
    }
  };


  /*
  struct TransfiniteInterpolationGeometry {
    double h;
    inline double get_c(int i) const
    {
      return h;
    }
  };

  struct GeneralGeometry {
    double *c;
    inline double get_c(int i) const
    {
      return c[i];
    }
  };
  */



}


#endif /* _GEOMETRY_H */
