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
  private:
    double h[dim];
    unsigned int n_nodes[dim];
    unsigned int n_nodes_total;
    dealii::Point<dim> lower_left;

  public:
    // FIXME: add default parameters which are [0,0,0,0....] and [1,1,1,1,...]
    CartesianGeometry(unsigned int n_nodes[dim],
                      const dealii::Point<dim> lower_left,
                      const dealii::Point<dim> upper_right);

    inline unsigned get_n_nodes(int d) const
    {
      return n_nodes[d];
    }

    inline unsigned get_n_nodes_total() const
    {
      return n_nodes_total;
    }

    inline double get_h(int d) const
    {
      return h[d];
    }


    inline double get_lower_left(int d) const
    {
      return lower_left(d);
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
  */


}

#endif /* _GEOMETRY_H */
