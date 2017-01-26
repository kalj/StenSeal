#include "stenseal/geometry.h"

namespace stenseal {

  template <int dim>
  CartesianGeometry<dim>::CartesianGeometry(unsigned int n_nodes[dim],
                                            const dealii::Point<dim> lower_left,
                                            const dealii::Point<dim> upper_right)
  {
    this->lower_left=lower_left;
    this->n_nodes_total=1;
    for(int d = 0; d < dim; ++d) {
      this->n_nodes[d] = n_nodes[d];
      this->n_nodes_total *=n_nodes[d];
      h[d] = (upper_right[d]-lower_left[d])/(n_nodes[d]-1);
    }
  }


  template class CartesianGeometry<1>;
  template class CartesianGeometry<2>;

}
