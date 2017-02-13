#include "stenseal/geometry.h"

namespace stenseal {

  template <int dim>
  CartesianGeometry<dim>::CartesianGeometry(const std::array<unsigned int,dim> n_nodes,
                                            const dealii::Point<dim> lower_left,
                                            const dealii::Point<dim> upper_right)
    : n_nodes(n_nodes), lower_left(lower_left)
  {
    n_nodes_total=1;
    for(int d = 0; d < dim; ++d) {
      n_nodes_total *=n_nodes[d];
      h[d] = (upper_right[d]-lower_left[d])/(n_nodes[d]-1);
    }
  }

  template <int dim>
  void CartesianGeometry<dim>::initialize_vector(dealii::Vector<double> &u,
                                                 dealii::Function<dim> &f) const
  {
    if(dim==1) {
      for(int i = 0; i < n_nodes[0]; ++i) {
        double x = i*h[0] + lower_left(0);
        dealii::Point<dim> p(x);
        u[i] = f.value(p);
      }
    }
    else if(dim==2){
      for(int i = 0; i < n_nodes[1]; ++i) {
        for(int j = 0; j < n_nodes[0]; ++j) {
          double x = j*h[0] + lower_left(0);
          double y = i*h[1] + lower_left(1);
          dealii::Point<dim> p(x,y);
          u[n_nodes[0]*i+j] = f.value(p);
        }
      }
    }
    else {
      AssertThrow(false,dealii::ExcNotImplemented());
    }
  }


  template class CartesianGeometry<1>;
  template class CartesianGeometry<2>;

}
