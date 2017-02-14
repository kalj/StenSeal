#include "stenseal/geometry.h"
#include <deal.II/base/function.h>

namespace stenseal {

  // cartesian geometry

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
  unsigned CartesianGeometry<dim>::get_n_nodes(int d) const
  {
    return n_nodes[d];
  }

  template <int dim>
  unsigned CartesianGeometry<dim>::get_n_nodes_total() const
  {
    return n_nodes_total;
  }

  template <int dim>
  double CartesianGeometry<dim>::get_mapped_h(int d) const
  {
    return h[d];
  }

  template <int dim>
  void CartesianGeometry<dim>::initialize_vector(dealii::Vector<double> &u,
                                                 const std::function<double (const dealii::Point<dim>&)> &f) const
  {
    initialize_vector(u, dealii::ScalarFunctionFromFunctionObject<dim>(f));
  }

  template <int dim>
  void CartesianGeometry<dim>::initialize_vector(dealii::Vector<double> &u,
                                                 const dealii::Function<dim> &f) const
  {
    if(dim==1) {
      for(int i = 0; i < n_nodes[0]; ++i) {
        const double x = i*h[0] + lower_left(0);
        const dealii::Point<dim> p(x);
        u[i] = f.value(p);
      }
    }
    else if(dim==2){
      for(int iy = 0; iy < n_nodes[1]; ++iy) {
        for(int ix = 0; ix < n_nodes[0]; ++ix) {
          const double x = ix*h[0] + lower_left(0);
          const double y = iy*h[1] + lower_left(1);
          const dealii::Point<dim> p(x,y);
          const int idx = iy*n_nodes[0]+ix;
          u[idx] = f.value(p);
        }
      }
    }
    else if(dim==3){
      for(int iz = 0; iz < n_nodes[2]; ++iz) {
        for(int iy = 0; iy < n_nodes[1]; ++iy) {
          for(int ix = 0; ix < n_nodes[0]; ++ix) {
            const double x = ix*h[0] + lower_left(0);
            const double y = iy*h[1] + lower_left(1);
            const double z = iz*h[2] + lower_left(2);
            const dealii::Point<dim> p(x,y,z);
            const int idx = (iz*n_nodes[1] + iy)*n_nodes[0]+ix;
            u[idx] = f.value(p);
          }
        }
      }
    }
    else {
      AssertThrow(false,dealii::ExcNotImplemented());
    }
  }

  // general geometry

  template <int dim>
  GeneralGeometry<dim>::GeneralGeometry(const std::array<unsigned int,dim> n_nodes,
                                        const std::vector<dealii::Point<dim>> points)
    : n_nodes(n_nodes), nodes(points)
  {
    n_nodes_total=1;
    for(int d = 0; d < dim; ++d) {
      n_nodes_total *=n_nodes[d];
      h_mapped[d] = 1.0/(n_nodes[0]-1);
    }
  }

  template <int dim>
  unsigned GeneralGeometry<dim>::get_n_nodes(int d) const
  {
    return n_nodes[d];
  }

  template <int dim>
  unsigned GeneralGeometry<dim>::get_n_nodes_total() const
  {
    return n_nodes_total;
  }

  template <int dim>
  double GeneralGeometry<dim>::get_mapped_h(int d) const
  {
    return h_mapped[d];
  }

  template <int dim>
  const std::vector<dealii::Point<dim>>& GeneralGeometry<dim>::get_node_points() const
  {
    return nodes;
  }

  template <int dim>
  void GeneralGeometry<dim>::initialize_vector(dealii::Vector<double> &u,
                                               const std::function<double (const dealii::Point<dim>&)> &f) const
  {
    initialize_vector(u, dealii::ScalarFunctionFromFunctionObject<dim>(f));
  }

  template <int dim>
  void GeneralGeometry<dim>::initialize_vector(dealii::Vector<double> &u,
                                               const dealii::Function<dim> &f) const
  {
    if(dim<4 && dim>0) {
      for(int i = 0; i < n_nodes_total; ++i) {
        u[i] = f.value(nodes[i]);
      }
    }
    else {
      AssertThrow(false,dealii::ExcNotImplemented());
    }
  }

  namespace internal
  {
    template <>
    dealii::Point<1> repeat_point(double d)
    {
      return dealii::Point<1>(d);
    }

    template <>
    dealii::Point<2> repeat_point(double d)
    {
      return dealii::Point<2>(d,d);
    }

    template <>
    dealii::Point<3> repeat_point(double d)
    {
      return dealii::Point<3>(d,d,d);
    }

  }

  template class CartesianGeometry<1>;
  template class CartesianGeometry<2>;
  template class CartesianGeometry<3>;

  template class GeneralGeometry<1>;
  template class GeneralGeometry<2>;
  template class GeneralGeometry<3>;

}
