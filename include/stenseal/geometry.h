/* -*- c-basic-offset:2; tab-width:2; indent-tabs-mode:nil -*-
 *
 *
 */

#ifndef _GEOMETRY_H
#define _GEOMETRY_H

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

namespace stenseal
{
  namespace internal
  {
    template <int dim>
    dealii::Point<dim> repeat_point(double d);
  }


  template <int dim>
  class CartesianGeometry
  {
  private:
    std::array<double,dim> h;
    std::array<unsigned int,dim> n_nodes;
    unsigned int n_nodes_total;
    dealii::Point<dim> lower_left;

  public:

    CartesianGeometry(const std::array<unsigned int,dim> n_nodes,
                      const dealii::Point<dim> lower_left=internal::repeat_point<dim>(0.0),
                      const dealii::Point<dim> upper_right=internal::repeat_point<dim>(1.0));

    unsigned get_n_nodes(int d) const
    {
      return n_nodes[d];
    }

    unsigned get_n_nodes_total() const
    {
      return n_nodes_total;
    }

    double get_h(int d) const
    {
      return h[d];
    }

    void initialize_vector(dealii::Vector<double> &u,
                           const dealii::Function<dim> &f) const;

    void initialize_vector(dealii::Vector<double> &u,
                           const std::function<double (const dealii::Point<dim> &)> &f) const;

  };

  template <int dim>
  class GeneralGeometry
  {
  private:
    const std::vector<dealii::Point<dim>> nodes;
    std::array<unsigned int,dim> n_nodes;
    unsigned int n_nodes_total;

  public:

    GeneralGeometry(const std::array<unsigned int,dim> n_nodes,
                    const std::vector<dealii::Point<dim>> points);

    unsigned get_n_nodes(int d) const
    {
      return n_nodes[d];
    }

    unsigned get_n_nodes_total() const
    {
      return n_nodes_total;
    }

    void initialize_vector(dealii::Vector<double> &u,
                           const dealii::Function<dim> &f) const;

    void initialize_vector(dealii::Vector<double> &u,
                           const std::function<double (const dealii::Point<dim> &)> &f) const;
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
