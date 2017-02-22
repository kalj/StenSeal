/* -*- c-basic-offset:2; tab-width:2; indent-tabs-mode:nil -*-
 *
 */

#ifndef _COMPACT_LAPLACE_H
#define _COMPACT_LAPLACE_H

#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>

#include "stenseal/metric.h"

namespace stenseal
{
  /**
   * Class representing an Compact Laplace operator on a `dim` dimensional grid
   * block with geometry defined by the type `Geometry`. The operator is defined
   * by the following stencils:
   *
   * - The compact first derivative Operator D1, including an interior Stencil
   *   and a boundary StencilArray
   * - The compact second derivative D2 defined by metric-coefficient
   *   StencilArrays for the interior and for the boundary (array of StencilArrays).
   */
  template <int dim, typename D2Operator, typename D1Operator, typename Geometry>
  class CompactLaplace
  {
  private:
    const D2Operator D2;
    const D1Operator D1;
    const Geometry geometry;
    const Metric<dim,Geometry> metric;
  public:
    /**
     * Constructor. Takes the 1D compact second-derivative SBP operator `d2`, the
     * 1D compact first-derivative SBP operator d1, and the geometry descriptor
     * `g` of this grid block.
     */
    CompactLaplace(const D2Operator d2, const D1Operator d1, const Geometry geom);

    /**
     * Apply this operator to the vector `src` writing the result into `dst`.
     */
    void apply(dealii::Vector<double> &dst, const dealii::Vector<double> &src) const;


    /**
     * Retuns operator in matrix-form
     */
    void matrix(dealii::SparseMatrix<double> &matrix_Laplace, dealii::SparsityPattern &sp_D2 ) const;
  };


  //---------------------------------------------------------------------------
  // implementations
  //---------------------------------------------------------------------------


  template <int dim, typename D2Operator, typename D1Operator, typename Geometry>
  CompactLaplace<dim,D2Operator,D1Operator,Geometry>
  ::CompactLaplace(const D2Operator d2, const D1Operator d1, const Geometry geom)
    : D2(d2), D1(d1), geometry(geom), metric(geom,d1)
  {}


  template <int dim, typename D2Operator, typename D1Operator, typename Geometry>
  void CompactLaplace<dim,D2Operator,D1Operator,Geometry>
  ::apply(dealii::Vector<double> &dst, const dealii::Vector<double> &src) const
  {
    if(dim == 1) {
      const unsigned int n = geometry.get_n_nodes(0);

      const auto &inv_jac = metric.inverse_jacobian();

      D2.apply(dst,src,inv_jac,n);

      // divide by h^2, and multiply by inverse jacobian
      const double h2 = geometry.get_mapped_h(0)*geometry.get_mapped_h(0);


      for(int i = 0; i < n; ++i) {
        dst[i] *= inv_jac.get(i)/h2;
      }
    }
    else {
      AssertThrow(false,dealii::ExcNotImplemented());
    }
  }

  template <int dim, typename D2Operator, typename D1Operator, typename Geometry>
  void CompactLaplace<dim,D2Operator,D1Operator,Geometry>
  ::matrix(dealii::SparseMatrix<double> &matrix_Laplace, dealii::SparsityPattern &sp_D2) const
  {
   if(dim == 1) {
    const unsigned int N = geometry.get_n_nodes(0);
    std::vector<unsigned int> row_lengths(N);

    const auto &inv_jac = metric.inverse_jacobian();
    const int width_i = D2Operator::width_i;
    const int width_b = D2Operator::width_b;
    const int height_b = D2Operator::height_b;
    const auto boundary = D2.get_boundary();
    const auto interior = D2.get_interior();

    for(int i = 0; i <height_b; ++i) {
      row_lengths[i] = width_b;
    }

    for(int i = height_b; i < N - height_b; i++){
      row_lengths[i] = width_i;
    }

    for(int i = N - height_b; i < N; i++ ){
      row_lengths[i] = width_b;
    }

    sp_D2.reinit(N,N,row_lengths);

    for(int i = 0; i <height_b; ++i) {
      for(int j = 0; j< width_b; j++){
        sp_D2.add(i,j);
      }
    }

    int offset = (width_i-1)/2;
    for(int i = height_b; i < N - height_b; i++ ){
      for(int j = 0; j < width_i; j++){
        sp_D2.add(i,j + i - offset);
      }
    }


    for(int i = N - height_b; i < N; i++ ){
      for(int j = 0; j< width_b; j++){
        sp_D2.add(i, N- width_b+j);
      }
    }

    sp_D2.compress();
    matrix_Laplace.reinit(sp_D2);

     // divide by h^2, and multiply by inverse jacobian
     // const double h2 = geometry.get_mapped_h(0)*geometry.get_mapped_h(0);

    for(int i = 0; i <height_b; ++i) {
      for(int j = 0; j< width_b; j++){
       matrix_Laplace.set(i,j, boundary[i][j].apply(metric.inverse_jacobian().template get_left_boundary_array<width_b>()));
      }
    }

   for(int i = height_b; i < N - height_b; i++ ){
    for(int j = 0; j < width_i; j++){
      matrix_Laplace.add(i,j + i - offset, interior[j].apply(metric.inverse_jacobian().template get_centered_array<width_i>(i)));
    }
   }

  for(int i = N - height_b; i < N; i++ ){
    for(int j = 0; j< width_b; j++){
     matrix_Laplace.set(i, N - j - 1, boundary[-i+N-1][j].apply(metric.inverse_jacobian().template get_right_boundary_array<width_b>()));
    }
   }

}
else {
  AssertThrow(false,dealii::ExcNotImplemented());
}
}
}

#endif /* _COMPACT_LAPLACE_H */
