/* -*- c-basic-offset:2; tab-width:2; indent-tabs-mode:nil -*-
 *
 *
 */

#ifndef _UPWIND_LAPLACE_H
#define _UPWIND_LAPLACE_H

#include <deal.II/lac/vector.h>
#include "stenseal/geometry.h"
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>

namespace stenseal
{

  /**
   * Class representing an Upwind Laplace operator defined by the 1D 1st
   * derivative "D-minus" operator DmT, on a `dim` dimensional grid block with
   * geometry defined by the type `Geometry`.
   */
  template <int dim, typename DmT, typename Geometry>
  class UpwindLaplace
  {
  private:
    const DmT Dm;
    const Geometry geometry;

  public:
    /**
     * Constructor. Takes the 1D D-minus Upwind SBP operator `dm`, and the
     * geometry descriptor `g` of this grid block.
     */
    UpwindLaplace(const DmT dm, const Geometry &g);

    /**
     * Apply this operator to the vector `src` writing the result into `dst`.
     */
    void apply(dealii::Vector<double> &dst, const dealii::Vector<double> &src) const;
    /**
     *Retuns operator in matrix-form
     */
    void matrix(dealii::SparseMatrix<double> &matrix_Laplace);
  };



  //---------------------------------------------------------------------------
  // implementations
  //---------------------------------------------------------------------------

  template <int dim, typename DmT, typename Geometry>
  UpwindLaplace<dim,DmT,Geometry>
  ::UpwindLaplace(const DmT dm, const Geometry &g)
    : Dm(dm), geometry(g)
  {}


  template <int dim, typename DmT, typename Geometry>
  void UpwindLaplace<dim,DmT,Geometry>
  ::apply(dealii::Vector<double> &dst,
          const dealii::Vector<double> &src) const
  {

    /*
      Ideas:
      - blocking
      - use several index variables stepping up/down to avoid complicated shift and index computations
    */

    if(dim==1) {

      const unsigned int n = geometry.get_n_nodes(0);

      // for now, use full temporary array
      // FIXME: this is stupid!! instead use blocking
      dealii::Vector<double> tmp (geometry.get_n_nodes_total());

      //-------------------------------------------------------------------------
      // apply Dm
      //-------------------------------------------------------------------------

      Dm.apply(tmp,src,n);

      // stupid inefficient way of multiplying with c and 1/h^2
      // FIXME: merge with above.
      const auto &coeff = geometry.get_metric_coefficient();
      for(int i = 0; i<n; ++i) {
        tmp[i] *= coeff.get(i);
      }

      //-------------------------------------------------------------------------
      // apply Dp
      //-------------------------------------------------------------------------
      Dm.apply_dp(dst,tmp,n);

      const double h2 = geometry.get_h(0)*geometry.get_h(0);
      dst /= h2;

    }
    else {
      AssertThrow(false,dealii::ExcNotImplemented());
    }

  }

  template <int dim, typename DmT, typename Geometry>
  void UpwindLaplace<dim,DmT,Geometry>
  ::matrix(dealii::SparseMatrix<double> &matrix_Laplace)
  {
    const unsigned int N = geometry.get_n_nodes(0);

    std::vector<unsigned int> row_lengths(N);

    const int height_r = Dm.height_r;
    const int height_l = Dm.height_l;
    const int width_l = Dm.width_l;
    const int width_r = Dm.width_r;
    const int width_i =  Dm.width_i;
    const auto lboundary = Dm.get_left_boundary_stencil();
    const auto rboundary = Dm.get_right_boundary_stencil();
    const auto interior = Dm.get_interior_stencil();


    //sprcity pattern Dm
    for(int i = 0; i <height_l; ++i) {
      row_lengths[i] = width_l;
    }

    for(int i = height_l; i < N - height_r; i++){
      row_lengths[i] = width_i;
    }

    for(int i = N - height_r; i < N; i++ ){
      row_lengths[i] = width_r;
    }

    dealii::SparsityPattern sp_Dm(N,N,row_lengths);
    for(int i = 0; i <height_r; ++i) {
      const auto offsetl = lboundary[i].get_offsets();
      for(int j = 0; j< width_l; j++){
        sp_Dm.add(i,i+offsetl[j]);
      }
    }

    const auto offseti = interior.get_offsets();
    for(int i = height_l; i < N-height_r; i++ ){
      for(int j = 0; j < width_i; j++){
        sp_Dm.add(i,i+offseti[j]);
      }
    }


    for(int i = N - height_r; i < N; i++ ){
      const auto offsetr = rboundary[i-N+height_l].get_offsets();
      for(int j = 0; j< width_r; j++){
        sp_Dm.add(i,i + offsetr[j]);
      }
    }


    sp_Dm.compress();

    dealii::SparseMatrix<double> matrixDm(sp_Dm);

    // Fill matrixDm Dm

    for(int i = 0; i < height_l; ++i) {
      const auto weights = lboundary[i].get_weights();
      const auto offsetl = lboundary[i].get_offsets();
      for(int j = 0; j< width_l; j++){
        matrixDm.set(i, i + offsetl[j],weights[j]);
      }
    }

    const auto weights = interior.get_weights();
    for(int i = height_l; i < N-height_r; i++ ){
      for(int j = 0; j < width_i; j++){
        matrixDm.set(i,i + offseti[j], weights[j]);
      }
    }

    for(int i = N - height_r; i < N; i++ ){
      const auto weights = rboundary[i-N+height_l].get_weights();
      const auto offsetr = rboundary[i-N+height_l].get_offsets();
      for(int j = 0; j< width_r; j++){
        matrixDm.set(i,i + offsetr[j], weights[j]);
      }
    }

    //matrixDm.print(std::cout);

    for(int i = 0; i <height_r; ++i) {
      row_lengths[i] = width_r;
    }

    for(int i = height_r; i < N - height_l; i++){
      row_lengths[i] = width_i;
    }

    for(int i = N - height_l; i < N; i++ ){
      row_lengths[i] = width_l;
    }

    dealii::SparsityPattern sp_Dp(N,N,row_lengths);


    for(int i = 0; i <height_r; ++i) {
      const auto offsetr = rboundary[height_r-i-1].get_offsets();
      for(int j = 0; j< width_r; j++){
        sp_Dp.add(i,i - offsetr[width_r-j-1]);
      }
    }

    const auto offsetip = interior.get_offsets();
    for(int i = height_r; i < N-height_l; i++ ){
      for(int j = 0; j < width_i; j++){
        sp_Dp.add(i,i - offsetip[width_i - j - 1]);
      }
    }


    for(int i = N - height_l; i < N; i++ ){
      const auto offsetr = lboundary[N-i-1].get_offsets();
      for(int j = 0; j< width_l; j++){
        sp_Dp.add(i,i - offsetr[width_l - j - 1]);
      }
    }


    sp_Dp.compress();
    dealii::SparseMatrix<double> matrixDp(sp_Dp);

    for(int i = 0; i < height_r; ++i) {
      const auto weights = rboundary[height_r - i - 1].get_weights();
      const auto offsetr = rboundary[height_r - i - 1].get_offsets();
      for(int j = 0; j< width_r; j++){
        matrixDp.set(i, i - offsetr[width_r-j-1], -weights[width_r-j-1]);
      }
    }

    const auto weights_i = interior.get_weights();
    for(int i = height_l; i < N-height_r; i++ ){
      for(int j = 0; j < width_i; j++){
        matrixDp.set(i,i - offsetip[width_i - j - 1], -weights_i[width_i - j - 1]);
      }
    }

    for(int i = N - height_l; i < N; i++ ){
      const auto weights = lboundary[N-i-1].get_weights();
      const auto offsetr = lboundary[N-i-1].get_offsets();
      for(int j = 0; j< width_l; j++){
        matrixDp.set(i,i - offsetr[width_l - j - 1], -weights[width_l - j - 1]);
      }
    }

    //dealii::SparsityPattern sp_Laplace(N,N,row_lengths);
    //sp_Laplace.compress();
    //dealii::SparseMatrix<double> matrix_Laplace(sp_Laplace);


    matrixDm.dealii::SparseMatrix<double>::mmult(matrix_Laplace,matrixDp);

  }

}

#endif /* _UPWIND_LAPLACE_H */
