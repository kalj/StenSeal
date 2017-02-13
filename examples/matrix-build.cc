
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>

#include "stenseal/geometry.h"
#include "stenseal/operator.h"
#include "stenseal/upwind_laplace.h"
#include "stenseal/operator_lib.h"

int main(int argc, char *argv[])
{

  int N = 10;

  //const auto Dm = stenseal::upwind_operator_2nd_order();
  const auto Dm = stenseal::upwind_operator_3rd_order();

  const int height_r = Dm.height_r;
  const int height_l = Dm.height_l;
  const int width_l = Dm.width_l;
  const int width_r = Dm.width_r;
  const int width_i =  Dm.width_i;
  const auto lboundary = Dm.get_left_boundary_stencil();
  const auto rboundary = Dm.get_right_boundary_stencil();
  const auto interior = Dm.get_interior_stencil();

  std::vector<unsigned int> row_lengths(N);

  //sprcity pattern
  for(int i = 0; i <height_r; ++i) {
    row_lengths[i] = width_l;
  }

  for(int i = height_r; i < N - height_l; i++){
    row_lengths[i] = width_i;
  }

  for(int i = N - height_l; i < N; i++ ){
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
 for(int i = height_l; i < N-height_l; i++ ){
  for(int j = 0; j < width_i; j++){
    sp_Dm.add(i,i+offseti[j]);
  }
}


for(int i = N - height_l; i < N; i++ ){
  const auto offsetr = rboundary[i-N+height_l].get_offsets();
  for(int j = 0; j< width_r; j++){
   sp_Dm.add(i,i + offsetr[j]);
 }
}


sp_Dm.compress();

dealii::SparseMatrix<double> matrix(sp_Dm);

// Fill matrix

for(int i = 0; i < height_r; ++i) {
 const auto weights = lboundary[i].get_weights();
 const auto offsetl = lboundary[i].get_offsets();
 for(int j = 0; j< width_l; j++){
  matrix.set(i, i + offsetl[j],weights[j]);
}
}

const auto weights = interior.get_weights();
for(int i = height_l; i < N-height_r; i++ ){
  for(int j = 0; j < width_i; j++){
    matrix.set(i,i + offseti[j], weights[j]);
  }
}

for(int i = N - height_r; i < N; i++ ){
  const auto weights = rboundary[i-N+height_l].get_weights();
  const auto offsetr = rboundary[i-N+height_l].get_offsets();
  for(int j = 0; j< width_r; j++){
   matrix.set(i,i + offsetr[j], weights[j]);
 }
}

matrix.print(std::cout);

return 0;

}