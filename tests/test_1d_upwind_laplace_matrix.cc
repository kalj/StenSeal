
#include <deal.II/lac/vector.h>
#include <deal.II/base/numbers.h>

#include "stenseal/operator.h"
#include "stenseal/upwind_laplace.h"
#include "stenseal/operator_lib.h"
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>

/**
 * Test second order upwind laplace
 */


template <typename OperatorType>
bool test_matrix(OperatorType Dm)
{
  const int dim = 1;
  const int n = 10;
  unsigned int n_nodes[dim] = { n };
  int n_nodes_tot = n_nodes[0];

  typedef stenseal::CartesianGeometry<dim> Geometry;
  Geometry geometry(n_nodes, dealii::Point<dim>(0.0),
                    dealii::Point<dim>(1.0));

  dealii::SparsityPattern sp_Laplace(n,n);
  sp_Laplace.compress();
  dealii::SparseMatrix<double> matrix_Laplace(sp_Laplace);


  stenseal::UpwindLaplace<dim,OperatorType,Geometry> op(Dm,geometry);

  op.matrix(matrix_Laplace);

  bool all_conv = true;
  matrix_Laplace.print(std::cout);
  //Add somthg that actually tests if the matrix is right :)
  return all_conv;
}

int main(int argc, char *argv[])
{
  bool all_conv = true;

  printf("Second order Upwind:\n");
  all_conv = test_matrix(stenseal::upwind_operator_2nd_order());
  printf("\n Kalles Second order Upwind:\n");
  all_conv = test_matrix(stenseal::upwind_operator_2nd_order_kalle()) && all_conv;

  printf("\n Third order Upwind:\n");
  all_conv = test_matrix(stenseal::upwind_operator_3rd_order()) && all_conv;

  printf("\n Fourth order Upwind:\n");
  all_conv = test_matrix(stenseal::upwind_operator_4th_order()) && all_conv;


  if(all_conv) {
    printf("All matrices are right\n");
    return 0;
  }
  else {
    printf("One or more matrices have the wrong coefficinets\n");
    return 1;
  }
}
