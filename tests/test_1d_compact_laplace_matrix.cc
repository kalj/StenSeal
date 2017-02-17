
#include <deal.II/lac/vector.h>
#include <deal.II/base/numbers.h>

#include "stenseal/operator.h"
#include "stenseal/compact_laplace.h"
#include "stenseal/operator_lib.h"
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>

template <typename OperatorTypeD2, typename OperatorTypeD1>
     bool test_matrix(std::pair<OperatorTypeD2,OperatorTypeD1> ops, int order)
     {
      const int dim = 1;
      const int n = 10;
      bool all_conv = false;
      std::array<unsigned int,dim> n_nodes{ n };
      int n_nodes_tot = n_nodes[0];


      typedef stenseal::CartesianGeometry<dim> Geometry;
      Geometry geometry(n_nodes);

      dealii::SparsityPattern sp_Laplace; //(n,n);
      //sp_Laplace.compress();
      dealii::SparseMatrix<double> matrix_Laplace;//(sp_Laplace);


      stenseal::CompactLaplace<dim,OperatorTypeD2,OperatorTypeD1,Geometry> op(ops.first,ops.second,geometry);

      op.matrix(matrix_Laplace, sp_Laplace);

      matrix_Laplace.print(std::cout);


}

int main(int argc, char *argv[])
{
  bool all_conv = true;

  printf("Second order Upwind: ");
  all_conv = test_matrix(stenseal::compact_operators_2nd_order(),2);

  if(all_conv) {
    printf("\nAll matrices are right\n");
    return 0;
  }
  else {
    printf("\nOne or more matrices have the wrong coefficinets\n");
    return 1;
  }
}