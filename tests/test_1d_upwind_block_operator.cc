
#include <deal.II/lac/vector.h>
#include <deal.II/base/numbers.h>

#include "stenseal/block_operator.h"

/**
   Test second order upwind laplace
*/


//=============================================================================
// Stencils
//=============================================================================

struct StencilDm {
  // boundary
  const static unsigned int boundary_width = 2;
  const static unsigned int boundary_height = 1;
  const static double boundary_stencil[];

  // interior
  const static unsigned interior_width = 3;
  const static int interior_start_idx = -1;
  const static double interior_stencil[];
};

const double StencilDm::boundary_stencil[] = { 0.5, -0.5 };
const double StencilDm::interior_stencil[] = { 0.5, 0.0, -0.5};


struct StencilDp {
  // boundary
  const static unsigned int boundary_width = 2;
  const static unsigned int boundary_height = 1;
  const static double boundary_stencil[];

  // interior
  const static unsigned interior_width = 3;
  const static int interior_start_idx = -1;
  const static double interior_stencil[];
};

const double StencilDp::boundary_stencil[] = { -0.5, 0.5 };
const double StencilDp::interior_stencil[] = {-0.5, 0.0, 0.5};


double compute_l2_norm(int n)
{
  const int dim = 1;

  int n_nodes[dim] = { n };
  int n_nodes_tot = n_nodes[0];
  double h = 1.0/(n_nodes_tot-1);
  const double PI = dealii::numbers::PI;

  stenseal::UpwindBlockOperator<dim,StencilDm,StencilDp> op(n_nodes);

  dealii::Vector<double> u(n_nodes_tot);

  for(int i = 0; i < n_nodes_tot; ++i) {
    u[i] = sin(PI*i*h);
  }

  dealii::Vector<double> v(n_nodes_tot);

  op.apply(v,u);

  double l2_norm = 0;
  for(int i = 0; i < n_nodes_tot; ++i) {
    double a = v[i]-(-(PI*h)*(PI*h)*sin(PI*i*h));
    l2_norm += a*a;
  }

  return l2_norm;
}


int main(int argc, char *argv[])
{
  const int n_tests = 7;
  double norms[n_tests];
  unsigned int size = 40;

  for(int i=0; i<n_tests; i++) {
    norms[i] = compute_l2_norm(size);
    size *= 2;
  }

  printf("     n        error     conv\n");
  printf(" ----------------------------\n");

  size=40;
  printf("%6d %12.4g        -\n",size,norms[0]);
  size *= 2;

  bool all_conv = true;
  for(int i=1; i<n_tests; i++) {
    double conv = norms[i-1] / norms[i];
    printf("%6d %12.4g %8.4g\n",size,norms[i],conv);
    size *= 2;

    all_conv = all_conv && conv > 3.999;
  }


  printf("\n");
  if(all_conv) {
    printf("Proper convergence order attained\n");
    return 0;
  }
  else {
    printf("Proper convergence NOT attained\n");
    return 1;
  }

}
