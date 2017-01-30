
#include <deal.II/lac/vector.h>
#include <deal.II/base/numbers.h>

#include "stenseal/operator.h"
#include "stenseal/block_operator.h"
#include "stenseal/operator_lib.h"

/**
   Test second order upwind laplace
*/

void compute_l2_norm(unsigned int n, double &l2_norm, double &l2_norm_interior)
{
  const int dim = 1;

  unsigned int n_nodes[dim] = { n };
  int n_nodes_tot = n_nodes[0];
  double h = 1.0/(n_nodes_tot-1);
  const double PI = dealii::numbers::PI;

  // FIXME: call without points once default values are supported
  typedef stenseal::CartesianGeometry<dim> Geometry;
  Geometry geometry(n_nodes, dealii::Point<dim>(0.0),
                    dealii::Point<dim>(1.0));

  // typedef stenseal::Operator<2,2,1> OperatorType;
  // const stenseal::Symbol usym;

  // define stencil
  // constexpr OperatorType Dm((-0.5)*usym[-1] + 0.5*usym[1],  // interior stencil
  //                           (-1.0)*usym[0] + 1.0*usym[1],   // left boundary
  //                           (-1.0)*usym[-1] + 1.0*usym[0]); // right boundary

  typedef stenseal::Operator<2,2,1,2,1> OperatorType;
  const stenseal::Symbol usym;

  constexpr OperatorType Dm = stenseal::upwind_operator_2nd_order();

  stenseal::UpwindBlockOperator<dim,OperatorType,Geometry> op(Dm,geometry);

  dealii::Vector<double> u(n_nodes_tot);

  for(int i = 0; i < n_nodes_tot; ++i) {
    u[i] = sin(PI*i*h);
  }

  dealii::Vector<double> v(n_nodes_tot);

  op.apply(v,u);


  // compute norms
  double sqsum = 0;
  double a;
  for(int i = 1; i < n_nodes_tot-1; ++i) {
    a = v[i] - (-PI*PI*sin(PI*i*h));
    sqsum += a*a;
  }

  l2_norm_interior = std::sqrt(h*sqsum);

  a = v[0] - (-PI*PI*sin(PI*0.0));
  sqsum += a*a;

  a = v[n_nodes_tot-1] - (-PI*PI*sin(PI*1.0));
  sqsum += a*a;

  l2_norm = std::sqrt(h*sqsum);
}


int main(int argc, char *argv[])
{
  const int n_tests = 7;
  double full_norms[n_tests];
  double interior_norms[n_tests];
  unsigned int size = 40;

  for(int i=0; i<n_tests; i++) {
    compute_l2_norm(size,full_norms[i],interior_norms[i]);
    size *= 2;
  }

  printf("     n   full error   factor  (p)   inter. error   factor   (p)\n");
  printf(" ---------------------------------------------------------------\n");

  size=40;
  printf("%6d %12.4g       - ( - )   %12.4g        - ( - )\n",size,full_norms[0],interior_norms[0]);
  size *= 2;

  bool all_conv = true;
  for(int i=1; i<n_tests; i++) {
    double full_conv = full_norms[i-1] / full_norms[i];
    double full_p = std::log2(full_conv);
    double interior_conv = interior_norms[i-1] / interior_norms[i];
    double interior_p = std::log2(interior_conv);
    printf("%6d %12.4g %7.3g (%.1f)   %12.4g %8.3g (%.1f)\n",size,full_norms[i],full_conv,full_p,
           interior_norms[i],interior_conv,interior_p);
    size *= 2;

    all_conv = all_conv && (interior_p > 1.99 && full_p > 1.499);
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
