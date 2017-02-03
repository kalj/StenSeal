
#include <deal.II/lac/vector.h>
#include <deal.II/base/numbers.h>

#include "stenseal/operator.h"
#include "stenseal/block_operator.h"
#include "stenseal/operator_lib.h"

#include <fstream>
/**
   Test first derivatives
*/

template <typename OperatorType>
void compute_l2_norm(OperatorType Dm, unsigned int n, double &l2_norm, double &l2_norm_interior,const char * fname)
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

  dealii::Vector<double> u(n_nodes_tot);

  for(int i = 0; i < n_nodes_tot; ++i) {
    u[i] = sin(PI*i*h);
  }

  dealii::Vector<double> v(n_nodes_tot);
  Dm.apply(v,u,n_nodes_tot);

  int height_r = OperatorType::height_r;
  int height_l = OperatorType::height_l;

  // compute norms
  double sqsum = 0;
  double a;

  for(int i = height_l; i < n_nodes_tot-(height_r); ++i) {
    a = 1/h*v[i] - (PI*cos(PI*i*h));
    sqsum += a*a;
  }

  l2_norm_interior = std::sqrt(h*sqsum);


  for(int i= 0; i < height_l; ++i) {
   a = 1/h*v[0] - (PI*cos(PI*i*h));
    sqsum += a*a;
  }

  for(int i = n_nodes_tot-(height_r); i < n_nodes_tot; ++i) {
    a = 1/h*v[n_nodes_tot-1] - (PI*cos(PI*i*h));
    sqsum += a*a;
  }

  l2_norm = std::sqrt(h*sqsum);

  std::ofstream ut(fname);
  v.print(ut);
}


template <typename OperatorType>
bool test_operator(OperatorType op, float interior_p_ref, float full_p_ref, const char *fname)
{
  const int n_tests = 7;
  double full_norms[n_tests];
  double interior_norms[n_tests];
  unsigned int size = 40;

  for(int i=0; i<n_tests; i++) {
    compute_l2_norm(op,size,full_norms[i],interior_norms[i],fname);
    size *= 2;
  }

  printf("     n   full error   factor  (p)   inter. error   factor   (p)\n");
  printf(" ---------------------------------------------------------------\n");

  size=40;
  printf("%6d %12.4g       - ( - )   %12.4g        - ( - )\n",size,full_norms[0],interior_norms[0]);
  size *= 2;

  double tol = 1e-12;
  bool all_conv = true;
  for(int i=1; i<n_tests; i++) {
    if(tol < interior_norms[i]){
      double full_conv = full_norms[i-1] / full_norms[i];
      double full_p = std::log2(full_conv);
      double interior_conv = interior_norms[i-1] / interior_norms[i];
      double interior_p = std::log2(interior_conv);
      printf("%6d %12.4g %7.3g (%.1f)   %12.4g %8.3g (%.1f)\n",size,full_norms[i],full_conv,full_p,
       interior_norms[i],interior_conv,interior_p);
      size *= 2;

      all_conv = all_conv && (interior_p > interior_p_ref && full_p > full_p_ref);
    }
  }
  printf("\n");
  if(all_conv) {
    printf("Proper convergence order attained\n");
  }
  else {
    printf("Proper convergence NOT attained\n");
  }
  return all_conv;
}

int main(int argc, char *argv[])
{
  printf("Second order Upwind:\n");
  bool all_conv = test_operator(stenseal::upwind_operator_2nd_order(),1.9,1.4,"1st_derivative_2nd_order.txt");

  printf("\n Kalles Second order Upwind:\n");
  all_conv = test_operator(stenseal::upwind_operator_2nd_order(),1.9,1.4,"1st_derivative_2nd_order_kalle.txt");

  printf("\n Third order Upwind:\n");
  all_conv = test_operator(stenseal::upwind_operator_3rd_order(),2.9,2.4,"1st_derivative_3rd_order.txt");

  printf("\n Fourth order Upwind:\n");
  all_conv = test_operator(stenseal::upwind_operator_4th_order(),3.8,2.4,"1st_derivative_4th_order.txt");

  return 0;
}
