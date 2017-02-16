
#include <deal.II/lac/vector.h>
#include <deal.II/base/numbers.h>

#include "stenseal/geometry.h"
#include "stenseal/operator.h"
#include "stenseal/compact_laplace.h"
#include "stenseal/metric_operator.h"
#include "stenseal/operator_lib.h"

#include <fstream>

template <typename OperatorTypeD2, typename OperatorTypeD1>
void compute_l2_norm(std::pair<OperatorTypeD2,OperatorTypeD1> ops, unsigned int n,
                     double &l2_norm, double &l2_norm_interior)
{
  const int dim = 1;

  std::array<unsigned int,dim> n_nodes{ n };
  int n_nodes_tot = n_nodes[0];
  double h = 1.0/(n_nodes_tot-1);
  const double PI = dealii::numbers::PI;

  // node positions
  std::vector<dealii::Point<dim>> nodes(n_nodes_tot);

  {
    // use cartesian geometry to set up initial points
    stenseal::CartesianGeometry<dim> g(n_nodes);

    dealii::Vector<double> xpos(n_nodes_tot);
    g.initialize_vector(xpos,[] (const dealii::Point<dim> &p)
                               { const double x= p(0);
                                 return (exp(x)-1)/(exp(1)-1);
                                });

    for(int i = 0; i < n_nodes[0]; ++i) {
      nodes[i](0) = xpos[i];
    }
  }

  // std::ofstream node_file("nodes.txt");
  // for(auto p : nodes)
  //   node_file << p << std::endl;

  typedef stenseal::GeneralGeometry<dim> Geometry;
  Geometry geometry(n_nodes,nodes);

  stenseal::CompactLaplace<dim,OperatorTypeD2,OperatorTypeD1,Geometry> op(ops.first,ops.second,geometry);

  dealii::Vector<double> u(n_nodes_tot);

  auto test_function =
    [=](const dealii::Point<dim> &p) { return sin(PI*p(0)); };
  auto test_function_2nd_derivative =
    [=](const dealii::Point<dim> &p) { return -PI*PI*sin(PI*p(0)); };

  geometry.initialize_vector(u,test_function);

  dealii::Vector<double> v(n_nodes_tot);

  op.apply(v,u);


  dealii::Vector<double> vref(n_nodes_tot);
  geometry.initialize_vector(vref,test_function_2nd_derivative);

  // std::ofstream u_file("initial_condition.txt");
  // u.print(u_file);

  // std::ofstream v_file("2nd_derivative.txt");
  // v.print(v_file);

  // std::ofstream vref_file("ref_2nd_derivative.txt");
  // vref.print(vref_file);

  // exclude points affected by boundary stencil
  const int height_bdry = 2;
  const int bdry_offset = height_bdry;

  // compute norms
  double sqsum = 0;
  double a;
  for(int i = bdry_offset; i < n_nodes_tot-bdry_offset; ++i) {
    a = v[i] - vref[i];
    sqsum += a*a;
  }

  l2_norm_interior = std::sqrt(h*sqsum);

  for(int i= 0; i < bdry_offset; ++i) {
    a = v[i] - vref[i];
    sqsum += a*a;
  }


  for(int i = n_nodes_tot-bdry_offset; i < n_nodes_tot; ++i) {
    a = v[i] - vref[i];
    sqsum += a*a;
  }

  l2_norm = std::sqrt(h*sqsum);
}


template <typename OperatorTypeD2, typename OperatorTypeD1>
bool test_operator(std::pair<OperatorTypeD2,OperatorTypeD1> ops, float interior_p_ref, float full_p_ref)
{
  const int n_tests = 7;
  double full_norms[n_tests];
  double interior_norms[n_tests];
  unsigned int size = 40;

  for(int i=0; i<n_tests; i++) {
    compute_l2_norm(ops,size,full_norms[i],interior_norms[i]);
    size *= 2;
  }

  printf("     n   full error   factor  (p)   inter. error   factor   (p)\n");
  printf(" ---------------------------------------------------------------\n");

  size=40;
  printf("%6d %12.4g       - ( - )   %12.4g        - ( - )\n",size,full_norms[0],interior_norms[0]);
  size *= 2;

  bool all_conv = true;
  double tol = 1e-08;
  for(int i=1; i<n_tests; i++) {
    double full_conv = full_norms[i-1] / full_norms[i];
    double full_p = std::log2(full_conv);
    double interior_conv = interior_norms[i-1] / interior_norms[i];
    double interior_p = std::log2(interior_conv);
    printf("%6d %12.4g %7.3g (%.1f)   %12.4g %8.3g (%.1f)\n",size,full_norms[i],full_conv,full_p,
           interior_norms[i],interior_conv,interior_p);

    size *= 2;
    all_conv = all_conv && (interior_p > interior_p_ref && full_p > full_p_ref);
    if(interior_norms[i] < tol) break;
  }
  printf("\n");
  return all_conv;
}

int main(int argc, char *argv[])
{
  bool all_conv = true;

  printf("Second order Compact:\n");

  const stenseal::Symbol sym;

  constexpr stenseal::Stencil<2> d1_interior((-0.5)*sym[-1] + 0.5*sym[1]);

  constexpr stenseal::StencilTensor2D<1,2> d1_boundary((-1.0)*sym[0] + 1.0*sym[1]);
  constexpr stenseal::StencilTensor2D<1,2> d1_boundary_r((-1.0)*sym[-1] + 1.0*sym[0]);


  constexpr stenseal::Operator<2,2,1,2,1> D1 (d1_interior,
                                              d1_boundary,
                                              d1_boundary_r);


  constexpr stenseal::StencilTensor2D<3,3> d2_interior((0.5)*sym[-1]  + (0.5)*sym[0] + 0.0*sym[1],
                                                       (-0.5)*sym[-1] + (-1.0)*sym[0]+ (-0.5)*sym[1],
                                                       0.0*sym[-1]    + (0.5)*sym[0] + (0.5)*sym[1]);


  constexpr stenseal::StencilTensor3D<2,3,3> d2_boundary( stenseal::StencilTensor2D<3,3>((2.0)*sym[0]   + (-1.0)*sym[1]  + (0.0)*sym[2],
                                                                                         (-3.0)*sym[-1] + (1.0)*sym[0]   + (0.0)*sym[1],
                                                                                         (1.0)*sym[-2]  + (0.0)*sym[-1]  + (0.0)*sym[0]),
                                                          stenseal::StencilTensor2D<3,3>((0.5)*sym[0]   + (0.5)*sym[1]   + (0.0)*sym[2],
                                                                                         (-0.5)*sym[-1] + (-1.0)*sym[0]  + (-0.5)*sym[1],
                                                                                         (0.0)*sym[-2]  + (0.5)*sym[-1] + (0.5)*sym[0]));


  constexpr stenseal::MetricOperator<3,3,2> D2 (d2_interior,
                                                d2_boundary);


  all_conv = test_operator(std::make_pair(D2,D1),1.9,1.4) && all_conv;


  if(all_conv) {
    printf("Proper convergence order attained\n");
    return 0;
  }
  else {
    printf("Proper convergence NOT attained\n");
    return 1;
  }
}
