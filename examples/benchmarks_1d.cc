
#include <deal.II/lac/vector.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/timer.h>

#include "stenseal/operator.h"
#include "stenseal/upwind_laplace.h"
#include "stenseal/operator_lib.h"

#define NREPS 1000

template <typename OperatorType>
void benchmark_operator(OperatorType Dm)
{
  dealii::Timer timer;
  const int dim = 1;
  typedef stenseal::CartesianGeometry<dim> Geometry;

  printf("%12s %17s\n","n_dofs","t_apply (ms)");

  for(unsigned int n=1<<19; n<(1<<23); n<<=1) {

    std::array<unsigned int,dim> n_nodes{ n };
    int n_nodes_tot = n_nodes[0];
    double h = 1.0/(n_nodes_tot-1);

    Geometry geometry(n_nodes, dealii::Point<dim>(0.0),
                      dealii::Point<dim>(1.0));

    stenseal::UpwindLaplace<dim,OperatorType,Geometry> op(Dm,geometry);

    dealii::Vector<double> u(n_nodes_tot);
    dealii::Vector<double> v(n_nodes_tot);

    for(auto &uval : u) {
      uval = (double)rand()/RAND_MAX;
    }

    // warm up
    for(int r = 0; r < 10; ++r) {
      op.apply(v,u);
      u.swap(v);
    }

    timer.restart();

    for(int r = 0; r < NREPS; ++r) {
      op.apply(v,u);
      u.swap(v);
    }
    timer.stop();
    double t = 1.0e3*timer.wall_time()/NREPS;

    printf("%12d %17.5g\n",n,t);
  }
}

int main(int argc, char *argv[])
{
  printf("Second order Upwind\n");
  benchmark_operator(stenseal::upwind_operator_2nd_order());
  printf("\n");

  printf("Kalles Second order Upwind\n");
  benchmark_operator(stenseal::upwind_operator_2nd_order_kalle());
  printf("\n");

  printf("Third order Upwind\n");
  benchmark_operator(stenseal::upwind_operator_3rd_order());
  printf("\n");

  printf("Fourth order Upwind\n");
  benchmark_operator(stenseal::upwind_operator_4th_order());
  printf("\n");

}
