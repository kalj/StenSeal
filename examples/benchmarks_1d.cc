
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/timer.h>

#include "stenseal/operator.h"
#include "stenseal/upwind_laplace.h"
#include "stenseal/compact_laplace.h"
#include "stenseal/operator_lib.h"

#define NREPS 1000

template <int dim, typename Geometry>
struct Wrapper;

template <int dim>
struct Wrapper<dim,stenseal::CartesianGeometry<dim>>
{
  static stenseal::CartesianGeometry<dim> make_geometry(std::array<unsigned int, dim> n_nodes)
  {
    return stenseal::CartesianGeometry<dim>(n_nodes);
  }
};

template <int dim>
struct Wrapper<dim,stenseal::GeneralGeometry<dim>>
{
  static stenseal::GeneralGeometry<dim> make_geometry(std::array<unsigned int, dim> n_nodes)
  {
    // use cartesian geometry to set up initial points
    stenseal::CartesianGeometry<dim> g(n_nodes);

    dealii::Vector<double> xpos(g.get_n_nodes_total());
    std::vector<dealii::Point<dim>> nodes(g.get_n_nodes_total());

    g.initialize_vector(xpos,[] (const dealii::Point<dim> &p)
                        { const double x= p(0);
                          return (exp(x)-1)/(exp(1)-1);
                        });

    for(int i = 0; i < n_nodes[0]; ++i) {
      nodes[i](0) = xpos[i];
    }


    return stenseal::GeneralGeometry<dim>(n_nodes,nodes);
  }
};

template <typename OpT>
double run_matrix_free_benchmark(const OpT &op, const unsigned int n_nodes_tot)
{
  dealii::Timer timer;

  dealii::Vector<double> u(n_nodes_tot);
  dealii::Vector<double> v(n_nodes_tot);

  for(auto &uval : u) {
    uval = 1.0;
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

  return 1.0e3*timer.wall_time()/NREPS;
}

// template <typename OpT>
// double run_matrix_based_benchmark(const OpT &op, const unsigned int n_nodes_tot)
// {
//   dealii::Timer timer;

//   dealii::Vector<double> u(n_nodes_tot);
//   dealii::Vector<double> v(n_nodes_tot);

//   dealii::SparseMatrix<double> spmat;
//   op.matrix(spmat);

//   for(auto &uval : u) {
//     uval = 1.0;
//   }

//   // warm up
//   for(int r = 0; r < 10; ++r) {
//     spmat.vmult(v,u);
//     u.swap(v);
//   }

//   timer.restart();

//   for(int r = 0; r < NREPS; ++r) {
//     spmat.vmult(v,u);
//     u.swap(v);
//   }
//   timer.stop();

//   return 1.0e3*timer.wall_time()/NREPS;

// }


template <bool use_matrix, typename OpT>
double run_benchmark(const OpT &op, const unsigned int n_nodes_tot)
{
  // if(!use_matrix)
    return run_matrix_free_benchmark(op,n_nodes_tot);
  // else
    // return run_matrix_based_benchmark(op,n_nodes_tot);
}

template <int dim, typename Geometry, bool use_matrix, typename OperatorType>
void benchmark_upwind_operator(const OperatorType &Dm,
                               const unsigned int minsize, const unsigned int maxsize)
{
  printf("%12s %17s\n","n_dofs","t_apply (ms)");

  for(unsigned int n=minsize; n<=maxsize; n*=2) {

    std::array<unsigned int,dim> n_nodes{ n };

    const Geometry geometry = Wrapper<dim,Geometry>::make_geometry(n_nodes);

    stenseal::UpwindLaplace<dim,OperatorType,Geometry> op(Dm,geometry);

    double t = run_benchmark<use_matrix>(op,geometry.get_n_nodes_total());

    printf("%12d %17.5g\n",n,t);
  }
}

template <int dim, typename Geometry, bool use_matrix, typename OperatorTypeD2, typename OperatorTypeD1>
void benchmark_compact_operator(const std::pair<OperatorTypeD2,OperatorTypeD1> &ops,
                                const unsigned int minsize, const unsigned int maxsize)
{
  printf("%12s %17s\n","n_dofs","t_apply (ms)");

  for(unsigned int n=minsize; n<=maxsize; n*=2) {

    std::array<unsigned int,dim> n_nodes{ n };

    const Geometry geometry = Wrapper<dim,Geometry>::make_geometry(n_nodes);

    stenseal::CompactLaplace<dim,OperatorTypeD2,OperatorTypeD1,Geometry> op(ops.first,ops.second,geometry);

    double t = run_benchmark<use_matrix>(op,geometry.get_n_nodes_total());

    printf("%12d %17.5g\n",n,t);
  }
}

template <int dim,typename Geometry, bool use_matrix>
void all_benchmarks()
{
  const unsigned int minsize = 1<<19; //512k
  const unsigned int maxsize = 1<<23; // 8M


  printf("Second order Upwind\n");
  benchmark_upwind_operator<dim,Geometry,use_matrix>(stenseal::upwind_operator_2nd_order(),
                                                     minsize,maxsize);
  printf("\n");

  printf("Kalles Second order Upwind\n");
  benchmark_upwind_operator<dim,Geometry,use_matrix>(stenseal::upwind_operator_2nd_order_kalle(),
                                                     minsize,maxsize);
  printf("\n");

  printf("Third order Upwind\n");
  benchmark_upwind_operator<dim,Geometry,use_matrix>(stenseal::upwind_operator_3rd_order(),
                                                     minsize,maxsize);
  printf("\n");

  printf("Fourth order Upwind\n");
  benchmark_upwind_operator<dim,Geometry,use_matrix>(stenseal::upwind_operator_4th_order(),
                                                     minsize,maxsize);
  printf("\n");

  if(!use_matrix) {
    printf("Second order Compact\n");
    benchmark_compact_operator<dim,Geometry,use_matrix>(stenseal::compact_operators_2nd_order(),
                                                        minsize,maxsize);
    printf("\n");

    printf("Fourth order Compact\n");
    benchmark_compact_operator<dim,Geometry,use_matrix>(stenseal::compact_operators_4th_order(),
                                                        minsize,maxsize);

  }
}

int main(int argc, char *argv[])
{
  const int dim = 1;
  printf("-------------------------------------------------------- \n");
  printf("------------------ Cartesian Geometry ------------------ \n");
  printf("-------------------------------------------------------- \n");

  printf("\n");
  printf(" == Matrix-free == \n");

  all_benchmarks<dim,stenseal::CartesianGeometry<dim>,false>();

  // printf("\n");
  // printf(" == Matrix-based == \n");

  // all_benchmarks<dim,stenseal::CartesianGeometry<dim>,true>();

  printf("\n");
  printf("\n");
  printf("-------------------------------------------------------- \n");
  printf("------------------- General Geometry ------------------- \n");
  printf("-------------------------------------------------------- \n");

  printf("\n");
  printf(" == Matrix-free == \n");

  all_benchmarks<dim,stenseal::GeneralGeometry<dim>,false>();

  // printf("\n");
  // printf(" == Matrix-based == \n");

  // all_benchmarks<dim,stenseal::GeneralGeometry<dim>,true>();

}
