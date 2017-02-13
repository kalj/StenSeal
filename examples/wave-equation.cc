
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>
#include <deal.II/base/function.h>

#include "stenseal/geometry.h"
#include "stenseal/operator.h"
#include "stenseal/upwind_laplace.h"


template <int dim>
class InitialValues : public dealii::Function<dim>
{

public:
  InitialValues () : dealii::Function<dim>() {}
  virtual ~InitialValues(){};

  virtual double value (const dealii::Point<dim>   &p,
                        const unsigned int  component = 0) const;
};

template <int dim>
double InitialValues<dim>::value (const dealii::Point<dim> &p,
                                   const unsigned int /*component*/) const
{

  double radius = 0.3;
  dealii::Point<dim> center;
  center[0] = 0.2;

  double res = 0;
  double dist = (p-center).norm();

  if( dist < radius)
    res = (1+cos(dealii::numbers::PI*dist/radius))/2;
  return res;
}

template <int dim, typename Geometry>
void initialize( dealii::Vector<double> &u, dealii::Function<dim> &f, const Geometry &g )
{
  if(dim==1) {
    for(int i = 0; i < g.get_n_nodes(0); ++i) {
      double x;
      x = i*g.get_h(0) + g.get_lower_left(0);
      dealii::Point<dim> p(x);
      u[i] = f.value(p);
    }
  } else if(dim==2){
    for(int i = 0; i < g.get_n_nodes(1); ++i) {
      for(int j = 0; j < g.get_n_nodes(0); ++j) {
        double x;
        double y;
        x = j*g.get_h(0) + g.get_lower_left(0);
        y = i*g.get_h(1) + g.get_lower_left(1);
        dealii::Point<dim> p(x,y);
        u[g.get_n_nodes(0)*i+j] = f.value(p);
      }
    }
  }
  else  {
    AssertThrow(false,dealii::ExcNotImplemented());
  }
}

int main(int argc, char *argv[])
{
  // const int dim = 1;

  // double xmin = 0;
  // double xmax = 1.0;

  // unsigned int n_nodes[dim] = {200};
  // dealii::Point<dim> lower_left(xmin);
  // dealii::Point<dim> upper_right(xmax);

  const int dim = 2;

  double xmin = 0;
  double xmax = 1.0;
  double ymin = 0;
  double ymax = 1.0;

  unsigned int n_nodes[dim] = { 10, 10 };
  dealii::Point<dim> lower_left(xmin,ymin);
  dealii::Point<dim> upper_right(xmax,ymax);

  typedef stenseal::CartesianGeometry<dim> Geometry;
  Geometry geometry(n_nodes,lower_left,upper_right);

  //=============================================================================
  // set up initial condition
  //=============================================================================
  dealii::Vector<double> u(geometry.get_n_nodes_total());
  InitialValues<dim> f;
  initialize<dim,Geometry>(u,f,geometry);

  u.print(std::cout);

  //=============================================================================
  // upwind 2nd order operators
  //=============================================================================
  typedef stenseal::Operator<2,2,1,2,1> OperatorType;
  const stenseal::Symbol usym;

  // now define operator like this:
  {
    constexpr stenseal::Stencil<2> interior((-1.0)*usym[-1] + 1.0*usym[0]);
    constexpr stenseal::Stencil<2> left_boundary((-1.0)*usym[0] + 1.0*usym[1]);
    constexpr stenseal::Stencil<2> right_boundary((-1.0)*usym[-1] + 1.0*usym[0]);

    constexpr stenseal::StencilTensor2D<1,2> left_boundary_block(left_boundary);
    constexpr stenseal::StencilTensor2D<1,2> right_boundary_block(right_boundary);
    constexpr OperatorType Dm(interior, left_boundary_block, right_boundary_block);
  }

  // or simpler:
  {
    constexpr stenseal::Stencil<2> interior((-1.0)*usym[-1] + 1.0*usym[0]);
    constexpr stenseal::Stencil<2> left_boundary((-1.0)*usym[0] + 1.0*usym[1]);
    constexpr stenseal::Stencil<2> right_boundary((-1.0)*usym[-1] + 1.0*usym[0]);
    constexpr OperatorType Dm(interior, left_boundary, right_boundary);
  }

  // or even:
  constexpr OperatorType Dm((-1.0)*usym[-1] + 1.0*usym[0],  // interior stencil
                            (-1.0)*usym[0] + 1.0*usym[1],   // left boundary
                            (-1.0)*usym[-1] + 1.0*usym[0]); // right boundary

  stenseal::UpwindLaplace<dim,OperatorType,Geometry> op(Dm,geometry);


  // TODO:
  //
  // Write code here for solving a wave equation in 2D
  //
  // We need (at least):
  //
  // - A time stepping scheme, either something simple, e.g. our own Euler fw,
  //   or (probably better) something good from dealii::TimeStepping (see
  //   [DEALSRC]/include/deal.II/base/time_stepping.h)
  // - Some sketch of how to apply the SAT, telling us what is needed for that
  // - Verification / outputting of solution
  //
  // In the long term:
  // - A more general way of setting up the domain
  // - A nice way of outputting the solution to .vtk
  //
  // Both of these should probably be done using dealii functionality

  return 0;
}
