
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

#include "stenseal/geometry.h"
#include "stenseal/operator.h"
#include "stenseal/block_operator.h"


int main(int argc, char *argv[])
{
  const int dim = 2;

  double xmin = 0;
  double xmax = 1.0;
  double ymin = 0;
  double ymax = 1.0;

  unsigned int n_nodes[dim] = { 200, 200 };
  dealii::Point<dim> lower_left(xmin,ymin);
  dealii::Point<dim> upper_right(xmax,ymax);

  typedef stenseal::CartesianGeometry<dim> Geometry;
  Geometry geometry(n_nodes,lower_left,upper_right);

  dealii::Vector<double> u;

  // TODO: initialize u

  //=============================================================================
  // upwind 2nd order operators
  //=============================================================================
  typedef stenseal::Operator<2,2,1> OperatorType;
  const stenseal::Symbol usym;

  // now define operator like this:
  {
    constexpr stenseal::Stencil<2> interior((-1.0)*usym[-1] + 1.0*usym[0]);
    constexpr stenseal::Stencil<2> left_boundary((-1.0)*usym[0] + 1.0*usym[1]);
    constexpr stenseal::Stencil<2> right_boundary((-1.0)*usym[-1] + 1.0*usym[0]);

    constexpr stenseal::BlockStencil<2,1> left_boundary_block(left_boundary);
    constexpr stenseal::BlockStencil<2,1> right_boundary_block(right_boundary);
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

  stenseal::UpwindBlockOperator<dim,OperatorType,Geometry> op(Dm,geometry);

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
  // - Some sensible initial condition
  // - Verification / outputting of solution
  //
  // In the long term:
  // - A more general way of setting up the domain
  // - A nice way of outputting the solution to .vtk
  //
  // Both of these should probably be done using dealii functionality

  return 0;
}
