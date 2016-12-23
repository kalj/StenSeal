
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

#include "stenseal/geometry.h"
#include "stenseal/block_operator.h"


//=============================================================================
// upwind 2nd order operators
//=============================================================================

struct StencilDm {
  // boundary
  const static unsigned int boundary_width = 2;
  const static unsigned int boundary_height = 1;
  const static double boundary_stencil[];

  // interior
  // const static unsigned interior_width = 3;
  const static unsigned interior_width = 2;
  const static int interior_start_idx = -1;
  const static double interior_stencil[];
};

const double StencilDm::boundary_stencil[] = { -1.0, 1.0 };
// const double StencilDm::interior_stencil[] = { -0.5, 0.0, 0.5};
const double StencilDm::interior_stencil[] = { -1.0, 1.0 };


struct StencilDp {
  // boundary
  const static unsigned int boundary_width = 2;
  const static unsigned int boundary_height = 1;
  const static double boundary_stencil[];

  // interior
  // const static unsigned interior_width = 3;
  // const static int interior_start_idx = -1;
  const static unsigned interior_width = 2;
  const static int interior_start_idx = 0;
  const static double interior_stencil[];
};

const double StencilDp::boundary_stencil[] = { -1.0, 1.0 };
// const double StencilDm::interior_stencil[] = { -0.5, 0.0, 0.5};
const double StencilDp::interior_stencil[] = {-1.0, 1.0};




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

  stenseal::CartesianGeometry<dim> geometry(n_nodes,lower_left,upper_right);

  dealii::Vector<double> u;

  // TODO: initialize u

  stenseal::UpwindBlockOperator<dim,StencilDm,StencilDp> op(geometry);

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
