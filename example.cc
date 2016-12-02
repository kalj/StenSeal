
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>

#include "stenseal/fd_element.h"

using namespace dealii;
using namespace stenseal;

int main(int argc, char *argv[])
{
  const int dim = 2;
  const int nsubdiv = 10;

  // Triangulation<dim> tria;
  // DoFHandler<dim> dof(tria);

  FDElement<dim> elem(nsubdiv);

  // dof.distribute_dofs(elem);

  // std::cout << "Number of unknowns on grid: " << dof.n_dofs() << std::endl;

  return 0;
}
