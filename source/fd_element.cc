#include "stenseal/fd_element.h"

using namespace stenseal;

template <int dim>
FDElement<dim>::FDElement(int n_subdiv)
  : n_subdivisions(n_subdiv)
{

}

// explicit instantiation
template class FDElement<2> ;
