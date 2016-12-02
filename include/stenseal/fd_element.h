#ifndef stenseal__fd_element_h
#define stenseal__fd_element_h

#include <deal.II/fe/fe.h>

namespace stenseal
{
  template <int dim>
  class FDElement
// : public dealii::FiniteElement<dim>
  {
  public:

    FDElement(int nsubdiv);


  private:

    /**
     * How many times is this block subdivided
     */
    int n_subdivisions;
  };
}

#endif /* stenseal__fd_element_h */
