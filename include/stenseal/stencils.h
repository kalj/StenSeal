#ifndef _STENCILS_H
#define _STENCILS_H


namespace stenseal
{
  namespace stencils {
    namespace second_order{
struct StencilDm {
  // boundary
  const static unsigned int boundary_width = 2;
  const static unsigned int boundary_height = 1;
  const static double boundary_stencil[];

  // interior
  const static unsigned interior_width = 3;
  const static int interior_start_idx = -1;
  const static double interior_stencil[];
};

const double StencilDm::boundary_stencil[] = { -1.0, 1.0 };
const double StencilDm::interior_stencil[] = { -0.5, 0.0, 0.5};


struct StencilDp {
  // boundary
  const static unsigned int boundary_width = 2;
  const static unsigned int boundary_height = 1;
  const static double boundary_stencil[];

  // interior
  const static unsigned interior_width = 3;
  const static int interior_start_idx = -1;
  const static double interior_stencil[];
};

const double StencilDp::boundary_stencil[] = { -1.0, 1.0 };
const double StencilDp::interior_stencil[] = {-0.5, 0.0, 0.5};
}
namespace third_order {

}

}
}

#endif /* _STENCILS_H */