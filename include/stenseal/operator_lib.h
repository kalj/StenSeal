  #ifndef _OPERATOR_LIB_H
  #define _OPERATOR_LIB_H

  #include "stenseal/stencil.h"
  #include "stenseal/operator.h"
  #include "stenseal/stencil_array.h"


namespace stenseal
{
  constexpr Operator<2,2,1,2,1> upwind_operator_2nd_order_kalle()
  {
    const stenseal::Symbol usym;

    constexpr stenseal::Stencil<2> interior((-1.0)*usym[-1] + 1.0*usym[0]);
    constexpr stenseal::Stencil<2> left_boundary((-1.0)*usym[0] + 1.0*usym[1]);
    constexpr stenseal::Stencil<2> right_boundary((-1.0)*usym[-1] + 1.0*usym[0]);

    constexpr stenseal::StencilArray<2,1> left_boundary_block(left_boundary);
    constexpr stenseal::StencilArray<2,1> right_boundary_block(right_boundary);

    return stenseal::Operator<2,2,1,2,1> (interior, left_boundary_block, right_boundary_block);
  }

  constexpr Operator<3,2,2,4,2> upwind_operator_2nd_order()
  {
    const stenseal::Symbol usym;

    constexpr stenseal::Stencil<3> interior((0.5)*usym[-2] + (-2.0)*usym[-1] + (1.5)*usym[0]);
    constexpr stenseal::StencilArray<2,2> left_boundary_block((-1.0)*usym[0] + 1.0*usym[1], (-1.0)*usym[-1] + 1.0*usym[0]);
    constexpr stenseal::StencilArray<4,2> right_boundary_block(0.4*usym[-2] + (-1.6)*usym[-1] + 1.0*usym[0] + 0.2*usym[1],
                                                               0.0*usym[-3] + 2.0*usym[-2] + (-5.0)*usym[-1] +3.0*usym[0]);

    return stenseal::Operator<3,2,2,4,2>(interior, left_boundary_block, right_boundary_block);
  }

}

  #endif /* _OPERATOR_LIB_H */