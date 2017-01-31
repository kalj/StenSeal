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

  constexpr Operator<4,4,3,5,3> upwind_operator_3rd_order()
  {
    const stenseal::Symbol usym;

    constexpr stenseal::Stencil<4> interior((1.0/6.0)*usym[-2] + (-1.0)*usym[-1] + (1.0/2.0)*usym[0] + (1.0/3.0)*usym[1]);
    constexpr stenseal::StencilArray<4,3> left_boundary_block((-11.0/9.0)*usym[0] + (13.0/9.0)*usym[1] + (-2.0/9.0)*usym[2] + 0*usym[3],
                                                   (-17.0/28.0)*usym[-1] + (3.0/14.0)*usym[0] + (11.0/28.0)*usym[1] + 0*usym[2],
                                                   (1.0/23.0)*usym[-2] + (-1)*usym[-1] + (11.0/23.0)*usym[0] + (8.0/23.0)*usym[1]);
    constexpr stenseal::StencilArray<5,3> right_boundary_block((4.0/23.0)*usym[-2] + (-24.0/23.0)*usym[-1] + (11.0/23.0)*usym[0] + (11.0/23.0)*usym[1] + (-2.0/23.0)*usym[2],
                                                           0*usym[-3] + (1.0/7.0)*usym[-2] + (-23.0/28.0)*usym[-1] + (3.0/14.0) *usym[0] + (13.0/28.0)*usym[1],
                                                           0*usym[-4] + 0*usym[-3] + (4.0/9.0)*usym[-2] + (-17.0/9.0)*usym[-1] + (13.0/9.0)*usym[0]);
    return stenseal::Operator<4,4,3,5,3>(interior, left_boundary_block, right_boundary_block);
  }

   constexpr Operator<5,5,4,7,4> upwind_operator_4th_order()
  {
    const stenseal::Symbol usym;

    constexpr stenseal::Stencil<5> interior((-1.0/12.0)*usym[-3] + (1.0/2.0)*usym[-2] + (-3.0/2.0)*usym[-1] + (5.0/6.0)*usym[0] + (1.0/4.0)*usym[1]);
    constexpr stenseal::StencilArray<5,4> left_boundary_block((-69.0/49.0)*usym[0] + (169.0/98.0)*usym[1] + (-11.0/49.0)*usym[2]  + (-9.0/98.0)*usym[3] + 0*usym[4],
                                                              (-205.0/366.0)*usym[-1] +  (11.0/61.0)*usym[0] + (39.0/122.0)*usym[1] + (11.0/183.0)*usym[2] + 0*usym[3],
                                                              (29.0/123.0)*usym[-2] + (-99.0/82.0)*usym[-1] + (29.0/41.0)*usym[0] + (65.0/246.0)*usym[1] + 0*usym[2],
                                                              (-3.0/298.0)*usym[-3] + (43.0/149.0)*usym[-2] + (-389.0/298.0)*usym[-1] + (117.0/149.0)*usym[0] + (36.0/149.0)*usym[1]);

    constexpr stenseal::StencilArray<7,4> right_boundary_block((-12.0/149.0)*usym[-3] + (72.0/149.0)*usym[-2] + (-216.0/149.0)*usym[-1] + (117.0/149.0)*usym[0] + (65.0/298.0)*usym[1] + (11.0/149.0)*usym[2] + (-9.0/298.0)*usym[3],
                                                                0*usym[-4] + (-4/41)*usym[-3] +  (24.0/41.0)*usym[-2] + (389.0/246.0)*usym[-1] + (29.0/41.0)*usym[0] + (39.0/82.0)*usym[1] + (-11.0/123.0)*usym[2],
                                                                0*usym[-5] +  0*usym[-4] + (-4.0/61.0)*usym[-3] + (43.0/183.0)*usym[-2] + (-99.0/122.0)*usym[-1] + (11.0/61.0)*usym[0] + (169.0/366.0)*usym[1],
                                                                0*usym[-6] + 0*usym[-5] + 0*usym[-4] + (-3.0/98.0)*usym[-3] + (29.0/49.0)*usym[-2] + (-205.0/98.0)*usym[-1] + (75.0/49.0)*usym[0]);
    return stenseal::Operator<5,5,4,7,4>(interior, left_boundary_block, right_boundary_block);
  }



}

  #endif /* _OPERATOR_LIB_H */