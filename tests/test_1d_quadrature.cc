
#include <deal.II/lac/vector.h>
#include <deal.II/base/numbers.h>

#include "stenseal/operator.h"
#include "stenseal/upwind_laplace.h"
#include "stenseal/operator_lib.h"

double qref_up_2[] = {  0.2500,    1.2500,    1.0000,    1.0000,  1.0000 , 1.0000,  1.0000,    1.0000,    1.000000000000000,   1.000000000000000,   1.000000000000000,   1.000000000000000,  1.0000,    1.0000,    1.2500,    0.2500};
double qref_up_3[] = {  0.375000000000000,   1.166666666666667,   0.958333333333333,   1.000000000000000,   1.000000000000000,   1.000000000000000,   1.000000000000000,   1.000000000000000,   1.000000000000000,   1.000000000000000,  1.000000000000000,   1.000000000000000,   1.000000000000000,   0.958333333333333,   1.166666666666667,   0.375000000000000 };
double qref_up_4[] = {  0.340277777777778,   1.270833333333333,   0.854166666666667,   1.034722222222222,   1.000000000000000,   1.000000000000000,   1.000000000000000,   1.000000000000000,   1.000000000000000,   1.000000000000000,  1.000000000000000,   1.000000000000000,   1.034722222222222,   0.854166666666667,   1.270833333333333,   0.340277777777778};
double qref_up_6[] = {  0.315115740740741,   1.394560185185185,   0.619212962962963,   1.248842592592593,   0.907523148148148,   1.014745370370370,   1.000000000000000,   1.000000000000000,   1.000000000000000,   1.000000000000000,  1.014745370370370,   0.907523148148148,   1.248842592592593,   0.619212962962963,   1.394560185185185,   0.315115740740741};

double qref_comp_2[] = {  0.500000000000000,   1.000000000000000,   1.000000000000000,   1.000000000000000,   1.000000000000000,   1.000000000000000,   1.000000000000000,   1.000000000000000,   1.000000000000000,   1.000000000000000,  1.000000000000000,   1.000000000000000,   1.000000000000000,   1.000000000000000,   1.000000000000000,   0.500000000000000};
double qref_comp_4[] = {  0.272727272727273,   1.621288167980855,   0.141741267470520,   1.744774462430583,   0.653104903834157,   1.066363925556612,   1.000000000000000,   1.000000000000000,   1.000000000000000,   1.000000000000000,  1.066363925556612,  0.653104903834157,  1.744774462430583,   0.141741267470520,   1.621288167980855,   0.272727272727273};
double qref_comp_6[] = {  0.295013975497606,   1.525036100088183,   0.259327876984127,   1.794691082451499,   0.417023533950617,   1.275002480158730,   0.924872960758377,    1.009031990110859,   1.009031990110859,   0.924872960758377,   1.275002480158730,   0.417023533950617,   1.794691082451499,   0.259327876984127,   1.525036100088183,   0.295013975497606};

template <typename OperatorType, typename QuadratureType>
bool test_quad_up(std::pair<OperatorType,QuadratureType> op, int order)
{
 bool all_conv = false;
 int n = 16;
 dealii::Vector<double> u(n);

for(int i = 0; i <n; ++i){
   u[i] = 1;
  }

dealii::Vector<double> v(n);
op.second.apply(v,u,n);

double diff =0;

if(order == 2){
 for(int i = 0; i <n; ++i){
   diff = diff + fabs( v[i] - qref_up_2[i] );
 }
}else if(order == 3){
 for(int i = 0; i <n; ++i){
   diff = diff + fabs( v[i] - qref_up_3[i] );
 }
}else if(order == 4){
 for(int i = 0; i <n; ++i){
   diff = diff + fabs( v[i] - qref_up_4[i] );
 }
}else if(order == 6){
  for(int i = 0; i <n; ++i){
    diff = diff + fabs( v[i] - qref_up_6[i] );
  }
}

 double tol = 1e-13;
 if( diff < tol && diff > -tol){
  all_conv = true;
  printf("OK\n");
}else{
  printf("NOT OK\n");
}

return all_conv;
}

template <typename OperatorType1, typename OperatorType2, typename QuadratureType>
bool test_quad_comp(std::tuple<OperatorType1, OperatorType2, QuadratureType> op, int order)
{
 bool all_conv = false;
 int n = 16;
 dealii::Vector<double> u(n);

for(int i = 0; i <n; ++i){
   u[i] = 1;
  }

dealii::Vector<double> v(n);
auto  quad = std::get<2>(op);
quad.apply(v,u,n);


double diff =0;

if(order == 2){
 for(int i = 0; i <n; ++i)
   diff = diff + fabs( v[i] - qref_comp_2[i] );
}else if(order == 4){
 for(int i = 0; i <n; ++i){
   diff = diff + fabs( v[i] - qref_comp_4[i] );
 }
 }else if(order == 6){
   for(int i = 0; i <n; ++i){
   diff = diff + fabs( v[i] - qref_comp_6[i] );
 }

}

 double tol = 1e-13;
 if( diff < tol && diff > -tol){
  all_conv = true;
  printf("OK\n");
}else{
  printf("NOT OK\n");
}

return all_conv;
}

int main(int argc, char *argv[])
{
  bool all_conv = true;

  printf("Second order Upwind: \n");
  all_conv = (test_quad_up(stenseal::upwind_operator_2nd_order(),2) && all_conv);

  printf("Third order Upwind: \n");
  all_conv = (test_quad_up(stenseal::upwind_operator_3rd_order(),3) && all_conv);

  printf("Fourth order Upwind: \n");
  all_conv = (test_quad_up(stenseal::upwind_operator_4th_order(),4) && all_conv);

  printf("Sixth order Upwind: \n");
  all_conv = (test_quad_up(stenseal::upwind_operator_6th_order(),6) && all_conv);

  printf("Second order Compact: \n");
  all_conv = (test_quad_comp(stenseal::compact_operators_2nd_order(),2) && all_conv);


  printf("Fourth order Compact: \n");
  all_conv = (test_quad_comp(stenseal::compact_operators_4th_order(),4) && all_conv);

  printf("Sixth order Compact: \n");
  all_conv = (test_quad_comp(stenseal::compact_operators_6th_order(),6) && all_conv);

  if(all_conv) {
    printf("\nAll quadratures are right\n");
    return 0;
  }
  else {
    printf("\nOne or more quadratures have the wrong coefficinets\n");
    return 1;
  }
}
