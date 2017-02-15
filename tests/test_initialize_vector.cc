
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>
#include <deal.II/base/function.h>

#include "stenseal/geometry.h"

#include "stenseal/utils.h"

#include <fstream>


double uref_cartesian_1d[] = {
 2.641215290314e-02,    1.269288388934e+00,    2.486291455265e+00,    1.269288388934e+00,
 2.641215290314e-02};

double uref_cartesian_2d[] = {
   1.591246589136e-02,    1.509919021644e-01,    3.182493105568e-02,    1.509918998654e-01,
   1.591246564901e-02,    1.509919021644e-01,    1.432748030873e+00,    3.019651681488e-01,
   1.432571281166e+00,    1.509732705833e-01,    3.182493105568e-02,    3.019651681488e-01,
   4.773739622000e-02,    1.510105314466e-01,    1.591246589136e-02,    1.509918998654e-01,
   1.432571281166e+00,    1.510105314466e-01,    3.535648615742e-04,    1.863388106479e-05,
   1.591246564901e-02,    1.509732705833e-01,    1.591246589136e-02,    1.863388106479e-05,
   4.846930553549e-10};

double uref_cartesian_3d[] = {
   3.638778236171e-05,    1.403631607185e-04,    6.014885512534e-06,    3.105685301287e-09,
   2.555827830023e-11,    3.452376812711e-04,    1.331728728121e-03,    5.706765938111e-05,
   2.946593876840e-08,    2.424901931663e-10,    3.640021236702e-05,    1.404810933253e-04,
   6.027315517843e-06,    3.120239341189e-09,    2.555846760896e-11,    1.605383214264e-07,
   1.283261801410e-06,    1.249753154813e-07,    1.417212716799e-10,    3.172177844467e-14,
   1.243055968361e-08,    1.179347463416e-07,    1.243009710500e-08,    1.455408742395e-11,
   1.896979840284e-16,    1.008932354567e-02,    3.891881654039e-02,    1.669717720522e-03,
   1.942550594658e-05,    1.963758436941e-06,    9.572495448696e-02,    3.692525691742e-01,
   1.584195857165e-02,    1.843042069766e-04,    1.863162129792e-05,    1.019004492191e-02,
   3.987443435123e-02,    1.770439096762e-03,    1.954343855340e-05,    1.963759970926e-06,
   9.674312087184e-04,    9.112218736980e-03,    9.575728638427e-04,    1.141658366052e-06,
   2.313878774096e-09,    1.007215314341e-04,    9.556184181258e-04,    1.007214032037e-04,
   1.179329044636e-07,    1.563892544574e-12,    3.107730453838e-02,    1.198803523848e-01,
   6.814205514886e-03,    1.591491086640e-02,    1.677161527663e-03,    2.948638706112e-01,
   1.137493824344e+00,    6.466199650490e-02,    1.509964825881e-01,    1.591246552953e-02,
   4.014395388814e-02,    2.059023263362e-01,    1.588085486465e-02,    1.592552682135e-02,
   1.677161665748e-03,    8.605836304392e-02,    8.162941428631e-01,    8.602995387573e-02,
   1.193558225739e-04,    1.965064366456e-06,    9.066649961151e-03,    8.602197708730e-02,
   9.066649591625e-03,    1.061619749794e-05,    1.636280241472e-10,    1.063407258334e-03,
   4.120645731916e-03,    1.608824551940e-02,    1.509733519628e-01,    1.591246540711e-02,
   1.009994179854e-02,    3.919630712741e-02,    1.526516382035e-01,    1.432395294182e+00,
   1.509732682880e-01,    1.013005660809e-02,    9.014261968333e-02,    2.515489486916e-02,
   1.509839679178e-01,    1.591246554520e-02,    8.602322038325e-02,    8.161586021817e-01,
   8.604081266017e-02,    2.774929988696e-04,    1.863289157467e-05,    9.066649504040e-03,
   8.602197532428e-02,    9.066649732868e-03,    1.061825443527e-05,    3.804314138910e-10,
   4.042588254089e-07,    3.523062658579e-06,    1.677228335323e-03,    1.591246543848e-02,
   1.677161514730e-03,    3.953418999570e-06,    3.454466319167e-05,    1.591321730085e-02,
   1.509732687235e-01,    1.591246540667e-02,    1.011256350653e-04,    9.591408734976e-04,
   1.777949711563e-03,    1.591258337109e-02,    1.677161516264e-03,    9.556182987297e-04,
   9.066653612910e-03,    9.575816578685e-04,    1.975049493547e-05,    1.963768792413e-06,
   1.007213777801e-04,    9.556178254467e-04,    1.007214033181e-04,    1.181749551398e-07,
   2.707712044131e-11};

double uref_general_1d[] = {
   2.641215290314e-02,    1.735637235532e+00,    2.482858325164e+00,    8.768967029247e-01,
   2.641215290314e-02};

double uref_general_2d[] = {
   1.591246589136e-02,    1.509919021644e-01,    3.182493105568e-02,    1.509918998654e-01,
   1.591246564901e-02,    1.509919021644e-01,    1.163814190849e+00,    3.709492941403e-01,
   1.298355206153e+00,    1.509732705833e-01,    3.182493105568e-02,    3.016976865743e-01,
   5.487805097361e-02,    3.777759014638e-02,    1.591246589136e-02,    1.509918998654e-01,
   1.226521011200e+00,    7.309228684390e-02,    2.055930438244e-03,    1.863388106479e-05,
   1.591246564901e-02,    1.509732705833e-01,    1.591246589136e-02,    1.863388106479e-05,
   4.846930553549e-10};

double uref_general_3d[] = {
   3.638778236171e-05,    1.403631607185e-04,    6.014885512534e-06,    3.105685301287e-09,
   2.555827830023e-11,    3.452376812711e-04,    1.331728728121e-03,    5.706765938111e-05,
   2.946593876840e-08,    2.424901931663e-10,    3.640021236702e-05,    1.404810933253e-04,
   6.027315517843e-06,    3.120239341189e-09,    2.555846760896e-11,    1.605383214264e-07,
   1.283261801410e-06,    1.249753154813e-07,    1.417212716799e-10,    3.172177844467e-14,
   1.243055968361e-08,    1.179347463416e-07,    1.243009710500e-08,    1.455408742395e-11,
   1.896979840284e-16,    1.008932354567e-02,    3.891881654039e-02,    1.669717720522e-03,
   1.942550594658e-05,    1.963758436941e-06,    9.572495448696e-02,    6.209123679078e-01,
   1.835262038094e-02,    5.559040069195e-05,    1.863162129792e-05,    1.019004492191e-02,
   2.091827790048e-02,    2.493731399295e-03,    1.640146841364e-06,    1.963759970926e-06,
   9.674312087184e-04,    2.406338913371e-02,    4.389149115333e-04,    1.072512138401e-05,
   2.313878774096e-09,    1.007215314341e-04,    9.556184181258e-04,    1.007214032037e-04,
   1.179329044636e-07,    1.563892544574e-12,    3.107730453838e-02,    1.198803523848e-01,
   6.814205514886e-03,    1.591491086640e-02,    1.677161527663e-03,    2.948638706112e-01,
   1.030789024577e+00,    1.380649120191e-01,    3.463038421383e-02,    1.591246552953e-02,
   4.014395388814e-02,    1.962219874283e-01,    8.043619339473e-02,    2.877116849638e-02,
   1.677161665748e-03,    8.605836304392e-02,    7.758560924391e-01,    2.415306615669e-02,
   1.903497547298e-04,    1.965064366456e-06,    9.066649961151e-03,    8.602197708730e-02,
   9.066649591625e-03,    1.061619749794e-05,    1.636280241472e-10,    1.063407258334e-03,
   4.120645731916e-03,    1.608824551940e-02,    1.509733519628e-01,    1.591246540711e-02,
   1.009994179854e-02,    1.418607516393e-02,    2.076581535410e-01,    1.160740477343e+00,
   1.509732682880e-01,    1.013005660809e-02,    7.771554535358e-02,    5.946616678222e-02,
   9.776832847948e-02,    1.591246554520e-02,    8.602322038325e-02,    3.576249042383e-01,
   7.857560704118e-02,    1.718099963295e-04,    1.863289157467e-05,    9.066649504040e-03,
   8.602197532428e-02,    9.066649732868e-03,    1.061825443527e-05,    3.804314138910e-10,
   4.042588254089e-07,    3.523062658579e-06,    1.677228335323e-03,    1.591246543848e-02,
   1.677161514730e-03,    3.953418999570e-06,    3.454466319167e-05,    1.591321730085e-02,
   1.509732687235e-01,    1.591246540667e-02,    1.011256350653e-04,    9.591408734976e-04,
   1.777949711563e-03,    1.591258337109e-02,    1.677161516264e-03,    9.556182987297e-04,
   9.066653612910e-03,    9.575816578685e-04,    1.975049493547e-05,    1.963768792413e-06,
   1.007213777801e-04,    9.556178254467e-04,    1.007214033181e-04,    1.181749551398e-07,
   2.707712044131e-11};


template <int dim>
class TestFunction : public dealii::Function<dim>
{
private:
  static const unsigned int n_source_centers = 3;
  static const dealii::Point<dim>   source_centers[n_source_centers];
  static const double       width;

public:
  TestFunction () : dealii::Function<dim>() {}

  virtual double value (const dealii::Point<dim>   &p,
                        const unsigned int  component = 0) const;
};

template <>
const dealii::Point<1>
TestFunction<1>::source_centers[TestFunction<1>::n_source_centers]
= { dealii::Point<1>(-1.0 / 3.0),
    dealii::Point<1>(0.0),
    dealii::Point<1>(+1.0 / 3.0)   };

template <>
const dealii::Point<2>
TestFunction<2>::source_centers[TestFunction<2>::n_source_centers]
= { dealii::Point<2>(-0.5, +0.5),
    dealii::Point<2>(-0.5, -0.5),
    dealii::Point<2>(+0.5, -0.5)   };

template <>
const dealii::Point<3>
TestFunction<3>::source_centers[TestFunction<3>::n_source_centers]
= { dealii::Point<3>(-0.5, +0.5, 0.25),
    dealii::Point<3>(-0.6, -0.5, -0.125),
    dealii::Point<3>(+0.5, -0.5, 0.5)   };

template <int dim>
const double
TestFunction<dim>::width = 1./3.;

template <int dim>
double TestFunction<dim>::value (const dealii::Point<dim>   &p,
                             const unsigned int) const
{
  const double pi = dealii::numbers::PI;
  double return_value = 0;
  for (unsigned int i=0; i<this->n_source_centers; ++i)
  {
    const dealii::Tensor<1,dim> x_minus_xi = p - this->source_centers[i];
    return_value += std::exp(-x_minus_xi.norm_square() /
                             (this->width * this->width));
  }

  const double scale = std::sqrt(2 * pi) * this->width;
  return return_value /(scale*scale);
}


template <int dim, typename Geometry>
bool evaluate_and_compare(const Geometry& geometry, const double *uref)
{
  const unsigned int ntot = geometry.get_n_nodes_total();
  dealii::Vector<double> u(ntot);

  constexpr auto lambda = [] (const dealii::Point<dim> &p) {
    const double radius = 0.6;
    dealii::Point<dim> center;
    center[0] = 0.2;

    const double dist = (p-center).norm();

    const double res = dist < radius? (1+cos(dealii::numbers::PI*dist/radius))/2 : 0.0;
    return res;
  };

  // dealii::ScalarFunctionFromFunctionObject<dim> f(lambda);
  TestFunction<dim> f;

  geometry.initialize_vector(u,f);

  double diff = 0.0;
  for(int i = 0; i < ntot; ++i) {
    double a = u[i] - uref[i];
    diff += a*a;
  }
  std::cout << "  Difference: "<< diff/ntot << std::endl << std::endl;

  return diff/ntot < 1e-24;
}

template <int dim>
bool test_cartesian_geometry( const double *uref)
{
  unsigned int n = 5;

  double xmin = -1.0;
  double xmax = 1.0;

  std::array<unsigned int,dim> n_nodes = stenseal::internal::repeat_value<dim>(n);
  dealii::Point<dim> lower_left = stenseal::internal::repeat_point<dim>(xmin);
  dealii::Point<dim> upper_right = stenseal::internal::repeat_point<dim>(xmax);

  typedef stenseal::CartesianGeometry<dim> Geometry;
  Geometry geometry(n_nodes,lower_left,upper_right);

  return evaluate_and_compare<dim>(geometry,uref);
}

template <int dim>
bool test_general_geometry(const double *uref)
{
  unsigned int n = 5;

  double xmin = -1.0;
  double xmax = 1.0;

  std::array<unsigned int,dim> n_nodes = stenseal::internal::repeat_value<dim>(n);
  dealii::Point<dim> lower_left = stenseal::internal::repeat_point<dim>(xmin);


  double h[dim];
  unsigned int n_nodes_tot = 1;
  for(int i = 0; i < dim; ++i) {
    h[i] = (xmax-xmin)/(n_nodes[i]-1);
    n_nodes_tot *= n_nodes[i];
  }

  std::vector<dealii::Point<dim>> nodes(n_nodes_tot);

  if(dim==1) {
    for(int i = 0; i < n_nodes[0]; ++i) {
      double x = lower_left(0)+h[0]*i;
      nodes[i](0) = x;
    }
  }
  else if(dim==2) {
    for(int j = 0; j < n_nodes[1]; ++j) {
      for(int i = 0; i < n_nodes[0]; ++i) {
        double x = lower_left(0)+h[0]*i;
        double y = lower_left(1)+h[1]*j;
        const unsigned int idx = j*n_nodes[0] + i;
        nodes[idx](0) = x;
        nodes[idx](1) = y;
      }
    }
  }
  else if(dim==3) {
    for(int k = 0; k < n_nodes[2]; ++k) {
      for(int j = 0; j < n_nodes[1]; ++j) {
        for(int i = 0; i < n_nodes[0]; ++i) {
          double x = lower_left(0)+h[0]*i;
          double y = lower_left(1)+h[1]*j;
          double z = lower_left(2)+h[2]*k;
          const unsigned int idx = (k*n_nodes[1]+j)*n_nodes[0] + i;
          nodes[idx](0) = x;
          nodes[idx](1) = y;
          nodes[idx](2) = z;
        }
      }
    }
  }


  if(dim==1) {
    for(int i = 1; i < n_nodes[0]-1; ++i) {
      nodes[i](0) += 0.3*h[0]*(2.0*(double)rand()/RAND_MAX - 1.0);
    }
  }
  else if(dim==2) {
    for(int j = 1; j < n_nodes[1]-1; ++j) {
      for(int i = 1; i < n_nodes[0]-1; ++i) {
        const unsigned int idx = j*n_nodes[0] + i;
        nodes[idx](0) += 0.3*h[0]*(2.0*(double)rand()/RAND_MAX - 1.0);
        nodes[idx](1) += 0.3*h[1]*(2.0*(double)rand()/RAND_MAX - 1.0);
      }
    }
  }
  else if(dim==3) {
    for(int k = 1; k < n_nodes[2]-1; ++k) {
      for(int j = 1; j < n_nodes[1]-1; ++j) {
        for(int i = 1; i < n_nodes[0]-1; ++i) {
          const unsigned int idx = (k*n_nodes[1]+j)*n_nodes[0] + i;
          nodes[idx](0) += 0.3*h[0]*(2.0*(double)rand()/RAND_MAX - 1.0);
          nodes[idx](1) += 0.3*h[1]*(2.0*(double)rand()/RAND_MAX - 1.0);
          nodes[idx](2) += 0.3*h[2]*(2.0*(double)rand()/RAND_MAX - 1.0);
        }
      }
    }
  }


  typedef stenseal::GeneralGeometry<dim> Geometry;
  Geometry geometry(n_nodes,nodes);

  return evaluate_and_compare<dim>(geometry,uref);
}

int main(int argc, char *argv[])
{
  srand(0);
  bool allpass=true;

  std::cout << "Tesing CartesianGeometry<1>" << std::endl;
  allpass = test_cartesian_geometry<1>(uref_cartesian_1d) && allpass;
  std::cout << "Tesing CartesianGeometry<2>" << std::endl;
  allpass = test_cartesian_geometry<2>(uref_cartesian_2d) && allpass;
  std::cout << "Tesing CartesianGeometry<3>" << std::endl;
  allpass = test_cartesian_geometry<3>(uref_cartesian_3d) && allpass;

  std::cout << "Tesing GeneralGeometry<1>" << std::endl;
  allpass = test_general_geometry<1>(uref_general_1d) && allpass;
  std::cout << "Tesing GeneralGeometry<2>" << std::endl;
  allpass = test_general_geometry<2>(uref_general_2d) && allpass;
  std::cout << "Tesing GeneralGeometry<3>" << std::endl;
  allpass = test_general_geometry<3>(uref_general_3d) && allpass;

  return allpass ? 0 : 1;
}
