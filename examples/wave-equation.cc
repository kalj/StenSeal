
    #include <deal.II/base/point.h>
    #include <deal.II/lac/vector.h>
    #include <deal.II/base/function.h>

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

    template <int dim>
    class InitialValues : public dealii::Function<dim>
    {

    public:
      InitialValues () : dealii::Function<dim>() {}
      virtual ~InitialValues(){};

      virtual double value (const dealii::Point<dim>   &p,
                            const unsigned int  component = 0) const;
    };

    template <int dim>
    double InitialValues<dim>::value (const dealii::Point<dim> &p,
                                       const unsigned int /*component*/) const
    {

      double radius = 0.3;
      dealii::Point<dim> center;
      center[0] = 0.2;

      double res = 0;
      double dist = (p-center).norm();

      if( dist < radius)
        res = (1+cos(dealii::numbers::PI*dist/radius))/2;
      return res;
    }

    template <int dim, typename Geometry>
    void initialize( dealii::Vector<double> &u, dealii::Function<dim> &f, const Geometry &g )
    {
    if(dim==1) {
      for(int i = 0; i < g.n_nodes[0]; ++i) {
        double x;
        x = i*g.h[0] + g.lower_left(0);
        dealii::Point<dim> p(x);
        u[i] = f.value(p);
      }
     } else if(dim==2){
       for(int i = 0; i < g.n_nodes[0]; ++i) {
         for(int j = 0; j < g.n_nodes[1]; ++j) {
            double x;
            double y;
            x = i*g.h[0] + g.lower_left(0);
            y = j*g.h[1] + g.lower_left(1);
            dealii::Point<dim> p(x,y);
            u[g.n_nodes[1]*i+j] = f.value(p);
          }
        }
     }
     else  {
          AssertThrow(false,dealii::ExcNotImplemented());
      }
    }

    int main(int argc, char *argv[])
    {
      const int dim = 2;

      double xmin = 0;
      double xmax = 1.0;
      double ymin = 0;
      double ymax = 1.0;

      unsigned int n_nodes[dim] = { 10, 10 };
      dealii::Point<dim> lower_left(xmin,ymin);
      dealii::Point<dim> upper_right(xmax,ymax);

      typedef stenseal::CartesianGeometry<dim> Geometry;
      Geometry geometry(n_nodes,lower_left,upper_right);

      dealii::Vector<double> u(geometry.n_nodes_total);  
      InitialValues<dim> f; 
      initialize<dim,Geometry>(u,f,geometry);

      stenseal::UpwindBlockOperator<dim,StencilDm,StencilDp,Geometry> op(geometry);
      u.print(std::cout);


      /*
      const int dim = 1;

      double xmin = 0;
      double xmax = 1.0;

      unsigned int n_nodes[dim] = {200};
      dealii::Point<dim> lower_left(xmin);
      dealii::Point<dim> upper_right(xmax);

      typedef stenseal::CartesianGeometry<dim> Geometry;
      Geometry geometry(n_nodes,lower_left,upper_right);

      dealii::Vector<double> u(geometry.n_nodes_total);  
      InitialValues<dim> f; 
      initialize<dim,Geometry>(u,f,geometry);

      stenseal::UpwindBlockOperator<dim,StencilDm,StencilDp,Geometry> op(geometry);
      u.print(std::cout);
      */
      
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
