#include "ConvectionDiffusion.h"
#include "NavierStokes.h"
#include "Darcy.h"
#include "TimeConvectionDiffusion.h"
#include "TimeNavierStokes.h"
#ifdef _MPI
#include "mpi.h"
#endif

#ifdef __2D__
constexpr int d = 2;
#else
constexpr int d = 3;
#endif

int main(int argc, char* argv[])
{
  bool i_am_root = true;
#ifdef _MPI
  MPI_Init(&argc, &argv);
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  i_am_root = (my_rank == 0);
#else
  (void) argc;
  (void) argv ;
#endif
  if(i_am_root)
  {
    auto cd_db = ConvectionDiffusion<d>::default_cd_database(true);
    cd_db.add("ANSATZ_ORDER", 1, "Polynomial degree of ansatz functions",
                   0, 1000);
    cd_db.write(d == 2 ? "cd2d.dat" : "cd3d.dat", true, false);
  }
  
  if(i_am_root)
  {
    auto nse_db = NavierStokes<d>::default_nse_database(true);
    nse_db.add("VELOCITY_SPACE", 2, "Polynomial order for the velocity space",
               0, 1000);
    nse_db.add("PRESSURE_SPACE", 1, "Polynomial order for the pressure space",
               0, 1000);
    nse_db.add("NSTYPE", 14, "indicate which kind of storage layout is used",
               {1, 2, 3, 4, 14});
    nse_db.write(d == 2 ? "nse2d.dat" : "nse3d.dat", true, false);
  }
  
  if(i_am_root)
  {
    auto darcy_db = Darcy<d>::default_darcy_database(true);
    darcy_db.add("VELOCITY_SPACE", 1000,
                 "Polynomial order for the velocity space", 0, 10000);
    darcy_db.add("PRESSURE_SPACE", -4711,
                 "Polynomial order for the pressure space", -4711, 10000);
    darcy_db.write(d == 2 ? "darcy2d.dat" : "darcy3d.dat", true, false);
  }
  
  if(i_am_root)
  {
    auto tcd_db = TimeConvectionDiffusion<d>::default_tcd_database(true);
    tcd_db.add("ANSATZ_ORDER", 1, "Polynomial degree of ansatz functions",
                   0, 1000);
    tcd_db.write(d == 2 ? "tcd2d.dat" : "tcd3d.dat", true, false);
  }
  
  if(i_am_root)
  {
    auto tnse_db = TimeNavierStokes<d>::default_tnse_database(true);
    tnse_db.add("VELOCITY_SPACE", 2, "Polynomial order for the velocity space",
               0, 1000);
    tnse_db.add("PRESSURE_SPACE", 1, "Polynomial order for the pressure space",
               0, 1000);
    tnse_db.add("NSTYPE", 14, "indicate which kind of storage layout is used",
               {1, 2, 3, 4, 14});
    tnse_db.write(d == 2 ? "nse2d.dat" : "nse3d.dat", true, false);
  }
#ifdef _MPI
  MPI_Finalize();
#endif
  return 0;
}
