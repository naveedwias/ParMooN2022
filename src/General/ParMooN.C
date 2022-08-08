#include "ParMooN.h"
#include "Database.h"
#include "MooNMD_Io.h"
#include "FEDatabase.h"
#include "Chrono.h"
#ifdef _MPI
#  include "mpi.h"
#endif

// to print the overall computation time at the end
static Chrono duration_all;

ParameterDatabase parmoon::parmoon_initialize(int argc, char* argv[])
{
  int my_rank = 0;
#ifdef _MPI
  // Construct and initialize the default MPI communicator.
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif
  TDatabase::create(argc > 1 ? argv[1] : nullptr);
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  if(argc > 1)
  {
    parmoon_db.read(argv[1]);
  }
  
  // open outfile, this is where all output is written to (additionally to
  // console)
  if(my_rank == 0)
  {
    Output::set_outfile(parmoon_db["outfile"], parmoon_db["script_mode"]);
  }
  Output::setVerbosity(parmoon_db["verbosity"]);
  
  // write all Parameters to the OUTFILE (not to console) for later reference
  if(argc > 0 && my_rank == 0)
  {
    parmoon_db.write(Output::get_outfile());
    TDatabase::WriteParamDB(argv[0]);
    TDatabase::WriteTimeDB();
  }
  
  duration_all.reset();
  duration_all.start();
  
  return parmoon_db;
}


void parmoon::parmoon_finalize()
{
  FEDatabase::destroy();
  TDatabase::destroy();
  duration_all.stop_and_print("entire ParMooN program");
  Output::print_warnings();
  Output::close_file();
#ifdef _MPI
  MPI_Finalize();
#endif
}
