// =======================================================================
// @(#)Database.C        1.37 06/27/00
// 
// Class:       TDatabase
// Purpose:     database of needed refinement, mapping and
//              shape descriptors
//              parameter database
//
// Author:      Volker Behns  29.07.97
//
// =======================================================================
#ifdef _MPI
#  include "mpi.h"
#endif

#include <Database.h>
#include <Mapper.h>
#include <Line.h>
#include <ParameterDatabase.h>
#include <MainUtilities.h>
#include "QuadratureFormulaDatabase.h"
// more explicit path due to name conflict with the 'triangle' library on mac
#include <../Geometry/Triangle.h> 
#include <Quadrangle.h>
#include <Parallelogram.h>
#include <Rectangle.h>
#include <RefNoRef.h>
#include <RefLineDesc.h>
#include <RefQuadRegDesc.h>
#include <RefQuadBis0Desc.h>
#include <RefQuadBis1Desc.h>
#include <RefQuad1Conf0Desc.h>
#include <RefQuad1Conf1Desc.h>
#include <RefQuad1Conf2Desc.h>
#include <RefQuad1Conf3Desc.h>
#include <RefQuad2Conf0Desc.h>
#include <RefQuad2Conf1Desc.h>
#include <RefQuad2Conf2Desc.h>
#include <RefQuad2Conf3Desc.h>
#include <RefQuadToTri0Desc.h>
#include <RefQuadToTri1Desc.h>
#include <RefTriRegDesc.h>
#include <RefTriBaryDesc.h>
#include <RefTriBis0Desc.h>
#include <RefTriBis1Desc.h>
#include <RefTriBis2Desc.h>
#include <RefTriBis01Desc.h>
#include <RefTriBis02Desc.h>
#include <RefTriBis10Desc.h>
#include <RefTriBis12Desc.h>
#include <RefTriBis20Desc.h>
#include <RefTriBis21Desc.h>
#include <It_Between.h>
#include <It_EQ.h>
#include <It_Finest.h>
#include <It_EQLevel.h>
#include <It_LELevel.h>
#include <It_OCAF.h>
#include <Utilities.h>
#include <MooNMD_Io.h>
#include <string.h>

#ifdef __3D__
  #include <Tetrahedron.h>
  #include <Hexahedron.h>
  #include <Brick.h>
  #include <RefTetraRegDesc.h>
  #include <RefTetraBaryDesc.h>
  #include <RefTetraReg0Desc.h>
  #include <RefTetraReg1Desc.h>
  #include <RefTetraReg2Desc.h>
  #include <RefTetraBis0Desc.h>
  #include <RefTetraBis1Desc.h>
  #include <RefTetraBis2Desc.h>
  #include <RefTetraBis3Desc.h>
  #include <RefTetraBis4Desc.h>
  #include <RefTetraBis5Desc.h>
  #include <RefTetraBis01Desc.h>
  #include <RefTetraBis02Desc.h>
  #include <RefTetraBis03Desc.h>
  #include <RefTetraBis04Desc.h>
  #include <RefTetraBis05Desc.h>
  #include <RefTetraBis10Desc.h>
  #include <RefTetraBis12Desc.h>
  #include <RefTetraBis13Desc.h>
  #include <RefTetraBis14Desc.h>
  #include <RefTetraBis15Desc.h>
  #include <RefTetraBis20Desc.h>
  #include <RefTetraBis21Desc.h>
  #include <RefTetraBis23Desc.h>
  #include <RefTetraBis24Desc.h>
  #include <RefTetraBis25Desc.h>
  #include <RefTetraBis30Desc.h>
  #include <RefTetraBis32Desc.h>
  #include <RefTetraBis34Desc.h>
  #include <RefTetraBis35Desc.h>
  #include <RefTetraBis40Desc.h>
  #include <RefTetraBis41Desc.h>
  #include <RefTetraBis43Desc.h>
  #include <RefTetraBis45Desc.h>
  #include <RefTetraBis51Desc.h>
  #include <RefTetraBis52Desc.h>
  #include <RefTetraBis53Desc.h>
  #include <RefTetraBis54Desc.h>
  #include <RefTetraQuad0Desc.h>
  #include <RefTetraQuad1Desc.h>
  #include <RefTetraQuad2Desc.h>
  #include <RefTetraQuad3Desc.h>
  #include <RefHexaRegDesc.h>
#endif


bool SearchParamInLine(std::string  line,
                       std::string  param_name,
                       std::string& param_value)
{
  bool is_found = false;
  // check if parameter is commented
  std::size_t param_pos = line.find(param_name);
  if((param_pos!=std::string::npos) && (param_pos<line.find("#")))
  {
    //check if value is given
    std::size_t end_param_pos = param_pos+param_name.length();
    std::size_t value_pos = line.find_first_not_of(" \n\t", end_param_pos);
    if(value_pos==std::string::npos)
    {
      int rank = 0;
#ifdef _MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
      if(rank==0)
      {
        ErrThrow("The parameter: ", param_name,
                  " is given but its value is missing.");
      }
    }
    else
    {
      param_value = line.substr(value_pos);
      is_found = true;
    }
  }
  return is_found;
}


template <typename T>
bool SetParameters(std::string line,
                   std::string param_name,
                   T&          param)
{
  std::string param_value;
  bool is_found = SearchParamInLine(line, param_name, param_value);
  if(is_found)
  {
    std::istringstream value(param_value);
    value>>param;
  }
  return is_found;
}


template<typename T>
bool SetParameters(std::string     line,
                   std::string     param_name,
                   std::vector<T>& param)
{
  std::string param_values;
  bool is_found = SearchParamInLine(line, param_name, param_values);
  if(is_found)
  {
    std::size_t value_pos;
    std::size_t end_value_pos = 0;;
    size_t nb_value = param.size();

    for(size_t i_val=0; i_val<nb_value; i_val++)
    {
      value_pos = param_values.find_first_not_of(" (,)\n\t", end_value_pos);
      if(value_pos==std::string::npos)
      {
        int rank = 0;
#ifdef _MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
        if(rank==0)
        {
          ErrThrow("The parameter: ", param_name, " seems to only have ", i_val,
                   " values, but ", nb_value, " are expected.");
        }
      }
      end_value_pos = param_values.find_first_of(" (,)\n\t", value_pos);
      std::istringstream value(param_values.substr(value_pos,
                                                   end_value_pos-value_pos));
      value>>param[i_val];
    }
  }
  return is_found;
}


void TParamDB::read_parameters(const char* ParamFile)
{
  int rank = 0;
#ifdef _MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  std::string line;
  int N_Param = 0;
  std::ifstream dat(ParamFile);

  if (!dat)
  {
    ErrThrow("cannot open '", ParamFile, "' for input");
  }

  while(std::getline(dat, line))
  {
    if(SetParameters(line, "VERSION:", VERSION)){N_Param++; continue;}
    if(SetParameters(line, "ANSATZ_ORDER:", ANSATZ_ORDER)){N_Param++; continue;}
    if(SetParameters(line, "VELOCITY_SPACE:", VELOCITY_SPACE)){N_Param++; continue;}
    if(SetParameters(line, "PRESSURE_SPACE:", PRESSURE_SPACE)){N_Param++; continue;}
    if(SetParameters(line, "NSTYPE:", NSTYPE)){N_Param++; continue;}
    if(SetParameters(line, "LAPLACETYPE:", LAPLACETYPE)){N_Param++; continue;}
    if(SetParameters(line, "SIGMA_PERM:", SIGMA_PERM)){N_Param++; continue;}
    if(SetParameters(line, "USE_ISOPARAMETRIC:", USE_ISOPARAMETRIC)){N_Param++; continue;}
    if(SetParameters(line, "UPWIND_ORDER:", UPWIND_ORDER)){N_Param++; continue;}
    if(SetParameters(line, "CELL_MEASURE:", CELL_MEASURE)){N_Param++; continue;}
    if(SetParameters(line, "SAMARSKI_DAMP:", UPWIND_FLUX_DAMP)){N_Param++; continue;}
    if(SetParameters(line, "UPWIND_APPLICATION:", UPWIND_APPLICATION)){N_Param++; continue;}
    if(SetParameters(line, "RE_NR:", RE_NR)){N_Param++; continue;}
    if(SetParameters(line, "FLOW_PROBLEM_TYPE:", FLOW_PROBLEM_TYPE)){N_Param++; continue;}
    if(SetParameters(line, "FACE_SIGMA:", FACE_SIGMA)){N_Param++; continue;}
    if(SetParameters(line, "WEAK_BC_SIGMA:", WEAK_BC_SIGMA)){N_Param++; continue;}
    if(SetParameters(line, "WEAK_BC:", WEAK_BC)){N_Param++; continue;}
    if(SetParameters(line, "PE_NR:", PE_NR)){N_Param++; continue;}
    if(SetParameters(line, "DELTA1:", DELTA1)){N_Param++; continue;}
    if(SetParameters(line, "DELTA0:", DELTA0)){N_Param++; continue;}
    if(SetParameters(line, "SDFEM_TYPE:", SDFEM_TYPE)){N_Param++; continue;}
    if(SetParameters(line, "FILTER_WIDTH_CONSTANT:", FILTER_WIDTH_CONSTANT)){N_Param++; continue;}
    if(SetParameters(line, "FILTER_WIDTH_POWER:", FILTER_WIDTH_POWER)){N_Param++; continue;}
    if(SetParameters(line, "GAUSSIAN_GAMMA:", GAUSSIAN_GAMMA)){N_Param++; continue;}
    if(SetParameters(line, "TURBULENT_VISCOSITY_TYPE:", TURBULENT_VISCOSITY_TYPE)){N_Param++; continue;}
    if(SetParameters(line, "TURBULENT_VISCOSITY_TENSOR:", TURBULENT_VISCOSITY_TENSOR)){N_Param++; continue;}
    if(SetParameters(line, "TURBULENT_VISCOSITY_CONSTANT:", TURBULENT_VISCOSITY_CONSTANT)){N_Param++; continue;}
    if(SetParameters(line, "TURBULENT_VISCOSITY_POWER:", TURBULENT_VISCOSITY_POWER)){N_Param++; continue;}
    if(SetParameters(line, "TURBULENT_VISCOSITY_SIGMA:", TURBULENT_VISCOSITY_SIGMA)){N_Param++; continue;}
    if(SetParameters(line, "ARTIFICIAL_VISCOSITY_CONSTANT:", ARTIFICIAL_VISCOSITY_CONSTANT)){N_Param++; continue;}
    if(SetParameters(line, "ARTIFICIAL_VISCOSITY_POWER:", ARTIFICIAL_VISCOSITY_POWER)){N_Param++; continue;}
    if(SetParameters(line, "FRICTION_CONSTANT:", FRICTION_CONSTANT)){N_Param++; continue;}
    if(SetParameters(line, "FRICTION_POWER:", FRICTION_POWER)){N_Param++; continue;}
    if(SetParameters(line, "FRICTION_U0:", FRICTION_U0)){N_Param++; continue;}
    if(SetParameters(line, "FRICTION_TYPE:", FRICTION_TYPE)){N_Param++; continue;}
    if(SetParameters(line, "PENETRATION_CONSTANT:", PENETRATION_CONSTANT)){N_Param++; continue;}
    if(SetParameters(line, "PENETRATION_POWER:", PENETRATION_POWER)){N_Param++; continue;}
    if(SetParameters(line, "OSEEN_ZERO_ORDER_COEFF:", OSEEN_ZERO_ORDER_COEFF)){N_Param++; continue;}
    if(SetParameters(line, "DIV_DIV_STAB_TYPE:", DIV_DIV_STAB_TYPE)){N_Param++; continue;}
    if(SetParameters(line, "DIV_DIV_STAB_C1:", DIV_DIV_STAB_C1)){N_Param++; continue;}
    if(SetParameters(line, "DIV_DIV_STAB_C2:", DIV_DIV_STAB_C2)){N_Param++; continue;}
    if(SetParameters(line, "LP_FULL_GRADIENT:", LP_FULL_GRADIENT)){N_Param++; continue;}
    if(SetParameters(line, "LP_STREAMLINE:", LP_STREAMLINE)){N_Param++; continue;}
    if(SetParameters(line, "LP_DIVERGENCE:", LP_DIVERGENCE)){N_Param++; continue;}
    if(SetParameters(line, "LP_PRESSURE:", LP_PRESSURE)){N_Param++; continue;}
    if(SetParameters(line, "LP_COEFF_TYPE:", LP_COEFF_TYPE)){N_Param++; continue;}
    if(SetParameters(line, "LP_FULL_GRADIENT_COEFF:", LP_FULL_GRADIENT_COEFF)){N_Param++; continue;}
    if(SetParameters(line, "LP_STREAMLINE_COEFF:", LP_STREAMLINE_COEFF)){N_Param++; continue;}
    if(SetParameters(line, "LP_PRESSURE_COEFF:", LP_PRESSURE_COEFF)){N_Param++; continue;}
    if(SetParameters(line, "LP_FULL_GRADIENT_EXPONENT:", LP_FULL_GRADIENT_EXPONENT)){N_Param++; continue;}
    if(SetParameters(line, "LP_STREAMLINE_EXPONENT:", LP_STREAMLINE_EXPONENT)){N_Param++; continue;}
    if(SetParameters(line, "LP_PRESSURE_EXPONENT:", LP_PRESSURE_EXPONENT)){N_Param++; continue;}
    if(SetParameters(line, "LP_ORDER_DIFFERENCE:", LP_ORDER_DIFFERENCE)){N_Param++; continue;}
    if(SetParameters(line, "LP_FULL_GRADIENT_ORDER_DIFFERENCE:", LP_FULL_GRADIENT_ORDER_DIFFERENCE)){N_Param++; continue;}
    if(SetParameters(line, "LP_STREAMLINE_ORDER_DIFFERENCE:", LP_STREAMLINE_ORDER_DIFFERENCE)){N_Param++; continue;}
    if(SetParameters(line, "LP_PRESSURE_ORDER_DIFFERENCE:", LP_PRESSURE_ORDER_DIFFERENCE)){N_Param++; continue;}
    if(SetParameters(line, "LP_CROSSWIND:", LP_CROSSWIND)){N_Param++; continue;}
    if(SetParameters(line, "LP_CROSSWIND_COEFF_TYPE:", LP_CROSSWIND_COEFF_TYPE)){N_Param++; continue;}
    if(SetParameters(line, "LP_CROSSWIND_COEFF:", LP_CROSSWIND_COEFF)){N_Param++; continue;}
    if(SetParameters(line, "LP_CROSSWIND_EXPONENT:", LP_CROSSWIND_EXPONENT)){N_Param++; continue;}
    if(SetParameters(line, "SOLD_TYPE:", SOLD_TYPE)){N_Param++; continue;}
    if(SetParameters(line, "SOLD_PARAMETER_TYPE:", SOLD_PARAMETER_TYPE)){N_Param++; continue;}
    if(SetParameters(line, "SOLD_CONST:", SOLD_CONST)){N_Param++; continue;}
    if(SetParameters(line, "SOLD_POWER:", SOLD_POWER)){N_Param++; continue;}
    if(SetParameters(line, "SOLD_S:", SOLD_S)){N_Param++; continue;}
    if(SetParameters(line, "SOLD_U0:", SOLD_U0)){N_Param++; continue;}
    if(SetParameters(line, "SOLD_PARAMETER_SCALING:", SOLD_PARAMETER_SCALING)){N_Param++; continue;}
    if(SetParameters(line, "SOLD_PARAMETER_SCALING_FACTOR:", SOLD_PARAMETER_SCALING_FACTOR)){N_Param++; continue;}
    if(SetParameters(line, "P0:", P0)){N_Param++; continue;}
    if(SetParameters(line, "P1:", P1)){N_Param++; continue;}
    if(SetParameters(line, "P2:", P2)){N_Param++; continue;}
    if(SetParameters(line, "P3:", P3)){N_Param++; continue;}
    if(SetParameters(line, "P4:", P4)){N_Param++; continue;}
    if(SetParameters(line, "P5:", P5)){N_Param++; continue;}
    if(SetParameters(line, "P6:", P6)){N_Param++; continue;}
    if(SetParameters(line, "P7:", P7)){N_Param++; continue;}
    if(SetParameters(line, "P8:", P8)){N_Param++; continue;}
    if(SetParameters(line, "P9:", P9)){N_Param++; continue;}
    if(SetParameters(line, "P10:", P10)){N_Param++; continue;}
    if(SetParameters(line, "P11:", P11)){N_Param++; continue;}
    if(SetParameters(line, "P12:", P12)){N_Param++; continue;}
    if(SetParameters(line, "P13:", P13)){N_Param++; continue;}
    if(SetParameters(line, "P14:", P14)){N_Param++; continue;}
    if(SetParameters(line, "P15:", P15)){N_Param++; continue;}
    if(SetParameters(line, "CIP_TYPE:", CIP_TYPE)){N_Param++; continue;}

    if(SetParameters(line, "VMS_ADAPT_LOWER:", VMS_ADAPT_LOWER)){N_Param++; continue;}
    if(SetParameters(line, "VMS_ADAPT_MIDDLE:", VMS_ADAPT_MIDDLE)){N_Param++; continue;}
    if(SetParameters(line, "VMS_ADAPT_UPPER:", VMS_ADAPT_UPPER)){N_Param++; continue;}
    if(SetParameters(line, "VMS_ADAPT_COMP:", VMS_ADAPT_COMP)){N_Param++; continue;}
    if(SetParameters(line, "SUPERCONVERGENCE_ORDER:", SUPERCONVERGENCE_ORDER)){N_Param++; continue;}
    if(SetParameters(line, "CYLINDER_22000_YPLUS_SIDES:", CYLINDER_22000_YPLUS_SIDES)){N_Param++; continue;}
    if(SetParameters(line, "CYLINDER_22000_YPLUS_FRONT:", CYLINDER_22000_YPLUS_FRONT)){N_Param++; continue;}
    if(SetParameters(line, "CYLINDER_22000_YPLUS_BACK:", CYLINDER_22000_YPLUS_BACK)){N_Param++; continue;}
    if(SetParameters(line, "INPUT_QUAD_RULE:", INPUT_QUAD_RULE)){N_Param++; continue;}
    if(SetParameters(line, "Par_P0:", Par_P0)){N_Param++; continue;}
    if(SetParameters(line, "Par_P1:", Par_P1)){N_Param++; continue;}
    if(SetParameters(line, "Par_P2:", Par_P2)){N_Param++; continue;}
    if(SetParameters(line, "Par_P3:", Par_P3)){N_Param++; continue;}
    if(SetParameters(line, "Par_P4:", Par_P4)){N_Param++; continue;}
    if(SetParameters(line, "Par_P5:", Par_P5)){N_Param++; continue;}
    if(SetParameters(line, "Par_P6:", Par_P6)){N_Param++; continue;}
    if(SetParameters(line, "Par_P7:", Par_P7)){N_Param++; continue;}
    if(SetParameters(line, "Par_P8:", Par_P8)){N_Param++; continue;}
    if(SetParameters(line, "Par_P9:", Par_P9)){N_Param++; continue;}

    // parameters for weakly imposing boundary/interface conditions
    if(SetParameters(line, "n_neumann_boundary:", n_neumann_boundary))
    {
      neumann_boundary_id.resize(n_neumann_boundary);
      neumann_boundary_value.resize(n_neumann_boundary);
      N_Param++;
      continue;
    }
    if(!neumann_boundary_id.empty())
      if(SetParameters(line, "neumann_boundary_id:", neumann_boundary_id)){N_Param++; continue;}
    if(!neumann_boundary_value.empty())
      if(SetParameters(line, "neumann_boundary_value:", neumann_boundary_value)){N_Param++; continue;}
    if(SetParameters(line, "n_unvn_boundary:", n_unvn_boundary))
    {
      unvn_boundary_id.resize(n_unvn_boundary);
      unvn_boundary_value.resize(n_unvn_boundary);
      N_Param++;
      continue;
    }
    if(!unvn_boundary_id.empty())
      if(SetParameters(line, "unvn_boundary_id:", unvn_boundary_id)){N_Param++; continue;}
    if(!unvn_boundary_value.empty())
      if(SetParameters(line, "unvn_boundary_value:", unvn_boundary_value)){N_Param++; continue;}
    if(SetParameters(line, "n_gradunv_boundary:", n_gradunv_boundary))
    {
      gradunv_boundary_id.resize(n_gradunv_boundary);
      gradunv_boundary_value.resize(n_gradunv_boundary);
      N_Param++;
      continue;
    }
    if(!gradunv_boundary_id.empty())
      if(SetParameters(line, "gradunv_boundary_id:", gradunv_boundary_id)){N_Param++; continue;}
    if(!gradunv_boundary_value.empty())
      if(SetParameters(line, "gradunv_boundary_value:", gradunv_boundary_value)){N_Param++; continue;}
    if(SetParameters(line, "n_u_v_boundary:", n_u_v_boundary))
    {
      u_v_boundary_id.resize(n_u_v_boundary);
      u_v_boundary_value.resize(n_u_v_boundary);
      N_Param++;
      continue;
    }
    if(!u_v_boundary_id.empty())
      if(SetParameters(line, "u_v_boundary_id:", u_v_boundary_id)){N_Param++; continue;}
    if(!u_v_boundary_value.empty())
      if(SetParameters(line, "u_v_boundary_value:", u_v_boundary_value)){N_Param++; continue;}
    if(SetParameters(line, "n_g_v_boundary:", n_g_v_boundary))
    {
      g_v_boundary_id.resize(n_g_v_boundary);
      g_v_boundary_value.resize(n_g_v_boundary);
      N_Param++;
      continue;
    }
    if(!g_v_boundary_id.empty())
      if(SetParameters(line, "g_v_boundary_id:", g_v_boundary_id)){N_Param++; continue;}
    if(!g_v_boundary_value.empty())
      if(SetParameters(line, "g_v_boundary_value:", g_v_boundary_value)){N_Param++; continue;}
    if(SetParameters(line, "n_p_v_n_boundary:", n_p_v_n_boundary))
    {
      p_v_n_boundary_id.resize(n_p_v_n_boundary);
      p_v_n_boundary_value.resize(n_p_v_n_boundary);
      N_Param++;
      continue;
    }
    if(!p_v_n_boundary_id.empty())
      if(SetParameters(line, "p_v_n_boundary_id:", p_v_n_boundary_id)){N_Param++; continue;}
    if(!p_v_n_boundary_value.empty())
      if(SetParameters(line, "p_v_n_boundary_value:", p_v_n_boundary_value)){N_Param++; continue;}

    // Nitsche Combi- weak Dirichlet
    if(SetParameters(line, "n_nitsche_boundary:", n_nitsche_boundary))
    {
      nitsche_boundary_id.resize(n_nitsche_boundary);
      nitsche_penalty.resize(n_nitsche_boundary);
      N_Param++;
      continue;
    }
    if(!nitsche_boundary_id.empty())
      if(SetParameters(line, "nitsche_boundary_id:", nitsche_boundary_id)){N_Param++; continue;}
    if(!nitsche_penalty.empty())
      if(SetParameters(line, "nitsche_penalty:", nitsche_penalty)){N_Param++; continue;}


    // ----------------------------------------------------------------


    if(SetParameters(line, "MapperType:", MapperType)){N_Param++; continue;}
  }


  if(rank==0)
    Output::info("read_parameters", "Parameter database (old) read with ",
                 N_Param, " parameters. Parameter file version ",
                 VERSION);
}


// Constructors
//TDatabase::TDatabase(const char *ParamFile)
void TDatabase::create(const char *ParamFile)
{
  // allocate databases
  ShapeDB = new const TShapeDesc*[N_SHAPES];
  RefDescDB = new const TRefDesc*[N_SHAPES + N_REFDESC];
  IteratorDB = new TIterator*[N_ITERATORS];
  ParamDB = new TParamDB;
  TimeDB = new TTimeDB;

  // initialize shape descriptors
  ShapeDB[S_Line] = new TLine();
  RefDescDB[S_Line] = new TRefNoRef(ShapeDB[S_Line]);

  ShapeDB[Triangle] = new TTriangle();
  RefDescDB[Triangle] = new TRefNoRef(ShapeDB[Triangle]);

  ShapeDB[Quadrangle] = new TQuadrangle();
  RefDescDB[Quadrangle] = new TRefNoRef(ShapeDB[Quadrangle]);

  ShapeDB[Parallelogram] = new TParallelogram();
  RefDescDB[Parallelogram] = new TRefNoRef(ShapeDB[Parallelogram]);

  ShapeDB[Rectangle] = new TRectangle();
  RefDescDB[Rectangle] = new TRefNoRef(ShapeDB[Rectangle]);

  #ifdef __3D__
    ShapeDB[Tetrahedron] = new TTetrahedron();
    RefDescDB[Tetrahedron] = new TRefNoRef(ShapeDB[Tetrahedron]);

    ShapeDB[Hexahedron] = new THexahedron();
    RefDescDB[Hexahedron] = new TRefNoRef(ShapeDB[Hexahedron]);

    ShapeDB[Brick] = new TBrick();
    RefDescDB[Brick] = new TRefNoRef(ShapeDB[Brick]);
  #endif

  // initialize refinement descriptors
  RefDescDB[N_SHAPES + LineReg] = new TRefLineDesc(ShapeDB[S_Line]);
  RefDescDB[N_SHAPES + TriReg]  = new TRefTriRegDesc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBary]  = new TRefTriBaryDesc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBis0] = new TRefTriBis0Desc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBis1] = new TRefTriBis1Desc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBis2] = new TRefTriBis2Desc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBis01]= new TRefTriBis01Desc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBis02]= new TRefTriBis02Desc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBis10]= new TRefTriBis10Desc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBis12]= new TRefTriBis12Desc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBis20]= new TRefTriBis20Desc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBis21]= new TRefTriBis21Desc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + QuadReg] = new TRefQuadRegDesc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES + ParallReg] = new TRefQuadRegDesc(ShapeDB[Parallelogram]);
  RefDescDB[N_SHAPES + RectReg] = new TRefQuadRegDesc(ShapeDB[Rectangle]);
  RefDescDB[N_SHAPES + QuadBis0] = new TRefQuadBis0Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES + QuadBis1] = new TRefQuadBis1Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES+Quad1Conf0] = new TRefQuad1Conf0Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES+Quad1Conf1] = new TRefQuad1Conf1Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES+Quad1Conf2] = new TRefQuad1Conf2Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES+Quad1Conf3] = new TRefQuad1Conf3Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES+Quad2Conf0] = new TRefQuad2Conf0Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES+Quad2Conf1] = new TRefQuad2Conf1Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES+Quad2Conf2] = new TRefQuad2Conf2Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES+Quad2Conf3] = new TRefQuad2Conf3Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES + QuadToTri0] = new
      TRefQuadToTri0Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES + QuadToTri1] = new
      TRefQuadToTri1Desc(ShapeDB[Quadrangle]);

  #ifdef __3D__
    RefDescDB[N_SHAPES + TetraReg] =
         new TRefTetraRegDesc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBary] =
         new TRefTetraBaryDesc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraReg0] =
         new TRefTetraReg0Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraReg1] =
         new TRefTetraReg1Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraReg2] =
         new TRefTetraReg2Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis0] = new TRefTetraBis0Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis1] = new TRefTetraBis1Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis2] = new TRefTetraBis2Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis3] = new TRefTetraBis3Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis4] = new TRefTetraBis4Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis5] = new TRefTetraBis5Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis01] = new TRefTetraBis01Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis02] = new TRefTetraBis02Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis03] = new TRefTetraBis03Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis04] = new TRefTetraBis04Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis05] = new TRefTetraBis05Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis10] = new TRefTetraBis10Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis12] = new TRefTetraBis12Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis13] = new TRefTetraBis13Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis14] = new TRefTetraBis14Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis15] = new TRefTetraBis15Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis20] = new TRefTetraBis20Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis21] = new TRefTetraBis21Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis23] = new TRefTetraBis23Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis24] = new TRefTetraBis24Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis25] = new TRefTetraBis25Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis30] = new TRefTetraBis30Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis32] = new TRefTetraBis32Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis34] = new TRefTetraBis34Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis35] = new TRefTetraBis35Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis40] = new TRefTetraBis40Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis41] = new TRefTetraBis41Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis43] = new TRefTetraBis43Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis45] = new TRefTetraBis45Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis51] = new TRefTetraBis51Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis52] = new TRefTetraBis52Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis53] = new TRefTetraBis53Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis54] = new TRefTetraBis54Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraQuad0] = new TRefTetraQuad0Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraQuad1] = new TRefTetraQuad1Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraQuad2] = new TRefTetraQuad2Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraQuad3] = new TRefTetraQuad3Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + HexaReg]  =
         new TRefHexaRegDesc(ShapeDB[Hexahedron]);
    RefDescDB[N_SHAPES + BrickReg]  =
         new TRefHexaRegDesc(ShapeDB[Brick]);
  #endif

  #ifdef __3D__
    //initialize mapper
    MapperDB = new const TMapper*[N_MAPPER];
    MapperDB[MapTriReg0] = new TMapper(MapTriReg0);
    MapperDB[MapTriReg1] = new TMapper(MapTriReg1);
    MapperDB[MapTriReg2] = new TMapper(MapTriReg2);

    MapperDB[MapTriBis00] = new TMapper(MapTriBis00);
    MapperDB[MapTriBis01] = new TMapper(MapTriBis01);
    MapperDB[MapTriBis02] = new TMapper(MapTriBis02);
    MapperDB[MapTriBis10] = new TMapper(MapTriBis10);
    MapperDB[MapTriBis11] = new TMapper(MapTriBis11);
    MapperDB[MapTriBis12] = new TMapper(MapTriBis12);
    MapperDB[MapTriBis20] = new TMapper(MapTriBis20);
    MapperDB[MapTriBis21] = new TMapper(MapTriBis21);
    MapperDB[MapTriBis22] = new TMapper(MapTriBis22);
    MapperDB[MapTriBis010] = new TMapper(MapTriBis010);
    MapperDB[MapTriBis011] = new TMapper(MapTriBis011);
    MapperDB[MapTriBis012] = new TMapper(MapTriBis012);
    MapperDB[MapTriBis020] = new TMapper(MapTriBis020);
    MapperDB[MapTriBis021] = new TMapper(MapTriBis021);
    MapperDB[MapTriBis022] = new TMapper(MapTriBis022);
    MapperDB[MapTriBis100] = new TMapper(MapTriBis100);
    MapperDB[MapTriBis101] = new TMapper(MapTriBis101);
    MapperDB[MapTriBis102] = new TMapper(MapTriBis102);
    MapperDB[MapTriBis120] = new TMapper(MapTriBis120);
    MapperDB[MapTriBis121] = new TMapper(MapTriBis121);
    MapperDB[MapTriBis122] = new TMapper(MapTriBis122);
    MapperDB[MapTriBis200] = new TMapper(MapTriBis200);
    MapperDB[MapTriBis201] = new TMapper(MapTriBis201);
    MapperDB[MapTriBis202] = new TMapper(MapTriBis202);
    MapperDB[MapTriBis210] = new TMapper(MapTriBis210);
    MapperDB[MapTriBis211] = new TMapper(MapTriBis211);
    MapperDB[MapTriBis212] = new TMapper(MapTriBis212);

    MapperDB[MapQuadReg0] = new TMapper(MapQuadReg0);
    MapperDB[MapQuadReg1] = new TMapper(MapQuadReg1);
    MapperDB[MapQuadReg2] = new TMapper(MapQuadReg2);
    MapperDB[MapQuadReg3] = new TMapper(MapQuadReg3);
  #endif

  // initialize iterators
  IteratorDB[It_EQ] = new TIt_EQ();
  IteratorDB[It_LE] = new TIt_LE();
  IteratorDB[It_Finest] = new TIt_Finest();
  IteratorDB[It_EQLevel] = new TIt_EQLevel();
  IteratorDB[It_LELevel] = new TIt_LELevel();
  IteratorDB[It_Between] = new TIt_Between();
  IteratorDB[It_OCAF] = new TIt_OCAF();
  
  // this is a temporary solution to make sure this is called at the beginning
  QuadratureFormulaDatabase::create();

  // Initialization of the default parameters
  SetDefaultParameters();
  if(ParamFile)
  {
    //read the param file and fil the old database
    Output::info<4>("READ-IN", "Constructing old database from file ",
                    ParamFile);
    read_parameters(ParamFile);
  }
}

const TShapeDesc **TDatabase::ShapeDB = nullptr;
const TRefDesc   **TDatabase::RefDescDB = nullptr;
const TMapper    **TDatabase::MapperDB = nullptr;
TIterator  **TDatabase::IteratorDB = nullptr;
TParamDB   *TDatabase::ParamDB = nullptr;
TTimeDB    *TDatabase::TimeDB = nullptr;

// Methods
void TDatabase::SetDefaultParameters()
{
  ParamDB->VERSION = 1;

  ParamDB->MapperType = 1;

  ParamDB->RE_NR=1.0;
  ParamDB->FLOW_PROBLEM_TYPE = 0;
  ParamDB->OSEEN_ZERO_ORDER_COEFF = 0.0;

  ParamDB->PE_NR=1.0;  

  ParamDB->ANSATZ_ORDER = 2;

  ParamDB->VELOCITY_SPACE = 22;
  ParamDB->PRESSURE_SPACE = -4711;

  ParamDB->UPWIND_ORDER = 0;
  ParamDB->UPWIND_FLUX_DAMP = 1;
  ParamDB->UPWIND_APPLICATION = 0;
  ParamDB->NSTYPE = 1;
    
  ParamDB->SIGMA_PERM = 1;

  ParamDB->LAPLACETYPE = 0;
  ParamDB->USE_ISOPARAMETRIC = 1;
  ParamDB->CELL_MEASURE = 0;

  ParamDB->FACE_SIGMA = 1;
  ParamDB->WEAK_BC_SIGMA = 1;
  ParamDB->WEAK_BC = 0;
 
  ParamDB->DELTA0=1.0;
  ParamDB->DELTA1=1.0;
  ParamDB->SDFEM_TYPE=2;
  ParamDB->CIP_TYPE=0;

  ParamDB->FILTER_WIDTH_CONSTANT = 2;
  ParamDB->FILTER_WIDTH_POWER = 1;
  ParamDB->GAUSSIAN_GAMMA = 6;
  
  ParamDB->TURBULENT_VISCOSITY_TYPE = 1;
  ParamDB->TURBULENT_VISCOSITY_TENSOR = 0;
  ParamDB->TURBULENT_VISCOSITY_CONSTANT = 0.01;
  ParamDB->TURBULENT_VISCOSITY_POWER = 1;
  ParamDB->TURBULENT_VISCOSITY_SIGMA = 6;

  ParamDB->ARTIFICIAL_VISCOSITY_CONSTANT = 1; // parameters for VMS
  ParamDB->ARTIFICIAL_VISCOSITY_POWER = 1;

  ParamDB->FRICTION_CONSTANT = 0.0;      // free slip 
  ParamDB->FRICTION_POWER = 0.0;         // free slip
  ParamDB->FRICTION_TYPE = 0;            // friction type
  ParamDB->FRICTION_U0 = 1.0;            // U_0
  ParamDB->PENETRATION_CONSTANT = 1e12;  // no penetration
  ParamDB->PENETRATION_POWER = -2;        // no penetration

  ParamDB->DIV_DIV_STAB_TYPE = 0;        // stabilization for div-div term 
  ParamDB->DIV_DIV_STAB_C1 = 2;
  ParamDB->DIV_DIV_STAB_C2 = 1;

  ParamDB->LP_FULL_GRADIENT = 1;
  ParamDB->LP_STREAMLINE = 0;
  ParamDB->LP_DIVERGENCE = 0;
  ParamDB->LP_PRESSURE = 0;
  ParamDB->LP_COEFF_TYPE = 0;

  ParamDB->LP_FULL_GRADIENT_COEFF = 1.0;
  ParamDB->LP_STREAMLINE_COEFF= 1.0;
  ParamDB->LP_PRESSURE_COEFF = 1.0;

  ParamDB->LP_FULL_GRADIENT_EXPONENT = 1.0;
  ParamDB->LP_STREAMLINE_EXPONENT = 1.0;
  ParamDB->LP_PRESSURE_EXPONENT = 1.0;

  ParamDB->LP_ORDER_DIFFERENCE = 1;
  ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE = -123;
  ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE = -123;
  ParamDB->LP_PRESSURE_ORDER_DIFFERENCE = -123;
  
  ParamDB->LP_CROSSWIND = 0;
  ParamDB->LP_CROSSWIND_COEFF_TYPE = 1;
  ParamDB->LP_CROSSWIND_COEFF = 1.0;
  ParamDB->LP_CROSSWIND_EXPONENT = 1.0;

  /** parameters for SOLD schemes */
  ParamDB->SOLD_TYPE = 0;
  ParamDB->SOLD_PARAMETER_TYPE = 11;
  ParamDB->SOLD_CONST = 1.0;
  ParamDB->SOLD_POWER = 1.0;
  ParamDB->SOLD_S = 1.0;
  ParamDB->SOLD_U0 = 1.0;
  ParamDB->SOLD_PARAMETER_SCALING = 0;
  ParamDB->SOLD_PARAMETER_SCALING_FACTOR = 1.0;

  /** the following parameters are for individual use */
  ParamDB->P2 = 1.0;


  /** parameters for weakly imposing boundary/interface conditions */
  ParamDB->n_neumann_boundary = 0.;
  ParamDB->neumann_boundary_id.clear();
  ParamDB->neumann_boundary_value.clear();
    
  ParamDB-> n_g_v_boundary = 0.;
  ParamDB-> g_v_boundary_id.clear();
  ParamDB->g_v_boundary_value.clear();
    
  ParamDB-> n_unvn_boundary = 0.;
  ParamDB-> unvn_boundary_id.clear();
  ParamDB->unvn_boundary_value.clear();

  ParamDB-> n_gradunv_boundary = 0.;
  ParamDB-> gradunv_boundary_id.clear();
  ParamDB->gradunv_boundary_value.clear();
    
  ParamDB-> n_u_v_boundary = 0.;
  ParamDB-> u_v_boundary_id.clear();
  ParamDB->u_v_boundary_value.clear();
    
  ParamDB-> n_p_v_n_boundary = 0.;
  ParamDB-> p_v_n_boundary_id.clear();
  ParamDB-> p_v_n_boundary_value.clear();

    
  // Nitsche Combination - Weak Dirichlet Boundary Conditions 
  ParamDB-> n_nitsche_boundary = 0.;
  ParamDB-> nitsche_boundary_id.clear();
  ParamDB-> nitsche_penalty.clear();

  // initialize TimeDB
  TimeDB->CURRENTTIME = 0;
  TimeDB->CURRENTTIMESTEPLENGTH = 1;
  TimeDB->TIMESTEPLENGTH = 1;
  TimeDB->SCALE_DIVERGENCE_CONSTRAINT = -1.0;

  // parameters for implicit Euler method
  TimeDB->THETA1 = 1;
  TimeDB->THETA2 = 0;
  TimeDB->THETA3 = 0;
  TimeDB->THETA4 = 1;
  TimeDB->TIME_DISC = 2;

  TimeDB->RB_TYPE = 3;
  TimeDB->RB_TYPE2 = -1;
  TimeDB->RB_SSC = 0;
  TimeDB->RB_SSC_TOL = 1;
  TimeDB->RB_SSC_ALPHA = 0.9;
  TimeDB->RB_SSC_ALPHA_MIN = 1.0;
  TimeDB->RB_SSC_ALPHA_MAX = 1.0;
  TimeDB->RB_SSC_MAX_ERROR = 1.0;

  TimeDB->RB_APPROX_J = 0;
  TimeDB->RB_APPROX_C = 0;

  ParamDB->INPUT_QUAD_RULE = 0;
  ParamDB->INTERNAL_PROBLEM_LINEAR = 0;
  ParamDB->INTERNAL_SLIP_WITH_FRICTION = 0;
  ParamDB->INTERNAL_SLIP_WITH_FRICTION_IDENTITY = 0;
  ParamDB->INTERNAL_SLIP_WEAK_FORM = 0;// 1=weak penetr.(original),0=strong penetr.(like 2D)
  ParamDB->INTERNAL_QUAD_HEXA = 0;
  ParamDB->INTERNAL_QUAD_TETRA = 0;
  ParamDB->INTERNAL_QUAD_QUAD = 0;
  ParamDB->INTERNAL_QUAD_TRIA = 0;
  ParamDB->INTERNAL_QUAD_RULE = 0;
  ParamDB->INTERNAL_LOCAL_DOF = 0;
  ParamDB->INTERNAL_PERIODIC_IDENTITY = 0;
  ParamDB->INTERNAL_PROBLEM_IDENTITY = 0;
  ParamDB->INTERNAL_LEVEL = 0;
  ParamDB->INTERNAL_POLYNOMIAL_DEGREE = 1;
  ParamDB->INTERNAL_COERCIVITY = -4711.0;

  ParamDB->INTERNAL_P1_Array = nullptr;

  // ******** parameters for VMS *********//
  // constants in AdaptProjectionSpace 
  ParamDB->VMS_ADAPT_LOWER = 0.25;
  ParamDB->VMS_ADAPT_MIDDLE = 1.0;
  ParamDB->VMS_ADAPT_UPPER = 3.0;
  ParamDB->VMS_ADAPT_COMP = 1;

  ParamDB->SUPERCONVERGENCE_ORDER = 0;

// parameters for turbulent flow around a squared cylinder
  ParamDB->CYLINDER_22000_YPLUS_SIDES = 1500;
  ParamDB->CYLINDER_22000_YPLUS_FRONT = 2300;
  ParamDB->CYLINDER_22000_YPLUS_BACK  = 1000;

  /** parameters for individual use */
  ParamDB->P0 = 0.;
  ParamDB->P1 = 0.;
  ParamDB->P2 = 0.;
  ParamDB->P3 = 0.;
  ParamDB->P4 = 0.;
  ParamDB->P5 = 0.;
  ParamDB->P6 = 0.;
  ParamDB->P7 = 0.;
  ParamDB->P8 = 0.;
  ParamDB->P9 = 0.;
  ParamDB->P10 = 0.;
  ParamDB->P11 = 0.;
  ParamDB->P12 = 0.;
  ParamDB->P13 = 0.;
  ParamDB->P14 = 0.;
  ParamDB->P15 = 0.;
  
  
  /** parameters for individual use in parallel computations */
  ParamDB->Par_P0 = 0;
  ParamDB->Par_P1 = 0;
  ParamDB->Par_P2 = 0;
  ParamDB->Par_P3 = 1;
  ParamDB->Par_P4 = 0;
  ParamDB->Par_P5 = 10;
  ParamDB->Par_P6 = 1;
  ParamDB->Par_P7 = 0;
  ParamDB->Par_P8 = 0;
  ParamDB->Par_P9 = 0;

  #ifdef _MPI
  ParamDB->Comm = MPI_COMM_WORLD;    
  #endif
  return;
}

void TDatabase::WriteParamDB(char *ExecutedFile)
{
  using namespace Output; // printToFile
  
  printToFile(">>> Start printing old database (ParamDB) <<<");
  printToFile("HOSTNAME: ", utilities::get_host_name(), " started on ", 
              utilities::get_date_and_time());
  printToFile("EXECUTED FILE: ", ExecutedFile);
  printToFile("VERSION: ", ParamDB->VERSION);
  printToFile("MapperType: ", ParamDB->MapperType);
  
  printToFile("ANSATZ_ORDER: ", ParamDB->ANSATZ_ORDER);
  
  printToFile("VELOCITY_SPACE: ", ParamDB->VELOCITY_SPACE);
  printToFile("PRESSURE_SPACE: ", ParamDB->PRESSURE_SPACE);

  printToFile("UPWIND_ORDER: ", ParamDB->UPWIND_ORDER);
  printToFile("UPWIND_FLUX_DAMP: ", ParamDB->UPWIND_FLUX_DAMP);
  printToFile("UPWIND_APPLICATION: ", ParamDB->UPWIND_APPLICATION);
  printToFile("NSTYPE: ", ParamDB->NSTYPE);
  printToFile("SIGMA_PERM: ", ParamDB->SIGMA_PERM);
  printToFile("LAPLACETYPE: ", ParamDB->LAPLACETYPE);
  printToFile("USE_ISOPARAMETRIC: ", ParamDB->USE_ISOPARAMETRIC);
  printToFile("CELL_MEASURE: ", ParamDB->CELL_MEASURE);
  printToFile("FACE_SIGMA: ", ParamDB->FACE_SIGMA);
  printToFile("WEAK_BC_SIGMA: ", ParamDB->WEAK_BC_SIGMA);
  printToFile("WEAK_BC: ", ParamDB->WEAK_BC);

  printToFile("RE_NR: ", ParamDB->RE_NR);
  printToFile("FLOW_PROBLEM_TYPE: ", ParamDB->FLOW_PROBLEM_TYPE);
  printToFile("PE_NR: ", ParamDB->PE_NR);  
  printToFile("DELTA0: ", ParamDB->DELTA0);
  printToFile("SDFEM_TYPE: ", ParamDB->SDFEM_TYPE);
  printToFile("DELTA1: ", ParamDB->DELTA1);
  printToFile("CIP_TYPE: ", ParamDB->CIP_TYPE);

  printToFile("SOLD_TYPE: ", ParamDB->SOLD_TYPE);
  printToFile("SOLD_PARAMETER_TYPE: ", ParamDB->SOLD_PARAMETER_TYPE);
  printToFile("SOLD_CONST: ", ParamDB->SOLD_CONST);
  printToFile("SOLD_POWER: ", ParamDB->SOLD_POWER);
  printToFile("SOLD_S: ", ParamDB->SOLD_S);
  printToFile("SOLD_U0: ", ParamDB->SOLD_U0);
  printToFile("SOLD_PARAMETER_SCALING: ", ParamDB->SOLD_PARAMETER_SCALING);
  printToFile("SOLD_PARAMETER_SCALING_FACTOR: ", ParamDB->SOLD_PARAMETER_SCALING_FACTOR);
  

  printToFile("FILTER_WIDTH_CONSTANT: ", ParamDB->FILTER_WIDTH_CONSTANT);
  printToFile("FILTER_WIDTH_POWER: ", ParamDB->FILTER_WIDTH_POWER);
  printToFile("GAUSSIAN_GAMMA: ", ParamDB->GAUSSIAN_GAMMA );

  printToFile("TURBULENT_VISCOSITY_TYPE: ", ParamDB->TURBULENT_VISCOSITY_TYPE);
  printToFile("TURBULENT_VISCOSITY_TENSOR: ", ParamDB->TURBULENT_VISCOSITY_TENSOR);
  printToFile("TURBULENT_VISCOSITY_CONSTANT: ", ParamDB->TURBULENT_VISCOSITY_CONSTANT);
  printToFile("TURBULENT_VISCOSITY_POWER: ", ParamDB->TURBULENT_VISCOSITY_POWER);
  printToFile("TURBULENT_VISCOSITY_SIGMA: ", ParamDB->TURBULENT_VISCOSITY_SIGMA);

  printToFile("ARTIFICIAL_VISCOSITY_CONSTANT: ", ParamDB->ARTIFICIAL_VISCOSITY_CONSTANT);
  printToFile("ARTIFICIAL_VISCOSITY_POWER: ", ParamDB->ARTIFICIAL_VISCOSITY_POWER);

  printToFile("FRICTION_CONSTANT: ", ParamDB->FRICTION_CONSTANT);
  printToFile("FRICTION_POWER: ", ParamDB->FRICTION_POWER);
  printToFile("FRICTION_TYPE: ", ParamDB->FRICTION_TYPE);
  printToFile("FRICTION_U0: ", ParamDB->FRICTION_U0);
  printToFile("PENETRATION_CONSTANT: ", ParamDB->PENETRATION_CONSTANT);
  printToFile("PENETRATION_POWER: ", ParamDB->PENETRATION_POWER);
  printToFile("DIV_DIV_STAB_TYPE: ", ParamDB->DIV_DIV_STAB_TYPE); 
  printToFile("DIV_DIV_STAB_C1: ", ParamDB->DIV_DIV_STAB_C1); 
  printToFile("DIV_DIV_STAB_C2: ", ParamDB->DIV_DIV_STAB_C2); 
  printToFile("OSEEN_ZERO_ORDER_COEFF: ", ParamDB->OSEEN_ZERO_ORDER_COEFF);

  printToFile("LP_FULL_GRADIENT: ", ParamDB->LP_FULL_GRADIENT);
  printToFile("LP_STREAMLINE: ", ParamDB->LP_STREAMLINE);
  printToFile("LP_DIVERGENCE: ", ParamDB->LP_DIVERGENCE);
  printToFile("LP_PRESSURE: ", ParamDB->LP_PRESSURE);
  printToFile("LP_COEFF_TYPE: ", ParamDB->LP_COEFF_TYPE);
  printToFile("LP_FULL_GRADIENT_COEFF: ", ParamDB->LP_FULL_GRADIENT_COEFF);
  printToFile("LP_STREAMLINE_COEFF: ", ParamDB->LP_STREAMLINE_COEFF);
  printToFile("LP_PRESSURE_COEFF: ", ParamDB->LP_PRESSURE_COEFF);
  printToFile("LP_FULL_GRADIENT_EXPONENT: ", ParamDB->LP_FULL_GRADIENT_EXPONENT);
  printToFile("LP_STREAMLINE_EXPONENT: ", ParamDB->LP_STREAMLINE_EXPONENT);
  printToFile("LP_PRESSURE_EXPONENT: ", ParamDB->LP_PRESSURE_EXPONENT);
  printToFile("LP_ORDER_DIFFERENCE: ", ParamDB->LP_ORDER_DIFFERENCE);
  printToFile("LP_FULL_GRADIENT_ORDER_DIFFERENCE: ", ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE);
  printToFile("LP_STREAMLINE_ORDER_DIFFERENCE: ", ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE);
  printToFile("LP_PRESSURE_ORDER_DIFFERENCE: ", ParamDB->LP_PRESSURE_ORDER_DIFFERENCE);

  printToFile("P0: ", ParamDB->P0);
  printToFile("P1: ", ParamDB->P1);
  printToFile("P2: ", ParamDB->P2);
  printToFile("P3: ", ParamDB->P3);
  printToFile("P4: ", ParamDB->P4);
  printToFile("P5: ", ParamDB->P5);
  printToFile("P6: ", ParamDB->P6);
  printToFile("P7: ", ParamDB->P7);
  printToFile("P8: ", ParamDB->P8);
  printToFile("P9: ", ParamDB->P9);
  printToFile("P10: ", ParamDB->P10);
  printToFile("P11: ", ParamDB->P11);
  printToFile("P12: ", ParamDB->P12);
  printToFile("P13: ", ParamDB->P13);
  printToFile("P14: ", ParamDB->P14);
  printToFile("P15: ", ParamDB->P15);

  printToFile("VORTICITY THICKNESS FOR MIXING LAYER (P8): ", ParamDB->P8);


  // ******** parameters for VMS *********//
  // constants in AdaptProjectionSpace 
  printToFile("VMS_ADAPT_LOWER: ", ParamDB->VMS_ADAPT_LOWER); 
  printToFile("VMS_ADAPT_MIDDLE: ", ParamDB->VMS_ADAPT_MIDDLE); 
  printToFile("VMS_ADAPT_UPPER: ", ParamDB->VMS_ADAPT_UPPER); 
  printToFile("VMS_ADAPT_COMP: ", ParamDB->VMS_ADAPT_COMP); 

  printToFile("SUPERCONVERGENCE_ORDER: ", ParamDB->SUPERCONVERGENCE_ORDER);
  printToFile("INPUT_QUAD_RULE: ", ParamDB->INPUT_QUAD_RULE);
}

void TDatabase::WriteTimeDB()
{
  using namespace Output; // printToFile
  
  printToFile("CURRENTTIME: ", TimeDB->CURRENTTIME);
  printToFile("CURRENTTIMESTEPLENGTH: ", TimeDB->CURRENTTIMESTEPLENGTH);
  printToFile("TIMESTEPLENGTH: ", TimeDB->TIMESTEPLENGTH);
  printToFile("SCALE_DIVERGENCE_CONSTRAINT: ", TimeDB->SCALE_DIVERGENCE_CONSTRAINT);

  printToFile("THETA1: ", TimeDB->THETA1);
  printToFile("THETA2: ", TimeDB->THETA2);
  printToFile("THETA3: ", TimeDB->THETA3);
  printToFile("THETA4: ", TimeDB->THETA4);

  printToFile("TIME_DISC: ", TimeDB->TIME_DISC);

  printToFile("RB_TYPE: ", TimeDB->RB_TYPE);
  printToFile("RB_TYPE2: ", TimeDB->RB_TYPE2);

  printToFile(">>> End printing old database (ParamDB) <<<");
  printToFile("");
    
//  printToFile(" n_neumann_boundary: ",ParamDB->n_neumann_boundary );
//  printToFile(" neumann_boundary_id: ", ParamDB->neumann_boundary_id);
//  printToFile("neumann_boundary_value: ", ParamDB->neumann_boundary_value);
}

//TDatabase::~TDatabase()
void TDatabase::destroy()
{
  // allocate databases
  delete ParamDB;
  delete TimeDB;

  // initialize shape descriptors
  delete ShapeDB[S_Line];
  delete ShapeDB[Triangle];
  delete ShapeDB[Quadrangle];
  delete ShapeDB[Parallelogram];
  delete ShapeDB[Rectangle];
  #ifdef __3D__
    delete ShapeDB[Tetrahedron];
    delete ShapeDB[Hexahedron];
    delete ShapeDB[Brick];
  #endif
  delete [] ShapeDB;
   
  delete RefDescDB[S_Line];
  delete RefDescDB[Triangle];
  delete RefDescDB[Quadrangle];
  delete RefDescDB[Parallelogram];
  delete RefDescDB[Rectangle];
  #ifdef __3D__
    delete RefDescDB[Tetrahedron];
    delete RefDescDB[Hexahedron];
    delete RefDescDB[Brick];
  #endif
  delete RefDescDB[N_SHAPES + LineReg];
  delete RefDescDB[N_SHAPES + TriReg];
  delete RefDescDB[N_SHAPES + TriBary];
  delete RefDescDB[N_SHAPES + TriBis0];
  delete RefDescDB[N_SHAPES + TriBis1];
  delete RefDescDB[N_SHAPES + TriBis2];
  delete RefDescDB[N_SHAPES + TriBis01];
  delete RefDescDB[N_SHAPES + TriBis02];
  delete RefDescDB[N_SHAPES + TriBis10];
  delete RefDescDB[N_SHAPES + TriBis12];
  delete RefDescDB[N_SHAPES + TriBis20];
  delete RefDescDB[N_SHAPES + TriBis21];
  delete RefDescDB[N_SHAPES + QuadReg];
  delete RefDescDB[N_SHAPES + ParallReg];
  delete RefDescDB[N_SHAPES + RectReg];
  delete RefDescDB[N_SHAPES + QuadBis0];
  delete RefDescDB[N_SHAPES + QuadBis1];
  delete RefDescDB[N_SHAPES+Quad1Conf0];
  delete RefDescDB[N_SHAPES+Quad1Conf1];
  delete RefDescDB[N_SHAPES+Quad1Conf2];
  delete RefDescDB[N_SHAPES+Quad1Conf3];
  delete RefDescDB[N_SHAPES+Quad2Conf0];
  delete RefDescDB[N_SHAPES+Quad2Conf1];
  delete RefDescDB[N_SHAPES+Quad2Conf2];
  delete RefDescDB[N_SHAPES+Quad2Conf3];
  delete RefDescDB[N_SHAPES+QuadToTri0];
  delete RefDescDB[N_SHAPES+QuadToTri1];

  #ifdef __3D__
    delete RefDescDB[N_SHAPES + TetraReg];
    delete RefDescDB[N_SHAPES + TetraBary];
    delete RefDescDB[N_SHAPES + TetraReg0];
    delete RefDescDB[N_SHAPES + TetraReg1];
    delete RefDescDB[N_SHAPES + TetraReg2];
    delete RefDescDB[N_SHAPES + TetraBis0];
    delete RefDescDB[N_SHAPES + TetraBis1];
    delete RefDescDB[N_SHAPES + TetraBis2];
    delete RefDescDB[N_SHAPES + TetraBis3];
    delete RefDescDB[N_SHAPES + TetraBis4];
    delete RefDescDB[N_SHAPES + TetraBis5];
    delete RefDescDB[N_SHAPES + TetraBis01];
    delete RefDescDB[N_SHAPES + TetraBis02];
    delete RefDescDB[N_SHAPES + TetraBis03];
    delete RefDescDB[N_SHAPES + TetraBis04];
    delete RefDescDB[N_SHAPES + TetraBis05];
    delete RefDescDB[N_SHAPES + TetraBis10];
    delete RefDescDB[N_SHAPES + TetraBis12];
    delete RefDescDB[N_SHAPES + TetraBis13];
    delete RefDescDB[N_SHAPES + TetraBis14];
    delete RefDescDB[N_SHAPES + TetraBis15];
    delete RefDescDB[N_SHAPES + TetraBis20];
    delete RefDescDB[N_SHAPES + TetraBis21];
    delete RefDescDB[N_SHAPES + TetraBis23];
    delete RefDescDB[N_SHAPES + TetraBis24];
    delete RefDescDB[N_SHAPES + TetraBis25];
    delete RefDescDB[N_SHAPES + TetraBis30];
    delete RefDescDB[N_SHAPES + TetraBis32];
    delete RefDescDB[N_SHAPES + TetraBis34];
    delete RefDescDB[N_SHAPES + TetraBis35];
    delete RefDescDB[N_SHAPES + TetraBis40];
    delete RefDescDB[N_SHAPES + TetraBis41];
    delete RefDescDB[N_SHAPES + TetraBis43];
    delete RefDescDB[N_SHAPES + TetraBis45];
    delete RefDescDB[N_SHAPES + TetraBis51];
    delete RefDescDB[N_SHAPES + TetraBis52];
    delete RefDescDB[N_SHAPES + TetraBis53];
    delete RefDescDB[N_SHAPES + TetraBis54];
    delete RefDescDB[N_SHAPES + TetraQuad0];
    delete RefDescDB[N_SHAPES + TetraQuad1];
    delete RefDescDB[N_SHAPES + TetraQuad2];
    delete RefDescDB[N_SHAPES + TetraQuad3];
    delete RefDescDB[N_SHAPES + HexaReg];
    delete RefDescDB[N_SHAPES + BrickReg];
  #endif
  delete [] RefDescDB;
  
  #ifdef __3D__
    //initialize mapper
    delete MapperDB[MapTriReg0];
    delete MapperDB[MapTriReg1];
    delete MapperDB[MapTriReg2];

    delete MapperDB[MapTriBis00];
    delete MapperDB[MapTriBis01];
    delete MapperDB[MapTriBis02];
    delete MapperDB[MapTriBis10];
    delete MapperDB[MapTriBis11];
    delete MapperDB[MapTriBis12];
    delete MapperDB[MapTriBis20];
    delete MapperDB[MapTriBis21];
    delete MapperDB[MapTriBis22];
    delete MapperDB[MapTriBis010];
    delete MapperDB[MapTriBis011];
    delete MapperDB[MapTriBis012];
    delete MapperDB[MapTriBis020];
    delete MapperDB[MapTriBis021];
    delete MapperDB[MapTriBis022];
    delete MapperDB[MapTriBis100];
    delete MapperDB[MapTriBis101];
    delete MapperDB[MapTriBis102];
    delete MapperDB[MapTriBis120];
    delete MapperDB[MapTriBis121];
    delete MapperDB[MapTriBis122];
    delete MapperDB[MapTriBis200];
    delete MapperDB[MapTriBis201];
    delete MapperDB[MapTriBis202];
    delete MapperDB[MapTriBis210];
    delete MapperDB[MapTriBis211];
    delete MapperDB[MapTriBis212];

    delete MapperDB[MapQuadReg0];
    delete MapperDB[MapQuadReg1];
    delete MapperDB[MapQuadReg2];
    delete MapperDB[MapQuadReg3];
  #endif
  delete [] MapperDB;

  // initialize iterators
  delete IteratorDB[It_EQ];
  delete IteratorDB[It_LE];
  delete IteratorDB[It_Finest];
  delete IteratorDB[It_EQLevel];
  delete IteratorDB[It_LELevel];
  delete IteratorDB[It_Between];
  delete IteratorDB[It_OCAF];
  delete [] IteratorDB;
  
  // this is a temporary solution to make sure this is called at the end
  QuadratureFormulaDatabase::destroy();
}

void TDatabase::read_parameters(const char* ParamFile)
{
  ParamDB->read_parameters(ParamFile);
  TimeDB->read_parameters(ParamFile);
}

void TTimeDB::read_parameters(const char* ParamFile)
{
  int rank = 0;
#ifdef _MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  std::string line;
  int N_Param = 0;
  std::ifstream dat(ParamFile);

  if (!dat)
  {
    ErrThrow("cannot open '", ParamFile, "' for input");
  }

  while(std::getline(dat, line))
  {
    // read in parameter for time discretization
    if(SetParameters(line, "STEPLENGTH:", TIMESTEPLENGTH)){N_Param++; continue;}
    if(SetParameters(line, "TIMESTEPLENGTH:", TIMESTEPLENGTH)){N_Param++; continue;}
    if(SetParameters(line, "SCALE_DIVERGENCE_CONSTRAINT:", SCALE_DIVERGENCE_CONSTRAINT)){N_Param++; continue;}
    if(SetParameters(line, "TIME_DISC:", TIME_DISC)){N_Param++; continue;}
    if(SetParameters(line, "RB_TYPE:", RB_TYPE)){N_Param++; continue;}
    if(SetParameters(line, "RB_TYPE2:", RB_TYPE2)){N_Param++; continue;}
    if(SetParameters(line, "STEPSIZECONTROL:", RB_SSC)){N_Param++; continue;}
    if(SetParameters(line, "RB_SSC_TOL:", RB_SSC_TOL)){N_Param++; continue;}
    if(SetParameters(line, "RB_SSC_ALPHA:", RB_SSC_ALPHA)){N_Param++; continue;}
    if(SetParameters(line, "RB_SSC_ALPHA_MAX:", RB_SSC_ALPHA_MAX)){N_Param++; continue;}
    if(SetParameters(line, "RB_SSC_ALPHA_MIN:", RB_SSC_ALPHA_MIN)){N_Param++; continue;}
    if(SetParameters(line, "RB_SSC_MAX_ERROR:", RB_SSC_MAX_ERROR)){N_Param++; continue;}
    if(SetParameters(line, "RB_APPROX_J:", RB_APPROX_J)){N_Param++; continue;}
    if(SetParameters(line, "RB_APPROX_C:", RB_APPROX_C)){N_Param++; continue;}
    if(SetParameters(line, "RB_APPROX_STEPS:", RB_APPROX_STEPS)){N_Param++; continue;}
  }


  if(rank==0)
    Output::info("read_parameters","Time database (old) read with ", N_Param,
                 " parameters.");
}


TParamDB::~TParamDB()
{
  Output::print<4>("deleted parameter database");
}


/* ************************************************************************** */
/// plays the role of the old CheckParameterConsistencyNSE
/// but with the new database.The remaining global parameters
/// can/must be removed here, progressively.
void check_parameters_consistency_NSE(ParameterDatabase& db)
{
#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#else
  int my_rank = 0;
#endif
  if(!db.contains("space_discretization_type"))
  {
    ErrThrow("check_parameters_consistency_NSE: you need to set the "
             "space_discretization_type");
  }

  if (db["space_discretization_type"].is("sdfem") && (TDatabase::ParamDB->NSTYPE==1))
  {
      //TDatabase::ParamDB->NSTYPE = 2;
      //Output::info("NSE Parameter Consistency","NSTYPE changed from 1 to 2 because of SDFEM discretization ");
    if(my_rank==0)
        Output::info("NSE Parameter Consistency","NSTYPE 1: only reduced SDFEM, only for 2D, fixed point, not skew !!!");
  }

  if ((db["space_discretization_type"].is("sdfem")) && (TDatabase::ParamDB->NSTYPE==3))
  {
    TDatabase::ParamDB->NSTYPE = 4;
    if(my_rank==0)
      Output::warn<1>("NSE Parameter Consistency","NSTYPE changed from 3 to 4 because of SDFEM discretization ");
  }

  // rotational form
  Parameter nonlin_form(db["nse_nonlinear_form"]);
  if(nonlin_form.is("rotational"))
  {
    if (TDatabase::ParamDB->NSTYPE<=2)
    {
      TDatabase::ParamDB->NSTYPE+=2;
      if(my_rank==0)
      {
        Output::warn<1>("NSE Parameter Consistency",
                        "NSTYPE changed to ", TDatabase::ParamDB->NSTYPE);
        Output::warn<1>("NSE Parameter Consistency",
                        " because of nse_nonlinear_form = ", nonlin_form);
      }
    }

    // change 'db["discretization_type]" for internal reasons
    if ( db["space_discretization_type"].is("galerkin") )
    {
      db["space_discretization_type"].set("smagorinsky");
      TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE = 0;
      if(my_rank==0)
        Output::warn<1>("NSE Parameter Consistency","discretization_type changed to 'smagorinsky' (4)"
            " for internal reasons, turbulent viscosity is switched off.");
    }

  }

  if (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE==5)
  {
    if (!(db["space_discretization_type"].is("vms_projection")||
        db["space_discretization_type"].is("vms_projection_expl")))
    {
      if(my_rank==0)
      {
        Output::warn<1>("NSE Parameter Consistency","TURBULENT_VISCOSITY_TYPE = 5 only defined for projection-based VMS methods");
        Output::warn<1>("NSE Parameter Consistency","Set different TURBULENT_VISCOSITY_TYPE !!!");
        ErrThrow("BOOM!");
      }
    }
  }

  // LOCAL_PROJECTION
  if (db["space_discretization_type"].is("local_projection"))
  {
    if (TDatabase::ParamDB->LP_FULL_GRADIENT)
    {
      if (TDatabase::ParamDB->LP_STREAMLINE)
      {
        TDatabase::ParamDB->LP_STREAMLINE = 0;
        if(my_rank==0)
        Output::warn<1>("NSE Parameter Consistency","LP_STREAMLINE changed to ", TDatabase::ParamDB->LP_STREAMLINE,
                      " due to LP_FULL_GRADIENT = ", TDatabase::ParamDB->LP_FULL_GRADIENT);
      }

      if (TDatabase::ParamDB->LP_DIVERGENCE)
      {
        TDatabase::ParamDB->LP_DIVERGENCE = 0;
        if(my_rank==0)
        Output::warn<1>("NSE Parameter Consistency","LP_DIVERGENCE changed to ", TDatabase::ParamDB->LP_DIVERGENCE,
                      " due to LP_FULL_GRADIENT = ", TDatabase::ParamDB->LP_FULL_GRADIENT);
      }
    } // end LP_FULL_GRADIENT

    if (TDatabase::ParamDB->LP_DIVERGENCE)
    {
      if (TDatabase::ParamDB->NSTYPE<=2)
      {
        TDatabase::ParamDB->NSTYPE+=2;
        if(my_rank==0)
        Output::warn<1>("NSE Parameter Consistency","NSTYPE changed to ", TDatabase::ParamDB->NSTYPE,
                      " LP_DIVERGENCE = ", TDatabase::ParamDB->LP_DIVERGENCE);
      }
    }

    if(TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE == -123)
      TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE = TDatabase::ParamDB->LP_ORDER_DIFFERENCE;

    if(TDatabase::ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE == -123)
      TDatabase::ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE = TDatabase::ParamDB->LP_ORDER_DIFFERENCE;

    if(TDatabase::ParamDB->LP_PRESSURE_ORDER_DIFFERENCE == -123)
      TDatabase::ParamDB->LP_PRESSURE_ORDER_DIFFERENCE = TDatabase::ParamDB->LP_ORDER_DIFFERENCE;

  } // end "discretization_type" = "LOCAL_PROJECTION"
  else
  {
    // switch off all local projection terms
    TDatabase::ParamDB->LP_FULL_GRADIENT = 0;
    TDatabase::ParamDB->LP_FULL_GRADIENT_COEFF = 0;
    TDatabase::ParamDB->LP_FULL_GRADIENT_EXPONENT = 1;

    TDatabase::ParamDB->LP_STREAMLINE = 0;
    TDatabase::ParamDB->LP_STREAMLINE_COEFF = 0;
    TDatabase::ParamDB->LP_STREAMLINE_EXPONENT = 1;

    TDatabase::ParamDB->LP_DIVERGENCE = 0;

    TDatabase::ParamDB->LP_PRESSURE = 0;
    TDatabase::ParamDB->LP_PRESSURE_COEFF = 0;
    TDatabase::ParamDB->LP_PRESSURE_EXPONENT = 1;

    TDatabase::ParamDB->LP_ORDER_DIFFERENCE = 1;
    TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE = 1;
    TDatabase::ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE = 1;
    TDatabase::ParamDB->LP_PRESSURE_ORDER_DIFFERENCE = 1;
  }

  if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE == STOKES)
    TDatabase::ParamDB->INTERNAL_PROBLEM_LINEAR = 1;
  if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE == OSEEN)
  {
    TDatabase::ParamDB->INTERNAL_PROBLEM_LINEAR = 1;
    switch (TDatabase::ParamDB->NSTYPE)
    {
      case 1:
        if(my_rank==0)
            Output::warn<1>("NSE Parameter Consistency","Galerkin discretization for Oseen because of NSTYPE ",
                      TDatabase::ParamDB->NSTYPE);
        db["space_discretization_type"].set("galerkin"); // DISCTYPE = 1
        break;
      case 2:
      case 4:
        db["space_discretization_type"].set("galerkin"); // DISCTYPE = 1
        break;
      case 14:
      if(my_rank==0)
          Output::warn<1>("NSE Parameter Consistency","SUPG/PSPG/grad-div discretization for Oseen because of NSTYPE ",
                      TDatabase::ParamDB->NSTYPE);
        db["space_discretization_type"].set("supg"); // DISCTYPE = 2
        break;
      default:
        if(my_rank==0)
          Output::warn<1>("NSE Parameter Consistency",
                          "No method for Oseen implemented for NSTYPE ",
                          TDatabase::ParamDB->NSTYPE);
        exit(4711);
    }
  }
}
