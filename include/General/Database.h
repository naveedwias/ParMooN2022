// =======================================================================
// @(#)Database.h        1.23 06/27/00
//
// Class:       TDatabase
// Purpose:     database of needed refinement, mapping and
//              shape descriptors as well as iterators
//              parameter database
//
// Author:      Volker Behns  24.07.97
//
// =======================================================================

//class TDatabase;

#ifdef _MPI
#  include "mpi.h"
#endif

#ifndef __DATABASE__
#define __DATABASE__

#include <Iterator.h>
#include <RefDesc.h>
#include <Mapper.h>

// forward declaration
class ParameterDatabase;

// typedef struct Foo { ... } Foo; --> see https://stackoverflow.com/a/612350
typedef struct TParamDB
{
  int VERSION;

  //======================================================================
  /** parameters for setting finite element spaces                      */
  //======================================================================
  int ANSATZ_ORDER;

  int VELOCITY_SPACE;
  int PRESSURE_SPACE;

  //======================================================================
  /** parameters for setting the discretization                         */
  //======================================================================
  // general
  int USE_ISOPARAMETRIC;
  int CELL_MEASURE;
  
  /** upwind methods */
  int UPWIND_ORDER;
  double UPWIND_FLUX_DAMP;
  int UPWIND_APPLICATION;

  /** stabilization with face integrals */
  int WEAK_BC;
  double FACE_SIGMA;
  double WEAK_BC_SIGMA;

  /** SUPG or SDFEM method */
  int    SDFEM_TYPE;
  double DELTA0;
  double DELTA1;
  int    CIP_TYPE;
  /** parameters for SOLD methods */
  int SOLD_TYPE;
  int SOLD_PARAMETER_TYPE;
  double SOLD_CONST;
  double SOLD_POWER;
  double SOLD_S;
  double SOLD_U0;
  int SOLD_PARAMETER_SCALING;
  double SOLD_PARAMETER_SCALING_FACTOR;

  //======================================================================
  /** parameter for local projection stabilization */
  //======================================================================
  int LP_FULL_GRADIENT;
  int LP_STREAMLINE;
  int LP_DIVERGENCE;
  int LP_PRESSURE;
  int LP_CROSSWIND;
  int LP_COEFF_TYPE;

  double LP_FULL_GRADIENT_COEFF;
  double LP_STREAMLINE_COEFF;
  double LP_PRESSURE_COEFF;

  double LP_FULL_GRADIENT_EXPONENT;
  double LP_STREAMLINE_EXPONENT;
  double LP_PRESSURE_EXPONENT;

  int LP_ORDER_DIFFERENCE;
  int LP_FULL_GRADIENT_ORDER_DIFFERENCE;
  int LP_STREAMLINE_ORDER_DIFFERENCE;
  int LP_PRESSURE_ORDER_DIFFERENCE;
  
  int LP_CROSSWIND_COEFF_TYPE;
  double LP_CROSSWIND_COEFF;
  double LP_CROSSWIND_EXPONENT; 

  /** parameter for superconvergence */
  int SUPERCONVERGENCE_ORDER;

  //======================================================================
  /** PARAMETERS FOR STOKES AND NAVIER-STOKES PROBLEMS                  */
  //======================================================================

  /** general parameters */
  double RE_NR; //FIXME This parameter must really be removed from global scope!
  int FLOW_PROBLEM_TYPE;
  int NSTYPE;
  int LAPLACETYPE;
  double OSEEN_ZERO_ORDER_COEFF;

  //======================================================================
  /** PARAMETERS FOR DARCY PROBLEM                  */
  //======================================================================
  double SIGMA_PERM;
  //======================================================================
    
  double PE_NR;
  
  /** parameters for LES */
  double FILTER_WIDTH_CONSTANT;
  double FILTER_WIDTH_POWER;
  double GAUSSIAN_GAMMA;

  /** parameters for turbulent viscosity */
  int TURBULENT_VISCOSITY_TYPE;
  int TURBULENT_VISCOSITY_TENSOR;
  double TURBULENT_VISCOSITY_CONSTANT;
  double TURBULENT_VISCOSITY_POWER;
  double TURBULENT_VISCOSITY_SIGMA;

  // constants in AdaptProjectionSpace
  double VMS_ADAPT_LOWER;
  double VMS_ADAPT_MIDDLE;
  double VMS_ADAPT_UPPER;
  int VMS_ADAPT_COMP;

  double ARTIFICIAL_VISCOSITY_CONSTANT;
  double ARTIFICIAL_VISCOSITY_POWER;

  //======================================================================
  /** parameters for slip with friction and penetration with resistance
      boundary conditions                                               */
  //======================================================================
  double FRICTION_CONSTANT;
  double FRICTION_POWER;
  int    FRICTION_TYPE;
  double FRICTION_U0;
  double PENETRATION_CONSTANT;
  double PENETRATION_POWER;

  //======================================================================
  /** parameters for div-div stabilization */
  //======================================================================
  int    DIV_DIV_STAB_TYPE;
  double DIV_DIV_STAB_C1;
  double DIV_DIV_STAB_C2;

  //======================================================================
  /** the following parameters are for individual use */
  //======================================================================
  double P0;
  double P1;
  double P2;
  double P3;
  double P4;
  double P5;
  double P6;
  double P7;
  double P8;
  double P9;
  double P10;
  double P11;
  double P12;
  double P13;
  double P14;
  double P15;
  //======================================================================
  /** PARAMETERS FOR APPLICATIONS */
  //======================================================================

  //======================================================================
  /** parameters for turbulent flow around a squared cylinder */
  //======================================================================
  double CYLINDER_22000_YPLUS_SIDES;
  double CYLINDER_22000_YPLUS_FRONT;
  double CYLINDER_22000_YPLUS_BACK;

  //======================================================================
  /** internal parameters
  cannot be set in the readin file
  are used as global variables */
  //======================================================================
  int    INTERNAL_PROBLEM_LINEAR;
  int    INTERNAL_SLIP_WITH_FRICTION;
  int    INTERNAL_SLIP_WITH_FRICTION_IDENTITY;
  int    INTERNAL_SLIP_WEAK_FORM;
  int    INPUT_QUAD_RULE;
  int    INTERNAL_QUAD_HEXA;
  int    INTERNAL_QUAD_TETRA;
  int    INTERNAL_QUAD_QUAD;
  int    INTERNAL_QUAD_TRIA;
  int    INTERNAL_QUAD_RULE;
  int    INTERNAL_LOCAL_DOF;
  int    INTERNAL_PERIODIC_IDENTITY;
  int    INTERNAL_PROBLEM_IDENTITY;
  int    INTERNAL_LEVEL;
  int    INTERNAL_POLYNOMIAL_DEGREE;
  int    INTERNAL_MESH_CELL_TYPE;
  double INTERNAL_VERTEX_X[8];
  double INTERNAL_VERTEX_Y[8];
  double INTERNAL_VERTEX_Z[8];
  double INTERNAL_HK_CONVECTION;
  int    INTERNAL_MOMENT;
  
  double INTERNAL_COERCIVITY;
  double *INTERNAL_P1_Array;
  int    INTERNAL_ARRAY_LENGTH;
  
  //======================================================================
  /** parameters for individual use in parallel computations */
  //======================================================================
  int Par_P0; //out rank
  int Par_P1; // 1 - root takes part in computation; 0 - not
  int Par_P2; // mesh partition type: 1 - dual; 0 - nodal
  int Par_P3; // 1 - use halocells; 0 - dont
  int Par_P4; // 1-redistribution of masters 0-otherwise
  int Par_P5;
  int Par_P6;
  int Par_P7;
  int Par_P8;
  int Par_P9;
  int MapperType;

  //======================================================================
  /** parameters for weakly imposing boundary/interface conditions for 2D Brinkman problems  */
  //======================================================================
    // number of boundary components with neumann conditions
    int n_neumann_boundary;
    // ID's of boundary components with neumann conditions
    std::vector<int> neumann_boundary_id;
    // factor for boundary integrals
    std::vector<double> neumann_boundary_value;
    
    int n_g_v_boundary;
    std::vector<int> g_v_boundary_id;
    std::vector<double> g_v_boundary_value;
    
    int n_unvn_boundary;
    std::vector<int> unvn_boundary_id;
    std::vector<double> unvn_boundary_value;
    
    int n_graduvn_boundary;
    std::vector<int> graduvn_boundary_id;
    std::vector<double> graduvn_boundary_value;
    
    int n_gradunv_boundary;
    std::vector<int> gradunv_boundary_id;
    std::vector<double> gradunv_boundary_value;
    
    int n_u_v_boundary;
    std::vector<int> u_v_boundary_id;
    std::vector<double> u_v_boundary_value;
    
    int n_p_v_n_boundary;
    std::vector<int> p_v_n_boundary_id;
    std::vector<double> p_v_n_boundary_value;
    

    
    // Nitsche Combination - Weak Dirichlet Boundary Conditions
    int n_nitsche_boundary;
    std::vector<int> nitsche_boundary_id;
    std::vector<double> nitsche_penalty;
    
    
  //======================================================================
    
  #ifdef _MPI
  /** MPI_Comm for which the computation is started (should not be changed during computation)*/
  MPI_Comm Comm;
 #endif
  
  TParamDB() = default;
  ~TParamDB();
  
  void read_parameters(const char* ParamFile);
} TParamDB;

typedef struct TTimeDB
{
  double CURRENTTIME;
  double CURRENTTIMESTEPLENGTH;

  double base_time_step_length;
  double TIMESTEPLENGTH;
  double SCALE_DIVERGENCE_CONSTRAINT;

  // control parameter
  double THETA1;
  double THETA2;
  double THETA3;
  double THETA4;

  int TIME_DISC;

  // parameters for Rosenbrock methods
  int RB_TYPE;
  int RB_TYPE2;
  int RB_SSC;
  double RB_SSC_TOL;
  double RB_SSC_ALPHA;
  double RB_SSC_ALPHA_MAX;
  double RB_SSC_ALPHA_MIN;
  double RB_SSC_MAX_ERROR;

  int RB_APPROX_J;
  int RB_APPROX_C;
  int RB_APPROX_STEPS;

  double RB_GAMMA_I;
  double RB_GAMMA_II;
  double RB_ALPHA_I;
  double RB_SIGMA_I;
  double RB_A_IJ[10];
  double RB_C_IJ[10];
  double RB_S_IJ[10];
  double RB_M_I;
  double RB_MS_I;

  void read_parameters(const char* ParamFile);
} TTimeDB;

/** database of needed refinement, mapping and
    shape descriptors as well as iterators */
namespace TDatabase
{
  /** database of shape descriptors */
  extern const TShapeDesc **ShapeDB;

  /** database of refinement descriptors */
  extern const TRefDesc **RefDescDB;

  /** database of mapper */
  extern const TMapper **MapperDB;

  /** database of iterators */
  extern TIterator **IteratorDB;

  /** general parameters */
  extern TParamDB *ParamDB;

  /** parameter for time discretization */
  extern TTimeDB *TimeDB;

  // Constructors
  /** @brief initialize the database.
   * 
   * If ParamFile is nullptr, all parameters are set to their default values
   * using the method 'SetDefaultParameters'. If a valid file is given as 
   * 'ParamFile', the respective default values are overwritten.
   * 
   * @note you have to call this once in basically every ParMooN program
   */
  void create(const char* ParamFile = nullptr);
  
  //~TDatabase();
  /** 
   * @brief free the memory of all databases in this namespace.
   * 
   * @note if you called create(), you also have to call this.
   */
  void destroy();

  // set default parameters
  void SetDefaultParameters();

  void WriteParamDB(char* ExecutedFile);

  void WriteTimeDB();
  
  void read_parameters(const char* ParamFile);
}


/** TODO: these checks still use intensively the global parameters
 * and are used in all the NSE programs. Finish to replace
 * all global dependencies by local ones, and find a relevant place
 * for these checks in the program. Their declaration and definition
 * in Database.C and .h are just temporary. */
void check_parameters_consistency_NSE(ParameterDatabase& db);

#endif
