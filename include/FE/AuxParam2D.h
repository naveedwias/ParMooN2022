// =======================================================================
// @(#)AuxParam2D.h        1.1 10/30/98
// 
// Class:       TAuxParam2D
// Purpose:     store parameter functions and FE functions
//
// Author:      Gunar Matthies (06.08.98)
//
// History:     start of implementation 06.08.98 (Gunar Matthies)
//
// =======================================================================

#ifndef __AUXPARAM2D__
#define __AUXPARAM2D__

class TAuxParam2D;
class TQuadFormula;

#include <Constants.h>
#include <FESpace2D.h>
#include <FEFunction2D.h>
#include <string>

/** store parameter functions and FE functions */
class TAuxParam2D
{
  public:
// =======================================================================
//  numbers of stored objects
// =======================================================================
    /** number of stored parameter function (ParamFct) */
    int N_ParamFct;

// =======================================================================
//  array of pointers to stored objects
// =======================================================================
    /** array of stored FEFunction2D */
    TFEFunction2D **FEFunctions2D;

    /** array of stored parameter function */
    ParamFct **ParameterFct;

// =======================================================================
//  information of FE values used by parameter functions
// =======================================================================
    /** number of FE values */
    int N_FEValues;

    /** index of FEFunction2D used for FE value i */
    int *FEValue_FctIndex;

    /** which multiindex is used for FE value i */
    MultiIndex2D *FEValue_MultiIndex;

// =======================================================================
//  information of parameter functions
// =======================================================================
    /** number of all parameters */
    int N_Parameters;

    /** index of first parameter produced by parameter function i */
    int *BeginParameter;

// =======================================================================
//  information of parameter functions
// =======================================================================
    /** storage for temporary FE values */
    double *Temp;

    const double **Values;
    double ***OrigValues;
    const int **Index;
    int *N_BaseFunct;

  public:
    /** constructor */
    TAuxParam2D(int n_paramfct, int n_fevalues, TFEFunction2D **fefunctions2d,
                ParamFct **parameterfct, int *fevalue_fctindex,
                MultiIndex2D *fevalue_multiindex, int n_parameters,
                int *beginparameter);

    /** @brief standard constructor
     * 
     * If you don't need values of a finite element function in your assembling,
     * choose this constructor. This is equivalent to calling 
     * TAuxParam2D(0, 0, 0, 0, nullptr, nullptr, nullptr, nullptr, nullptr, 0, nullptr);
     */
    TAuxParam2D();

    
    /** destructor */
    ~TAuxParam2D();

    /** return all parameters at all quadrature points */
    void GetParameters(const TQuadFormula& qf, int cellnum,
                       double **Parameters);

    int GetN_Parameters() const
    { return N_Parameters; }
    int GetN_ParamFct() const
    { return N_ParamFct; }
    ParamFct * get_ParameterFct(int i) const
    { return ParameterFct[i]; }
    
    int get_BeginParameter(int i) const
    { return BeginParameter[i]; }

    int get_N_FEValues() const
    { return N_FEValues; }
    
    TFEFunction2D **get_FEFunctions2D() const
    { return FEFunctions2D; }
    
    int get_FEValue_FctIndex(int i) const
    { return FEValue_FctIndex[i]; }
    
    MultiIndex2D get_FEValue_MultiIndex(int i) const
    { return FEValue_MultiIndex[i]; }
};

// standard function to use for Navier-Stokes
void Velocity_Fct(double *inputList, double *outputValues);


#endif
