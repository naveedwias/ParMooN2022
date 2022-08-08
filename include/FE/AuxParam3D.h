// =======================================================================
// %W% %G%
// 
// Class:       TAuxParam3D
// Purpose:     store parameter functions and FE functions
//
// Author:      Gunar Matthies (06.08.98)
//
// History:     start of implementation 06.08.98 (Gunar Matthies)
//
// =======================================================================

#ifndef __AUXPARAM3D__
#define __AUXPARAM3D__

class TFEFunction3D;
class TQuadFormula;
#include <Constants.h>
#include <FESpace3D.h>
#include <FEFunction3D.h>
#include "MooNMD_Io.h"

/** store parameter functions and FE functions */
class TAuxParam3D
{
  protected:
// =======================================================================
//  numbers of stored objects
// =======================================================================
    /** number of stored FESpace3D */
    int N_FESpace3D;

    /** number of stored FEFunction3D */
    int N_FEFunction3D;

    /** number of stored parameter function (ParamFct) */
    int N_ParamFct;

// =======================================================================
//  array of pointers to stored objects
// =======================================================================
    /** array of stored FESpace3D */
    const TFESpace3D **FESpaces3D;

    /** array of stored FEFunction3D */
    TFEFunction3D **FEFunctions3D;

    /** array of stored parameter function */
    ParamFct **ParameterFct;

// =======================================================================
//  information of FE values used by parameter functions
// =======================================================================
    /** number of FE values */
    int N_FEValues;

    /** index of FEFunction3D used for FE value i */
    int *FEValue_FctIndex;

    /** which multiindex is used for FE value i */
    MultiIndex3D *FEValue_MultiIndex;

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

    double **Values;
    double ***OrigValues;
    const int **Index;
    int *N_BaseFunct;

  public:
    /** constructor */
    TAuxParam3D(int n_fespace3d, int n_fefunction3d, int n_paramfct,
              int n_fevalues,
              const TFESpace3D **fespaces3d, TFEFunction3D **fefunctions3d,
              ParamFct **parameterfct,
              int *fevalue_fctindex, MultiIndex3D *fevalue_multiindex,
              int n_parameters, int *beginparameter);
    /** @brief standard constructor
     * 
     * If you don't need values of a finite element function in your assembling,
     * choose this constructor. This is equivalent to calling 
     * TAuxParam3D(0, 0, 0, 0, nullptr, nullptr, nullptr, nullptr, nullptr, 0, nullptr);
     */
    TAuxParam3D();

    TAuxParam3D(const std::string& name, TFEFunction3D **fefunctions3d);
    

    /** destructor */
    ~TAuxParam3D();

    /** return all parameters at all quadrature points */
    void GetParameters(const TQuadFormula& qf, int cellnum,
                       double **Parameters);

    //! Getter marked as deprecated, because it does not follow the new naming conventions.
    [[deprecated]] int GetN_Parameters()
    { return N_Parameters; }

    //Getter methods which follow new naming conventions.

	int* getBeginParameter() const {
		return BeginParameter;
	}

    int getBeginParameter(int i) const
    {
    	if (i >= N_ParamFct)
    	{
    		Output::warn("TAuxParam3D::getBeginParameter",
                     "Array out of bound: i >= N_ParamFct.");
    	}
    	return BeginParameter[i];
    }

	TFEFunction3D** getFeFunctions3D() const {
		return FEFunctions3D;
	}

	int* getFeValueFctIndex() const {
		return FEValue_FctIndex;
	}

    int getFeValueFctIndex(int i) const
    {
    	if (i >= N_FEValues)
    	{
    		Output::warn("TAuxParam3D::getFeValueFctIndex",
                     "Array out of bound: i >= N_FEValues.");
    	}
    	return FEValue_FctIndex[i];
    }

	MultiIndex3D* getFeValueMultiIndex() const {
		return FEValue_MultiIndex;
	}

    MultiIndex3D getFeValueMultiIndex(int i) const
    {
    	if (i >= N_FEValues)
    	{
    		Output::warn("TAuxParam3D::getFeValueMultiIndex",
                     "Array out of bound: i >= N_FEValues.");
    	}
    	return FEValue_MultiIndex[i];
    }


	const int** getIndex() const {
		return Index;
	}

	int* getNBaseFunct() const {
		return N_BaseFunct;
	}

	int getNFeFunction3D() const {
		return N_FEFunction3D;
	}

	int getNFeSpace3D() const {
		return N_FESpace3D;
	}

	int getNFeValues() const {
		return N_FEValues;
	}

	int getNParameters() const {
		return N_Parameters;
	}

	int getNParamFct() const {
		return N_ParamFct;
	}

	double*** getOrigValues() const {
		return OrigValues;
	}

	ParamFct** getParameterFct() const {
		return ParameterFct;
	}

    ParamFct* getParameterFct(int i) const
    {
    	if (i >= N_ParamFct)
    	{
    		Output::warn("TAuxParam3D::getParameterFct",
                     "Array out of bound: i >= N_ParamFct.");
    	}
    	return ParameterFct[i];
    }

	double* getTemp() const {
		return Temp;
	}

	double** getValues() const {
		return Values;
	}

};

// standard function to use for Navier-Stokes
void Velocity_Fct(const double *inputList, double *outputValues);

#endif // __AUXPARAM3D__
