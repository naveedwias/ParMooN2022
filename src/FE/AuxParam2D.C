// =======================================================================
// @(#)AuxParam2D.C        1.2 09/17/99
// 
// Class:       TAuxParam2D
// Purpose:     store parameter functions and FE functions
//
// Author:      Gunar Matthies (06.08.98)
//
// History:     start of implementation 06.08.98 (Gunar Matthies)
//
// =======================================================================

#include <AuxParam2D.h>
#include "FEDatabase.h"
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <stdio.h>

TAuxParam2D::TAuxParam2D(int n_paramfct, int n_fevalues,
                         TFEFunction2D **fefunctions2d, ParamFct **parameterfct,
                         int *fevalue_fctindex,
                         MultiIndex2D *fevalue_multiindex, int n_parameters,
                         int *beginparameter)
{
  N_ParamFct = n_paramfct;
  N_FEValues = n_fevalues;

  FEFunctions2D = fefunctions2d;
  ParameterFct = parameterfct;

  FEValue_FctIndex = fevalue_fctindex;
  FEValue_MultiIndex = fevalue_multiindex;

  N_Parameters = n_parameters;
  BeginParameter = beginparameter;

  Temp = new double[2 + N_FEValues];

  Values = new const double* [N_FEValues];
  OrigValues = new double** [N_FEValues];
  Index = new const int* [N_FEValues];
  N_BaseFunct = new int[N_FEValues];
}

TAuxParam2D::TAuxParam2D() 
 : TAuxParam2D(0, 0, nullptr, nullptr, nullptr, nullptr, 0, nullptr)
{
}



/** return all parameters at all quadrature points */
void TAuxParam2D::GetParameters(const TQuadFormula& qf, int cellnum,
                                double **Parameters)
{
  // collect information
  for(int j=0;j<N_FEValues;j++)
  {
    const TFEFunction2D *fefunction = FEFunctions2D[FEValue_FctIndex[j]];
    Values[j] = fefunction->GetValues();
    //  if (N_FEValues==8)
    //  Output::print("aac ", (int) fefunction, " ", Values[j][0]);

    auto fespace = fefunction->GetFESpace2D();
    auto  element = fespace->get_fe(cellnum);
    N_BaseFunct[j]=element.GetN_DOF();
    OrigValues[j] = FEDatabase::GetOrigElementValues(*element.GetBaseFunct(),
                                                     FEValue_MultiIndex[j]);
    Index[j] = fespace->GetGlobalDOF(cellnum);
  } // endfor j

  int N_Points = qf.GetN_QuadPoints();
  // loop over all quadrature points
  for(int i=0;i<N_Points;i++)
  {
    double * param = Parameters[i];

    auto p = qf.get_point(i);
    Temp[0] = p.x;
    Temp[1] = p.y;

    // loop to calculate all FE values
    for(int k=2,j=0;j<N_FEValues;j++,k++)
    {
      double s = 0;
      int n = N_BaseFunct[j];
      const double *CurrValues = Values[j];
      const double *CurrOrigValues = OrigValues[j][i];
      const int *CurrIndex = Index[j];
      for(int l=0;l<n;l++)
        s += CurrValues[CurrIndex[l]]*CurrOrigValues[l];
      Temp[k] = s;
    }  // endfor j

    // loop to calculate all parameters
    for(int j=0;j<N_ParamFct;j++)
    {
      double * currparam = param + BeginParameter[j];
      ParameterFct[j](Temp, currparam);
    } // endfor j
  } // endfor i
}

/** destructor */
TAuxParam2D::~TAuxParam2D()
{
  delete [] Temp;
  delete [] Values;
  delete [] OrigValues;
  delete [] Index;
  delete [] N_BaseFunct;
}
