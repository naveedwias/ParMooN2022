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

#include <AuxParam3D.h>
#include "FEDatabase.h"
#include <stdlib.h>
#include "BaseCell.h"

/** constructor */
TAuxParam3D::TAuxParam3D(
        int n_fespace3d, int n_fefunction3d, int n_paramfct,
        int n_fevalues,
        const TFESpace3D **fespaces3d, TFEFunction3D **fefunctions3d,
        ParamFct **parameterfct,
        int *fevalue_fctindex, MultiIndex3D *fevalue_multiindex,
        int n_parameters, int *beginparameter)
{
  N_FESpace3D = n_fespace3d;
  N_FEFunction3D = n_fefunction3d;
  N_ParamFct = n_paramfct;
  N_FEValues = n_fevalues;

  FESpaces3D = fespaces3d;
  FEFunctions3D = fefunctions3d;
  ParameterFct = parameterfct;

  FEValue_FctIndex = fevalue_fctindex;
  FEValue_MultiIndex = fevalue_multiindex;

  N_Parameters = n_parameters;
  BeginParameter = beginparameter;

  Temp = new double[3 + N_FEValues];

  Values = new double* [N_FEValues];
  OrigValues = new double** [N_FEValues];
  Index = new const int* [N_FEValues];
  N_BaseFunct = new int[N_FEValues];
}

TAuxParam3D::TAuxParam3D()
 : TAuxParam3D(0, 0, 0, 0, nullptr, nullptr, nullptr, nullptr, nullptr, 0,
               nullptr)
{
}


/** set aux parameters giving a keyword*/
TAuxParam3D::TAuxParam3D(const std::string& name, TFEFunction3D **fefunctions3d)
{
  if(name=="Velocity") // for Navier-Stokes
  {
    N_FESpace3D = 0;
    FESpaces3D = nullptr;
    
    N_FEFunction3D = 3;
    FEFunctions3D = fefunctions3d;
    
    N_ParamFct = 1;
    N_FEValues = 3;
    
    //see TNSE3D_FixPo.C
    ParameterFct=new ParamFct*[1];
    ParameterFct[0] = &Velocity_Fct; // see below in this file
    FEValue_FctIndex =new int[3];
    FEValue_FctIndex[0] = 0;
    FEValue_FctIndex[1] = 1;
    FEValue_FctIndex[2] = 2;
    
    FEValue_MultiIndex=new MultiIndex3D[3];
    FEValue_MultiIndex[0]= MultiIndex3D::D000;
    FEValue_MultiIndex[1]= MultiIndex3D::D000;
    FEValue_MultiIndex[2]= MultiIndex3D::D000;
    
    N_Parameters = 3;
    BeginParameter =new int[1];
    BeginParameter[0]=0;
    
  }
  else
  {
    ErrThrow(" AuxParam3D:: Constructor: ERROR, name ", name,
             " for initialization not imlpemented ");
  }
  
  Temp = new double[3 + N_FEValues];

  Values = new double* [N_FEValues];
  OrigValues = new double** [N_FEValues];
  Index = new const int* [N_FEValues];
  N_BaseFunct = new int[N_FEValues];
}



/** return all parameters at all quadrature points */
void TAuxParam3D::GetParameters(const TQuadFormula& qf, int cellnum,
                                double **Parameters)
{
  int i, j, k, l, n;
  double *param, *currparam, s;
  TFEFunction3D *fefunction;

  double *CurrValues, *CurrOrigValues;
  const int *CurrIndex;

   // collect information
  for(j=0;j<N_FEValues;j++)
  {
    fefunction = FEFunctions3D[FEValue_FctIndex[j]];
    Values[j] = fefunction->GetValues();

    auto fespace = fefunction->GetFESpace3D();
    auto fe = fespace->get_fe(cellnum);

    N_BaseFunct[j]=fe.GetN_DOF();
    OrigValues[j] = FEDatabase::GetOrigElementValues(*fe.GetBaseFunct(),
                                                     FEValue_MultiIndex[j]);

    Index[j] = fespace->GetGlobalDOF(cellnum);
  } // endfor j

  int N_Points = qf.GetN_QuadPoints();

  // loop over all quadrature points
  for(i=0;i<N_Points;i++)
  {
    // parameters in this quadrature point  
    param = Parameters[i];
    // first three parameters are the coordinates
    auto p = qf.get_point(i);
    Temp[0] = p.x;
    Temp[1] = p.y;
    Temp[2] = p.z;

    // loop to calculate all FE values
    for(k=3,j=0;j<N_FEValues;j++,k++)
    {
      s = 0;
      n = N_BaseFunct[j];
      CurrValues = Values[j];
      CurrOrigValues = OrigValues[j][i];
      CurrIndex = Index[j];
      for(l=0;l<n;l++)
        s += CurrValues[CurrIndex[l]]*CurrOrigValues[l];
      Temp[k] = s;
    }  // endfor j

    // loop to calculate all parameters
    for(j=0;j<N_ParamFct;j++)
    {
      currparam = param + BeginParameter[j];
      ParameterFct[j](Temp, currparam);
    } // endfor j
  } // endfor i
}

/** destructor */
TAuxParam3D::~TAuxParam3D()
{
  delete [] Temp;
  delete [] Values;
  delete [] OrigValues;
  delete [] Index;
  delete [] N_BaseFunct;
}


void Velocity_Fct(const double *inputList, double *outputValues)
{
  outputValues[0] = inputList[3];                // u1old
  outputValues[1] = inputList[4];                // u2old
  outputValues[2] = inputList[5];                // u3old
}
