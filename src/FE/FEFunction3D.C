#include <Database.h>
#include <Joint.h>
#include <BoundFace.h>
#include "FEDatabase.h"
#include <FEFunction3D.h>
#include <AllRefTrans3D.h>
#include <GridCell.h>
#include <ConvDiff.h>
#include "AuxParam3D.h"
#include "QuadratureFormulaDatabase.h"

#include <stdlib.h>
#include <InterfaceJoint.h>
#include <BdPlane.h>
#include <LinAlg.h>
#include <cstring> // memset

#ifdef _MPI
#include <ParFECommunicator3D.h>
#endif

#include <algorithm>
#include <cmath>


void OnlyDirichlet(double, double, double, BoundCond &cond)
{
	cond = DIRICHLET;
}

TFEFunction3D::TFEFunction3D() :
    Name("dummy_fe_fct_3d")
{
  FESpace3D=nullptr;
  Values=nullptr;
}

/** constructor with vector initialization */
TFEFunction3D::TFEFunction3D(std::shared_ptr<const TFESpace3D> fespace3D,
                             const std::string& name, double *values)
: Name(name)
{
  Output::print<5>("Constructor of TFEFunction3D");
  FESpace3D=fespace3D;
  Values=values;
}

TFEFunction3D& TFEFunction3D::operator=(const TFEFunction3D& other)
{
  this->Name        = other.Name;
  this->FESpace3D   = other.FESpace3D;
  this->Values      = other.Values;

  return *this;
}



/** calculate errors to given function */
void TFEFunction3D::GetErrors(DoubleFunct3D *Exact, int N_Derivatives,
                              MultiIndex3D *NeededDerivatives,
                              int N_Errors, ErrorMethod *ErrorMeth, 
                              CoeffFct3D Coeff, 
                              TAuxParam3D *Aux,
                              int n_fespaces, const TFESpace3D **fespaces,
                              double *errors,
                              std::function<bool(const TBaseCell*, int)>funct)
  const
{
  int i,j,k,l;
  int N_Cells, N_Points, N_Parameters, N_;
  const double *weights;
  const double *X, *Y, *Z;
  std::vector<double> AbsDetjk(MaxN_QuadPoints_3D, 1.);
  double *Param[MaxN_QuadPoints_3D], *aux;
  double *Derivatives[MaxN_QuadPoints_3D];
  double *ExactVal[MaxN_QuadPoints_3D];
  double *AuxArray[MaxN_QuadPoints_3D];
  double **OrigFEValues, *Orig, value;
  double FEFunctValues[MaxN_BaseFunctions3D];
  std::vector<double> LocError(N_Errors, 0.0);
  double hK;
  bool *SecondDer;

  
#ifdef _MPI
   int ID, rank;
   
   MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);    
#endif

  SecondDer = new bool[n_fespaces];
  for(i=0;i<n_fespaces;i++)
    SecondDer[i] = false;

  N_Parameters = Aux->GetN_Parameters();
  
  if(N_Parameters==0)
   aux = nullptr;
  else
   aux = new double [MaxN_QuadPoints_3D*N_Parameters];
  
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Param[j] = aux + j*N_Parameters;

  aux = new double [MaxN_QuadPoints_3D*N_Derivatives];
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Derivatives[j] = aux + j*N_Derivatives;
  
  aux = new double [MaxN_QuadPoints_3D * 5];
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    ExactVal[j] = aux + j*5;

  // 20 <= number of term
  aux = new double [MaxN_QuadPoints_3D*20]; 
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    AuxArray[j] = aux + j*20;

  for(i=0;i<N_Errors;i++)
    errors[i] = 0.0;

// ########################################################################
// loop over all cells
// ########################################################################
  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTetra); // dummy type
  TQuadFormula qf_orig(qf_ref);
  std::vector<const FiniteElement*> LocalUsedElements(n_fespaces, nullptr);
  auto Coll = fespaces[0]->GetCollection(); // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
  
  for(i=0;i<N_Cells;i++)
  {
    auto cell = Coll->GetCell(i);
    if(funct(cell, i))
      continue;

#ifdef _MPI
    ID = cell->GetSubDomainNo();
    if(rank!=ID) continue;
#endif
    
    hK = cell->GetDiameter();

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    for(j=0;j<n_fespaces;j++)
    {
      LocalUsedElements[j] = &fespaces[j]->get_fe(i);
    }
    // ####################################################################
    // calculate values on original element
    // ####################################################################
    FEDatabase::GetOrig(LocalUsedElements, Coll, cell, SecondDer, qf_ref,
                        qf_orig);
    qf_orig.GetFormulaData(N_Points, weights, X, Y, Z);

    if(N_Parameters>0)
      Aux->GetParameters(qf_orig, i, Param);
    

    // calculate all needed derivatives of this FE function
      auto fe = FESpace3D->get_fe(i);
    N_ = fe.GetN_DOF();

    auto DOF = FESpace3D->GetGlobalDOF(i);
    for(l=0;l<N_;l++)
      FEFunctValues[l] = Values[DOF[l]];

    for(k=0;k<N_Derivatives;k++)
    {
      OrigFEValues = FEDatabase::GetOrigElementValues(*fe.GetBaseFunct(),
                                                      NeededDerivatives[k]);
      for(j=0;j<N_Points;j++)
      {
        Orig = OrigFEValues[j];
        value = 0;
        for(l=0;l<N_;l++)
        {
          value += FEFunctValues[l] * Orig[l];
        } // endfor l
        Derivatives[j][k] = value;
      } // endfor j
    } // endfor k
 
    for(j=0;j<N_Points;j++)
    {
      auto p = qf_orig.get_point(j);
      Exact(p.x, p.y, p.z, ExactVal[j]);
    }

    if(Coeff)
      Coeff(N_Points, X, Y, Z, Param, AuxArray);

    ErrorMeth(N_Points, {{X, Y, Z}}, AbsDetjk.data(), weights, hK, Derivatives,
              ExactVal, AuxArray, LocError.data());


    // ########################################################################
    // DG-error
    // ########################################################################
    // The DG Error is not implemented completely. Up to now only the volume
    // part is implemented. See FEFunction2D::GetErrors for an idea of an
    // implementation
    if(FESpace3D->is_discontinuous() &&
        ErrorMeth == &conv_diff_l2_h1_linf_error<3>)
    {
      LocError[3] = std::nan("");
    }


    if(!(ErrorMeth == &conv_diff_l2_h1_linf_error<3>))
    {
      for(int j = 0; j < N_Errors; j++)
        errors[j] += LocError[j];
    }
    else
    {
      for(int j = 0; j < N_Errors-1; j++)
        errors[j] += LocError[j];
      // L_infty error
      if(errors[N_Errors-1] <  LocError[N_Errors-1])
        errors[N_Errors-1] = LocError[N_Errors-1];
    }
  } // endfor i

#ifndef _MPI // std::sqrt(errors[j]) in the main programm after collecting error from all subdomains
  if(!(ErrorMeth == &conv_diff_l2_h1_linf_error<3>))
  {
    for(j=0;j<N_Errors;j++)
      errors[j] = std::sqrt(errors[j]);
  }
  else
  {
    for(j=0;j<N_Errors-1;j++)
      errors[j] = std::sqrt(errors[j]);
  }
#endif
  
  delete [] AuxArray[0];
  delete [] SecondDer;
  delete [] ExactVal[0];
  delete [] Derivatives[0];
  
  if(Param[0])
   delete [] Param[0];
  
} // TFEFunction3D::GetErrors

void TFEFunction3D::GetErrorsForVectorValuedFunction(
    DoubleFunct3D * const * const Exact, ErrorMethod * const ErrMeth,
  double * const errors)
{
  // write zeros into the array "errors"
  std::fill(errors, errors+3, 0.0);
  auto Coll = FESpace3D->GetCollection(); 
  int N_Cells = Coll->GetN_Cells(); // number of cells
  
  // second derivatives are not needed
  bool Needs2ndDer = false;
  TQuadFormula qf_orig(QuadratureFormula_type::BaryCenterTetra); // dummy type
  // for loop over cells
  for(int i =0; i < N_Cells; i++)
  {
    auto cell = Coll->GetCell(i);
    auto fe = FESpace3D->get_fe(i);
    // number of basis functions
    int N_BaseFuncts = fe.GetN_DOF();
    auto bf = fe.GetBaseFunct();
    const int * DOF = FESpace3D->GetGlobalDOF(i);
    ReferenceTransformation_type RefTrans = fe.GetRefTransID();
    // quadrature formula
    const TQuadFormula *qf = QuadratureFormulaDatabase::qf_from_degree(
        TDatabase::ParamDB->INPUT_QUAD_RULE, bf->GetRefElement());
    
    TRefTrans3D *F_K = FEDatabase::GetRefTrans3D(RefTrans);
    F_K->SetCell(cell);
    // fill arrays, so call of GetOrigElementValues is now possible
    FEDatabase::GetOrigValues(RefTrans, {bf}, Coll, cell, *qf, &Needs2ndDer);
    F_K->GetOrigFromRef(*qf, qf_orig);
    
    int N_Points = qf_orig.GetN_QuadPoints();
    // fill ExactVal
    double ** ExactVal = new double*[N_Points];
    for(int j=0;j<N_Points;j++)
    {
      auto p = qf_orig.get_point(j);
      // 5 values for each of the three components:
      // (value, 3 derivatives, Laplace)
      ExactVal[j] = new double[15]; 
      Exact[0](p.x, p.y, p.z, ExactVal[j]     ); // x-component
      Exact[1](p.x, p.y, p.z, ExactVal[j] +  5); // y-component
      Exact[2](p.x, p.y, p.z, ExactVal[j] + 10); // z-component
    }
    // will store values of this FEFunction and its first derivatives at all 
    // quadrature points
    double ** Derivatives = new double*[N_Points];
    for(int j=0;j<N_Points;j++)
    {
      // the 12 means values and 3 first derivatives for all three components
      Derivatives[j] = new double[12];
    }
    // Get the function values of this FE-function at the local dofs.
    // some local dofs get a negative sign according to global orientation of 
    // the normals
    double * FEFunctValues = new double[N_BaseFuncts];
    for(int l=0;l<N_BaseFuncts;l++)
    {
      FEFunctValues[l] = Values[DOF[l]];
    }
    
    // these arrays were created in GetOrigValues called earlier
    double **AllOrigValues[4];
    AllOrigValues[0] = FEDatabase::GetOrigElementValues(*bf,
                                                        MultiIndex3D::D000);
    AllOrigValues[1] = FEDatabase::GetOrigElementValues(*bf,
                                                        MultiIndex3D::D100);
    AllOrigValues[2] = FEDatabase::GetOrigElementValues(*bf,
                                                        MultiIndex3D::D010);
    AllOrigValues[3] = FEDatabase::GetOrigElementValues(*bf,
                                                        MultiIndex3D::D001);
    
    // loop over all needed derivatives
    for(int k = 0; k < 4; k++)
    {
      // loop over all quadrature points
      for(int j = 0; j < N_Points; j++) 
      {
        double value_x = 0;
        double value_y = 0;
        double value_z = 0;
        // loop over all basis functions
        for(int l=0;l<N_BaseFuncts;l++)
        {
          value_x += FEFunctValues[l] * AllOrigValues[k][j][l                ];
          value_y += FEFunctValues[l] * AllOrigValues[k][j][l+   N_BaseFuncts];
          value_z += FEFunctValues[l] * AllOrigValues[k][j][l+ 2*N_BaseFuncts];
        }
        Derivatives[j][k    ] = value_x;
        Derivatives[j][k + 4] = value_y;
        Derivatives[j][k + 8] = value_z;
      }
    }
    // cell diameter, we set it one here since it is not needed in 
    // ErrorMeth=L2DivH1Errors
    double hK = 1;
    double LocError[3]; // L^2 error in value, divergence and first derivative
    std::vector<double> AbsDetjk(N_Points, 1.);
    std::vector<double> weights(N_Points);
    std::vector<double> XYZ(3*N_Points);
    for(int j = 0; j < N_Points; j++) 
    {
      auto p = qf_orig.get_point(j);
      XYZ[j] = p.x;
      XYZ[j+N_Points] = p.y;
      XYZ[j+2*N_Points] = p.z;
      weights[j] = qf_orig.get_weight(j);
    }
    ErrMeth(N_Points,
            {{XYZ.data(), XYZ.data() + N_Points, XYZ.data() + 2*N_Points}},
            AbsDetjk.data(), weights.data(), hK, Derivatives, ExactVal, nullptr,
            LocError);
    for(int j=0;j<3;j++) 
    {
      errors[j] += LocError[j];
    }
    // delete everything which was created with "new" within this loop
    // otherwise one would get (many) memory leaks
    for (int j=0; j<N_Points; j++)
    {
      delete [] ExactVal[j];    ExactVal[j] = nullptr;
      delete [] Derivatives[j]; Derivatives[j] = nullptr;
    }
    delete [] ExactVal;      ExactVal = nullptr;
    delete [] Derivatives;   Derivatives = nullptr;
    delete [] FEFunctValues; FEFunctValues = nullptr;
    
  } // end loop over all cells
  
  for(int j=0;j<3;j++)  
    errors[j] = std::sqrt(errors[j]);
}


void TFEFunction3D::GetMeshCellParams(DoubleFunct3D *Exact, int N_Derivatives,
                              MultiIndex3D *NeededDerivatives,
                              int N_Errors, ErrorMethod *ErrorMeth, 
                              CoeffFct3D Coeff, 
                              TAuxParam3D *Aux,
                              int n_fespaces, const TFESpace3D **fespaces,
                              double *errors, double *cell_parameters)
{
  int i,j,k,l;
  int N_Cells, N_Points, N_Parameters, N_;
  const double *weights;
  const double *X, *Y, *Z;
  std::vector<double> AbsDetjk(MaxN_QuadPoints_3D, 1.);
  double *Param[MaxN_QuadPoints_3D], *aux;
  double *Derivatives[MaxN_QuadPoints_3D];
  double *ExactVal[MaxN_QuadPoints_3D];
  double *AuxArray[MaxN_QuadPoints_3D];
  double **OrigFEValues, *Orig, value;
  double FEFunctValues[MaxN_BaseFunctions3D];
  double LocError[4];
  double hK;
  bool *SecondDer;

  SecondDer = new bool[n_fespaces];
  for(i=0;i<n_fespaces;i++)
    SecondDer[i] = false;

  N_Parameters = Aux->GetN_Parameters();
  aux = new double [MaxN_QuadPoints_3D*N_Parameters];
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Param[j] = aux + j*N_Parameters;

  aux = new double [MaxN_QuadPoints_3D*N_Derivatives];
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Derivatives[j] = aux + j*N_Derivatives;
  
  aux = new double [MaxN_QuadPoints_3D * 5];
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    ExactVal[j] = aux + j*5;

  // 20 <= number of term
  aux = new double [MaxN_QuadPoints_3D*20]; 
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    AuxArray[j] = aux + j*20;

  for(i=0;i<N_Errors;i++)
    errors[i] = 0.0;

// ########################################################################
// loop over all cells
// ########################################################################
  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTetra); // dummy type
  TQuadFormula qf_orig(qf_ref);
  std::vector<const FiniteElement*> LocalUsedElements(n_fespaces, nullptr);
  auto Coll = fespaces[0]->GetCollection(); // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    auto cell = Coll->GetCell(i);

    hK = cell->GetDiameter();

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    for(j=0;j<n_fespaces;j++)
    {
      LocalUsedElements[j] = &fespaces[j]->get_fe(i);
    }

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    FEDatabase::GetOrig(LocalUsedElements,  Coll, cell, SecondDer, qf_ref,
                        qf_orig);
    qf_orig.GetFormulaData(N_Points, weights, X, Y, Z);
    
    if(N_Parameters>0)
      Aux->GetParameters(qf_orig, i, Param); 

    // calculate all needed derivatives of this FE function
    auto fe = FESpace3D->get_fe(i);
    N_ = fe.GetN_DOF();

    auto DOF = FESpace3D->GetGlobalDOF(i);
    for(l=0;l<N_;l++)
      FEFunctValues[l] = Values[DOF[l]];

    for(k=0;k<N_Derivatives;k++)
    {
      OrigFEValues = FEDatabase::GetOrigElementValues(*fe.GetBaseFunct(),
                                                      NeededDerivatives[k]);
      for(j=0;j<N_Points;j++)
      {
        Orig = OrigFEValues[j];
        value = 0;
        for(l=0;l<N_;l++)
        {
          value += FEFunctValues[l] * Orig[l];
        } // endfor l
        Derivatives[j][k] = value;
      } // endfor j
    } // endfor k

    for(j=0;j<N_Points;j++)
      Exact(X[j], Y[j], Z[j], ExactVal[j]);

    if(Coeff)
      Coeff(N_Points, X, Y, Z, Param, AuxArray);

    ErrorMeth(N_Points, {{X, Y, Z}}, AbsDetjk.data(), weights, hK, Derivatives,
              ExactVal, AuxArray, LocError);

    for(j=0;j<N_Errors;j++)
    {
      errors[j] += LocError[j];
      cell_parameters[i + j *N_Cells] = LocError[j];
    }
  } // endfor i

  for(j=0;j<N_Errors;j++)
    errors[j] = std::sqrt(errors[j]);

  delete AuxArray[0];
  delete SecondDer;
  delete ExactVal[0];
  delete Derivatives[0];
  delete Param[0];
  
} // TFEFunction3D::GetErrors



void TFEFunction3D::getValuesFromOriginalOnes(double* retValues, int valuesOffset,
                                              int N_BaseFunct, int cellNumber,
                                              double* uorig, double* uxorig,
                                              double* uyorig, double* uzorig) const
{
  std::memset(retValues, 0, 4*sizeof(double));
  
  const int* Numbers = FESpace3D->GetGlobalDOF(cellNumber);
  
  for(int j = 0; j < N_BaseFunct; j++)
  {
    double val = Values[valuesOffset + Numbers[j]];
    retValues[0] += uorig[j]*val;
    retValues[1] += uxorig[j]*val;
    retValues[2] += uyorig[j]*val;
    retValues[3] += uzorig[j]*val;
  }
}


void TFEFunction3D::getOrigValuesForCell(int cellNumber, double x, double y, double z, double* uorig, double* uxorig, double* uyorig, double* uzorig) const
{
  const FiniteElement FE_Obj = FESpace3D->get_fe(cellNumber);
  ReferenceTransformation_type refTrans = FE_Obj.GetRefTransID();
  TBaseCell* cellPtr = FESpace3D->GetCollection()->GetCell(cellNumber);
  
  bool isIsoparametric = TDatabase::ParamDB->USE_ISOPARAMETRIC
                         && cellPtr->has_isoparametric_joint();
  if(isIsoparametric)
  {
    BFRefElements refElement = FE_Obj.GetBaseFunct()->GetRefElement();
    
    switch(refElement)
    {
      case BFRefElements::BFUnitHexahedron:
        refTrans = ReferenceTransformation_type::HexaIsoparametric;
      break;
  
      case BFRefElements::BFUnitTetrahedron:
        refTrans = ReferenceTransformation_type::TetraIsoparametric;
      break;
      
      default:
        ;
    }
  }
  
  FEDatabase::SetCellForRefTrans(cellPtr, refTrans);
  
  double xi, eta, zeta;
  FEDatabase::GetRefFromOrig(refTrans, x, y, z, xi, eta, zeta);
  
  const BaseFunctions* bf = FE_Obj.GetBaseFunct();
  int N_BaseFunct = bf->GetDimension();
  
  double* uref = new double[N_BaseFunct];
  double* uxiref = new double[N_BaseFunct];
  double* uetaref = new double[N_BaseFunct];
  double* uzetaref = new double[N_BaseFunct];

  bf->GetDerivatives(MultiIndex3D::D000, xi, eta, zeta, uref);
  bf->GetDerivatives(MultiIndex3D::D100, xi, eta, zeta, uxiref);
  bf->GetDerivatives(MultiIndex3D::D010, xi, eta, zeta, uetaref);
  bf->GetDerivatives(MultiIndex3D::D001, xi, eta, zeta, uzetaref);
  
  FEDatabase::GetOrigValues(refTrans, xi, eta, zeta, bf,
                            FESpace3D->GetCollection(), (TGridCell *)cellPtr,
                            uref, uxiref, uetaref, uzetaref,
                            uorig, uxorig, uyorig, uzorig);
  
  delete[] uref;
  delete[] uxiref;
  delete[] uetaref;
  delete[] uzetaref;
}


bool TFEFunction3D::FindGradient(double x, double y, double z,
                                 std::vector<double>& values) const
{

  if(values.size() != 4)
  {
    ErrThrow("TFEFunction3D::FindGradient expects vector of size 4 !=", values.size());
  }
  
  int N_Found = 0;
  
  double valuesFromCells[4];
  
  std::memset(values.data(), 0, sizeof(double)*values.size());

  const TCollection* Coll = FESpace3D->GetCollection();
  int N_Cells = Coll->GetN_Cells();

  for(int i=0;i<N_Cells;i++)
  {
    TBaseCell* cell = Coll->GetCell(i);
    
    if(cell->PointInCell(parmoon::Point(x,y,z)))
    {
      N_Found++;    
      FindGradientLocal(i, x, y, z, valuesFromCells);
      
      for(int j = 0; j < 4; j++)
        values[j] += valuesFromCells[j];
    }
  }
  
  if(N_Found>0)
  {
    for(int j = 0; j < 4; j++)
      values[j] /= N_Found;
    return true;
  }
  return false;
} 

void TFEFunction3D::FindGradientLocal(int cell_no, 
                                      double x, double y, double z, 
                                      double *values) const
{
  const FiniteElement& fe = FESpace3D->get_fe(cell_no);
  const BaseFunctions* bf = fe.GetBaseFunct();
  int N_BaseFunct = bf->GetDimension();
  
  double* uorig = new double[N_BaseFunct];
  double* uxorig = new double[N_BaseFunct];
  double* uyorig = new double[N_BaseFunct];
  double* uzorig = new double[N_BaseFunct];
  
  getOrigValuesForCell(cell_no, x, y, z, uorig, uxorig, uyorig, uzorig);
  
  getValuesFromOriginalOnes(values, 0, N_BaseFunct, cell_no, uorig, uxorig, uyorig, uzorig);
  
  delete [] uorig;
  delete [] uxorig;
  delete [] uyorig;
  delete [] uzorig;
}

void TFEFunction3D::FindValueLocal(const TBaseCell *cell, int cell_no, 
                                   double x, double y, double z, 
                                   double *values) const
{
  double xi, eta, zeta;
  ReferenceTransformation_type RefTrans;
  int N_BaseFunct;
  double *uorig, *uxorig, *uyorig, *uzorig, *uref, *uxiref, *uetaref, *uzetaref;
  
  auto Coll = FESpace3D->GetCollection();
  auto fe = FESpace3D->get_fe(cell_no);
  RefTrans = fe.GetRefTransID();

  // set cell for reference transformation
  FEDatabase::SetCellForRefTrans(cell, RefTrans);
  
  // find local coordinates of the given point
  FEDatabase::GetRefFromOrig(RefTrans, x, y, z, xi, eta, zeta);
  // cout << " xi: " << xi << endl;
  // cout << "eta: " << eta << endl;
  // cout << "zeta: " << zeta << endl;
  
  // get base function object
  auto bf = fe.GetBaseFunct();
  N_BaseFunct = bf->GetDimension();
  int BaseVectDim = bf->GetBaseVectDim(); // either 1 or 3
  
  uorig = new double[N_BaseFunct*BaseVectDim];
  uxorig = new double[N_BaseFunct*BaseVectDim];
  uyorig = new double[N_BaseFunct*BaseVectDim];
  uzorig = new double[N_BaseFunct*BaseVectDim];
  
  uref = new double[N_BaseFunct*BaseVectDim];
  uxiref = new double[N_BaseFunct*BaseVectDim];
  uetaref = new double[N_BaseFunct*BaseVectDim];
  uzetaref = new double[N_BaseFunct*BaseVectDim];
  
  bf->GetDerivatives(MultiIndex3D::D000, xi, eta, zeta, uref);
  bf->GetDerivatives(MultiIndex3D::D100, xi, eta, zeta, uxiref);
  bf->GetDerivatives(MultiIndex3D::D010, xi, eta, zeta, uetaref);
  bf->GetDerivatives(MultiIndex3D::D001, xi, eta, zeta, uzetaref);
  
  FEDatabase::GetOrigValues(RefTrans, xi, eta, zeta, bf, Coll, cell,
                            uref, uxiref, uetaref, uzetaref,
                            uorig, uxorig, uyorig, uzorig);
  
  const int *Numbers = FESpace3D->GetGlobalDOF(cell_no);
  for (int i = 0; i < BaseVectDim; i++)
  {
    double u = 0;
    for(int j = 0; j < N_BaseFunct; j++)
    {
      double val = Values[Numbers[j]];
      u  +=  uorig[j + i*N_BaseFunct]*val;
    }
    values[i] = u;
  }

  
  delete [] uorig;
  delete [] uxorig;
  delete [] uyorig;
  delete [] uzorig;
  delete [] uref;
  delete [] uxiref;
  delete [] uetaref;
  delete [] uzetaref;
  
}

/** calculate the interpolation of an exact function */
void TFEFunction3D::Interpolate(DoubleFunct3D *Exact)
{
  int i,j;
  int N_Cells;
  int N_DOFs, N_LocalDOFs;
  int N_Points;
  const double *xi, *eta, *zeta;
  ReferenceTransformation_type F_K = ReferenceTransformation_type::TetraAffin; //avoid uninit warning
  TRefTrans3D *rt;
  double X[MaxN_PointsForNodal3D], Y[MaxN_PointsForNodal3D];
  double Z[MaxN_PointsForNodal3D];
  double PointValues[MaxN_PointsForNodal3D];
  double FunctionalValues[MaxN_PointsForNodal3D];
  double FctVal[5];
  int PolynomialDegree, ApproxOrder;
  bool IsIsoparametric;
  BFRefElements RefElement;
  ReferenceTransformation_type RefTrans;

  // begin code
  
  auto Coll = FESpace3D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  N_DOFs = FESpace3D->get_n_dof();

  memset(Values, 0, sizeof(double)*N_DOFs);

  for(i=0;i<N_Cells;i++)
  {
    auto cell = Coll->GetCell(i);
    auto Element = FESpace3D->get_fe(i);
    auto nf = Element.GetNodalFunctional();
    nf->GetPointsForAll(N_Points, xi, eta, zeta);
    N_LocalDOFs = Element.GetN_DOF();

    PolynomialDegree = Element.GetBaseFunct()->GetPolynomialDegree();
    ApproxOrder = Element.GetBaseFunct()->GetAccuracy();

    RefElement = Element.GetBaseFunct()->GetRefElement();
    auto QuadFormula = QuadratureFormulaDatabase::qf_from_degree(
        3*PolynomialDegree, RefElement);
    RefTrans = Element.GetRefTransID();

    IsIsoparametric = TDatabase::ParamDB->USE_ISOPARAMETRIC
                      && cell->has_isoparametric_joint();
    if(IsIsoparametric)
    {
      switch(RefElement)
      {
        case BFRefElements::BFUnitHexahedron:
          RefTrans = ReferenceTransformation_type::HexaIsoparametric;
        break;
  
        case BFRefElements::BFUnitTetrahedron:
          RefTrans = ReferenceTransformation_type::TetraIsoparametric;
        break;
        
        default:
          ;
      }
    } // endif IsIsoparametric
    // cout << "IsIsoparametric: " << IsIsoparametric << endl;
  
    switch(RefTrans)
    {
      case ReferenceTransformation_type::HexaAffin:
        rt = FEDatabase::GetRefTrans3D(ReferenceTransformation_type::HexaAffin);
        ((THexaAffin *)rt)->SetCell(cell);
        F_K = ReferenceTransformation_type::HexaAffin;
        break;
      case ReferenceTransformation_type::HexaTrilinear:
        rt = FEDatabase::GetRefTrans3D(
          ReferenceTransformation_type::HexaTrilinear);
        ((THexaTrilinear *)rt)->SetCell(cell);
        F_K = ReferenceTransformation_type::HexaTrilinear;
        break;
      case ReferenceTransformation_type::HexaIsoparametric:
        rt = FEDatabase::GetRefTrans3D(
          ReferenceTransformation_type::HexaIsoparametric);
        ((THexaIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((THexaIsoparametric *)rt)->SetQuadFormula(QuadFormula->get_type());
        ((THexaIsoparametric *)rt)->SetCell(cell);
        F_K = ReferenceTransformation_type::HexaIsoparametric;
        break;
      case ReferenceTransformation_type::TetraAffin:
        rt = FEDatabase::GetRefTrans3D(
          ReferenceTransformation_type::TetraAffin);
        ((TTetraAffin *)rt)->SetCell(cell);
        F_K = ReferenceTransformation_type::TetraAffin;
        break;
      case ReferenceTransformation_type::TetraIsoparametric:
        rt = FEDatabase::GetRefTrans3D(
          ReferenceTransformation_type::TetraIsoparametric);
        ((TTetraIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((TTetraIsoparametric *)rt)->SetQuadFormula(QuadFormula->get_type());
        ((TTetraIsoparametric *)rt)->SetCell(cell);
        F_K = ReferenceTransformation_type::TetraIsoparametric;
        break;
        
      default:
        ;
    }
    FEDatabase::GetOrigFromRef(F_K, N_Points, xi, eta, zeta, X, Y, Z);

    // cout << "----------------" << endl;
    for(j=0;j<N_Points;j++)
    {
      // cout << j << " ";
      // cout << "ref: " << xi[j] << " " << eta[j] << " " << zeta[j] << endl;
      // cout << "Ori: " << X[j] << " " << Y[j] << " " << Z[j] << endl;
      Exact(X[j], Y[j], Z[j], FctVal);
      // cout << FctVal[0] << endl;
      PointValues[j] = FctVal[0];
    }

    nf->GetAllFunctionals(Coll, (TGridCell *)cell, PointValues,FunctionalValues);

    auto DOF = FESpace3D->GetGlobalDOF(i);

    for(j=0;j<N_LocalDOFs;j++)
    {
      Values[DOF[j]] = FunctionalValues[j];
    }
  }
}

void TFEFunction3D::Interpolate_vector_valued_function(
    const std::vector<DoubleFunct3D*>& Exact)
{
  if(Exact.size() != 3)
  {
    ErrThrow("TFEFunction3D::Interpolate_vector_valued_function You have to "
             "provide three functions describing the exact solution ",
             Exact.size());
  }

  // reset function values
  int N_DOFs = this->FESpace3D->get_n_dof();
  memset(this->Values, 0, sizeof(double)*N_DOFs);


  // loop over all cells
  auto Coll = this->FESpace3D->GetCollection();
  int N_Cells = Coll->GetN_Cells();
  for(int i = 0; i < N_Cells; i++)
  {
    auto cell = Coll->GetCell(i);
    auto Element = FESpace3D->get_fe(i);
    auto nf = Element.GetNodalFunctional();
    int n_functionals = nf->n_functionals();
    // number of points needed to evaluate nodal functionals
    int N_Points = nf->n_points();
    // coordinates of points on reference element to evaluate nodal
    // functionals
    const double *xi, *eta, *zeta;
    nf->GetPointsForAll(N_Points, xi, eta, zeta);
    int N_LocalDOFs = Element.GetN_DOF();

    ReferenceTransformation_type RefTrans = Element.GetRefTransID();

    if(TDatabase::ParamDB->USE_ISOPARAMETRIC && cell->has_isoparametric_joint())
    {
      ErrThrow("Interpolate_vector_valued_function with isoparametric elements"
               " not yet implemented");
    }

    FEDatabase::SetCellForRefTrans(cell, RefTrans);

    // coordinates of points needed to evaluate the nodal functionals on
    // cell in grid.
    std::vector<double> X(N_Points), Y(N_Points), Z(N_Points);
    FEDatabase::GetOrigFromRef(RefTrans, N_Points, xi, eta, zeta,
                               X.data(), Y.data(), Z.data());

    std::vector<double> PointValues(Exact.size()*N_Points); // Exact.size() == 3
    for(int j = 0; j < N_Points; j++)
    {
/*       cout << j << " ";
       cout << "ref: " << xi[j] << " " << eta[j] << " " << zeta[j] << endl;
       cout << "Ori: " << X[j] << " " << Y[j] << " " << Z[j] << endl;*/
      for(unsigned int dim = 0; dim < Exact.size(); dim++) // Exact.size() == 3
      {
        double FctVal[5];
        Exact[dim](X[j], Y[j], Z[j], FctVal);

        // cout << FctVal[0] << endl;
        PointValues[j + dim*N_Points] = FctVal[0];
      }
    }

    std::vector<double> FunctionalValues(n_functionals);
    nf->GetAllFunctionals(Coll, (TGridCell *)cell, PointValues.data(),
                          FunctionalValues.data());

    const int *DOF = this->GetFESpace3D()->GetGlobalDOF(i);

    for(int j = 0; j < N_LocalDOFs; j++)
    {
      this->Values[DOF[j]] = FunctionalValues[j];
    }
  }
}


void TFEFunction3D::Interpolate(TFEFunction3D *funct)
{
  int i,j, N_Cells;
  int N_LocalDOFs;
  int N_Points;

  const double *xi, *eta, *zeta;
  double X[MaxN_PointsForNodal3D], Y[MaxN_PointsForNodal3D], Z[MaxN_PointsForNodal3D];
  double PointValues[MaxN_PointsForNodal3D];
  double FunctionalValues[MaxN_PointsForNodal3D];
  std::vector<double> values(4);


  auto Coll = FESpace3D->GetCollection();
  N_Cells = Coll->GetN_Cells();


  for(i=0;i<N_Cells;i++)
  {
    auto cell = Coll->GetCell(i);
    
    auto Element = FESpace3D->get_fe(i);
    auto nf = Element.GetNodalFunctional();
    nf->GetPointsForAll(N_Points, xi, eta, zeta);
    N_LocalDOFs = Element.GetN_DOF();

    
    ReferenceTransformation_type RefTrans = Element.GetRefTransID();

    
    bool IsIsoparametric = TDatabase::ParamDB->USE_ISOPARAMETRIC
                           && cell->has_isoparametric_joint();
    if(IsIsoparametric)
    {
      BFRefElements RefElement = Element.GetBaseFunct()->GetRefElement();
      
      switch(RefElement)
      {
        case BFRefElements::BFUnitHexahedron:
          RefTrans = ReferenceTransformation_type::HexaIsoparametric;
        break;
  
        case BFRefElements::BFUnitTetrahedron:
          RefTrans = ReferenceTransformation_type::TetraIsoparametric;
        break;
        
        default:
          ErrThrow("unexpected reference element ", RefElement);
          break;
      }
    }
    
    
    FEDatabase::SetCellForRefTrans(cell, RefTrans);
    FEDatabase::GetOrigFromRef(RefTrans, N_Points, xi, eta, zeta, X, Y, Z);
    
    for(j=0;j<N_Points;j++)
    {
      funct->FindGradient(X[j], Y[j], Z[j], values);
      PointValues[j] = values[0]; 
    }
    
    nf->GetAllFunctionals(Coll, (TGridCell *)cell, PointValues, FunctionalValues);

    auto DOF = FESpace3D->GetGlobalDOF(i);

    for(j=0;j<N_LocalDOFs;j++)
    {
      Values[DOF[j]] = FunctionalValues[j];
    }
  } //for i
}

void TFEFunction3D::add(AnalyticFunction f)
{
  int N_Points;
  const double *xi, *eta, *zeta;
  // begin code
 
  auto Coll = FESpace3D->GetCollection();
  int N_Cells = Coll->GetN_Cells();
  int N_DOFs = FESpace3D->get_n_dof();
  std::vector<int> IntIndex(N_DOFs, 0);
 
  for(int i = 0; i < N_Cells; i++)
  {
    TBaseCell * cell = Coll->GetCell(i);
    const FiniteElement& Element = FESpace3D->get_fe(i);
    auto nf = Element.GetNodalFunctional();
    nf->GetPointsForAll(N_Points, xi, eta, zeta);
    int N_LocalDOFs = Element.GetN_DOF();
 
    if(Element.GetBaseFunct()->GetBaseVectDim() != 1)
      ErrThrow("TFEFunction3D::add not implemented for vector-valued basis "
               "functions");
 
    BFRefElements RefElement = Element.GetBaseFunct()->GetRefElement();
    
    // else BFRefElements::BFUnitTetrahedron
    ReferenceTransformation_type RefTrans = Element.GetRefTransID();
    bool IsIsoparametric = TDatabase::ParamDB->USE_ISOPARAMETRIC
                           && cell->has_isoparametric_joint();
    if(IsIsoparametric)
    {
      RefTrans = (RefElement==BFRefElements::BFUnitHexahedron)
                 ? ReferenceTransformation_type::HexaIsoparametric
                 : ReferenceTransformation_type::TetraIsoparametric;
    }
    FEDatabase::SetCellForRefTrans(cell, RefTrans);
    std::vector<double> X(N_Points, 0.);
    std::vector<double> Y(N_Points, 0.);
    std::vector<double> Z(N_Points, 0.);
    FEDatabase::GetOrigFromRef(RefTrans, N_Points, xi, eta, zeta,
                               X.data(), Y.data(), Z.data());
    std::vector<double> PointValues(N_Points, 0.);
 
    for(int j = 0; j < N_Points; j++)
    {
      PointValues[j] = f(cell, i, {{X[j], Y[j], Z[j]}});
    }
    std::vector<double> FunctionalValues(N_LocalDOFs, 0.);
    nf->GetAllFunctionals(Coll, (TGridCell *)cell, PointValues.data(),
        FunctionalValues.data());
 
    const int * DOF = FESpace3D->GetGlobalDOF(i);
 
    for(int j = 0; j < N_LocalDOFs; j++)
    {
      if(IntIndex[DOF[j]] == 0)
        Values[DOF[j]] += FunctionalValues[j];
      IntIndex[DOF[j]] ++;
    }
  }
  for(int i = 0; i < N_DOFs; i++)
  {
    if(IntIndex[i] == 0)
    {
      ErrThrow("TFEFunction3D::add: unable to set dof ", i);
    }
  }
}

void TFEFunction3D::add_constant(double b)
{
  int length = FESpace3D->get_n_dof();
  // vector of the same length as this TFEFunction3D. It represents a function
  // which has the constant value 'mean' for all nodal functionals. The last
  // step in this projection will be to add this vector from the vector of
  // this TFEFunction3D
  // for standard P_k or Q_k finite elements this is a constant function
  std::vector<double> interpol(length);

  auto coll = FESpace3D->GetCollection();
  const int n_cells = coll->GetN_Cells();
  for(int i = 0; i < n_cells; i++)
  {
    auto cell = coll->GetCell(i); // current cell
    // finite element on the current cell
    auto fe = FESpace3D->get_fe(i);
    const int n_loc_dof = fe.GetN_DOF(); // number of local dofs
    const int * DOF = FESpace3D->GetGlobalDOF(i);
    auto nf = fe.GetNodalFunctional();
    int n_points; // number of evaluation points to compute nodal functionals
    const double *xi, *eta, *zeta; // coordinates of evaluation points in reference cell
    nf->GetPointsForAll(n_points, xi, eta, zeta);
    double *point_values = new double[n_points];
    for(int j = 0; j < n_points; j++)
      point_values[j] = b;
    // evaluate nodal functionals
    double *functional_values = new double[n_loc_dof];
    nf->GetAllFunctionals(coll, cell, point_values, functional_values);
    for(int j = 0; j < n_loc_dof; j++)
      interpol[DOF[j]] = functional_values[j];

    delete [] point_values;
    delete [] functional_values;
  }

  // the vector 'interpol' is now complete
  // substract it from the vector of this TFEFunction3D
  for(int i = 0; i < length; i++)
    Values[i] += interpol[i];
}

/** compute integral and measure */
void TFEFunction3D::compute_integral_and_measure(double& integral,
                                                 double& measure) const
{
  auto coll = FESpace3D->GetCollection();

  integral = 0.0; // variable to store integral value of this TFEFunction3D
  measure = 0.0; // variable to store the measure of the domain

  // loop over all cells, find out integral value of this FEFunction3D and the`
  // measure of its domain
  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTetra); // dummy type
  TQuadFormula qf_orig(qf_ref);
  const int n_cells = coll->GetN_Cells();
  for(int i = 0; i < n_cells; i++)
  {
    TBaseCell *cell = coll->GetCell(i); // current cell
#ifdef _MPI // skip halo cells
    if (cell->IsHaloCell())
    {
      continue;
    }
#endif
    // calculate values on original element (i.e. prepare reference
    // transformation)
    bool SecondDer = false;
    int n_points = 0;
    auto fe = FESpace3D->get_fe(i);
    FEDatabase::GetOrig({&fe}, coll, cell, &SecondDer, qf_ref, qf_orig);
    n_points = qf_orig.GetN_QuadPoints();

    // finite element on the current cell
    const int n_loc_dof = fe.GetN_DOF(); // number of local dofs
    const int * DOF = FESpace3D->GetGlobalDOF(i);

    // id of the local basis functions
    auto bf = fe.GetBaseFunct();
    // transformed values of basis functions
    double **orig_values = FEDatabase::GetOrigElementValues(*bf,
                                                            MultiIndex3D::D000);
    // local integration (loop over all quadrature points)
    for(int j = 0; j < n_points; j++)
    {
      // local transformed values on this quadrature point
      double * orig = orig_values[j];
      double value = 0; // value of this TFEFunction3D at this quadrature point
      for(int l = 0; l < n_loc_dof; l++)
      {
        // entry in the vector of this TFEFunction3D times basis function
        value += Values[DOF[l]] * orig[l];
      } // endfor l

      const double w = qf_orig.get_weight(j);
      integral += w * value;
      measure += w;
    } // endfor j
  }

#ifdef _MPI //communicate the results and add up over all processes
  double sendbuf[2] = {0.0, 0.0};
  double recvbuf[2] = {0.0, 0.0};
  sendbuf[0] = integral; //partial values
  sendbuf[1] = measure;
  MPI_Allreduce(sendbuf, recvbuf, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  integral = recvbuf[0]; //fill in the summed up values
  measure = recvbuf[1];
#endif
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double TFEFunction3D::get_L2_norm() const
{
  auto collection = FESpace3D->GetCollection();
  auto n_cells = collection->GetN_Cells();
  double l2_norm = 0.;
  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTetra); // dummy type
  TQuadFormula qf_orig(qf_ref);
  for(int i_cell = 0; i_cell < n_cells; ++i_cell)
  {
    auto cell = collection->GetCell(i_cell);
#ifdef _MPI // skip halo cells
    if (cell->IsHaloCell())
    {
      continue;
    }
#endif
    auto fe = FESpace3D->get_fe(i_cell);
    auto basis_functions = fe.GetBaseFunct();
    auto n_local_basis_functions = basis_functions->GetDimension();
    bool SecondDer = false;

    // Compute transformation of basis functions (and their derivatives) to 
    // original cell and the evaluation for all quadrature points 
    std::vector<const FiniteElement*> used_fe(1, &fe);
    FEDatabase::GetOrig(used_fe, collection, cell, &SecondDer, qf_ref, qf_orig);
    
    int N_QuadPoints = qf_orig.GetN_QuadPoints();

    // calculate all needed derivatives of this FE function
    std::vector<double> FEFunctValues(n_local_basis_functions);
    const int *DOF = FESpace3D->GetGlobalDOF(i_cell);
    for(int l = 0; l < n_local_basis_functions; l++)
    {
      FEFunctValues[l] = Values[DOF[l]];
    }

    // create a pointer to the QuadPoint-values of the basis on original element and its derivatives
    double **OrigFEValues = FEDatabase::GetOrigElementValues(
        *basis_functions, MultiIndex3D::D000);
    for(int j = 0; j < N_QuadPoints; j++)
    {
      double *Orig = OrigFEValues[j];
      double value = 0.;

      for(int l = 0; l < n_local_basis_functions; l++)
      {
        value += FEFunctValues[l] * Orig[l];
      }
      l2_norm += value * value * qf_orig.get_weight(j);
    }
  }
#ifdef _MPI //communicate the results and add up over all processes
  double sendbuf[1] = {l2_norm};
  double recvbuf[1] = {0.0};
  MPI_Allreduce(sendbuf, recvbuf, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  l2_norm = recvbuf[0]; //fill in the summed up values
#endif
  return std::sqrt(l2_norm);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
std::pair<double,double> TFEFunction3D::get_L2_and_H1_norm() const
{
  auto collection = FESpace3D->GetCollection();
  auto n_cells = collection->GetN_Cells();
  double l2_norm = 0.;
  double h1_norm = 0.;
  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTetra); // dummy type
  TQuadFormula qf_orig(qf_ref);
  for(int i_cell = 0; i_cell < n_cells; ++i_cell)
  {
    auto cell = collection->GetCell(i_cell);
#ifdef _MPI // skip halo cells
    if (cell->IsHaloCell())
    {
      continue;
    }
#endif
    auto fe = FESpace3D->get_fe(i_cell);
    auto basis_functions = fe.GetBaseFunct();
    auto n_local_basis_functions = basis_functions->GetDimension();
    bool SecondDer = false;

    // Compute transformation of basis functions (and their derivatives) to 
    // original cell and the evaluation for all quadrature points 
    std::vector<const FiniteElement*> used_fe(1, &fe);
    FEDatabase::GetOrig(used_fe, collection, cell, &SecondDer, qf_ref, qf_orig);
    
    int N_QuadPoints = qf_orig.GetN_QuadPoints();

    // calculate all needed derivatives of this FE function
    std::vector<double> FEFunctValues(n_local_basis_functions);
    const int *DOF = FESpace3D->GetGlobalDOF(i_cell);
    for(int l = 0; l < n_local_basis_functions; l++)
    {
      FEFunctValues[l] = Values[DOF[l]];
    }

    // create a pointer to the QuadPoint-values of the basis on original element and its derivatives
    double **OrigFEValues = FEDatabase::GetOrigElementValues(
        *basis_functions, MultiIndex3D::D000);
    double **OrigFED100 = FEDatabase::GetOrigElementValues(
        *basis_functions, MultiIndex3D::D100);
    double **OrigFED010 = FEDatabase::GetOrigElementValues(
        *basis_functions, MultiIndex3D::D010);
    double **OrigFED001 = FEDatabase::GetOrigElementValues(
        *basis_functions, MultiIndex3D::D001);
    for(int j = 0; j < N_QuadPoints; j++)
    {
      double *Orig = OrigFEValues[j];
      double *Orig100 = OrigFED100[j];
      double *Orig010 = OrigFED010[j];
      double *Orig001 = OrigFED001[j];
      double value = 0.;
      double d100 = 0.;
      double d010 = 0.;
      double d001 = 0.;

      for(int l = 0; l < n_local_basis_functions; l++)
      {
        value += FEFunctValues[l] * Orig[l];
        d100 += FEFunctValues[l] * Orig100[l];
        d010 += FEFunctValues[l] * Orig010[l];
        d001 += FEFunctValues[l] * Orig001[l];
      }
      l2_norm += value * value * qf_orig.get_weight(j);
      h1_norm += (d100*d100 + d010*d010 + d001*d001) * qf_orig.get_weight(j);
    }
  }
#ifdef _MPI //communicate the results and add up over all processes
  double sendbuf[2] = {l2_norm, h1_norm};
  double recvbuf[2] = {0.0, 0.0};
  MPI_Allreduce(sendbuf, recvbuf, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  l2_norm = recvbuf[0]; //fill in the summed up values
  h1_norm = recvbuf[1];
#endif
  return {std::sqrt(l2_norm), std::sqrt(h1_norm)};
}

/** project function into the space L20 (having zero mean value, or in general a mean value) */
void TFEFunction3D::project_into_L20(double a)
{
  // compute current integral and measure of the domain:
  double integral, measure;
  this->compute_integral_and_measure(integral, measure);
  double new_mean = (integral - a)/measure;

  // substract the new mean to this FEFunction3D
  this->add_constant(-new_mean);
}

/**Set Dirichlet values according to boundary conditions*/
void TFEFunction3D::SetDirichletBC(BoundCondFunct3D *BoundaryCondition,
                                   BoundValueFunct3D *BoudaryValue)
{
  int i,j, m;
  int N_Cells, N_Points;
  const double *xi, *eta, *zeta;
  ReferenceTransformation_type F_K = ReferenceTransformation_type::TetraAffin;
  TRefTrans3D *rt;
  double X[MaxN_PointsForNodal3D], Y[MaxN_PointsForNodal3D];
  double Z[MaxN_PointsForNodal3D];
  double PointValues[MaxN_PointsForNodal3D];
  double FunctionalValues[MaxN_PointsForNodal3D];  
  int ApproxOrder;
  bool IsIsoparametric;
  TJoint *joint;
  BFRefElements RefElement;
  ReferenceTransformation_type RefTrans;
  int N_EdgeDOF, N_Joints;
  bool InnerBoundary, OuterBoundary;
  BoundCond Cond0;
  const int *TmpFV, *TmpLen;
  int MaxLen;
  double t0, xf, yf, zf;
  const double *t, *s;
  double LinComb[4] = {0,0,0,0};
  int *EdgeDOF;
  
  
  auto Coll = FESpace3D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  
  for(i=0;i<N_Cells;i++)
  {
    auto cell = Coll->GetCell(i);
    auto Element = FESpace3D->get_fe(i);
    auto nf = Element.GetNodalFunctional();
    nf->GetPointsForAll(N_Points, xi, eta, zeta);
    
    auto DOF = FESpace3D->GetGlobalDOF(i);
    
    ApproxOrder = Element.GetBaseFunct()->GetAccuracy();

    RefElement = Element.GetBaseFunct()->GetRefElement();
    int degree = Element.GetBaseFunct()->GetPolynomialDegree();
    auto QuadFormula = QuadratureFormulaDatabase::qf_from_degree(3*degree,
                                                                 RefElement);
    RefTrans = Element.GetRefTransID();
    
    IsIsoparametric = TDatabase::ParamDB->USE_ISOPARAMETRIC
                      && cell->has_isoparametric_joint();
    if(IsIsoparametric)
    {      
      switch(RefElement)
      {
        case BFRefElements::BFUnitHexahedron:
          RefTrans = ReferenceTransformation_type::HexaIsoparametric;
        break;
  
        case BFRefElements::BFUnitTetrahedron:
          RefTrans = ReferenceTransformation_type::TetraIsoparametric;
        break;
        
        default:
          ;
      }
    } 
    
  
    switch(RefTrans)
    {
      case ReferenceTransformation_type::HexaAffin:
        rt = FEDatabase::GetRefTrans3D(ReferenceTransformation_type::HexaAffin);
        ((THexaAffin *)rt)->SetCell(cell);
        F_K = ReferenceTransformation_type::HexaAffin;
        // cout << "HexaAffin: " << endl;
        break;
      case ReferenceTransformation_type::HexaTrilinear:
        rt = FEDatabase::GetRefTrans3D(
          ReferenceTransformation_type::HexaTrilinear);
        ((THexaTrilinear *)rt)->SetCell(cell);
        F_K = ReferenceTransformation_type::HexaTrilinear;
        // cout << "HexaTrilinear: " << endl;
        break;
      case ReferenceTransformation_type::HexaIsoparametric:
        rt = FEDatabase::GetRefTrans3D(
          ReferenceTransformation_type::HexaIsoparametric);
        ((THexaIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((THexaIsoparametric *)rt)->SetQuadFormula(QuadFormula->get_type());
        ((THexaIsoparametric *)rt)->SetCell(cell);
        F_K = ReferenceTransformation_type::HexaIsoparametric;
        // cout << "HexaIsoparametric: " << endl;
        break;
      case ReferenceTransformation_type::TetraAffin:
        rt = FEDatabase::GetRefTrans3D(
          ReferenceTransformation_type::TetraAffin);
        ((TTetraAffin *)rt)->SetCell(cell);
        F_K = ReferenceTransformation_type::TetraAffin;
        // cout << "TetraAffin: " << endl;
        break;
      case ReferenceTransformation_type::TetraIsoparametric:
        rt = FEDatabase::GetRefTrans3D(
          ReferenceTransformation_type::TetraIsoparametric);
        ((TTetraIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((TTetraIsoparametric *)rt)->SetQuadFormula(QuadFormula->get_type());
        ((TTetraIsoparametric *)rt)->SetCell(cell);
        F_K = ReferenceTransformation_type::TetraIsoparametric;
        // cout << "TetraIsoparametric: " << endl;
        break;
        
      default:
        ;
    }
    
    auto FEDesc_Obj = Element.GetFEDesc();
    N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();
    N_Joints = cell->GetN_Faces();
    
    for(m=0;m<N_Joints;m++)
    {
      joint = cell->GetJoint(m);
      InnerBoundary = false;
      OuterBoundary = false;
      
      if(joint->GetType() == BoundaryFace ||
           joint->GetType() == IsoBoundFace)
          OuterBoundary = true;
          
      if(InnerBoundary || OuterBoundary)
      {
        cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
        
        t0 = 1./TmpLen[m];
        
        xf =0; yf =0; zf =0;
        for(j=0;j<TmpLen[m];j++)
        {
          cell->GetVertex(TmpFV[m*MaxLen+j])->GetCoords(X[j], Y[j], Z[j]);
          
          xf += t0*X[j];
          yf += t0*Y[j];
          zf += t0*Z[j];                    
        }

        auto boundface              = (TBoundFace*)joint;
        const TBoundComp* BoundComp = boundface->GetBoundComp();
        int comp                    = BoundComp->get_physical_id();
        BoundaryCondition(comp, xf, yf, zf, Cond0);

        switch(Cond0)
        {
          case DIRICHLET:
            if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
            {
              /// cout << "USE_ISOPARAMETRIC: " << endl;
              nf->GetPointsForFace(N_Points, t, s);
              for(j=0;j<N_Points;j++)
              {
                switch(TmpLen[m])
                {
                  case 4:
                    LinComb[0] = (1-t[j])*(1-s[j]);
                    LinComb[1] = t[j]*(1-s[j]);
                    LinComb[2] = t[j]*s[j];
                    LinComb[3] = (1-t[j])*s[j];
                    
                    xf = LinComb[0]*X[0] + LinComb[1]*X[1]
                        +LinComb[2]*X[2] + LinComb[3]*X[3];
                    yf = LinComb[0]*Y[0] + LinComb[1]*Y[1]
                        +LinComb[2]*Y[2] + LinComb[3]*Y[3];
                    zf = LinComb[0]*Z[0] + LinComb[1]*Z[1]
                        +LinComb[2]*Z[2] + LinComb[3]*Z[3];
                    break;
                  case 3:
                    LinComb[0] = 1-t[j]-s[j];
                    LinComb[1] = t[j];
                    LinComb[2] = s[j];
                    
                    xf = LinComb[0]*X[0] + LinComb[1]*X[1]
                        +LinComb[2]*X[2];
                    yf = LinComb[0]*Y[0] + LinComb[1]*Y[1]
                        +LinComb[2]*Y[2];
                    zf = LinComb[0]*Z[0] + LinComb[1]*Z[1]
                        +LinComb[2]*Z[2];
                    break;
                }
                BoudaryValue(comp, xf, yf, zf, PointValues[j]);
              }
            }
            else
            {
              /// cout << "Non_ISOPARAMETRIC: " << endl;
              nf->GetPointsForFace(m, N_Points, xi, eta, zeta);
              FEDatabase::GetOrigFromRef(F_K, N_Points, xi, eta, zeta, X, Y, Z);

            
              for(j=0;j<N_Points;j++)             
                BoudaryValue(comp, X[j], Y[j], Z[j], PointValues[j]);
            }
            nf->GetFaceFunctionals(Coll, cell, m, PointValues, FunctionalValues);
            
            EdgeDOF = FEDesc_Obj->GetJointDOF(m);
            N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();
              
            for(j=0;j<N_EdgeDOF;j++)
            {
               Values[DOF[EdgeDOF[j]]] = FunctionalValues[j];
               // cout << i << setw(20) << Values[DOF[EdgeDOF[j]]] << endl;
             }
             break;
	  default:
	    
	  break;	     
	     
         }         
       }///endif
    }///for m<N_Joints
  }///endfor i<N_Cells
}

// this function computes the values of the basis functions on a grid on the 
// reference cell. When using these values it is important to go through them in
// exactly the same order as here.
// nppd = number of points per dim
void update_data(std::vector<double>& bf_values, int nppd,
                 BaseFunction_type& current_type, const BaseFunctions* bf)
{
  if(current_type == bf->GetID())
    return;
  current_type = bf->GetID();
  int tetra = (bf->GetRefElement() == BFRefElements::BFUnitTetrahedron);
  auto reference_element = bf->GetRefElement();
  int n_base_functs = bf->GetDimension();
  int n_points = tetra ? ((nppd+1) * (nppd+2) * (nppd+3) / 6)
                       : (nppd+1) * (nppd+1) * (nppd+1);
  bf_values.resize(n_points * n_base_functs);
  if(reference_element == BFRefElements::BFUnitHexahedron)
  {
    int index = 0;
    for(int ix = 0; ix <= nppd; ++ix)
    {
      double x = -1. + ix*2./nppd;
      for(int iy = 0; iy <= nppd; ++iy)
      {
        double y = -1. + iy*2./nppd;
        for(int iz = 0; iz <= nppd; ++iz)
        {
          double z = -1. + iz*2./nppd;
          bf->GetDerivatives(MultiIndex3D::D000, x, y, z,
                             &bf_values[index*n_base_functs]);
          index++;
        }
      }
    }
  }
  else if(reference_element == BFRefElements::BFUnitTetrahedron)
  {
    int index = 0;
    for(int ix = 0; ix <= nppd; ++ix)
    {
      double x = (double)ix / nppd;
      for(int iy = 0; iy+ix <= nppd; ++iy)
      {
        double y = (double)iy / nppd;
        for(int iz = 0; iz+iy+ix <= nppd; ++iz)
        {
          double z = (double)iz / nppd;
          bf->GetDerivatives(MultiIndex3D::D000, x, y, z,
                             &bf_values[index*n_base_functs]);
          index++;
        }
      }
    }
  }
  else
  {
    ErrThrow("Reference element ", reference_element,
             " not supported here.");
  }
}



// Compute minimal and maximal value of nodal functionals of a Pk element in
// a given cell evaluated at an FE-function, where the Pk element are of the
// same order as the given FE-function. Attention: This works only correctly if
// the reference transformation preserves point values.
std::pair<double, double> compute_cell_min_max_nodal_fctn(const
    TFEFunction3D& fe_function, int cell_i)
{
  auto fespace = fe_function.GetFESpace();
  std::pair<double, double> minmax(1e100, -1e100);
  // Here the nodal functionals of the Pk elements are used that
  // are of the same order to compute the minimum and maximum value in each
  // cell. Due to the Kronecker-Delta property of the Pk elements finding the
  // minimum and maximum corresponds to finding the minimal and maximal value
  // of the FE-vector in the Pk basis. This is in particular easy if the
  // FESpace itself is Pk space. Otherwise, the local FE-vector has to be
  // transformed.
  auto is_pk_space = !fespace->is_discontinuous();
  auto values = fe_function.GetValues();
  auto global_dofs = fespace->GetGlobalDOF(cell_i);
  if (is_pk_space)
  {
    auto n_loc_dof = fespace->get_n_local_dof(cell_i);
    for (unsigned int dof_j = 0; dof_j < n_loc_dof; ++dof_j)
    {
      auto val = values[global_dofs[dof_j]];
      minmax.first = std::min(minmax.first, val);
      minmax.second = std::max(minmax.second, val);
    }
  }
  else
  {
    // The idea here is to locally map the polynomial basis functions of the
    // discontinuous elements onto the basis functions of a sufficiently large
    // continuous finite element and determine the min (max) value in this basis
    // by evaluating the min (max) of the vector.
    // For this a transformation matrix on the reference cell is computed. The
    // FE-vector is then computed by multiplying the FE-vector to the matrix.
    // Actually the values have to be transformed from the reference cell to the
    // original cell, but for value-preserving transformation this is not
    // needed. Therefore, this method works not correctly for instance for Hdiv
    // elements.

    // Get the matrix. This information can be gathered from
    // FEDatabase::GetProlongationMatrix3D
    auto fe_original = fespace->get_fe(cell_i);
    auto fe_type_original = fe_original.GetID();
    auto degree = fe_original.GetBaseFunct()->GetPolynomialDegree();
    auto ele_for_shape = TFESpace3D::get_element_for_shape(degree, 3);
    auto shape = fespace->GetCollection()->GetCell(cell_i)->GetType();
    auto element = ele_for_shape[shape];

    auto matrix = FEDatabase::GetProlongationMatrix3D(fe_type_original,
        NoRef, element, 0);

    // Apply matrix vector multiplication on reference cell. Transformation back
    // to the original cell is not used here, see also above.
    auto n_loc_dof = fespace->get_n_local_dof(cell_i);
    for (unsigned int dof_i = 0; dof_i < n_loc_dof; ++dof_i)
    {
      double val = 0;
      for (unsigned int dof_j = 0; dof_j < n_loc_dof; ++dof_j)
      {
        val += matrix[dof_j + dof_i * MaxN_BaseFunctions3D]
          * values[global_dofs[dof_j]];
      }
      minmax.first = std::min(val, minmax.first);
      minmax.second = std::max(val, minmax.second);
    }
  }
  if (minmax.first > minmax.second)
  {
    Output::warn<1>("TFEFunction3D::compute_cell_wise_min_max_function_eval",
        "Method was not successful on cell ", cell_i);
  }
  return minmax;
}


// Compute minimal and maximal value of uh using point evaluations
std::pair<double, double>
compute_cell_min_max_function_eval(const TFEFunction3D& fe_function,
    std::vector<double>& bf_values, int n_points_per_dim, BaseFunction_type&
    current_type, const BaseFunctions* bf, int cell_i)
{
  int n_base_functs = bf->GetDimension();
  update_data(bf_values, n_points_per_dim, current_type, bf);
  int tetra = (bf->GetRefElement() == BFRefElements::BFUnitTetrahedron);
  int index = 0;
  auto global_dofs = fe_function.GetFESpace()->GetGlobalDOF(cell_i);
  auto values = fe_function.GetValues();
  std::pair<double, double> minmax(1e100, -1e100);
  for(int ix = 0; ix <= n_points_per_dim; ++ix)
  {
    for(int iy = 0; iy+ix*tetra <= n_points_per_dim; ++iy)
    {
      for(int iz = 0; iz+(iy+ix)*tetra <= n_points_per_dim; ++iz)
      {
        double val = 0.;
        for(int i = 0; i < n_base_functs; ++i)
        {
          int global_dof = global_dofs[i];
          val += bf_values[index*n_base_functs + i]*values[global_dof];
        }
        index++;
        minmax.first = std::min(minmax.first, val);
        minmax.second = std::max(minmax.second, val);
      }
    }
  }
  if (minmax.first > minmax.second)
  {
    Output::warn<1>("TFEFunction3D::compute_cell_wise_min_max_function_eval",
        "Method was not successful on cell ", cell_i);
  }
  return minmax;
}



std::pair<double, double> TFEFunction3D::compute_cell_min_max(int cell_nr,
    std::vector<double>& bf_values, BaseFunction_type& current_type, bool
    use_pk_nodal_fctn) const
{
  // Gather some information used in the case use_pk_nodal_fctn == false
  constexpr int n_points_per_dim = 21;

  // Find local minimum and maximum over all cells and update the global minimum
  // and maximum accordingly
  std::pair<double, double> minmax;
  if (use_pk_nodal_fctn)
  {
    minmax = compute_cell_min_max_nodal_fctn(*this, cell_nr);
  }
  else
  {
    auto bf = FESpace3D->get_fe(cell_nr).GetBaseFunct();
    minmax = compute_cell_min_max_function_eval(*this, bf_values,
        n_points_per_dim, current_type, bf, cell_nr);
  }
  return minmax;
}

void TFEFunction3D::MinMax(double & min, double & max, const bool&
    use_pk_nodal_fctn) const
{

  //double t = GetTime();
  if(this->GetFESpace()->GetBaseVectDim() != 1)
  {
    ErrThrow("computing the minimum and maximum of a vector valued function is "
             "not implemented.");
  }
  // Find local minimum and maximum over all cells and update the global minimum
  // and maximum accordingly
  max = -1e100, min = 1e100;

  // Find local minimum and maximum over all cells and update the global minimum
  // and maximum accordingly
  auto n_cells = GetFESpace()->GetCollection()->GetN_Cells();
  std::vector<double> bf_values; // values of basis functions on reference cell
  BaseFunction_type current_type = BF_C_L_P0_1D; // dummy type
  for (int cell_i = 0; cell_i < n_cells; ++cell_i)
  {
#ifdef _MPI
    auto coll = GetFESpace()->GetCollection();
    if(coll->GetCell(cell_i)->IsHaloCell())
    {
      continue;
    }
#endif
    auto minmax = this->compute_cell_min_max(cell_i, bf_values, current_type,
        use_pk_nodal_fctn);
    // Update min and max accordingly
    min = std::min(min, minmax.first);
    max = std::max(max, minmax.second);
  }
#ifdef _MPI
  // reduce min and max in root
  double sendbuf_min = min;
  double sendbuf_max = max;
  MPI_Allreduce(&sendbuf_min, &min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&sendbuf_max, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
  //t = GetTime() - t;
  //Output::print<3>("time for min/max computation: ", t, " seconds");
}

void TFEFunction3D::PrintMinMax(const std::string& name, const bool&
    use_pk_nodal_fctn) const
{
#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#else
  int my_rank = 0;
#endif

  double min, max;
  this->MinMax(min, max, use_pk_nodal_fctn);

  if(my_rank ==0) //only root has results - only root prints
  {
    if( min <= max )
    {
      if(name.empty())
        Output::stat("MinMax", this->Name, " min ", std::setprecision(12), min,
                     ", max ", max);
      else
        Output::stat("MinMax", name, " min ", std::setprecision(12), min,
                     ", max ", max);
    }
    else
    {
        Output::warn("TFEFunction3D::MinMax was not successful!");
    }
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double TFEFunction3D::compute_mean_oscillation(double u_min, double u_max, const
    bool& use_pk_nodal_fctn) const
{
  if(this->GetFESpace()->GetBaseVectDim() != 1)
  {
    ErrThrow("computing the minimum and maximum of a vector valued function is "
        "not implemented.");
  }

  // Gather some information used in the case use_nodal_fctn == false
  constexpr int n_points_per_dim = 21;
  BaseFunction_type current_type = BF_C_L_P0_1D; // dummy type
  std::vector<double> bf_values; // values of basis functions on reference cell

  // Find local minimum and maximum over all cells and update the mean
  // oscillation accordingly
  double mean_osc = 0;
  auto n_cells = GetFESpace()->GetCollection()->GetN_Cells();
  for (int cell_i = 0; cell_i < n_cells; ++cell_i)
  {
#ifdef _MPI
    auto coll = GetFESpace()->GetCollection();
    if(coll->GetCell(cell_i)->IsHaloCell())
    {
      continue;
    }
#endif
    std::pair<double, double> minmax;
    if (use_pk_nodal_fctn)
    {
      minmax = compute_cell_min_max_nodal_fctn(*this, cell_i);
    }
    else
    {
      auto bf = FESpace3D->get_fe(cell_i).GetBaseFunct();
       minmax = compute_cell_min_max_function_eval(*this, bf_values,
          n_points_per_dim, current_type, bf, cell_i);
    }
    // Update min and max accordingly
    mean_osc += std::max(0., minmax.second - u_max) + std::max(0., u_min
        - minmax.first);
  }
#ifdef _MPI
  // reduce min and max in root
  double sendbuf = mean_osc;
  MPI_Allreduce(&sendbuf, &mean_osc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  auto coll = GetFESpace()->GetCollection();
  int n_own_cells = coll->GetN_OwnCells();
  MPI_Allreduce(&n_own_cells, &n_cells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
  mean_osc /= n_cells;
  return mean_osc;
}





