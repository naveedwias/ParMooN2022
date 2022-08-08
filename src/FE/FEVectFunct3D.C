#ifdef _MPI
# include "mpi.h"
#endif

#include <FEVectFunct3D.h>
#include "FEDatabase.h"
#include <HexaAffin.h>
#include <HexaTrilinear.h>
#include <HexaIsoparametric.h>
#include <TetraAffin.h>
#include <TetraIsoparametric.h>
#include <Database.h>
#include <BoundFace.h>
#include "BaseCell.h"
#include <GridCell.h>
#include "AuxParam3D.h"
#include "QuadratureFormulaDatabase.h"

#include <MooNMD_Io.h>
#include <cmath>
#include <cstring> // memset

/// Default constructor. Constructs an empty object.
TFEVectFunct3D::TFEVectFunct3D()
{
  N_Components = 0;
}



/** constructor with vector initialization */
TFEVectFunct3D::TFEVectFunct3D(std::shared_ptr<const TFESpace3D> fespace3D,
                               const std::string& name, double *values,
                               int n_components)
: TFEFunction3D(fespace3D, name, values)
{
  N_Components = n_components;
}

TFEVectFunct3D& TFEVectFunct3D::operator=( const TFEVectFunct3D & other)
{
  //call base class copy assignment
  TFEFunction3D::operator=(other);

  this->N_Components  = other.N_Components;

  return *this;
}

/** compute integral and measure */
double TFEVectFunct3D::compute_flux(int surface_id) const
{
  const int dim = 3;
  double flux = 0.;

  auto coll = FESpace3D->GetCollection();
  std::array<std::unique_ptr<TFEFunction3D>, dim> fe_fct_component;
  for(int icoor=0; icoor<dim; icoor++)
  {
    fe_fct_component[icoor] = this->GetComponent(icoor);
  }

  for(int i=0; i< coll->GetN_Cells(); i++)
  {
    auto cell = coll->GetCell(i); //boundaryCells[i];

#ifdef _MPI
    // only compute flux on own cells
    if( cell->IsHaloCell() )
    {
      continue;
    }
#endif

    const int* DOF = FESpace3D->GetGlobalDOF(cell->GetCellIndex());

    for(size_t joint_id=0; joint_id< (size_t) cell->GetN_Faces(); joint_id++)
    {
      TJoint* joint = cell->GetJoint(joint_id);

      if(  joint->GetType() == BoundaryFace
        || joint->GetType() == IsoBoundFace )
      {
        // convert the joint to an object of BoundFace type
        TBoundFace* boundface = (TBoundFace*)joint;

        // check if the face is on the desired component
        if( boundface->GetBoundComp()->get_physical_id() == surface_id )
        {
          // ===================
          // get quadrature data
          // ===================
          // set quadrature formula and compute quadrature info
          auto fe = FESpace3D->get_fe(cell->GetCellIndex());
          int fe_degree = fe.GetBaseFunct()->GetPolynomialDegree();
          const Shapes* face_types;
          cell->GetShapeDesc()->GetFaceType(face_types);
          // get a quadrature formula good enough for the velocity FE space
          const TQuadFormula* qf2 = QuadratureFormulaDatabase::qf_from_degree(
              2*fe_degree, face_types[joint_id]);
          int N_Points = qf2->GetN_QuadPoints();

          // ===============================================================
          // generate data on reference mesh cell for the 2d face of 3d cell
          // ===============================================================
          // values of base functions in all quadrature points on face
          double** JointValues = FEDatabase::GetJointDerivatives3D(
            *fe.GetBaseFunct(), *qf2, joint_id, MultiIndex3D::D000);

          fe.GetBaseFunct()->ChangeBF(FESpace3D->GetCollection(), cell,
                                      N_Points, JointValues);

          // compute normal vector
          std::vector<double> normal;
          double transformationDeterminant;

          cell->computeNormalAndTransformationData(joint_id,
                                                   normal,
                                                  transformationDeterminant);

          // note: the normal computed above is not always directed outward
          //       (for boundary cells)
          normal.resize(dim);
          double x, y, z;
          ///@attention we assume that the bound.face is planar
          auto p0 = qf2->get_point(0);
          boundface->GetXYZofTS(p0.x, p0.y, x, y, z);
          boundface->get_normal_vector(x,y,z,normal[0],normal[1],normal[2]);
          Output::print<5>(" ** computed normal vector on (", x, ",", y, ",", z,
                           ") => n = (", normal[0], ",", normal[1],
                           ",", normal[2], ")");

          // compute \int_F u.n = sum_{gauss pt} w_k \sum_j u_j.n phi_j(x_k)
          double value = 0;
          for(int l=0; l < N_Points; l++)
          {
            double u_n = 0;

            // compute u.n on l-th Gauss point
            for(int k=0; k<fe.GetN_DOF(); k++)
            {
              int global_dof_from_local = DOF[k];
              for(int icoor=0; icoor<dim; icoor++)
              {
                double* u_icoor_values = fe_fct_component[icoor]->GetValues();
                double u_icoor_on_x_k = u_icoor_values[ global_dof_from_local ];
                u_n += JointValues[l][k]* (u_icoor_on_x_k*normal[icoor]) ;
              }
            }
            value += qf2->get_weight(l) * transformationDeterminant * u_n;
          }
          flux += value;
        }
      }
    }
  }

#ifdef _MPI
  // sum the contribution of each process to the total flux
  MPI_Allreduce(MPI_IN_PLACE, &flux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  return flux;
}

/** calculate errors to given vector function */
void TFEVectFunct3D::GetDeformationTensorErrors( 
  DoubleFunct3D *Exact, DoubleFunct3D *Exact1,
  DoubleFunct3D *Exact2,
  int N_Derivatives,
  MultiIndex3D *NeededDerivatives,
  int N_Errors, TFEFunction3D::ErrorMethod *ErrorMeth, 
  const CoeffFct3D& Coeff, 
  TAuxParam3D *Aux,
  int n_fespaces, TFESpace3D **fespaces,
  double *errors)
{
  int i,j,k,l;
  int N_Cells, N_Points, N_Parameters, N_;
  TFESpace3D *fespace;
  const double *weights;
  const double *X, *Y, *Z;
  std::vector<double> AbsDetjk(MaxN_QuadPoints_3D, 1.);
  double *Param[MaxN_QuadPoints_3D], *aux, *aux1, *aux2, *aux3;
  double *Derivatives[3*MaxN_QuadPoints_3D];
  double *ExactVal[3*MaxN_QuadPoints_3D];
  double *AuxArray[MaxN_QuadPoints_3D];
  double **OrigFEValues, *Orig, value, value1, value2;
  double FEFunctValues[MaxN_BaseFunctions3D];
  double FEFunctValues1[MaxN_BaseFunctions3D];
  double FEFunctValues2[MaxN_BaseFunctions3D];
  double LocError[4], *Values0,*Values1, *Values2;
  double hK;
  bool *SecondDer;

  SecondDer = new bool[n_fespaces];
  for(i=0;i<n_fespaces;i++)
    SecondDer[i] = false;

  N_Parameters = Aux->GetN_Parameters();
  aux1 = new double [MaxN_QuadPoints_3D*N_Parameters];
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Param[j] = aux1 + j*N_Parameters;

  aux2 = new double [3*MaxN_QuadPoints_3D*N_Derivatives];
  for(j=0;j<3*MaxN_QuadPoints_3D;j++)
    Derivatives[j] = aux2 + j*N_Derivatives;
  
  aux3 = new double [3*MaxN_QuadPoints_3D * 4];
  for(j=0;j<3*MaxN_QuadPoints_3D;j++)
    ExactVal[j] = aux3 + j*4;

  // 20 <= number of term
  aux = new double [MaxN_QuadPoints_3D*20]; 
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    AuxArray[j] = aux + j*20;

  int Length = FESpace3D->get_n_dof();
  fespace = fespaces[0];
  Values0 = Values;
  Values1 = Values+Length;
  Values2 = Values1+Length;

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
    auto fe = fespace->get_fe(i);
    auto& bf = *fe.GetBaseFunct();
    N_ = fe.GetN_DOF();

    auto DOF = fespace->GetGlobalDOF(i);
    for(l=0;l<N_;l++)
    {
      FEFunctValues[l] = Values0[DOF[l]];
      FEFunctValues1[l] = Values1[DOF[l]];
      FEFunctValues2[l] = Values2[DOF[l]];
    }

    // for all needed derivatives
    for(k=0;k<N_Derivatives;k++)
    {
      OrigFEValues = FEDatabase::GetOrigElementValues(bf, NeededDerivatives[k]);
      // for all quadrature points
      for(j=0;j<N_Points;j++)
      {
        Orig = OrigFEValues[j];
        value = 0;
        value1 = 0;
        value2 = 0;
        for(l=0;l<N_;l++)
        {
          value += FEFunctValues[l] * Orig[l];
          value1 += FEFunctValues1[l] * Orig[l];
          value2 += FEFunctValues2[l] * Orig[l];
        } // endfor l
        Derivatives[j][k] = value;
        Derivatives[j+N_Points][k] = value1;
        Derivatives[j+2*N_Points][k] = value2;
      } // endfor j
    } // endfor k

    // exact value for first component
    for(j=0;j<N_Points;j++)
    {
      auto p = qf_orig.get_point(j);
      Exact(p.x, p.y, p.z, ExactVal[j]);
      Exact1(p.x, p.y, p.z, ExactVal[j+N_Points]);
      Exact2(p.x, p.y, p.z, ExactVal[j+2*N_Points]);
    }

    if(Coeff)
      Coeff(N_Points, X, Y, Z, Param, AuxArray);      

    ErrorMeth(N_Points, {{X, Y, Z}}, AbsDetjk.data(), weights, hK, Derivatives, 
              ExactVal, AuxArray, LocError);

    for(j=0;j<N_Errors;j++)
      errors[j] += LocError[j];

  } // endfor i

  for(j=0;j<N_Errors;j++)
  {
    if (errors[j]>0)
      errors[j] = std::sqrt(errors[j]);
  }

  delete aux;
  delete aux1;
  delete aux2;
  delete aux3;
  delete SecondDer;
} // TFEFunction3D::GetDeformationTensorErrors


double TFEVectFunct3D::GetL2NormDivergenceError(DoubleFunct3D* Exact_u1,
                                                DoubleFunct3D* Exact_u2,
                                                DoubleFunct3D* Exact_u3)
{
  std::vector<double> values(1);
  std::array<double, 5> exact_sol; // avoid re-allocations
  auto f = [Exact_u1, Exact_u2, Exact_u3, &exact_sol]
           (std::vector<double>& v, std::array<double, 15> e)
           {
             double local_divergence_fe = e[6] + e[10] + e[14];
             Exact_u1(e[0], e[1], e[2], exact_sol.data());
             double local_divergence_exact = exact_sol[1];
             Exact_u2(e[0], e[1], e[2], exact_sol.data());
             local_divergence_exact += exact_sol[2];
             Exact_u3(e[0], e[1], e[2], exact_sol.data());
             local_divergence_exact += exact_sol[3];
             v[0] += local_divergence_fe - local_divergence_exact;
           };
  this->get_functional_value(values, f);
  return values[0];
} // TFEVectFunct3D::GetL2NormDivergenceError

//================================================================
std::pair<double, double> TFEVectFunct3D::get_L2_norm_divergence_curl() const
{
  std::vector<double> values(2);
  auto f = [](std::vector<double>& v, std::array<double, 15> e)
           {
             v[0] += e[6] + e[10] + e[14]; // divergence
             double curl_u_1 = e[11] - e[13];  // (u3y - u2z);
             double curl_u_2 = e[12] - e[8];   // (u1z - u3x); 
             double curl_u_3 = e[7]  - e[9];   // (u2x - u1y);
             double norm_curl_squared = curl_u_1*curl_u_1 + 
                                        curl_u_2*curl_u_2 +
                                        curl_u_3*curl_u_3;
             v[1] += std::sqrt(norm_curl_squared); // curl
           };
  this->get_functional_value(values, f);
  return {values[0], values[1]};

} // TFEVectFunct3D::get_L2_norm_divergence_curl

void TFEVectFunct3D::get_functional_value(std::vector<double>& values,
  const std::function<void(std::vector<double>&, std::array<double, 15>)>&
    functional)
const
{
  bool SecondDer[1] = { false };
  double FEValues0[MaxN_BaseFunctions3D];
  double FEValues1[MaxN_BaseFunctions3D];
  double FEValues2[MaxN_BaseFunctions3D];
  
  auto n_values = values.size();
  std::fill(values.begin(), values.end(), 0.0);

  int Length = FESpace3D->get_n_dof();
  double *Values0 = Values;
  double *Values1 = Values+  Length;
  double *Values2 = Values+2*Length;

  auto Coll = FESpace3D->GetCollection();
  auto N_Cells = Coll->GetN_Cells();
  
  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTetra); // dummy type
  TQuadFormula qf_orig(qf_ref);
  for(int i = 0; i < N_Cells; i++)
  {
    auto cell = Coll->GetCell(i);
#ifdef _MPI
    if(cell->IsHaloCell())
    {
      continue;
    }
#endif
    auto& fe = FESpace3D->get_fe(i);
    FEDatabase::GetOrig({&fe}, Coll, cell, SecondDer, qf_ref, qf_orig);
    int N_Points = qf_orig.GetN_QuadPoints();

    // calculate all needed derivatives of this FE function
    auto& BaseFunct = *fe.GetBaseFunct();
    auto N_Bf = fe.GetN_DOF();

    auto DOF = FESpace3D->GetGlobalDOF(i);
    
    for(int j = 0; j < N_Bf; j++)
    {
      int k = DOF[j];
      FEValues0[j] = Values0[k];
      FEValues1[j] = Values1[k];
      FEValues2[j] = Values2[k];
    }

    auto OrigFEValues  = FEDatabase::GetOrigElementValues(BaseFunct,
                                                          MultiIndex3D::D000);
    auto OrigFEValuesX = FEDatabase::GetOrigElementValues(BaseFunct,
                                                          MultiIndex3D::D100);
    auto OrigFEValuesY = FEDatabase::GetOrigElementValues(BaseFunct,
                                                          MultiIndex3D::D010);
    auto OrigFEValuesZ = FEDatabase::GetOrigElementValues(BaseFunct,
                                                          MultiIndex3D::D001);
    

    // for all quadrature points
    for(int j = 0; j < N_Points; j++)
    {
      auto p = qf_orig.get_point(j);
      double * Orig = OrigFEValues[j];
      double * OrigX = OrigFEValuesX[j];
      double * OrigY = OrigFEValuesY[j];
      double * OrigZ = OrigFEValuesZ[j];
      std::vector<double> local_values(n_values, 0.0);
      std::array<double, 15> evaluations{{p.x, p.y, p.z,
                                          0., 0., 0., 0., 0., 0.,
                                          0., 0., 0., 0., 0., 0.}};
      for(int l = 0; l < N_Bf; l++)
      {
        evaluations[3]  += Orig[l]  * FEValues0[l];
        evaluations[4]  += Orig[l]  * FEValues1[l];
        evaluations[5]  += Orig[l]  * FEValues2[l];
        evaluations[6]  += OrigX[l] * FEValues0[l];
        evaluations[7]  += OrigX[l] * FEValues1[l];
        evaluations[8]  += OrigX[l] * FEValues2[l];
        evaluations[9]  += OrigY[l] * FEValues0[l];
        evaluations[10] += OrigY[l] * FEValues1[l];
        evaluations[11] += OrigY[l] * FEValues2[l];
        evaluations[12] += OrigZ[l] * FEValues0[l];
        evaluations[13] += OrigZ[l] * FEValues1[l];
        evaluations[14] += OrigZ[l] * FEValues2[l];
      } // endfor l
      functional(local_values, evaluations);
      double local_weight = qf_orig.get_weight(j);
      for(auto k = 0ul; k < n_values; ++k)
      {
        values[k] += local_weight * (local_values[k] * local_values[k]);
      }
    } // endfor j
  } // endfor i
#ifdef _MPI
  // sum the contribution of each process to the total error
  MPI_Allreduce(MPI_IN_PLACE, values.data(), n_values, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
#endif
  for(auto k = 0ul; k < n_values; ++k)
  {
    values[k] = std::sqrt(values[k]);
  }
} // TFEVectFunct3D::get_functional_value

void TFEVectFunct3D::FindValueLocal(const TBaseCell* cell, int cell_no, 
				    double x, double y, double z, 
				    double* values) const
{
 this->TFEFunction3D::FindValueLocal(cell, cell_no, x, y, z, values);
 auto u2 = this->GetComponent(1);
 u2->FindValueLocal(cell, cell_no, x, y, z, values+1);
 auto u3 = this->GetComponent(2);
 u3->FindValueLocal(cell, cell_no, x, y, z, values+2);
}


void TFEVectFunct3D::FindVectGradientLocal(int cell_no, double x,
                                           double y, double z,
                                           double* val1, double* val2, double* val3) const
{
  FiniteElement FE_Obj = FESpace3D->get_fe(cell_no);
  
  const BaseFunctions* bf = FE_Obj.GetBaseFunct();
  int N_BaseFunct = bf->GetDimension();
  
  double *uorig = new double[N_BaseFunct];
  double *uxorig = new double[N_BaseFunct];
  double *uyorig = new double[N_BaseFunct];
  double *uzorig = new double[N_BaseFunct];
  
  getOrigValuesForCell(cell_no, x, y, z, uorig, uxorig, uyorig, uzorig);
  
  getValuesFromOriginalOnes(val1, 0, N_BaseFunct, cell_no, uorig, uxorig, uyorig, uzorig);
  getValuesFromOriginalOnes(val2, FESpace3D->get_n_dof(), N_BaseFunct,
                            cell_no, uorig, uxorig, uyorig, uzorig);
  getValuesFromOriginalOnes(val3, 2*FESpace3D->get_n_dof(), N_BaseFunct,
                            cell_no, uorig, uxorig, uyorig, uzorig);
  
  delete [] uorig;
  delete [] uxorig;
  delete [] uyorig;
  delete [] uzorig;
}


bool TFEVectFunct3D::FindVectGradient(double x, double y, double z,
                                      std::vector<double>& val1,
                                      std::vector<double>& val2,
                                      std::vector<double>& val3) const
{
  double valuesFromCells1[4];
  double valuesFromCells2[4];
  double valuesFromCells3[4];
  
  int N_Found = 0;
  
  for(int i = 0; i < 4; i++)
  {
    val1.at(i) = 0.0;
    val2.at(i) = 0.0;
    val3.at(i) = 0.0;
  }
  
  const TCollection* Coll = FESpace3D->GetCollection();
  int N_Cells = Coll->GetN_Cells();
  
  for(int i = 0; i < N_Cells; i++)
  {
    TBaseCell* cell = Coll->GetCell(i);
    
    if(cell->PointInCell(parmoon::Point(x,y,z)))
    {
      N_Found++;
      
      FindVectGradientLocal(i, x, y, z, valuesFromCells1, valuesFromCells2, valuesFromCells3);
      
      for(int j = 0; j < 3; j++)
      {
        val1[j] += valuesFromCells1[j];
        val2[j] += valuesFromCells2[j];
        val3[j] += valuesFromCells3[j];
      }
    }
  }
  
  
  if(N_Found>0)
  {
    for(int i = 0; i < 4; i++)
    {
      val1[i] /= (double)N_Found;
      val2[i] /= (double)N_Found;
      val3[i] /= (double)N_Found;
    }
    return true;
  }
  else 
    return false;
}

void TFEVectFunct3D::Interpolate(TFEVectFunct3D *OldVectFunct)
{
  int N_Points;

  const double *xi, *eta, *zeta;
  double X[MaxN_PointsForNodal3D], Y[MaxN_PointsForNodal3D], Z[MaxN_PointsForNodal3D];
  double PointValuesx[MaxN_PointsForNodal3D];
  double PointValuesy[MaxN_PointsForNodal3D];
  double PointValuesz[MaxN_PointsForNodal3D];
  double FunctionalValues[MaxN_PointsForNodal3D];
  std::vector<double> valx(4);
  std::vector<double> valy(4);
  std::vector<double> valz(4);

  
  const TCollection* Coll = FESpace3D->GetCollection();
  int N_Cells = Coll->GetN_Cells();
  int N_DOFs = FESpace3D->get_n_dof();


  for(int i = 0; i < N_Cells; i++)
  {
    TBaseCell* cell = Coll->GetCell(i);
    
    const FiniteElement& Element = FESpace3D->get_fe(i);
    const NodalFunctional* nf = Element.GetNodalFunctional();
    nf->GetPointsForAll(N_Points, xi, eta, zeta);
    int N_LocalDOFs = Element.GetN_DOF();

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
          ;
      }
    }
    
    
    FEDatabase::SetCellForRefTrans(cell, RefTrans);
    FEDatabase::GetOrigFromRef(RefTrans, N_Points, xi, eta, zeta, X, Y, Z);
    
    for(int j = 0; j < N_Points; j++)
    {
      OldVectFunct->FindVectGradient(X[j], Y[j], Z[j], valx, valy, valz);
      PointValuesx[j] = valx[0]; 
      PointValuesy[j] = valy[0]; 
      PointValuesz[j] = valz[0]; 
    }

    const int* DOF = FESpace3D->GetGlobalDOF(i);

    nf->GetAllFunctionals(Coll, (TGridCell *)cell, PointValuesx, FunctionalValues);
    for(int j = 0; j < N_LocalDOFs; j++)
      Values[DOF[j]] = FunctionalValues[j];
    
    nf->GetAllFunctionals(Coll, (TGridCell *)cell, PointValuesy, FunctionalValues);
    for(int j = 0; j < N_LocalDOFs; j++)
      Values[N_DOFs + DOF[j]] = FunctionalValues[j];
    
    nf->GetAllFunctionals(Coll, (TGridCell *)cell, PointValuesz, FunctionalValues);
    for(int j = 0; j < N_LocalDOFs; j++)
      Values[N_DOFs*2 + DOF[j]] = FunctionalValues[j];
    
  } //for i
}


