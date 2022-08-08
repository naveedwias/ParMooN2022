/**
 * @brief A test program to test some H(div) finite elements
 *
 * If the nodal functionals are \f$N_i\f$ and the basis functions \f$b_i\f$ 
 * then, this function checks that
 * \f[ N_i(b_j) = \delta_{ij} \f]
 * with the Kronecker \f$\delta\f$.
 */
#include <cmath>
#include <Domain.h>
#include <Database.h>
#include "FEDatabase.h"
#include "BaseCell.h"
#include "FiniteElement.h"
#include "ParMooN.h"
#include <list>

// =======================================================================
// main program
// =======================================================================

bool CheckFEOnCell(const FiniteElement &hdiv_fe, const TBaseCell & cell, 
                   TCollection & coll)
{
  // get the points on the reference element, at which we have to evaluate the 
  // basis functions in order to evaluate the nodal functional
  int N_Points;
  const double *xi, *eta;
  hdiv_fe.GetNodalFunctional()->GetPointsForAll(N_Points, xi, eta);
  
  // the basis functions object
  const BaseFunctions* bf = hdiv_fe.GetBaseFunct();
  
  // evaluate the basis functions at these points (xi[i],eta[i])
  // dimension of the basis function (usually 1, for H(div) elements it is 2)
  int baseVectDim = bf->GetBaseVectDim();
  // number of basis functions
  int nDof = hdiv_fe.GetN_DOF();
  // number of basis functions, this is the length of the array needed to 
  // evaluate the basis functions (therefore the factor baseVectDim)
  int nBaseFunc = nDof * baseVectDim;
  
  // the id of the reference transformation
  ReferenceTransformation_type refTransID = hdiv_fe.GetRefTransID();
  FEDatabase::SetCellForRefTrans(&cell, refTransID);
  
  std::vector<std::vector<double>> uorig(
    N_Points, std::vector<double>(baseVectDim*nDof));
  std::vector<std::vector<double>> AllPointValues(
    N_Points, std::vector<double>(nBaseFunc));
  for(int k = 0; k < N_Points; k++)
  {
    bf->GetDerivatives(MultiIndex2D::D00, xi[k], eta[k],
                       AllPointValues[k].data());
    // Piola transform
    FEDatabase::GetOrigValues(refTransID, xi[k], eta[k], bf, &coll,
                              (const TGridCell *) &cell,
                              AllPointValues[k].data(), nullptr, nullptr,
                              uorig[k].data(), nullptr, nullptr);
  }
  
  std::vector<double> PointValues(N_Points * baseVectDim);
  std::vector<double> FunctionalValues(nDof);
  bool ret = true;
  for(int k = 0; k < nDof; k++)
  {
    for(int l = 0; l < N_Points; l++)
    {
      for(int i = 0; i < baseVectDim; ++i)
      {
        PointValues[l + i * N_Points] = uorig[l][k + i*nDof];
      }
    }
    
    hdiv_fe.GetNodalFunctional()->GetAllFunctionals(&coll, &cell,
                                                    PointValues.data(),
                                                    FunctionalValues.data());
    
    for(int i = 0; i < nDof; i++)
    {
      if( std::abs(FunctionalValues[i]) < 1e-10 )
      {
        FunctionalValues[i] = 0;
      }
      //Output::print(k, " ", i, " ", FunctionalValues[i]);
      if( i == k && std::abs(FunctionalValues[i]-1) > 1e-8 )
      {
        Output::print("basis function: ", k, " nodal functional: ", i, " ", 
                      FunctionalValues[i]);
        ret = false;
      }
      if( i != k && std::abs(FunctionalValues[i]-0) > 1e-8 )
      {
        Output::print("basis function: ", k, " nodal functional: ", i, " ",
                      FunctionalValues[i]);
        ret = false;
      }
    }
    if(ret == false)
      break;
  }
  
  return ret;
}




// check the properties on a given mesh and given elements
bool check(TDomain & domain, const std::list<FE_type>& elements)
{
  // call FiniteElement::CheckNFandBF on all these finite elements
  for(auto e : elements)
  {
    Output::print("starting with ", e, " on the reference cell");
    const FiniteElement hdiv_fe(e);
    // this checks only on the reference cell and does not return true or false
    hdiv_fe.CheckNFandBF();
  }
  
#ifdef __2D__
  unsigned int nRefinements = 3;
  // refine grid up to the coarsest level
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  
  // get the collection of cells on the finest mesh
  TCollection * coll = domain.GetCollection(It_Finest, 0);
  
  
  unsigned int nCells = coll->GetN_Cells();
  // compute global normals on all cells
  for(unsigned int c = 0; c < nCells; ++c)
  {
    coll->GetCell(c)->SetNormalOrientation();
  }
  
  for(auto e : elements)
  {
    Output::print("starting with ", e, " on a mesh");
    const FiniteElement hdiv_fe(e);
    for(unsigned int c = 0; c < nCells; ++c)
    {
      bool succesful = CheckFEOnCell(hdiv_fe, *coll->GetCell(c), *coll);
      if(!succesful) // this test failed
      {
        Output::print("test failed with element ", e, " on cell ", c);
        return false;
      }
    }
  }
  delete coll;
#else // 3D
  (void)domain; // avoid compiler warning
#endif // 2D
  
  return true;
}


int main(int, char**)
{
  auto db = parmoon::parmoon_initialize();
#ifdef __2D__
  db.add("refinement_n_initial_steps", (size_t) 1,"");

  Output::print("************************************************************");
  Output::print("\ntest with quads");
  {
    // the domain is initialised with default description and default
    // initial mesh
    db.add("boundary_file", "Default_UnitSquare", "");
    db.add("geo_file", "UnitSquare", "", {"UnitSquare", "TwoTriangles"});
    TDomain domain(db);
    
    std::list<FE_type> elements = { N_RT0_2D_Q_M, N_RT1_2D_Q_M, N_RT2_2D_Q_M, 
                                 N_RT3_2D_Q_M, N_BDM1_2D_Q_M, N_BDM2_2D_Q_M,
                                 N_BDM3_2D_Q_M };
    
    if(!check(domain, elements))
      return 1;
  }
  
  Output::print("************************************************************");
  Output::print("test with triangles");
  {
    // the domain is initialised with default description and default
    // initial mesh
    db["geo_file"]= "TwoTriangles";
    TDomain domain(db);
    
    std::list<FE_type> elements = { N_RT0_2D_T_A, N_RT1_2D_T_A, N_RT2_2D_T_A,
                                    N_RT3_2D_T_A, N_BDM1_2D_T_A, N_BDM2_2D_T_A,
                                    N_BDM3_2D_T_A };
    
    if(!check(domain, elements))
      return 1;
  }
  
#else // 2D -> 3D
  Output::print("************************************************************");
  Output::print("\ntest with hexahedra");
  {
    // the domain is initialised with default description and default
    // initial mesh
    db.add("boundary_file", "Default_UnitCube", "");
    db.add("geo_file", "Default_UnitCube_Hexa", ""); // Default_UnitCube_Tetra
    TDomain domain(db);
    
    std::list<FE_type> elements = { N_RT0_3D_H_A, N_RT1_3D_H_A };
    
    if(!check(domain, elements))
      return 1;
  }
#endif // 3D
  parmoon::parmoon_finalize();
}
