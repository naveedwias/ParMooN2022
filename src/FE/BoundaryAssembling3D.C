// ============================================================================ //
// @(#)BoundaryAssembling3D.C        28.10.16  revision started 10.10.18                                 //
//                                                                              //
// Functions for (external and internal) boundary integral                      //
//                                                                              //
// ============================================================================ //

#include "FEDatabase.h"
#include <FEVectFunct3D.h>
#include <BoundaryAssembling3D.h>
#include <Collection.h>
#include <BoundFace.h>
#include <BaseCell.h>
#include "Database.h"
#include "QuadratureFormulaDatabase.h"
#include <cmath>

#if _MPI
#include <mpi.h>
#endif

using namespace parmoon;

void BoundaryAssembling3D::rhs_directional_do_nothing(BlockVector &rhs,
         const TFESpace3D *U_Space, BoundValueFunct3D *given_boundary_data,
         const TFEVectFunct3D &velocity,
         std::vector<TBoundFace*>& boundaryFaceList, int componentID,
         double mult)
{
  if (boundaryFaceList.empty())
  {
    auto coll = U_Space->GetCollection();
    coll->get_face_list_on_component(componentID, boundaryFaceList);
  }

  // n has to be a vector because of the call to
  // computeNormalAndTransformationData below
  std::vector<double> n(3);
  std::vector<double> u(3);

  for (size_t m = 0; m < boundaryFaceList.size(); m++)
  {
    TBoundFace *boundface = boundaryFaceList[m];
    const TBoundComp* BoundComp = boundface->GetBoundComp();
    int comp = BoundComp->get_physical_id();
    TBaseCell *cell = ((TJoint *)boundface)->GetNeighbour(0);

    int icell = cell->GetCellIndex();
    // mapping from local (cell) DOF to global DOF
    const int *DOF = U_Space->GetGlobalDOF(icell);

    int joint_id = boundface->get_index_in_neighbour(cell);

    // get all data necessary for computing the integral:
    // quadrature weights, points, functions values, normal, determinant

    std::vector<double> qWeights, qPointsT, qPointsS;
    std::vector<std::vector<double>> basisFunctionsValues;

    this->getQuadratureData(U_Space, cell, joint_id, qWeights,
          qPointsT, qPointsS, basisFunctionsValues);

    double transformationDeterminant;
    cell->computeNormalAndTransformationData(joint_id, n,
               transformationDeterminant);

    double x, y, z;
    double value;

    boundface->GetXYZofTS(qPointsT[0], qPointsS[0], x, y, z);
    boundface->get_normal_vector(x, y, z, n[0], n[1], n[2]);

    // loop over Gauss points
    for (size_t l = 0; l < qWeights.size(); l++)
    {
      boundface->GetXYZofTS(qPointsT[l], qPointsS[l], x, y, z);

      if (given_boundary_data != nullptr)
      {
        given_boundary_data(comp, x, y, z, value);
      }
      else
      {
        value = 1.0;
      }

      velocity.FindValueLocal(cell, icell, x, y, z, u.data());

      double u_dot_n = u[0] * n[0] + u[1] * n[1] + u[2] * n[2];

      if (u_dot_n >= 0.0)
      {
        continue;
      }

      double scale_factor = 0.5 * u_dot_n * mult * value
        * qWeights[l] * transformationDeterminant;

      for (size_t k = 0; k < basisFunctionsValues[l].size(); k++)
      {
        int global_dof_from_local = DOF[k];

        // if the DOF is Dirichlet, continue
        if (global_dof_from_local < U_Space->get_n_active())
        {
          double v = scale_factor * basisFunctionsValues[l][k];

          rhs.block(0)[global_dof_from_local] += v * u[0];
          rhs.block(1)[global_dof_from_local] += v * u[1];
          rhs.block(2)[global_dof_from_local] += v * u[2];
        }
      }
    }
  }
}

void BoundaryAssembling3D::rhs_directional_do_nothing_smoothstep(BlockVector &rhs,
         const TFESpace3D *U_Space,
         const TFEVectFunct3D &velocity,
         const TFEVectFunct3D &dt_velocity,
         std::vector<TBoundFace*>& boundaryFaceList,
         int componentID,
         double nu_D0, double beta, double delta)
{
  if (boundaryFaceList.empty())
  {
    auto coll = U_Space->GetCollection();
    coll->get_face_list_on_component(componentID, boundaryFaceList);
  }

  // n has to be a vector because of the call to
  // computeNormalAndTransformationData below
  std::vector<double> n(3);
  std::vector<double> u(3);
  std::vector<double> dtu(3);

  for (size_t m = 0; m < boundaryFaceList.size(); m++)
  {
    TBoundFace *boundface = boundaryFaceList[m];
    TBaseCell *cell = ((TJoint *)boundface)->GetNeighbour(0);

    int icell = cell->GetCellIndex();
    // mapping from local (cell) DOF to global DOF
    const int *DOF = U_Space->GetGlobalDOF(icell);

    int joint_id = boundface->get_index_in_neighbour(cell);

    // get all data necessary for computing the integral:
    // quadrature weights, points, functions values, normal, determinant

    std::vector<double> qWeights, qPointsT, qPointsS;
    std::vector<std::vector<double>> basisFunctionsValues;

    this->getQuadratureData(U_Space, cell, joint_id, qWeights,
          qPointsT, qPointsS, basisFunctionsValues);

    double transformationDeterminant;
    cell->computeNormalAndTransformationData(joint_id, n,
               transformationDeterminant);

    double x, y, z;

    boundface->GetXYZofTS(qPointsT[0], qPointsS[0], x, y, z);
    boundface->get_normal_vector(x, y, z, n[0], n[1], n[2]);

    // loop over Gauss points
    for (size_t l = 0; l < qWeights.size(); l++)
    {
      boundface->GetXYZofTS(qPointsT[l], qPointsS[l], x, y, z);

      velocity.FindValueLocal(cell, icell, x, y, z, u.data());
      dt_velocity.FindValueLocal(cell, icell, x, y, z, dtu.data());

      double u_dot_n = u[0] * n[0] + u[1] * n[1] + u[2] * n[2];
      double u_dot_u = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
      double Theta = 0.0;

      if (delta > 0.0)
      {
        Theta = 0.5 * (1.0 - std::tanh(u_dot_n / delta));
      }
      else if (u_dot_n < 0.0)
      {
        Theta = 1.0;
      }

      double u_term_mul = 0.5 * beta * Theta * u_dot_n;
      double n_term_mul = 0.5 * beta * Theta * u_dot_u;
      double dtu_term_mul = -nu_D0;

      double w = qWeights[l] * transformationDeterminant;

      for (size_t k = 0; k < basisFunctionsValues[l].size(); k++)
      {
        int global_dof_from_local = DOF[k];

        // if the DOF is Dirichlet, continue
        if (global_dof_from_local < U_Space->get_n_active())
        {
          double v = w * basisFunctionsValues[l][k];
          double vu = u_term_mul * v;
          double vn = n_term_mul * v;
          double vdtu = dtu_term_mul * v;

          rhs.block(0)[global_dof_from_local] += vu * u[0] + vn * n[0] + vdtu * dtu[0];
          rhs.block(1)[global_dof_from_local] += vu * u[1] + vn * n[1] + vdtu * dtu[1];
          rhs.block(2)[global_dof_from_local] += vu * u[2] + vn * n[2] + vdtu * dtu[2];
        }
      }
    }
  }
}

void BoundaryAssembling3D::rhs_directional_do_nothing_smoothstep(BlockVector &rhs,
         const TFESpace3D *U_Space,
         const TFEVectFunct3D &velocity,
         std::vector<TBoundFace*>& boundaryFaceList,
         int componentID,
         double beta, double delta)
{
  if (boundaryFaceList.empty())
  {
    auto coll = U_Space->GetCollection();
    coll->get_face_list_on_component(componentID, boundaryFaceList);
  }

  // n has to be a vector because of the call to
  // computeNormalAndTransformationData below
  std::vector<double> n(3);
  std::vector<double> u(3);

  for (size_t m = 0; m < boundaryFaceList.size(); m++)
  {
    TBoundFace *boundface = boundaryFaceList[m];
    TBaseCell *cell = ((TJoint *)boundface)->GetNeighbour(0);

    int icell = cell->GetCellIndex();
    // mapping from local (cell) DOF to global DOF
    const int *DOF = U_Space->GetGlobalDOF(icell);

    int joint_id = boundface->get_index_in_neighbour(cell);

    // get all data necessary for computing the integral:
    // quadrature weights, points, functions values, normal, determinant

    std::vector<double> qWeights, qPointsT, qPointsS;
    std::vector<std::vector<double>> basisFunctionsValues;

    this->getQuadratureData(U_Space, cell, joint_id, qWeights,
          qPointsT, qPointsS, basisFunctionsValues);

    double transformationDeterminant;
    cell->computeNormalAndTransformationData(joint_id, n,
               transformationDeterminant);

    double x, y, z;

    boundface->GetXYZofTS(qPointsT[0], qPointsS[0], x, y, z);
    boundface->get_normal_vector(x, y, z, n[0], n[1], n[2]);

    // loop over Gauss points
    for (size_t l = 0; l < qWeights.size(); l++)
    {
      boundface->GetXYZofTS(qPointsT[l], qPointsS[l], x, y, z);

      velocity.FindValueLocal(cell, icell, x, y, z, u.data());

      double u_dot_n = u[0] * n[0] + u[1] * n[1] + u[2] * n[2];
      double u_dot_u = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
      double Theta = 0.0;

      if (delta > 0.0)
      {
        Theta = 0.5 * (1.0 - std::tanh(u_dot_n / delta));
      }
      else if (u_dot_n < 0.0)
      {
        Theta = 1.0;
      }

      double u_term_mul = 0.5 * beta * Theta * u_dot_n;
      double n_term_mul = 0.5 * beta * Theta * u_dot_u;

      double w = qWeights[l] * transformationDeterminant;

      for (size_t k = 0; k < basisFunctionsValues[l].size(); k++)
      {
        int global_dof_from_local = DOF[k];

        // if the DOF is Dirichlet, continue
        if (global_dof_from_local < U_Space->get_n_active())
        {
          double v = w * basisFunctionsValues[l][k];
          double vu = u_term_mul * v;
          double vn = n_term_mul * v;

          rhs.block(0)[global_dof_from_local] += vu * u[0] + vn * n[0];
          rhs.block(1)[global_dof_from_local] += vu * u[1] + vn * n[1];
          rhs.block(2)[global_dof_from_local] += vu * u[2] + vn * n[2];
        }
      }
    }
  }
}

//=================================================================================
// int_{Gamma} mult * given_boundary_data(x,y,z) * < v, normal >
void BoundaryAssembling3D::rhs_g_v_n(BlockVector &rhs,
    const TFESpace3D *U_Space, BoundValueFunct3D *given_boundary_data,
    std::vector<TBoundFace*>& boundaryFaceList, int componentID, double mult)
{
  if (boundaryFaceList.empty())
  {
    auto coll = U_Space->GetCollection();
    coll->get_face_list_on_component(componentID, boundaryFaceList);
  }

  for (size_t m = 0; m < boundaryFaceList.size(); m++)
  {
    TBoundFace *boundface = boundaryFaceList[m];
    const TBoundComp* BoundComp = boundface->GetBoundComp();
    int comp = BoundComp->get_physical_id();
    TBaseCell *cell = ((TJoint *)boundface)->GetNeighbour(0);

    // mapping from local (cell) DOF to global DOF
    const int *DOF = U_Space->GetGlobalDOF(cell->GetCellIndex());

    int joint_id = boundface->get_index_in_neighbour(cell);

    // get all data necessary for computing the integral:
    // quadrature weights, points, functions values, normal, determinant

    std::vector<double> qWeights, qPointsT, qPointsS;
    std::vector<std::vector<double>> basisFunctionsValues;

    this->getQuadratureData(U_Space, cell, joint_id, qWeights,
			    qPointsT, qPointsS, basisFunctionsValues);

    std::vector<double> n;

    double transformationDeterminant;
    cell->computeNormalAndTransformationData(joint_id, n,
					     transformationDeterminant);

    double x, y, z;
    double value;

    n.resize(3);
    boundface->GetXYZofTS(qPointsT[0], qPointsS[0], x, y, z);
    boundface->get_normal_vector(x, y, z, n[0], n[1], n[2]);

    // loop over Gauss points
    for (size_t l = 0; l < qWeights.size(); l++)
    {
      if (given_boundary_data != nullptr)
      {
        boundface->GetBoundComp()->GetXYZofTS(qPointsT[l], qPointsS[l], x, y, z);
        given_boundary_data(comp, x, y, z, value);
      }
      else
      {
        value = 1.0;
      }

      double scale_factor = mult * value
        * qWeights[l] * transformationDeterminant;

      for (size_t k = 0; k < basisFunctionsValues[l].size(); k++)
      {
        int global_dof_from_local = DOF[k];

        // if the DOF is Dirichlet, continue
        if (global_dof_from_local < U_Space->get_n_active())
        {
          double v1 = basisFunctionsValues[l][k]; // value of test function (vtest = vx = vy =vz)
          double v2 = v1;
          double v3 = v1;

          rhs.block(0)[global_dof_from_local] += scale_factor * v1 * n[0];
          rhs.block(1)[global_dof_from_local] += scale_factor * v2 * n[1];
          rhs.block(2)[global_dof_from_local] += scale_factor * v3 * n[2];
        }
      }
    }
  }
}

//=================================================================================
void BoundaryAssembling3D::rhs_g_v_n(BlockVector &rhs,
				     const TFESpace3D *U_Space,
				     BoundValueFunct3D *given_boundary_data,
				     std::vector<TBaseCell*> &,
				     int componentID,
				     double mult)
{  
  // now we loop always over the whole collection: inefficient!
  auto coll = U_Space->GetCollection();

  for(int i = 0; i < coll->GetN_Cells(); i++)
  {
    TBaseCell* cell = coll->GetCell(i); //boundaryCells[i];

    // mapping from local (cell) DOF to global DOF
    const int *DOF = U_Space->GetGlobalDOF(cell->GetCellIndex());
    const auto& bf_u = *(U_Space->get_fe(cell->GetCellIndex()).GetBaseFunct());
    
    for(size_t joint_id = 0; joint_id < (size_t) cell->GetN_Faces(); joint_id++)
    {
      TJoint* joint = cell->GetJoint(joint_id);

      if (joint->GetType() == BoundaryFace || joint->GetType() == IsoBoundFace) {

        // convert the joint to an object of BoundFace type
        TBoundFace *boundface = (TBoundFace *)joint;
        const TBoundComp* BoundComp = boundface->GetBoundComp();
        int comp  = BoundComp->get_physical_id();

        if (boundface->GetBoundComp()->get_physical_id() == componentID)
        {
	  // --------------------------------------------------------------
          // get all data necessary for computing the integral:
          // quadrature weights, points, functions values, normal, determinant
          std:: vector<double> qWeights, qPointsT, qPointsS;
          std::vector< std::vector<double> > basisFunctionsValues;
          getFaceQuadratureData(*coll, bf_u, cell, joint_id, qWeights, qPointsT,
                                qPointsS, basisFunctionsValues);
	  // --------------------------------------------------------------
	  

	  // --------------------------------------------------------------
	  // get normal and face area
          std::vector<double> normal;
          double transformationDeterminant;
          cell->computeNormalAndTransformationData(joint_id, normal,
              transformationDeterminant);

          double x, y, z;
          double value;
          normal.clear();
          normal.resize(3);
          boundface->GetXYZofTS(qPointsT[0], qPointsS[0], x, y, z);
          boundface->get_normal_vector(x, y, z, normal[0], normal[1], normal[2]);
	  // --------------------------------------------------------------


	  // --------------------------------------------------------------
	  // ***** COMPUTE INTEGRAL *****
          // loop over Gauss points
          for (size_t l = 0; l < qWeights.size(); l++)
          {
            
            if(given_boundary_data != nullptr)
	    {
	      boundface->GetBoundComp()->GetXYZofTS(qPointsT[l], qPointsS[l], x, y, z);
              given_boundary_data(comp, x, y, z, value);
	    }
            else
              value = 1.;

            double scale_factor = mult * qWeights[l] * transformationDeterminant;

            for (size_t k = 0; k < basisFunctionsValues[l].size(); k++)
            {
              int global_dof_from_local = DOF[k];

              if (global_dof_from_local < U_Space->get_n_active())
              {
                double v1 = basisFunctionsValues[l][k]; // value of test function (vtest = vx = vy =vz)
                double v2 = v1;
                double v3 = v1;
		rhs.block(0)[global_dof_from_local] += scale_factor * value * (v1*normal[0]);
                rhs.block(1)[global_dof_from_local] += scale_factor * value * (v2*normal[1]);
                rhs.block(2)[global_dof_from_local] += scale_factor * value * (v3*normal[2]);
              }
            }
          }
        }
      }
    }
  }
}


// ===========================================================================
// int_{Gamma} mult * < u, v >
void BoundaryAssembling3D::matrix_u_v(BlockFEMatrix &M,
				      const TFESpace3D *U_Space,
				      std::vector<TBoundFace*>& boundaryFaceList,
				      int componentID,
				      double mult,
				      bool rescale_by_h)
{
  auto coll = U_Space->GetCollection();
  if (boundaryFaceList.empty())
  {
    coll->get_face_list_on_component(componentID, boundaryFaceList);
  }

  // loop over the faces (joints)
  for (size_t m = 0; m < boundaryFaceList.size(); m++)
  {
    TBoundFace *boundface = boundaryFaceList[m];
    TBaseCell *cell = ( (TJoint *)boundface)->GetNeighbour(0);

    // mapping from local (cell) DOF to global DOF
    const int *DOF = U_Space->GetGlobalDOF(cell->GetCellIndex());
    const auto& bf_u = *(U_Space->get_fe(cell->GetCellIndex()).GetBaseFunct());

    int joint_id = boundface->get_index_in_neighbour(cell);

    // --------------------------------------------------------------
    // get all data necessary for computing the integral:
    // quadrature weights, points, functions values, normal, determinant
    std:: vector<double> qWeights, qPointsT, qPointsS;
    std::vector< std::vector<double> > basisFunctionsValues;
    getFaceQuadratureData(*coll, bf_u, cell, joint_id, qWeights, qPointsT,
                          qPointsS, basisFunctionsValues);
    // --------------------------------------------------------------

    // --------------------------------------------------------------
    // get normal and face area
    std::vector<double> normal;
    double transformationDeterminant;
    cell->computeNormalAndTransformationData(joint_id, normal,
					     transformationDeterminant);

    // --------------------------------------------------------------

    std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();

    // rescale local integral by mesh-size (important for Nitsche boundary)
    double h = 1.;
    if (rescale_by_h)
      h = cell->Get_hK(0);
    
    for (size_t l = 0; l < qWeights.size(); l++)
    {
      double scale_factor = mult * qWeights[l] * transformationDeterminant / h;

      // test 
      for (size_t k1 = 0; k1 < basisFunctionsValues[l].size(); k1++)
      {
        if ( DOF[k1] < U_Space->get_n_active())
        {
          double v1 = basisFunctionsValues[l][k1]; 
          double v2 = v1;
          double v3 = v1;

          // ansatz
          for (size_t k2 = 0; k2 < basisFunctionsValues[l].size(); k2++)
          {
            double u1 = basisFunctionsValues[l][k2];
            double u2 = u1;
            double u3 = u2;

            // add for all three components
            blocks[0]->add(DOF[k1], DOF[k2], scale_factor * u1 * v1 ); // A11
            blocks[5]->add(DOF[k1], DOF[k2], scale_factor * u2 * v2 ); // A22
            blocks[10]->add(DOF[k1], DOF[k2], scale_factor * u3 * v3 ); // A33
          }
        }
      }
    }
  }
}


// ===========================================================================
// int_{Gamma} p v n
void BoundaryAssembling3D::matrix_p_v_n(BlockFEMatrix &M,
					const TFESpace3D *U_Space,
					const TFESpace3D *P_Space,
					std::vector<TBoundFace*>& boundaryFaceList,
					int componentID,
					double mult)
{
  auto coll = U_Space->GetCollection();
  if (boundaryFaceList.empty())
  {
    coll->get_face_list_on_component(componentID, boundaryFaceList);
  }

  for(size_t m = 0; m < boundaryFaceList.size(); m++)
  {
    TBoundFace *boundface = boundaryFaceList[m];
    TBaseCell *cell = ( (TJoint *)boundface)->GetNeighbour(0);

    // mapping from local (cell) DOF to global DOF
    const int *DOF_u = U_Space->GetGlobalDOF(cell->GetCellIndex());
    const int *DOF_p = P_Space->GetGlobalDOF(cell->GetCellIndex());
    const auto& bf_u = *(U_Space->get_fe(cell->GetCellIndex()).GetBaseFunct());

    int joint_id = boundface->get_index_in_neighbour(cell);
    // --------------------------------------------------------------
    // get all data necessary for computing the integral:
    // quadrature weights, points, functions values, normal, determinant
    std:: vector<double> qWeights, qPointsT, qPointsS;
    std::vector< std::vector<double> > basisFunctionsValues_u, basisFunctionsValues_p;

    int quad_degree = U_Space->getFEDegree(cell) * P_Space->getFEDegree(cell);
    getFaceQuadratureData(*coll, bf_u, cell, joint_id, qWeights, qPointsT,
                          qPointsS, basisFunctionsValues_u, quad_degree);

    // we need to get the formula in order to use the same quadrature points for P
    // (in the case that fe spaces are different)
    const TQuadFormula * faceQuadFormula = getFaceQuadratureFormula(
        cell, joint_id, quad_degree);
    
    const auto& bf_p = *(P_Space->get_fe(cell->GetCellIndex()).GetBaseFunct());
    getFaceQuadratureValue(*coll, bf_p, cell, joint_id, faceQuadFormula,
                           basisFunctionsValues_p);
    // --------------------------------------------------------------

    std::vector<double> normal;
    double transformationDeterminant;
    cell->computeNormalAndTransformationData(joint_id, normal,
              transformationDeterminant);

    double x, y, z;
    normal.clear();
    normal.resize(3);
    boundface->GetXYZofTS(qPointsT[0], qPointsS[0], x, y, z);
    boundface->get_normal_vector(x, y, z, normal[0], normal[1], normal[2]);
    double n1 = normal[0];
    double n2 = normal[1];
    double n3 = normal[2];

    std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();

    for (size_t l = 0; l < qWeights.size(); l++)
    {
      double scale_factor  = mult * qWeights[l] * transformationDeterminant;

      // test
      for (size_t k1 = 0; k1 < basisFunctionsValues_u[l].size(); k1++)
      {
        if (DOF_u[k1] < U_Space->get_n_active())
        {
          double v1 = basisFunctionsValues_u[l][k1]; 
          double v2 = v1;
          double v3 = v1;

          // ansatz
          for (size_t k2 = 0; k2 < basisFunctionsValues_p[l].size(); k2++)
          {
            double p = basisFunctionsValues_p[l][k2]; // value of ansatz function

            // add for all three components
            blocks[3]->add(DOF_u[k1],  DOF_p[k2], scale_factor * p * v1 * n1 ); // B1T
            blocks[7]->add(DOF_u[k1],  DOF_p[k2], scale_factor * p * v2 * n2 ); // B2T
            blocks[11]->add(DOF_u[k1], DOF_p[k2], scale_factor * p * v3 * n3 ); // B3T
          }
        }
      }
    }
  }
}


// TO CHECK 
// ===========================================================================
// int_{Gamma} < q, u.n >
void BoundaryAssembling3D::matrix_q_u_n(BlockFEMatrix &M,
    const TFESpace3D *U_Space,
    const TFESpace3D *P_Space,
    std::vector<TBoundFace*>& boundaryFaceList,
    int componentID,
    double mult)
{
  std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();
  /**
   * @todo: check if the matrix structure is correct:
   * we need 6 square matrices with the same FE spaces
   */

  auto coll = U_Space->GetCollection();
  if (boundaryFaceList.empty())
  {
    coll->get_face_list_on_component(componentID, boundaryFaceList);
  }

  for(size_t m = 0; m < boundaryFaceList.size(); m++)
  {
    TBoundFace *boundface = boundaryFaceList[m];
    TBaseCell *cell = ( (TJoint *)boundface)->GetNeighbour(0);

    // mapping from local (cell) DOF to global DOF
    const int *DOF_u = U_Space->GetGlobalDOF(cell->GetCellIndex());
    // mapping from local (cell) DOF to global DOF
    const int *DOF_p = P_Space->GetGlobalDOF(cell->GetCellIndex());
    
    const auto& bf_u = *(U_Space->get_fe(cell->GetCellIndex()).GetBaseFunct());

    int joint_id = boundface->get_index_in_neighbour(cell);

    // get all data necessary for computing the integral:
    // quadrature weights, points, functions values, normal, determinant
    std:: vector<double> qWeights_u, qPointsT_u, qPointsS_u, qWeights_p, qPointsT_p, qPointsS_p; // S,T for parametrization in 3D
    std::vector< std::vector<double> > basisFunctionsValues_u, basisFunctionsValues_p;

    int quad_degree = U_Space->getFEDegree(cell) * P_Space->getFEDegree(cell);

    getFaceQuadratureData(*coll, bf_u, cell, joint_id, qWeights_u, qPointsT_u,
                          qPointsS_u, basisFunctionsValues_u, quad_degree);
    // we need to get the formula in order to use the same quadrature points for P
    // (in the case that fe spaces are different)
    const TQuadFormula * faceQuadFormula = getFaceQuadratureFormula(
        cell, joint_id, quad_degree);

    const auto& bf_p = *(P_Space->get_fe(cell->GetCellIndex()).GetBaseFunct());
    getFaceQuadratureValue(*coll, bf_p, cell, joint_id, faceQuadFormula,
                           basisFunctionsValues_p);



  /* OLD 15.10.18 LB:
    this->getQuadratureData(U_Space, cell,joint_id,
        qWeights_u, qPointsT_u, qPointsS_u, basisFunctionsValues_u);
    this->getQuadratureData(P_Space, cell, joint_id,
        qWeights_p, qPointsT_p, qPointsS_p, basisFunctionsValues_p);
*/

    std::vector<double> normal;
    double transformationDeterminant;
    this->computeNormalAndTransformationData(cell, joint_id, normal,
        transformationDeterminant);

    double x, y, z;
    normal.clear();
    normal.resize(3);
    boundface->GetXYZofTS(qPointsT_u[0], qPointsS_u[0], x, y, z);
    boundface->get_normal_vector(x, y, z, normal[0], normal[1], normal[2]);

    // rescale local integral by mesh-size (important for Nitsche boundary)
    for (size_t l = 0; l < qWeights_u.size(); l++)
    {

      /*cout << "l: "<<l<<endl;
      cout << "basisFunctionsValues_u[l].size(): "<< basisFunctionsValues_u[l].size() <<endl;
      cout << "basisFunctionsValues_p[l].size(): "<< basisFunctionsValues_p[l].size() <<endl;*/


      double scale_factor  = mult * qWeights_u[l] * transformationDeterminant;

      for (size_t k1 = 0; k1 < basisFunctionsValues_u[l].size(); k1++)
      {
        int global_dof_from_local_u = DOF_u[k1]; // ansatz-DOF
        if (global_dof_from_local_u < U_Space->get_n_active())
        {
          double u1 = basisFunctionsValues_u[l][k1]; // value of ansatz function
          double u2 = u1; // value of ansatz function
          double u3 = u1; // value of ansatz function

          for (size_t k2 = 0; k2 < basisFunctionsValues_p[l].size(); k2++)
          {
            int global_dof_from_local_p = DOF_p[k2]; // test-DOF

            if (global_dof_from_local_p < P_Space->get_n_active())
            {
              double q = basisFunctionsValues_p[l][k2]; // value of test function (vtest = vx = vy =vz)

              double  n1 = normal[0];
              double  n2 = normal[1];
              double  n3 = normal[2];

              // add for all three components
              blocks[12]->add(DOF_p[k2], DOF_u[k1], scale_factor * q * u1 * n1 ); // B1  //NEW 14.10.18
              blocks[13]->add(DOF_p[k2], DOF_u[k1], scale_factor * q * u2 * n2 ); // B2
              blocks[14]->add(DOF_p[k2], DOF_u[k1], scale_factor * q * u3 * n3 ); // B3
            }
          }
        }
      }
    }
  }
}


// ===========================================================================
// int_{Gamma}  < q, u.n >
void BoundaryAssembling3D::rhs_q_uD_n(BlockVector &rhs,
    const TFESpace3D *U_Space,
    const TFESpace3D *P_Space,
    BoundValueFunct3D *given_boundary_data0,
    BoundValueFunct3D *given_boundary_data1,
    BoundValueFunct3D *given_boundary_data2,
				      TFEVectFunct3D *given_fe_fct,
    std::vector<TBoundFace*>& boundaryFaceList,
    int componentID,
    double mult)
{
  /**
   * @todo: check if the matrix structure is correct:
   * we need 4 square matrices with the same FE spaces
   */

  if (boundaryFaceList.empty())
  {
    auto coll = U_Space->GetCollection();
    coll->get_face_list_on_component(componentID, boundaryFaceList);
  }

  for(size_t m = 0; m < boundaryFaceList.size(); m++)
  {
    TBoundFace *boundface = boundaryFaceList[m];
    const TBoundComp* BoundComp = boundface->GetBoundComp();
    int comp  = BoundComp->get_physical_id();
    TBaseCell *cell = ( (TJoint *)boundface)->GetNeighbour(0);

    // mapping from local (cell) DOF to global DOF
    const int *DOF_u = U_Space->GetGlobalDOF(cell->GetCellIndex());
    // mapping from local (cell) DOF to global DOF
    const int *DOF_p = P_Space->GetGlobalDOF(cell->GetCellIndex());

    int joint_id = boundface->get_index_in_neighbour(cell);

    // get all data necessary for computing the integral:
    // quadrature weights, points, functions values, normal, determinant
    std:: vector<double> qWeights_u, qPointsT_u, qPointsS_u, qWeights_p,qPointsT_p,qPointsS_p; // S,T for parametrization in 3D
    std::vector< std::vector<double> > basisFunctionsValues_u, basisFunctionsValues_p;
    this->getQuadratureData(U_Space, cell, joint_id,
        qWeights_u, qPointsT_u, qPointsS_u, basisFunctionsValues_u);
    this->getQuadratureData(P_Space, cell, joint_id,
        qWeights_p, qPointsT_p, qPointsS_p, basisFunctionsValues_p);

    std::vector<double> normal;
    double transformationDeterminant;
    this->computeNormalAndTransformationData(cell, joint_id, normal,
        transformationDeterminant);

    double x, y, z;
    normal.clear();
    normal.resize(3);
    boundface->GetXYZofTS(qPointsT_u[0], qPointsS_u[0], x, y, z);
    boundface->get_normal_vector(x, y, z, normal[0], normal[1], normal[2]);


    std::vector<double> uDirichlet(3);

    // rescale local integral by mesh-size (important for Nitsche boundary)
    for (size_t l = 0; l < qWeights_u.size(); l++)
    {
      boundface->GetBoundComp()->GetXYZofTS(qPointsT_u[l], qPointsS_u[l], x, y, z);

      // get the boundary values of rhs
      if (given_fe_fct == nullptr)
      { 
	if(given_boundary_data0 != nullptr)  given_boundary_data0(comp, x, y, z, uDirichlet[0]);
	else uDirichlet[0] = 0.; // 1.; //
	
	if(given_boundary_data1 != nullptr) given_boundary_data1(comp, x, y, z, uDirichlet[1]);
	else uDirichlet[1] = 0.; // 1.; //
	
	if(given_boundary_data2 != nullptr) given_boundary_data2(comp, x, y, z, uDirichlet[2]);
	else uDirichlet[2] = 0.; //1.; //
	
      } else {

  for(size_t icoor=0; icoor<3; icoor++)
  {
    uDirichlet[icoor] = 0.;
    for (size_t k = 0; k < basisFunctionsValues_p[l].size(); k++)
    {
      double v1 = basisFunctionsValues_u[l][k];
      int global_dof_from_local = DOF_u[k];
      auto given_fe_fct_component = given_fe_fct->GetComponent(icoor);
      double *u_icoor_values = given_fe_fct_component->GetValues();
      uDirichlet[icoor] += u_icoor_values[ global_dof_from_local ] * v1;
    }
  }
      }
      double scale_factor  = mult * qWeights_u[l] * transformationDeterminant;

      for (size_t k1 = 0; k1 < basisFunctionsValues_p[l].size(); k1++)
      {
        int global_dof_from_local_p = DOF_p[k1]; // Test-DOF

        if (global_dof_from_local_p < P_Space->get_n_active())
        {
          double q = basisFunctionsValues_p[l][k1]; // value of test function (vtest = vx = vy =vz)

          double  n1 = normal[0];
          double  n2 = normal[1];
          double  n3 = normal[2];

          // add for all three components
          rhs.block(3)[global_dof_from_local_p] += scale_factor * q * uDirichlet[0] * n1;
          rhs.block(3)[global_dof_from_local_p] += scale_factor * q * uDirichlet[1] * n2;
          rhs.block(3)[global_dof_from_local_p] += scale_factor * q * uDirichlet[2] * n3;
          //todo: check if scaling with transformation size is necessary (see Boundary_Assembling_2D.C)
        }
      }
    }
  }
}



// ===========================================================================
// int_{Gamma} mult * ( (grad u).n, v )
void BoundaryAssembling3D::matrix_gradu_n_v(BlockFEMatrix &M,
    const TFESpace3D *U_Space,
    std::vector<TBoundFace*>& boundaryFaceList,
    int componentID,
    double mult)
{
  std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();

  if (boundaryFaceList.empty())
  {
    auto coll = U_Space->GetCollection();
    coll->get_face_list_on_component(componentID, boundaryFaceList);
  }

  for (size_t m = 0; m < boundaryFaceList.size(); m++)
  {
    TBoundFace *boundface = boundaryFaceList[m];
    TBaseCell *cell = ( (TJoint *)boundface)->GetNeighbour(0);

    // mapping from local (cell) DOF to global DOF
    const int *DOF = U_Space->GetGlobalDOF(cell->GetCellIndex());

    int joint_id = boundface->get_index_in_neighbour(cell);

    // get all data necessary for computing the integral:
    // quadrature weights, points, functions values, normal, determinant
    std:: vector<double> qWeights, qPointsT, qPointsS; // S,T for parametrization in 3D
    std::vector< std::vector<double> > basisFunctionsValues, basisFunctionsValues_derivative_x, basisFunctionsValues_derivative_y, basisFunctionsValues_derivative_z;
    this->getQuadratureDataIncludingFirstDerivatives(U_Space, cell, joint_id,
        qWeights, qPointsT, qPointsS, basisFunctionsValues,
        basisFunctionsValues_derivative_x, basisFunctionsValues_derivative_y, basisFunctionsValues_derivative_z );

    std::vector<double> normal;
    double transformationDeterminant;

    this->computeNormalAndTransformationData(cell, joint_id, normal,
        transformationDeterminant);
    //cout << "transformationDeterminant: "<< transformationDeterminant << endl;

    double x, y, z;
    normal.clear();
    normal.resize(3);
    boundface->GetXYZofTS(qPointsT[0], qPointsS[0], x, y, z);
    boundface->get_normal_vector(x, y, z, normal[0], normal[1], normal[2]);

    /*cout<< "Nitsche Joint_ID: "<< joint_id  << endl;
    cout<< "cell->GetCellIndex(): "<< cell->GetCellIndex() << endl;
    cout << "n1, n2, n3: "<< normal[0] << ", "<< normal[1] << ", " << normal[2] <<endl;*/

    // rescale local integral by mesh-size (important for Nitsche boundary)
    for (size_t l = 0; l < qWeights.size(); l++)
    {
      double scale_factor = mult * qWeights[l] * transformationDeterminant;

      for (size_t k1 = 0; k1 < basisFunctionsValues[l].size(); k1++)
      {
        int global_dof_from_local = DOF[k1]; // Test-DOF

        if (global_dof_from_local < U_Space->get_n_active())
        {
          double v1 = basisFunctionsValues[l][k1]; // value of test function (vtest = vx = vy =vz)
          double v2 = v1;
          double v3 = v1;

          /*cout << "l: " << l << "k1: " << k1 << endl;
	    cout << "basisFunctionsValues[l][k1]: " << basisFunctionsValues[l][k1] << endl;*/

          for (size_t k2 = 0; k2 < basisFunctionsValues[l].size(); k2++)
          {
            int global_dof_from_local = DOF[k2]; // Ansatz-DOF

            if (global_dof_from_local < U_Space->get_n_active())
            {
              double  n1 = normal[0];
              double  n2 = normal[1];
              double  n3 = normal[2];

              double u_dx = basisFunctionsValues_derivative_x[l][k2];
              double u_dy = basisFunctionsValues_derivative_y[l][k2];
              double u_dz = basisFunctionsValues_derivative_z[l][k2];

              /*cout << "l: " << l << "k2: " << k2 << endl;
              cout << "basisFunctionsValues_derivative_x[l][k2]: " << basisFunctionsValues_derivative_x[l][k2] << endl;
              cout << "basisFunctionsValues_derivative_y[l][k2]: " << basisFunctionsValues_derivative_y[l][k2] << endl;
              cout << "basisFunctionsValues_derivative_z[l][k2]: " << basisFunctionsValues_derivative_z[l][k2] << endl;*/

              // (see the note about blocks at the beginning of the function)
              blocks[0]->add(DOF[k1], DOF[k2],  scale_factor * v1 * u_dx * n1 ); // A11
              blocks[0]->add(DOF[k1], DOF[k2],  scale_factor * v1 * u_dy * n2 ); // A11
              blocks[0]->add(DOF[k1], DOF[k2],  scale_factor * v1 * u_dz * n3 ); // A11

              blocks[5]->add(DOF[k1], DOF[k2],  scale_factor * v2 * u_dx * n1 ); // A22
              blocks[5]->add(DOF[k1], DOF[k2],  scale_factor * v2 * u_dy * n2 ); // A22
              blocks[5]->add(DOF[k1], DOF[k2],  scale_factor * v2 * u_dz * n3 ); // A22

              blocks[10]->add(DOF[k1], DOF[k2], scale_factor * v3 * u_dx * n1 ); // A33
              blocks[10]->add(DOF[k1], DOF[k2], scale_factor * v3 * u_dy * n2 ); // A33
              blocks[10]->add(DOF[k1], DOF[k2], scale_factor * v3 * u_dz * n3 ); // A33
            }
          }
        }
      }
    }
  }
}

// ===========================================================================
// int_{Gamma} mult * < (grad v).n, u >
void BoundaryAssembling3D::matrix_gradv_n_u(BlockFEMatrix &M,
    const TFESpace3D *U_Space,
    std::vector<TBoundFace*>& boundaryFaceList,
    int componentID,
    double mult)

{
  std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();

  if (boundaryFaceList.empty())
  {
    auto coll = U_Space->GetCollection();
    coll->get_face_list_on_component(componentID, boundaryFaceList);
  }

  // loop over the faces (joints)
  for(size_t m = 0; m < boundaryFaceList.size(); m++)
  {
    TBoundFace *boundface = boundaryFaceList[m];
    TBaseCell *cell = ( (TJoint *)boundface)->GetNeighbour(0);

    // mapping from local (cell) DOF to global DOF
    const int *DOF = U_Space->GetGlobalDOF(cell->GetCellIndex());

    int joint_id = boundface->get_index_in_neighbour(cell);

    // get all data necessary for computing the integral:
    // quadrature weights, points, functions values, normal, determinant
    std:: vector<double> qWeights, qPointsT, qPointsS; // S,T for parametrization in 3D
    std::vector< std::vector<double> > basisFunctionsValues, basisFunctionsValues_derivative_x, basisFunctionsValues_derivative_y, basisFunctionsValues_derivative_z ;
    this->getQuadratureDataIncludingFirstDerivatives(U_Space, cell,joint_id,
        qWeights, qPointsT, qPointsS,
        basisFunctionsValues,
        basisFunctionsValues_derivative_x,
        basisFunctionsValues_derivative_y,
        basisFunctionsValues_derivative_z );
    std::vector<double> normal;
    double transformationDeterminant;
    this->computeNormalAndTransformationData(cell, joint_id, normal,
        transformationDeterminant);

    double x, y, z;
    normal.clear();
    normal.resize(3);
    boundface->GetXYZofTS(qPointsT[0], qPointsS[0], x, y, z);
    boundface->get_normal_vector(x, y, z, normal[0], normal[1], normal[2]);

    // rescale local integral by mesh-size (important for Nitsche boundary)
    for (size_t l = 0; l < qWeights.size(); l++)
    {
      double scale_factor = mult * qWeights[l] * transformationDeterminant;

      for (size_t k1 = 0; k1 < basisFunctionsValues[l].size(); k1++)
      {
        int global_dof_from_local = DOF[k1]; // Test-DOF

        if (global_dof_from_local < U_Space->get_n_active())
        {
          double v_dx = basisFunctionsValues_derivative_x[l][k1];
          double v_dy = basisFunctionsValues_derivative_y[l][k1];
          double v_dz = basisFunctionsValues_derivative_z[l][k1];

          for (size_t k2 = 0; k2 < basisFunctionsValues[l].size(); k2++)
          {
            int global_dof_from_local = DOF[k2]; // Ansatz-DOF

            if(global_dof_from_local < U_Space->get_n_active())
            {
              double u1 = basisFunctionsValues[l][k2]; // value of test function (vtest = vx = vy =vz)
              double u2 = u1;
              double u3 = u1;

              double  n1 = normal[0];
              double  n2 = normal[1];
              double  n3 = normal[2];

              // (see the note about blocks at the beginning of the function)
              blocks[0]->add(DOF[k1], DOF[k2], scale_factor * u1 * v_dx * n1 ); // A11
              blocks[0]->add(DOF[k1], DOF[k2], scale_factor * u1 * v_dy * n2 ); // A11
              blocks[0]->add(DOF[k1], DOF[k2], scale_factor * u1 * v_dz * n3 ); // A11

              blocks[5]->add(DOF[k1], DOF[k2], scale_factor * u2 * v_dx * n1 ); // A22
              blocks[5]->add(DOF[k1], DOF[k2], scale_factor * u2 * v_dy * n2 ); // A22
              blocks[5]->add(DOF[k1], DOF[k2], scale_factor * u2 * v_dz * n3 ); // A22

              blocks[10]->add(DOF[k1], DOF[k2], scale_factor * u3 * v_dx * n1 ); // A33
              blocks[10]->add(DOF[k1], DOF[k2], scale_factor * u3 * v_dy * n2 ); // A33
              blocks[10]->add(DOF[k1], DOF[k2], scale_factor * u3 * v_dz * n3 ); // A33
            }
          }
        }
      }
    }
  }
}


// ===========================================================================
// int_{Gamma} mult * < (grad v).n, uD >
void BoundaryAssembling3D::rhs_gradv_n_uD(BlockVector &rhs,
					  const TFESpace3D *U_Space,
					  BoundValueFunct3D *given_boundary_data0,
					  BoundValueFunct3D *given_boundary_data1,
					  BoundValueFunct3D *given_boundary_data2,
					  TFEVectFunct3D *given_fe_fct,
					  std::vector<TBoundFace*>& boundaryFaceList,
					  int componentID,
					  double mult)
{
  if (boundaryFaceList.empty())
  {
    auto coll = U_Space->GetCollection();
    coll->get_face_list_on_component(componentID, boundaryFaceList);
  }

  // loop over the faces (joints)
  for(size_t m = 0; m < boundaryFaceList.size(); m++)
  {
    TBoundFace *boundface = boundaryFaceList[m];
    const TBoundComp* BoundComp = boundface->GetBoundComp();
    int comp  = BoundComp->get_physical_id();
    TBaseCell *cell = ( (TJoint *)boundface)->GetNeighbour(0);

    // mapping from local (cell) DOF to global DOF
    const int *DOF = U_Space->GetGlobalDOF(cell->GetCellIndex());

    int joint_id = boundface->get_index_in_neighbour(cell);

  //  cout << "joint_id from get index in neighbour: "<< joint_id << endl;

    // get all data necessary for computing the integral:
    // quadrature weights, points, functions values, normal, determinant
    std:: vector<double> qWeights, qPointsT, qPointsS; // S,T for parametrization in 3D
    std::vector< std::vector<double> > basisFunctionsValues, basisFunctionsValues_derivative_x, basisFunctionsValues_derivative_y, basisFunctionsValues_derivative_z ;
    this->getQuadratureDataIncludingFirstDerivatives(U_Space, cell, joint_id,
        qWeights, qPointsT, qPointsS,
        basisFunctionsValues, basisFunctionsValues_derivative_x, basisFunctionsValues_derivative_y, basisFunctionsValues_derivative_z );
    std::vector<double> normal;
    double transformationDeterminant;
    this->computeNormalAndTransformationData(cell, joint_id, normal,
        transformationDeterminant);

    double x, y, z;
    normal.clear();
    normal.resize(3);
    boundface->GetXYZofTS(qPointsT[0], qPointsS[0], x, y, z);
    boundface->get_normal_vector(x, y, z, normal[0], normal[1], normal[2]);



    std::vector<double> uDirichlet(3);

    // loop over Gauss points
    for (size_t l = 0; l < qWeights.size(); l++)
    {
      // get the boundary values of rhs
      boundface->GetBoundComp()->GetXYZofTS(qPointsT[l], qPointsS[l], x, y, z);

      if (given_fe_fct == nullptr)
      {
	if(given_boundary_data0 != nullptr)  given_boundary_data0(comp, x, y, z, uDirichlet[0]);
	else  uDirichlet[0] = 0.; //1.; //
	
	if(given_boundary_data1 != nullptr)  given_boundary_data1(comp, x, y, z, uDirichlet[1]);
	else uDirichlet[1] = 0.; // 1.; //
	
	if(given_boundary_data2 != nullptr) given_boundary_data2(comp, x, y, z, uDirichlet[2]);
	else  uDirichlet[2] = 0.; //0.; //
      } else {

  for(size_t icoor=0; icoor<3; icoor++)
  {
    uDirichlet[icoor] = 0.;
    for (size_t k = 0; k < basisFunctionsValues[l].size(); k++)
    {
      double v1 = basisFunctionsValues[l][k];
      int global_dof_from_local = DOF[k];
      auto given_fe_fct_component = given_fe_fct->GetComponent(icoor);
      double *u_icoor_values = given_fe_fct_component->GetValues();
      uDirichlet[icoor] += u_icoor_values[ global_dof_from_local ] * v1;
    }
  }

      }
      // rescale local integral by mesh-size (important for Nitsche boundary)
      double scale_factor = mult * qWeights[l] * transformationDeterminant;

      for (size_t k1 = 0; k1 < basisFunctionsValues[l].size(); k1++)
      {
        int global_dof_from_local = DOF[k1]; // Test-DOF

        if (global_dof_from_local < U_Space->get_n_active())
        {
          double v_dx = basisFunctionsValues_derivative_x[l][k1];
          double v_dy = basisFunctionsValues_derivative_y[l][k1];
          double v_dz = basisFunctionsValues_derivative_z[l][k1];

          double  n1 = normal[0];
          double  n2 = normal[1];
          double  n3 = normal[2];

          rhs.block(0)[global_dof_from_local] += scale_factor * uDirichlet[0] * (v_dx * n1 + v_dy * n2 + v_dz * n3);
          rhs.block(1)[global_dof_from_local] += scale_factor * uDirichlet[1] * (v_dx * n1 + v_dy * n2 + v_dz * n3);
          rhs.block(2)[global_dof_from_local] += scale_factor * uDirichlet[2] * (v_dx * n1 + v_dy * n2 + v_dz * n3);
        }
      }
    }
  }
}





// ===========================================================================

// int_{Gamma} mult * given_boundary_data(x,y,z) * v
void BoundaryAssembling3D::rhs_uD_v(BlockVector &rhs,
				    const TFESpace3D *U_Space,
				    BoundValueFunct3D *given_boundary_data0,
				    BoundValueFunct3D *given_boundary_data1,
				    BoundValueFunct3D *given_boundary_data2,
				    TFEVectFunct3D *given_fe_fct,
				    std::vector<TBoundFace*>& boundaryFaceList,
				    int componentID,
				    double mult,
				    bool rescale_by_h)
{
  if (boundaryFaceList.empty())
  {
    auto coll = U_Space->GetCollection();
    coll->get_face_list_on_component(componentID, boundaryFaceList);
  }

  for(size_t m = 0; m < boundaryFaceList.size(); m++)
  {
    TBoundFace *boundface = boundaryFaceList[m];
    const TBoundComp* BoundComp = boundface->GetBoundComp();
    int comp  = BoundComp->get_physical_id();
    TBaseCell *cell = ( (TJoint *)boundface)->GetNeighbour(0);

    // mapping from local (cell) DOF to global DOF
    const int *DOF = U_Space->GetGlobalDOF(cell->GetCellIndex());

    int joint_id = boundface->get_index_in_neighbour(cell);

    // get all data necessary for computing the integral:
    // quadrature weights, points, functions values, normal, determinant
    std:: vector<double> qWeights, qPointsT, qPointsS;
    std::vector< std::vector<double> > basisFunctionsValues;
    this->getQuadratureData(U_Space, cell, joint_id, qWeights,
        qPointsT, qPointsS, basisFunctionsValues);

    std::vector<double> normal;
    double transformationDeterminant;
    this->computeNormalAndTransformationData(cell, joint_id, normal,
        transformationDeterminant);

    double x, y, z;
    normal.clear();
    normal.resize(3);
    boundface->GetXYZofTS(qPointsT[0], qPointsS[0], x, y, z);
    boundface->get_normal_vector(x, y, z, normal[0], normal[1], normal[2]);

    std::vector<double> uDirichlet(3);

    // loop over Gauss points
    for (size_t l = 0; l < qWeights.size(); l++)
    {

      
      // get the boundary values of rhs
      boundface->GetXYZofTS(qPointsT[l], qPointsS[l], x, y, z);

      if (given_fe_fct == nullptr) {
	
	if (given_boundary_data0 != nullptr) given_boundary_data0(comp, x, y, z, uDirichlet[0]);
	else uDirichlet[0] = 0.;

	if(given_boundary_data1 != nullptr) given_boundary_data1(comp, x, y, z, uDirichlet[1]);
	else uDirichlet[1] = 0.;

	if(given_boundary_data2 != nullptr) given_boundary_data2(comp, x, y, z, uDirichlet[2]);
	else uDirichlet[2] = 0.; //1.; //
	
      } else {
	
  for(size_t icoor=0; icoor<3; icoor++)
  {
    uDirichlet[icoor] = 0.;

    for (size_t k = 0; k < basisFunctionsValues[l].size(); k++)
    {
      double v1 = basisFunctionsValues[l][k];
      int global_dof_from_local = DOF[k];
      auto given_fe_fct_component = given_fe_fct->GetComponent(icoor);
      double *u_icoor_values = given_fe_fct_component->GetValues();
      uDirichlet[icoor] += u_icoor_values[ global_dof_from_local ] * v1;
    }
  }
      }

      Output::print<4>(" Q=(", qPointsT[l], " ", qPointsS[l], ") X=(" , x, " ", y, " ",z, 
		    ") UD=( ", uDirichlet[0], " ", uDirichlet[1], " ", uDirichlet[2],") ");

      
      
      double scale_factor = mult * qWeights[l] * transformationDeterminant;

      for (size_t k = 0; k < basisFunctionsValues[l].size(); k++)
      {
        int global_dof_from_local = DOF[k];

        if (global_dof_from_local < U_Space->get_n_active())
        {
          double v1 = basisFunctionsValues[l][k]; // value of test function (vtest = vx = vy =vz)
          double v2 = v1;
          double v3 = v1;

          // add for all three components
          if (!rescale_by_h)
          {
            rhs.block(0)[global_dof_from_local] += scale_factor * uDirichlet[0] * v1;
            rhs.block(1)[global_dof_from_local] += scale_factor * uDirichlet[1] * v2;
            rhs.block(2)[global_dof_from_local] += scale_factor * uDirichlet[2] * v3;
          }
          else
          {
            double h = cell->Get_hK(0);

            rhs.block(0)[global_dof_from_local] += ( scale_factor * uDirichlet[0] * v1) /h;
            rhs.block(1)[global_dof_from_local] += ( scale_factor * uDirichlet[1] * v2) /h;
            rhs.block(2)[global_dof_from_local] += ( scale_factor * uDirichlet[2] * v3) /h;
          }
        }
      }
    }
  }
}


void TransformReferenceTS(TBaseCell* cell, int joint_index,
  std::vector<double> &qPointsT, std::vector<double> &qPointsS)
{
  const int *faceVertexMap, *faceVertexMapLength;
  int maxNVerticesPerFace;

  const TJoint* joint = cell->GetJoint(joint_index);

  if (joint->GetType() != BoundaryFace && joint->GetType() != IsoBoundFace)
  {
    ErrThrow("Called TransformReferenceTS on a non-boundary joint!");
  }

  const auto bd_face = reinterpret_cast<const TBoundFace*>(joint);
  const auto bd_comp = bd_face->GetBoundComp();

  cell->GetShapeDesc()->GetFaceVertex(faceVertexMap, faceVertexMapLength,
    maxNVerticesPerFace);

  size_t nFaceVertices = faceVertexMapLength[joint_index];
  std::vector<Point> faceVertices(nFaceVertices, Point(3u));

  for (size_t i = 0; i < nFaceVertices; i++)
  {
    double _x,_y,_z;
    cell->GetVertex(faceVertexMap[joint_index * maxNVerticesPerFace + i])
      ->GetCoords(_x,_y,_z);

    faceVertices[i].x = _x;
    faceVertices[i].y = _y;
    faceVertices[i].z = _z;
  }

  switch (faceVertices.size())
  {
    case 3:
    {
      Point vt = faceVertices[1] - faceVertices[0];
      Point vs = faceVertices[2] - faceVertices[0];

      size_t num_points = qPointsT.size();
      if (qPointsS.size() != num_points)
      {
        ErrThrow("Point array mismatch!");
      }

      for (size_t j = 0; j < num_points; j++)
      {
        Point p = faceVertices[0] + vt * qPointsT[j] + vs * qPointsS[j];

        bd_comp->GetTSofXYZ(p.x, p.y, p.z, qPointsT[j], qPointsS[j]);
      }

      break;
    }

    case 4:
    {
      Point p01 = faceVertices[1] - faceVertices[0];
      Point p10 = faceVertices[2] - faceVertices[0];
      Point p11 = faceVertices[3] - faceVertices[0];

      size_t num_points = qPointsT.size();
      if (qPointsS.size() != num_points)
      {
        ErrThrow("Point array mismatch!");
      }

      for (size_t j = 0; j < num_points; j++)
      {
        double t = qPointsT[j];
        double s = qPointsS[j];

        Point p = faceVertices[0]
        + (t * s) * p11
        + (t * (1 - s)) * p01
        + ((1 - t) * s) * p10;

        bd_comp->GetTSofXYZ(p.x, p.y, p.z, qPointsT[j], qPointsS[j]);
      }

      break;
    }

    default:
      ErrThrow("Unknown cell type in TransformReferenceTS!");
  }
}


// ===========================================================================
void BoundaryAssembling3D::getQuadratureDataIncludingFirstDerivatives(
    const TFESpace3D *fespace,TBaseCell *cell, int m,
    std::vector<double>& qWeights,std::vector<double>& qPointsT,
    std::vector<double>& qPointsS,
    std::vector< std::vector<double> >& basisFunctionsValues,
    std::vector< std::vector<double> >& basisFunctionsValues_derivative_x,
    std::vector< std::vector<double> >& basisFunctionsValues_derivative_y,
    std::vector< std::vector<double> >& basisFunctionsValues_derivative_z)
{
  // set quadrature formula and compute quadrature info
  auto fe = fespace->get_fe(cell->GetCellIndex());
  int fe_degree = fe.GetBaseFunct()->GetPolynomialDegree();
  const Shapes *face_types;
  cell->GetShapeDesc()->GetFaceType(face_types);
  // get a quadrature formula good enough for the velocity FE space
  const TQuadFormula *qf2 = QuadratureFormulaDatabase::qf_from_degree(
      2*fe_degree, face_types[m]);
  int N_Points = qf2->GetN_QuadPoints();
  // ====================================
  qWeights.resize(N_Points);
  qPointsT.resize(N_Points);
  qPointsS.resize(N_Points);

  for (size_t k = 0; k < (size_t) N_Points; k++)
  {
    auto p = qf2->get_point(k);
    qWeights[k] = qf2->get_weight(k);
    qPointsT[k] = p.x;
    qPointsS[k] = p.y;
  }

  TransformReferenceTS(cell, m, qPointsT, qPointsS);

  // ====================================
  // values of base functions in all quadrature points on face
  double **JointValues = FEDatabase::GetJointDerivatives3D(
      *fe.GetBaseFunct(), *qf2, m, MultiIndex3D::D000);

  double **JointValues_derivative_xi = FEDatabase::GetJointDerivatives3D(
      *fe.GetBaseFunct(), *qf2, m, MultiIndex3D::D100);

  double **JointValues_derivative_eta = FEDatabase::GetJointDerivatives3D(
      *fe.GetBaseFunct(), *qf2, m, MultiIndex3D::D010);

  double **JointValues_derivative_rho = FEDatabase::GetJointDerivatives3D(
      *fe.GetBaseFunct(), *qf2, m, MultiIndex3D::D001);

  fe.GetBaseFunct()->ChangeBF(fespace->GetCollection(), cell, qWeights.size(), JointValues);
  fe.GetBaseFunct()->ChangeBF(fespace->GetCollection(), cell, qWeights.size(), JointValues_derivative_xi);
  fe.GetBaseFunct()->ChangeBF(fespace->GetCollection(), cell, qWeights.size(), JointValues_derivative_eta);
  fe.GetBaseFunct()->ChangeBF(fespace->GetCollection(), cell, qWeights.size(), JointValues_derivative_rho);

  // convert the double** to a vector
  basisFunctionsValues.resize(qWeights.size());
  basisFunctionsValues_derivative_x.resize(qWeights.size());
  basisFunctionsValues_derivative_y.resize(qWeights.size());
  basisFunctionsValues_derivative_z.resize(qWeights.size());

  for (unsigned int l = 0; l < qWeights.size(); l++)
  {
    auto n_base_functions = fe.GetN_DOF();
    basisFunctionsValues[l].resize(n_base_functions);
    basisFunctionsValues_derivative_x[l].resize(n_base_functions);
    basisFunctionsValues_derivative_y[l].resize(n_base_functions);
    basisFunctionsValues_derivative_z[l].resize(n_base_functions);

    for (unsigned int k = 0; k < basisFunctionsValues[l].size(); k++)
    {
      basisFunctionsValues[l][k] = JointValues[l][k];
      basisFunctionsValues_derivative_x[l][k]=JointValues_derivative_xi[l][k];
      basisFunctionsValues_derivative_y[l][k]=JointValues_derivative_eta[l][k];
      basisFunctionsValues_derivative_z[l][k]=JointValues_derivative_rho[l][k];
    }
  }
}

// ===========================================================================
void BoundaryAssembling3D::computeNormalAndTransformationData(TBaseCell *cell, int m,
    std::vector<double>& normal,
    double &transformationDeterminant)
{
  const int *faceVertexMap, *faceVertexMapLength;
  int maxNVerticesPerFace;
  // For the current cell, get information of faces and local vertices
  // faceVertexMap should be seen as an array of arrays, e.g.
  // faceVertexMap = { {a,b,c},{b,c,d},{a,c,d},{a,b,d}}
  // where faceVertexMap[i] contains the id of vertices defining face i
  // faceVertexMapLength is an array specifying the length of each list
  // note: in the case that faces of an element have differennt number of
  // vertices (e.g. a pyramid), the faceVertexMap lists have all lenght equal to
  // maxNVerticesPerFace, and these are filled with 0 for the faces with less vertices
  cell->GetShapeDesc()->GetFaceVertex(faceVertexMap, faceVertexMapLength, maxNVerticesPerFace);
  // simplify: number of vertices on face m (m=joint_id)
  size_t nFaceVertices = faceVertexMapLength[m];
  std::vector< Point > faceVertices(nFaceVertices, Point((unsigned int) 3));

  for (size_t l1 = 0; l1 < nFaceVertices; l1++)
  {
    double _x,_y,_z;
    cell->GetVertex(faceVertexMap[ m * maxNVerticesPerFace + l1 ])->GetCoords(_x,_y,_z);
    faceVertices[l1].x = _x;
    faceVertices[l1].y = _y;
    faceVertices[l1].z = _z;
  }
  normal.clear();
  normal.resize(3);
  double xc1, yc1, zc1, xc2, yc2, zc2;
  double areaT, areaT1, areaT2, area_parallelogramm, area_parallelogramm_1, area_parallelogramm_2;

  switch (faceVertices.size())
  {
  case 3:
    xc1 = faceVertices[1].x - faceVertices[0].x;
    xc2 = faceVertices[2].x - faceVertices[0].x;

    yc1 = faceVertices[1].y - faceVertices[0].y;
    yc2 = faceVertices[2].y - faceVertices[0].y;

    zc1 = faceVertices[1].z - faceVertices[0].z;
    zc2 = faceVertices[2].z - faceVertices[0].z;

    // normal vector
    normal[0] = yc1*zc2 - zc1*yc2;
    normal[1] = zc1*xc2 - xc1*zc2;
    normal[2] = xc1*yc2 - yc1*xc2;
   /* // compute the 2 vectors that span the plane containing the current face
    xc1 = faceVertices[1].x() - faceVertices[0].x();
    xc2 = faceVertices[2].x() - faceVertices[0].x();

    yc1 = faceVertices[1].y() - faceVertices[0].y();
    yc2 = faceVertices[2].y() - faceVertices[0].y();

    zc1 = faceVertices[1].z() - faceVertices[0].z();
    zc2 = faceVertices[2].z() - faceVertices[0].z();

    // plane spanned by vectors v1=(xc1, yc1, zc1) and v2=(xc2, yc2, zc2)
    // Area of the triangle: 0.5*||v1 x v2||
    // normed Normal vector = (v1 x v2)/||v1 x v2||
    // Area of reference triangle (0,0)-(0,1)-(1,0): 1/2*g*h=0.5
 */   // Determinant of tranform.: A(triangle)/A(ref triangle) = ||v1 x v2||

    /* Old 15.10.18
    normal[0] = yc2*zc1 - zc2*yc1;// yc1*zc2 - zc1*yc2;
    normal[1] = zc2*xc1 - xc2*zc1; //zc1*xc2 - xc1*zc2;
    normal[2] = xc2*yc1 - yc2*xc1; //xc1*yc2 - yc1*xc2;
*/
    // determinant of reference trafo in order to get a normed normal vector
    area_parallelogramm = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
    normal[0] /= area_parallelogramm;
    normal[1] /= area_parallelogramm;
    normal[2] /= area_parallelogramm;
    areaT = area_parallelogramm / 2.0;
    transformationDeterminant = areaT / 0.5;
    break;
  case 4:
    // We consider a quadrilateral (P0,P1,P2,P3) as composed by 2 triangles
    // T1: P0,P1,P2
    // T2: P2,P3,P0
    // and we do the same as above (twice)
    // normed normal: ( (P1-P0) x (P2-P0) ) / || (P1-P0) x (P2-P0) ||
    // area: || (P1-P0) x (P2-P0) || / 2 + || (P3-P2) x (P0-P2) || / 2
    // area reference element [-1,1]x[-1,1]: 4
    // first triangle
    xc1 = faceVertices[1].x - faceVertices[0].x;
    xc2 = faceVertices[3].x - faceVertices[0].x; //faceVertices[2].x() - faceVertices[0].x();

    yc1 = faceVertices[1].y - faceVertices[0].y;
    yc2 = faceVertices[3].y - faceVertices[0].y; //faceVertices[2].y() - faceVertices[0].y();

    zc1 = faceVertices[1].z - faceVertices[0].z;
    zc2 = faceVertices[3].z - faceVertices[0].z; //faceVertices[2].z() - faceVertices[0].z();

    // normal vector (the same (except for length) for T1 and T2)
    normal[0] = yc2*zc1 - zc2*yc1;// yc1*zc2 - zc1*yc2;
    normal[1] = zc2*xc1 - xc2*zc1; //zc1*xc2 - xc1*zc2;
    normal[2] = xc2*yc1 - yc2*xc1; //xc1*yc2 - yc1*xc2;

    // determinant of reference transformation in order to get a normed normal vector
    area_parallelogramm_1 = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);

    normal[0] /= area_parallelogramm_1;
    normal[1] /= area_parallelogramm_1;
    normal[2] /= area_parallelogramm_1;

    areaT1 = area_parallelogramm_1 / 2.0;
    // second triangle
    xc1 = faceVertices[3].x - faceVertices[2].x;
    xc2 = faceVertices[1].x - faceVertices[2].x; //faceVertices[0].x() - faceVertices[2].x();

    yc1 = faceVertices[3].y - faceVertices[2].y;
    yc2 = faceVertices[1].y - faceVertices[2].y; //faceVertices[0].y() - faceVertices[2].y();

    zc1 = faceVertices[3].z - faceVertices[2].z;
    zc2 = faceVertices[1].z - faceVertices[2].z; //faceVertices[0].z() - faceVertices[2].z();


    // normal vector (the same (except for length) for T1 and T2)
    normal[0] = yc2*zc1 - zc2*yc1;// yc1*zc2 - zc1*yc2;
    normal[1] = zc2*xc1 - xc2*zc1; //zc1*xc2 - xc1*zc2;
    normal[2] = xc2*yc1 - yc2*xc1; //xc1*yc2 - yc1*xc2;

    // determinant of reference trasformation in order to get a normed normal vector
    area_parallelogramm_2 = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
    normal[0] /= area_parallelogramm_2;
    normal[1] /= area_parallelogramm_2;
    normal[2] /= area_parallelogramm_2;

    areaT2 = area_parallelogramm_2 / 2.0;

    // note: the reference element is [-1,1] x [-1,1]
    transformationDeterminant = (areaT1 + areaT2)/4.;

    break;
  default:
    ErrThrow("Unknown cell type in BoundaryAssembling3D::computeNormalAndTransformationData()");

  } // tria or quads
}

//// ===========================================================================
/////@todo This function should be a member of Boundface.C
//void BoundaryAssembling3D::compute_h(TBaseCell *cell,
//                                     int m,
//                                     double &h)
//{
//    const int *faceVertexMap, *faceVertexMapLength;
//    int maxNVerticesPerFace;
//    // For the current cell, get information of faces and local vertices
//    // faceVertexMap should be seen as an array of arrays, e.g.
//    // faceVertexMap = { {a,b,c},{b,c,d},{a,c,d},{a,b,d}}
//    // where faceVertexMap[i] contains the id of vertices defining face i
//    // faceVertexMapLength is an array specifying the length of each list
//    // note: in the case that faces of an element have differennt number of
//    // vertices (e.g. a pyramid), the faceVertexMap lists have all lenght equal to
//    // maxNVerticesPerFace, and these are filled with 0 for the faces with less vertices
//    cell->GetShapeDesc()->GetFaceVertex(faceVertexMap,faceVertexMapLength,maxNVerticesPerFace);
//    // simplify: number of vertices on face m (m=joint_id)
//    size_t nFaceVertices = faceVertexMapLength[ m ];
//    std::vector< Point > faceVertices(nFaceVertices,Point((unsigned int) 3));
//    for (size_t l1=0; l1<nFaceVertices; l1++)
//    {
//        double _x,_y,_z;
//        cell->GetVertex(faceVertexMap[ m*maxNVerticesPerFace+l1 ])->GetCoords(_x,_y,_z);
//        faceVertices[l1].x() = _x;
//        faceVertices[l1].y() = _y;
//        faceVertices[l1].z() = _z;
//    }
//    
//    std::vector<double> edge_lengths;
//    double min_edge_length, max_edge_length;
//    double edge_length1, edge_length2, edge_length3, edge_length4;
//    switch(faceVertices.size()) {
//       case 3:
//            // we consider a triangle (P0,P1,P2)
//            // and compute the edge lengths ||P1-P0||_2, ||P2-P0||_2, ||P1-P2||_2
//            edge_lengths[0] = std::sqrt  ((faceVertices[1].x - faceVertices[0].x)*(faceVertices[1].x - faceVertices[0].x)+
//                                     (faceVertices[1].y - faceVertices[0].y)*(faceVertices[1].y - faceVertices[0].y)+
//                                     (faceVertices[1].z - faceVertices[0].z)*(faceVertices[1].z - faceVertices[0].z)  );
//            
//            edge_lengths[1]= std::sqrt  ((faceVertices[2].x - faceVertices[0].x)*(faceVertices[2].x - faceVertices[0].x)+
//                                    (faceVertices[2].y - faceVertices[0].y)*(faceVertices[2].y - faceVertices[0].y)+
//                                    (faceVertices[2].z - faceVertices[0].z)*(faceVertices[2].z - faceVertices[0].z)  );
//            
//            edge_lengths[2] = std::sqrt  ((faceVertices[1].x - faceVertices[2].x)*(faceVertices[1].x - faceVertices[2].x)+
//                                     (faceVertices[1].y - faceVertices[2].y)*(faceVertices[1].y - faceVertices[2].y)+
//                                     (faceVertices[1].z - faceVertices[2].z)*(faceVertices[1].z - faceVertices[2].z)  );
//            
//            // compute the maximum and minimum edge length
//            min_edge_length = *min_element(edge_lengths.begin(), edge_lengths.end());
//            max_edge_length = *max_element(edge_lengths.begin(), edge_lengths.end());
//            
//            // set h to be the maximum edge length
//            h=max_edge_length;
//            
//            break;
//            
//        case 4:
//            // we consider a quadrilateral (P0,P1,P2,P3)
//            // and compute the edge lengths ||P1-P0||_2, ||P2-P1||_2, ||P3-P2||_2, ||P3-P0||_2
//            edge_lengths[0] = std::sqrt  ((faceVertices[1].x - faceVertices[0].x)*(faceVertices[1].x - faceVertices[0].x)+
//                                     (faceVertices[1].y - faceVertices[0].y)*(faceVertices[1].y - faceVertices[0].y)+
//                                     (faceVertices[1].z - faceVertices[0].z)*(faceVertices[1].z - faceVertices[0].z)  );
//            
//            edge_lengths[1] = std::sqrt  ((faceVertices[2].x - faceVertices[1].x)*(faceVertices[2].x - faceVertices[1].x)+
//                                     (faceVertices[2].y - faceVertices[1].y)*(faceVertices[2].y - faceVertices[1].y)+
//                                     (faceVertices[2].z - faceVertices[1].z)*(faceVertices[2].z - faceVertices[1].z)  );
//            
//            edge_lengths[2] = std::sqrt  ((faceVertices[3].x - faceVertices[2].x)*(faceVertices[3].x - faceVertices[2].x)+
//                                     (faceVertices[3].y - faceVertices[2].y)*(faceVertices[3].y - faceVertices[2].y)+
//                                     (faceVertices[3].z - faceVertices[2].z)*(faceVertices[3].z - faceVertices[2].z)  );
//            
//            edge_lengths[3] = std::sqrt  ((faceVertices[3].x - faceVertices[0].x)*(faceVertices[3].x - faceVertices[0].x)+
//                                     (faceVertices[3].y - faceVertices[0].y)*(faceVertices[3].y - faceVertices[0].y)+
//                                     (faceVertices[3].z - faceVertices[0].z)*(faceVertices[3].z - faceVertices[0].z)  );
//            
//            // compute the maximum and minimum edge length
//            min_edge_length = *min_element(edge_lengths.begin(), edge_lengths.end());
//            max_edge_length = *max_element(edge_lengths.begin(), edge_lengths.end());
//            
//            // set h to be the maximum edge length
//            h=max_edge_length;
//            
//            break;
//            
//    } // triangles or quadrilaterals
//    
//}



// ===========================================================================

void BoundaryAssembling3D::nitsche_bc(BlockFEMatrix &s_matrix, BlockVector &s_rhs,
				      const TFESpace3D * v_space, const TFESpace3D *p_space,
				      BoundValueFunct3D * U1, BoundValueFunct3D *U2, BoundValueFunct3D *U3,
				      TFEVectFunct3D* U_FE,
				      std::vector<TBoundFace*>& boundaryFaceList,
				      int bd_comp, double gamma, double mu,
				      int sym_u, int sym_p)
{

  if (TDatabase::ParamDB->NSTYPE != 14)
  {
    Output::print("WARNING BoundaryAssembling2D::nitsche_bc(..): The NSTYPE is not equal to 14. This might result in errors related to bad access.");
  }

  //============================== PENALTY TERMS ===================================
  // gamma/h (u,v)
  // rescale local integral by edge values
  matrix_u_v(s_matrix, v_space, boundaryFaceList, bd_comp, gamma * mu, true);

  // gamma/h (uD,v) [rhs]
  // rescale local integral by edge values
  rhs_uD_v(s_rhs, v_space, U1, U2, U3, U_FE,boundaryFaceList, bd_comp, gamma * mu, true);

  /* // todo
    // sigma * L_0^2 * gamma/h (u.n,v.n)
  matrix_u_n_v_n(s_matrix, v_space,  bd_comp, gamma * sigma * L_0 * L_0, true); // true = rescale local integral by edge values

  // sigma * L_0^2 * gamma/h (uD.n,v.n)
  rhs_uD_n_v_n(s_rhs, v_space, U1, U2, U3, bd_comp, gamma * sigma * L_0 * L_0, true);   // true = rescale local integral by edge values
  */
  
  if ((mu == 0.) & (gamma != 0))
    ErrThrow("WARNING/ERROR: Penalty terms of Nitsche method not yet implemented in 3D for mueff = 0 (Darcy Limit)");

  //=========================== PENALTY-FREE TERMS =================================
  // - (mu grad(u)n,v)
  matrix_gradu_n_v(s_matrix, v_space, boundaryFaceList, bd_comp, (-1.) * mu);  // OK

  // - sign_u * (u, mu grad(v)n) [sign_u=1: symmetric, -1: skew-symmetric]
  matrix_gradv_n_u(s_matrix, v_space, boundaryFaceList, bd_comp, (-1.) * sym_u * mu);

  // - sign_u * (uD,mu grad(v)n) [rhs]
  rhs_gradv_n_uD(s_rhs, v_space, U1, U2, U3, U_FE, boundaryFaceList, bd_comp, (-1) * sym_u * mu );

  // (pn,v)
  matrix_p_v_n(s_matrix, v_space, p_space, boundaryFaceList, bd_comp, 1.);

  // sign_p * (u,qn)
  matrix_q_u_n(s_matrix, v_space, p_space, boundaryFaceList, bd_comp, sym_p);

  // sign_p * (uD,qn) [rhs]
  rhs_q_uD_n(s_rhs, v_space, p_space, U1, U2, U3, U_FE, boundaryFaceList, bd_comp, sym_p);
}


void BoundaryAssembling3D::nitsche_bc_nonlinear_iteration(BlockFEMatrix &s_matrix,
							  const TFESpace3D * v_space, const TFESpace3D * /*p_space*/,
							  std::vector<TBoundFace*>& boundaryFaceList,
							  int bd_comp, double gamma, double mu,
							  int sym_u, int /*sym_p*/)
{

  if (TDatabase::ParamDB->NSTYPE != 14)
  {
    Output::print("WARNING BoundaryAssembling2D::nitsche_bc(..): The NSTYPE is not equal to 14. This might result in errors related to bad access.");
  }

  //============================== PENALTY TERMS ===================================
  // gamma/h (u,v)
  // rescale local integral by edge values
  matrix_u_v(s_matrix, v_space, boundaryFaceList, bd_comp, gamma * mu, true);


  // - (mu grad(u)n,v)
  matrix_gradu_n_v(s_matrix, v_space, boundaryFaceList, bd_comp, (-1.) * mu);  // OK

  // - sign_u * (u, mu grad(v)n) [sign_u=1: symmetric, -1: skew-symmetric]
  matrix_gradv_n_u(s_matrix, v_space, boundaryFaceList, bd_comp, (-1.) * sym_u * mu);

}


void BoundaryAssembling3D::nitsche_bc_matrix(BlockFEMatrix &s_matrix, 
					     const TFESpace3D * v_space, const TFESpace3D *p_space,
					     std::vector<TBoundFace*>& boundaryFaceList,
					     int bd_comp, double gamma, double mu,
					     int sym_u, int sym_p)
{

  if (TDatabase::ParamDB->NSTYPE != 14)
  {
    Output::print("WARNING BoundaryAssembling2D::nitsche_bc(..): The NSTYPE is not equal to 14. This might result in errors related to bad access.");
  }

  //============================== PENALTY TERMS ===================================
  // gamma/h (u,v)
  // rescale local integral by edge values
  matrix_u_v(s_matrix, v_space, boundaryFaceList, bd_comp, gamma * mu, true);

if ((mu == 0.) & (gamma != 0))
  ErrThrow("WARNING/ERROR: Penalty terms of Nitsche method not yet implemented in 3D for mueff = 0 (Darcy Limit)");

  //=========================== PENALTY-FREE TERMS =================================
  // - (mu grad(u)n,v)
  matrix_gradu_n_v(s_matrix, v_space, boundaryFaceList, bd_comp, (-1.) * mu);  // OK

  // - sign_u * (u, mu grad(v)n) [sign_u=1: symmetric, -1: skew-symmetric]
  matrix_gradv_n_u(s_matrix, v_space, boundaryFaceList, bd_comp, (-1.) * sym_u * mu);

  // (pn,v)
  matrix_p_v_n(s_matrix, v_space, p_space, boundaryFaceList, bd_comp, 1.);

  // sign_p * (u,qn)
  matrix_q_u_n(s_matrix, v_space, p_space, boundaryFaceList, bd_comp, sym_p);

}



// ===========================================================================
void BoundaryAssembling3D::getQuadratureData(const TFESpace3D *fespace, TBaseCell *cell, int m,
    std::vector<double>& qWeights, std::vector<double>& qPointsT,
    std::vector<double>& qPointsS,
    std::vector< std::vector<double> >& basisFunctionsValues)
{
  // set quadrature formula and compute quadrature info
  auto fe = fespace->get_fe(cell->GetCellIndex());
  int fe_degree = fe.GetBaseFunct()->GetPolynomialDegree();

  const Shapes *face_types;
  cell->GetShapeDesc()->GetFaceType(face_types);

  // get a quadrature formula good enough for the velocity FE space
  const TQuadFormula *qf2 = QuadratureFormulaDatabase::qf_from_degree(
      2 * fe_degree, face_types[m]);

  int N_Points = qf2->GetN_QuadPoints();

  // ====================================
  // values of base functions in all quadrature points on face

  double **JointValues = FEDatabase::GetJointDerivatives3D(
      *fe.GetBaseFunct(), *qf2, m, MultiIndex3D::D000);

  fe.GetBaseFunct()->ChangeBF(fespace->GetCollection(), cell, N_Points,
                              JointValues);

  // ====================================
  // convert the double* to vectors

  qWeights.resize(N_Points);
  qPointsT.resize(N_Points);
  qPointsS.resize(N_Points);

  for (size_t k = 0; k < (size_t) N_Points; k++)
  {
    auto p = qf2->get_point(k);
    qWeights[k] = qf2->get_weight(k);

    qPointsT[k] = p.x;
    qPointsS[k] = p.y;
  }

  TransformReferenceTS(cell, m, qPointsT, qPointsS);

  basisFunctionsValues.resize(qWeights.size());

  for (unsigned int l = 0; l < qWeights.size(); l++)
  {
    basisFunctionsValues[l].resize(fe.GetN_DOF());

    for (unsigned int k = 0; k < basisFunctionsValues[l].size(); k++)
    {
      basisFunctionsValues[l][k] = JointValues[l][k];
    }
  }
  // ====================================
}

const TQuadFormula* BoundaryAssembling3D::getFaceQuadratureFormula(
  TBaseCell *cell, int m, int _degree) const
{
  const Shapes *face_types;
  cell->GetShapeDesc()->GetFaceType(face_types);

  auto FaceQuadFormula = QuadratureFormulaDatabase::qf_from_degree(
      _degree, face_types[m]);

  return FaceQuadFormula;
}

// ===========================================================================
void BoundaryAssembling3D::getFaceQuadratureData(
    const TCollection& coll, const BaseFunctions& bf, TBaseCell *cell, int m,
    std::vector<double>& qWeights, std::vector<double>& qPointsT,
    std::vector<double>& qPointsS,
    std::vector< std::vector<double> >& basisFunctionsValues, int _degree) const
{
  // set quadrature formula and compute quadrature info
  const TQuadFormula *qf2 = getFaceQuadratureFormula(cell, m, _degree);
  int N_Points = qf2->GetN_QuadPoints();

  // ====================================
  // values of base functions in all quadrature points on face
  double **JointValues = FEDatabase::GetJointDerivatives3D(
      bf, *qf2, m, MultiIndex3D::D000);

  bf.ChangeBF(&coll, cell, N_Points, JointValues);

  // ====================================
  // convert the double* to vectors

  qWeights.resize(N_Points);
  qPointsT.resize(N_Points);
  qPointsS.resize(N_Points);

  for (size_t k = 0; k < (size_t)N_Points; k++)
  {
    auto p = qf2->get_point(k);
    qWeights[k] = qf2->get_weight(k);
    qPointsT[k] = p.x;
    qPointsS[k] = p.y;
  }

  TransformReferenceTS(cell, m, qPointsT, qPointsS);

  basisFunctionsValues.resize(qWeights.size());
  for (unsigned int l = 0; l < qWeights.size(); l++)
  {
    basisFunctionsValues[l].resize(bf.GetDimension());

    for (unsigned int k = 0; k < basisFunctionsValues[l].size(); k++)
    {
      basisFunctionsValues[l][k] = JointValues[l][k];
    }
  }

  // ====================================
}


// ===========================================================================
// analogous to getFaceQuadratureData, but taking as input the quad. points
void BoundaryAssembling3D::getFaceQuadratureValue(
    const TCollection& coll, const BaseFunctions& bf, TBaseCell *cell, int m,
    const TQuadFormula * qf2,
    std::vector< std::vector<double> >& basisFunctionsValues) const
{
  int N_Points = qf2->GetN_QuadPoints();

  // ====================================
  // generate data on reference mesh cell for the 2d face of 3d cell
  // set quadrature formula and compute quadrature info
  // values of base functions in all quadrature points on face
  double **JointValues = FEDatabase::GetJointDerivatives3D(
      bf, *qf2, m, MultiIndex3D::D000);

  bf.ChangeBF(&coll, cell, N_Points, JointValues);

  basisFunctionsValues.resize(N_Points);
  for(int l=0; l<N_Points; l++)
  {
    basisFunctionsValues[l].resize(bf.GetDimension());
    for(unsigned int k=0; k<basisFunctionsValues[l].size(); k++)
    {
      basisFunctionsValues[l][k]=JointValues[l][k];
    }
  }
  // ====================================
  
}

