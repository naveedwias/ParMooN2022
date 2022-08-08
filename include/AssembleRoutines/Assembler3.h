/** ***************************************************************************
 *
 * @name   2D_Assembler
 * @brief  Assemble
 *
 *
 *
 * @author     Alfonso Caiazzo and Laura Blank
 * @date       16.08.2016
 *
 ******************************************************************************/

#pragma once

#include <Constants.h>
#include <BlockVector.h>
#include <BlockFEMatrix.h>
#include <Example2D.h>
#include "LocalAssembling.h"

/**
   @brief Todo...
 */

class Assembler3{
    
protected:
    
public:
    
    /** @brief constructor (taken from LocalAssembling)
     *
     * In the current version, the Assembler object is created using
     * a constructor analogous to LocalAssembling.
     * This is just a temporary solution and should probably done differently
     */
  Assembler3(LocalAssembling2D_type type, TFEFunction2D **fefunctions2d,
	     CoeffFct2D *coeffs);
  
  ///@brief The variational form to be assembled
  LocalAssembling2D la;
  ///@brief The collection on which we want to assemble
  TCollection *Coll;
  
  ///@brief Vector for handling hanging nodes
  std::vector< std::vector<double> > hangingEntries, hangingRhs;
  
  ///@brief blocks to be assembled
  std::vector< TSquareMatrix2D* > square_matrices;
  std::vector< TMatrix2D* > rectangular_matrices;
  std::vector< double* > rhs_blocks;
  
  /** a function from a finite element space */
  void Assemble2D(BlockFEMatrix &M,BlockVector &b_rhs,
		  std::vector <const TFESpace2D*>& fespaces,
		  std::vector <const TFESpace2D*>& ferhs,
		  const Example2D& example,
		  int AssemblePhaseID=-1
		  );
  
  /**
     @brief Set the variables of the class (hanging nodes, allocate pointers)
     
     This is a temporary solution.
  */
  void init(BlockFEMatrix &M,BlockVector &b_rhs,
	    std::vector<const TFESpace2D*>& fespaces,
	    std::vector<const TFESpace2D*>& ferhs);

  /**
     @brief loop over all cells and assemble local terms
   */
  void loop_over_cells(std::vector <const TFESpace2D*>& fespaces,
		       int AssemblePhaseID,
		       double ***LocMatrices,
		       double *righthand,
		       int N_AllMatrices, double **Matrices,
		       double **Param,
		       double **AuxArray,double **LocRhs,
		       std::vector<const TFESpace2D*>& ferhs,
		       const Example2D& example);
  
  void imposeBoundaryConditions(const TFESpace2D *fespace, int j,
				const Example2D& example,
				TBaseCell* cell,
				int i,
				double *RHS,
				int *DOF);

  
  void handleHangingNodes(std::vector<const TFESpace2D*>& ferhs);

 
  void addLocalToGlobalRhs(TBaseCell *cell, int i,
			   std::vector<const TFESpace2D*>& ferhs,
			   double *righthand);
  
  ~Assembler3();
  
  
  
  
  
};
