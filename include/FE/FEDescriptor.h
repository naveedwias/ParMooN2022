#ifndef __FEDESC3D__
#define __FEDESC3D__

#include "Enumerations_fe.h"

/** @brief store a finite element descriptor for a 3D element */
class FEDescriptor
{
  protected:
    /** @brief identifier for this fe descriptor */
    FEDescriptor_type type;
    
    /** @brief number of degrees of freedom */
    int N_DOF;

    /** @brief number of degrees of freedom on closure of each joint */
    int N_JointDOF;

    /** @brief local numbers of all degrees of freedom on the joints */
    int **JointDOF;

    /** @brief number of inner degrees of freedom */
    int N_InnerDOF;

    /** @brief local numbers of all inner degrees of freedom */
    int *InnerDOF;

    /** @brief number of degrees of freedom on cell boundary */
    int N_OuterDOF;

    /** @brief local numbers of all degrees of freedom on cell boundary */
    int *OuterDOF;

    /// @brief members needed for mpi
    //@{
    /** @brief flag, true if all data needed for mpi is set. */
    bool EdgeVertData_Filled;

    /** @brief number of degrees of freedom on closure of each edge */
    int N_EdgeDOF;

    /** @brief local numbers of all degrees of freedom on the edge */
    int **EdgeDOF;

    /** @brief number of degrees of freedom on each vertex */
    int N_VertDOF;

    /** @brief local numbers of all degrees of freedom on the vertices */
    int *VertDOF;
    //@}

  public:
    
    explicit FEDescriptor ( FEDescriptor_type id);
    
    /** constructor, setting all data with dof on cell boundary */
    FEDescriptor (char *description, int n_dof, int n_jointdof,
              int **jointdof, int n_innerdof, int *innerdof,
              int n_outerdof = 0, int *outerdof = nullptr);

    /// @brief members needed for mpi
    //@{
    /** @brief return number of degrees of freedom per closure of each edge */
    int GetN_EdgeDOF() const
    { return N_EdgeDOF; }

    /** @brief return local numbers of degrees of freedom on each edge */
    int **GetEdgeDOF() const
    { return EdgeDOF; }

    /** @brief return local numbers of degrees of freedom on edge i */
    int *GetEdgeDOF(int i) const
    { return EdgeDOF[i]; } 

    /** @brief return number of degrees of freedom per closure of each vertex */
    int GetN_VertDOF() const
    { return N_VertDOF; }

    /** @brief return local numbers of degrees of freedom on vertex i */
    int  GetVertDOF(int i) const
    { return VertDOF[i]; }

    bool IsEdgeVertData_Filled() const
    {return EdgeVertData_Filled; }
    //@}
    
    /** @brief return number of degrees of freedom */
    int GetN_DOF() const
    { return N_DOF; }

    /** @brief return number of degrees of freedom per closure of each joint */
    int GetN_JointDOF() const
    { return N_JointDOF; }

    /** @brief return number of inner degrees of freedom */
    int GetN_InnerDOF() const
    { return N_InnerDOF; }

    /** @brief return local numbers of inner degrees of freedom */
    int *GetInnerDOF() const
    { return InnerDOF; }

    /** @brief return number of degrees of freedom on cell boundary */
    int GetN_OuterDOF() const
    { return N_OuterDOF; }

    /** @brief return local numbers of degrees of freedom on cell boundary */
    int *GetOuterDOF() const
    { return OuterDOF; }

    /** @brief return total number and local numbers of degrees of freedom
        on cell boundary */ 
    void GetOuterDOF(int &n_outerdof, int* &outerdof) const
    { n_outerdof = N_OuterDOF; outerdof = OuterDOF; }

    /** @brief return local numbers of degrees of freedom on each joint */
    int **GetJointDOF() const
    { return JointDOF; }

    /** @brief return local numbers of degrees of freedom on joint i */
    int *GetJointDOF(int i) const
    { return JointDOF[i]; }
    
    /** @brief return ID for this fe descriptor */
    FEDescriptor_type GetID() const
    { return type; }

    /** @brief return face on which the i-th local degree of freedom is.
    If i is not a dof on a face, return -1.
    If i is a dof on two faces (e.g. on a vertex), one of these two faces is 
    returned. Don't use this function in this case.
    */
    int GetJointOfThisDOF(int localDOF) const;

};

#endif
