#ifndef __FESPACE3D__
#define __FESPACE3D__

#include "FESpace.h"
#include <vector>

#ifdef _MPI
#include <memory>
class TParFECommunicator3D;
class TParFEMapper3D;
#endif

/** class for all 3D finite element spaces */
class TFESpace3D : public TFESpace
{
  protected:
    /**
     *  @brief Boundary condition used to create this space.
     *
     *  In order to reproduce the information, which boundary condition was used
     *  in the creation of the space, store a function pointer to it.
     */
   BoundCondFunct3D* boundCondition_;

#ifdef  _MPI
   /// There belongs a ParFECommunicator to this space. Store it!
   std::shared_ptr<TParFECommunicator3D> comm_;

   ///The communicator needs a mapper TODO Move mapper into comm as member!
   std::shared_ptr<TParFEMapper3D> mapper_;

#endif

    /** @brief construct space
     * 
     * @warning this method is a huge black box. Noone is able to understand 
     * this code, it needs refactoring. We rely on this heavily.
     */
    void ConstructSpace();

    /** @brief find an upper bound for the number of boundary dofs for each type
     * 
     * For example, the first index in the returned array corresponds to
     * Dirichlet boundary dofs. For the order of the boundary types, see the
     * enumeration `BoundCond` in Constants.h.
     * 
     * This computes only an upper bound because it counts local degrees of 
     * freedom. Essentially this loops over all boundary faces and adds the 
     * number of dofs on that face to the respective entry in the array. Since
     * some degrees of freedom are counted multiple times, the returned numbers
     * are greater or equal to the actual numbers.
     */
    std::array<int, TFESpace::N_DiffBoundNodeTypes>
      find_boundary_upper_bounds() const;

  public:
    /** constructor for building a space with elements of order k */
    TFESpace3D(const TCollection *coll, const std::string& name,
               BoundCondFunct3D *BoundaryCondition, int k);

    /** constructor for building a space with the given elements */
    TFESpace3D(const TCollection *coll, const std::string& name,
               BoundCondFunct3D *BoundaryCondition, FE_type *fes);

    TFESpace3D(const TCollection *coll, const std::string& name,
               BoundCondFunct3D *BoundaryCondition, SpaceType type, int ord);

    /** destructor */
    ~TFESpace3D() = default;

    /** @return The boundary condition function pointer. */
    BoundCondFunct3D* get_boundary_condition() const
    { return boundCondition_; }

    /** @brief return the boundary condition on a given boundary joint. */
    BoundCond get_boundary_condition(const BoundaryJoint& bd_joint) const;

    /** return position of one given DOF */
    void GetDOFPosition(int dof, double &x, double &y, double &z) const;

    /**
     * velocity space, if there is a element that only has dirichlet dof's??
     */
    bool CheckMesh() const;

#ifdef  _MPI
    const TParFECommunicator3D& get_communicator() const
    { return *comm_; }

    TParFECommunicator3D& get_communicator()
    { return *comm_; }

    const TParFEMapper3D& get_mapper() const
    { return *mapper_; }

    TParFEMapper3D& get_mapper()
    { return *mapper_; }
#endif // _MPI

};

#endif // __FESPACE3D__
