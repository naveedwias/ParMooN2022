#ifndef __FESPACE2D__
#define __FESPACE2D__

#include "FESpace.h"
#include <vector>
#include <map>

/** @brief class for all 2D finite element spaces */
class TFESpace2D : public TFESpace
{
  protected:
    /** @brief boundary condition used to create this space */
    BoundCondFunct2D *BoundCondition;

    /** @brief construct space */
    void ConstructSpace();

    /** @brief find an upper bound for the number of boundary dofs for each type
     * 
     * For example, the first index in the returned array corresponds to
     * Dirichlet boundary dofs. For the order of the boundary types, see the
     * enumeration `BoundCond` in Constants.h.
     * 
     * This computes only an upper bound because it counts local degrees of 
     * freedom. Essentially this loops over all boundary edges and adds the 
     * number of dofs on that edge to the respective entry in the array. Since
     * some degrees of freedom are counted multiple times, the returned numbers
     * are greater or equal to the actual numbers.
     */
    std::array<int, TFESpace::N_DiffBoundNodeTypes>
      find_boundary_upper_bounds() const;

  public:

    /** @brief constructor for building a space with elements of order k */
    TFESpace2D(const TCollection *coll, const std::string& name,
               BoundCondFunct2D *BoundaryCondition, int k);

    TFESpace2D(const TCollection *coll, const std::string& name,
               BoundCondFunct2D *BoundaryCondition, SpaceType type, int k);

    /** @brief constructor for building a space with the given elements */
    TFESpace2D(const TCollection *coll, const std::string& name,
               BoundCondFunct2D *BoundaryCondition, FE_type *fes);

    TFESpace2D(const  TFESpace2D&)=delete;
    TFESpace2D& operator=(TFESpace2D) = delete;
    /** destructor */
    ~TFESpace2D() = default;

    /** @brief return position of one given DOF */
    void GetDOFPosition(int dof, double &x, double &y) const;

    /** @brief return boundary condition */
    BoundCondFunct2D *get_boundary_condition() const
    { return BoundCondition; }
    
    /** @brief return the boundary condition on a given boundary joint. */
    BoundCond get_boundary_condition(const BoundaryJoint& bd_joint) const;
    
    friend  bool operator== (const TFESpace2D &lhs, const TFESpace2D &rhs);
    friend  bool operator!= (const TFESpace2D &lhs, const TFESpace2D &rhs);
};


#endif // __FESPACE2D__
