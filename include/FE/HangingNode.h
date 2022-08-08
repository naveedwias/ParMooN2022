#ifndef __HANGINGNODE__
#define __HANGINGNODE__

#include "HNDesc.h"
#include <vector>

/** represent a hanging node */
class THangingNode
{
  protected:
    /** type of hanging node */
    const THNDesc * hn_descriptor;

    /** numbers of degrees of freedom in coupling */
    std::vector<int> DOF;
    
    /** hanging nodes mostly come with other hanging nodes on the same joint,
     *  these are the dofs of these other hanging nodes for this hanging node
     * 
     * If there is only one hanging dof on a joint (e.g. P1/Q1) this vector is
     * empty.
     */
    std::vector<int> partner_dofs;

  public:
    /** constructor, filling all data */
    THangingNode(const THNDesc *hn_desc, const int *all_dofs,
                 const int *coupling, const std::vector<int>& partners);

    /** destructor */
    ~THangingNode() = default;

    // Methods
    /** return number of nodes in coupling */
    HNDesc GetType() const
    { return hn_descriptor->get_type(); }

    /** @brief return number of nodes in coupling */
    int GetN_Nodes() const
    { return hn_descriptor->GetN_Nodes(); }

    /** @brief return coefficients of degrees of freedom in coupling */
    const double *GetCoeff() const
    { return hn_descriptor->GetCoeff(); }

    /** return numbers of degrees of freedom in coupling */
    const int *GetDOF() const
    { return DOF.data(); }

    /** @brief find out if `p_dof` is one of the other hanging nodes on the 
     * joint where this hanging nodes is located at. */
    bool has_partner(int p_dof) const;
    
    /** @brief return vector of dofs of other hanging nodes on this joint */
    std::vector<int> get_partners() const
    { return partner_dofs; }

    int *GetDOF()
    { return DOF.data(); }

    std::vector<int>& get_partners()
    { return partner_dofs; }
};

#endif
