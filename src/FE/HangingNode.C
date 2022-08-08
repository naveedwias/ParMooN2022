#include "HangingNode.h"
#include "MooNMD_Io.h"
#include <algorithm>

THangingNode::THangingNode(const THNDesc *hn_desc, const int *all_dofs,
                           const int *coupling,
                           const std::vector<int>& partners)
: hn_descriptor(hn_desc), DOF(), partner_dofs(partners)
{
  if(!hn_desc)
    ErrThrow("please provide a valid THNDesc pointer.");
  int n_nodes = hn_desc->GetN_Nodes();
  DOF.resize(n_nodes);
  for(int i = 0; i < n_nodes; i++)
    DOF[i] = all_dofs[coupling[i]];
}

bool THangingNode::has_partner(int p_dof) const
{
  auto end = partner_dofs.end();
  return std::find(partner_dofs.begin(), end, p_dof) != end;
}
