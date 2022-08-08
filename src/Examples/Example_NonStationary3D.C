#include <Example_NonStationary3D.h>

Example_NonStationary3D::Example_NonStationary3D(const ParameterDatabase & db)
 : Example3D(db)
 , timeDependentRhs()
 , timeDependentCoeffs()
 , initialCondtion()
{

}

Example_NonStationary3D::Example_NonStationary3D(
  const std::vector<DoubleFunct3D*>& exact,
  const std::vector<BoundCondFunct3D*>& bc,
  const std::vector<BoundValueFunct3D*>& bd, 
  const CoeffFct3D& coeffs, bool timedependentrhs, bool timedependentcoeffs,
  const std::vector<DoubleFunct3D*>& init_cond)
: Example3D(exact, bc, bd, coeffs)
  , timeDependentRhs(timedependentrhs)
  , timeDependentCoeffs(timedependentcoeffs)
  , initialCondtion(init_cond)
{

}
