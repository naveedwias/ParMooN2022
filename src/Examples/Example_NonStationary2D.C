#include <Example_NonStationary2D.h>

Example_NonStationary2D::Example_NonStationary2D(const ParameterDatabase & db) 
 : Example2D(db), timeDependentRhs(true), 
   timeDependentCoeffs(true), initialCondition()
{
}

Example_NonStationary2D::Example_NonStationary2D(
  const std::vector<DoubleFunct2D*>& exact,
  const std::vector<BoundCondFunct2D*>& bc,
  const std::vector<BoundValueFunct2D*>& bd,
  const CoeffFct2D& coeffs, bool timedependentrhs, bool timedependentcoeffs,
  std::vector<DoubleFunct2D*> init_cond)
 :Example2D(exact, bc, bd, coeffs)
 , timeDependentRhs(timedependentrhs) 
 , timeDependentCoeffs(timedependentcoeffs)
 , initialCondition(init_cond)
{
  
}

