#include "FEFunction1D.h"
#include "BaseCell.h"
#include "MooNMD_Io.h"

TFEFunction1D::TFEFunction1D(std::shared_ptr<TFESpace1D> fespace1D,
                             const std::string& name, double *values,
                             int length)
: Name(name), FESpace1D(fespace1D), Values(values), Length(length)
{
}

TFEFunction1D::~TFEFunction1D()
{
}

