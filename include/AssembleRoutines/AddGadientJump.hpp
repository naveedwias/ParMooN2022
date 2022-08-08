#ifndef ADDGADIENTJUMP_H
#define ADDGADIENTJUMP_H

#include "FESpace2D.h"
#include "LocalAssembling.h"
#include "SquareMatrix2D.h"

void ComputeAddStab(const TFESpace2D** space, LocalAssembling2D& la, 
                    TSquareMatrix2D **sqmatrices, TSquareMatrix2D **A22);
#endif // ADDGADIENTJUMP_H
