
#include <MacroCell.h>

// Constructors
TMacroCell::TMacroCell(const TRefDesc *refdesc, int reflevel) :
              TGridCell(refdesc, reflevel)
{
  SubGridID = 0;
}

// Methods
