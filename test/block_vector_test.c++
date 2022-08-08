#include <vector>
#include <BlockVector.h>

#include <MooNMD_Io.h>

void normTest();

int main(int, char**)
{
  Output::print("Testing BlockVector.");
  
  normTest();
  
  return 0;
}

void normTest()
{
  // Test for BlockVector::norm in sequential case
  #ifndef _MPI
  Output::print("Testing BlockVector::norm in sequential case.");
  
  std::vector<unsigned int> lengths = {1, 2, 1, 2};
  std::vector<unsigned int> selectedBlocks = {1, 2};
  BlockVector blockVector(lengths);
  
  blockVector[0] = 13;
  blockVector[1] = 3;
  blockVector[2] = 4;
  blockVector[3] = 12;
  blockVector[4] = 13;
  blockVector[5] = 13;
  
  const int norm = 26;
  const int normOfSelectedBlocks = 13;
  
  if( blockVector.norm() != norm )
    ErrThrow("Norm of blockVector is ", blockVector.norm(), ", but it should be ", norm, ".");
  
  if( blockVector.norm(selectedBlocks) != normOfSelectedBlocks )
    ErrThrow("Norm of selected blocks is ", blockVector.norm(selectedBlocks), ", but it should be ", normOfSelectedBlocks, ".");
  #endif
}
