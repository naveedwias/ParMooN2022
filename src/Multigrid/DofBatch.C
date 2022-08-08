/*
 * DofBatch.C
 *
 *  @date Mar 31, 2015
 *  @author Clemens Bartsch
 *
 *  Implementation of class DofBatch declared in include/FE/DofBatch.h.
 *
 */

#include <DofBatch.h>
#include <algorithm>
#include <vector>

/*!
 * @brief Construct a DofBatch with known size.
 *
 * @param[in] nEntries The number of dofs which will be added later. Defaults to 1.
 */
DofBatch::DofBatch(const size_t nDofs) :
dofList_()
{
	//Reserve space for the entries to be added later.
	dofList_.reserve(nDofs);
}

/*!
 * @brief Add a dof to the batch.
 *
 * @param[in] newDof The number of the dof to be added to the batch.
 */
void DofBatch::addDof(const int newDof){
	dofList_.push_back(newDof);
}

/*!
 * @brief Return the number of the dof at place 'index'.
 * @return The MooNMD internal number of the dof.
 *
 * @param[in] index The place of the dof in the batch.
 */
int DofBatch::getDof(const size_t index) const{
	return dofList_.at(index);
}

/*!
 * @brief Removes duplicates and sorts the dofList_.
 *
 * This must be called when the procedure of determining which dofs to add and adding them is finished.
 */
void DofBatch::tidyUp(){
	// Sort the entries.
	std::sort(dofList_.begin(), dofList_.end());

	// Throw out duplicates, by means of std methods.
	intVec::iterator it = std::unique(dofList_.begin(), dofList_.end());
	dofList_.resize( std::distance(dofList_.begin(),it) );

	// Adjust capacity.
	dofList_.shrink_to_fit();
}

/*!
 * @brief Return the number of dofs stored in the batch.
 *
 * @return The number of dofs in the batch.
 */
size_t DofBatch::getSize() const{
	return dofList_.size();;
}
