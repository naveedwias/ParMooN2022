/*
 * DofPatch.C
 *
 *  @date Mar 31, 2015
 *  @author Clemens Bartsch
 *
 *  Implementation of class DofPatch declared in include/FE/DofPatch.h.
 *
 */

#include <DofPatch.h>
#include <algorithm>
#include <vector>

/*!
 * @brief Construct a DofPatch with known size.
 *
 * @param[in] nEntries The number of dofs which will be added later. Defaults to 1.
 */
DofPatch::DofPatch(const size_t nDofs) :
dofList_()
{
	//Reserve space for the entries to be added later.
	dofList_.reserve(nDofs);
}

/*!
 * @brief Add a dof to the patch.
 *
 * @param[in] newDof The number of the dof to be added to the patch.
 */
void DofPatch::addDof(const int newDof){
	dofList_.push_back(newDof);
}

/*!
 * @brief Return the number of the dof at place 'index'.
 * @return The MooNMD internal number of the dof.
 *
 * @param[in] index The place of the dof in the patch.
 */
int DofPatch::getDof(const size_t index) const{
	return dofList_.at(index);
}

/*!
 * @brief Removes duplicates and sorts the dofList_.
 *
 * This must be called when the procedure of determining which dofs to add and adding them is finished.
 */
void DofPatch::tidyUp(){
	// Sort the entries.
	std::sort(dofList_.begin(), dofList_.end());

	// Throw out duplicates, by means of std methods.
	intVec::iterator it = std::unique(dofList_.begin(), dofList_.end());
	dofList_.resize( std::distance(dofList_.begin(),it) );

	// Adjust capacity.
	dofList_.shrink_to_fit();
}

/*!
 * @brief Return the number of dofs stored in the patch.
 *
 * @return The number of dofs in the patch.
 */
size_t DofPatch::getSize() const{
	return dofList_.size();;
}
