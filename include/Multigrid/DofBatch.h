/*
 * DofBatch.h
 *
 *  @date Mar 31, 2015
 *  @author Clemens Bartsch
 *
 *  Declaration of a class which holds degrees of freedom (pressure or velocity) to be used in a Vanka smoothing step.
 */

#ifndef DOFBATCH_H_
#define DOFBATCH_H_

// Standard includes.
#include<vector>
#include<stdlib.h>


/*!
 * @brief A DofBatch is a collection of degrees of freedom (pressure or velocity), identified by integer numbers.
 * Basically this class just wraps a vector<int> to use it more conveniently in classes derived from VankaSmoother.
 */
class DofBatch {

typedef std::vector<int> intVec;

public:

	//! Construct a dof batch, reserving space for a vector<int> of size nEntries.
	DofBatch(const size_t nDofs=1);

	//! Add a further dof to the end of the dofList.
	void addDof(const int newDof);

	//! Delegate an iterator to the class std::vector.
	intVec::const_iterator begin() const{ return dofList_.begin(); }
	intVec::const_iterator end() const { return dofList_.end(); }

	//! Return the number of the dof at place index.
	int getDof(const size_t index) const;

	//! Tidy up - removes duplicates, sorts the dofList_ and shrinks it to fit its contents.
	void tidyUp();

	//! Return the number of dofs stored in the batch.
	size_t getSize() const;


private:
	//! @brief In MoonMD dofs are stored by an int number. Store the dofs belonging to the batch in an std::vector<int>.
	intVec dofList_;

};


#endif /* DOFBATCH_H_ */
