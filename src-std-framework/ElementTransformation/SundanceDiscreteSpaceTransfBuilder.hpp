/*
 * SundanceDiscreteSpaceTransfBuilder.hpp
 *
 *  Created on: Mar 21, 2010
 *      Author: benk
 */

#ifndef SUNDANCEDISCRETESPACETRANSFBUILDER_HPP_
#define SUNDANCEDISCRETESPACETRANSFBUILDER_HPP_

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceDOFMapBase.hpp"
#include "SundanceHNDoFMapBase.hpp"
#include "SundanceTransformationHN.hpp"
#include "SundanceTransformationBase.hpp"

#include "TSFVectorType.hpp"
#include "TSFVectorDecl.hpp"

namespace Sundance {

using namespace Teuchos;
using namespace std;
using namespace TSFExtended;

/** This builds a transformation for one discrete space <br>
 *  The class also calls the transformation method of the created transformation object
 *  */
class DiscreteSpaceTransfBuilder {
public:

	/** Empty Ctor, in this case there will be no transformation*/
	DiscreteSpaceTransfBuilder();

	/** */
	DiscreteSpaceTransfBuilder(const Mesh& mesh, const BasisArray& basis,
			const RCP<DOFMapBase>& map );

	/** Dtor */
	virtual ~DiscreteSpaceTransfBuilder() {;}

	/** If there is a valid transformation then returns true, else false */
	const inline bool validTransformation() const { return hasTransformation_; }

	/** Function to make the local transformation for the discrete space
	 * @param dofs [in]
	 * @param funcID [in] (the number of the type of the function, not how many of them are)
	 * @param nrDoFs [in]
	 * @param ghostView [in] the array where we have to get values from
	 * @param localValues [out] */
	void getDoFsWithTransformation(const Array<int>& dofs,
			                             const int funcID ,
			                             const int nrDoFs ,
			                 			 const int cellDim ,
			                 			 const Array<int>& cellLIDs,
			                             const RCP<GhostView<double> >& ghostView ,
			                             Array<double>& localValues) const;

protected:

	/** for verbosity */
	int verb() const {return verb_;}

private:

	int verb_;

	/** Number of different function defined on the DoFMap*/
	int basisSize_;

	/** true if has a valid transformation */
	mutable bool hasTransformation_;

	/** true if has a valid transformation */
	mutable bool oneTransformation_;

	/** The transformation for the Discrete space */
	mutable RCP<TransformationBase> singleTransformation_;

	/** One transformation per space */
	//mutable Array< RCP<TransformationBase> > transformations_;
};

}

#endif /* SUNDANCEDISCRETESPACETRANSFBUILDER_HPP_ */
