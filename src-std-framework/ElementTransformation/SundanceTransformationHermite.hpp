/*
 * SundanceTransformationHermite.hpp
 *
 *  Created on: Mar 16, 2010
 *      Author: benk
 */

#ifndef SUNDANCETRANSFORMATIONHERMITE_HPP_
#define SUNDANCETRANSFORMATIONHERMITE_HPP_

#include "SundanceTransformationBase.hpp"

namespace Sundance {

class TransformationHermite: public Sundance::TransformationBase {
public:

	TransformationHermite();

	virtual ~TransformationHermite();

	/** */
	virtual void preApply(
			const int multiplicity ,
			const int vectorPerCell,
			const CellJacobianBatch& JTrans,
		    const CellJacobianBatch& JVol,
			const Array<int>& facetIndex,
			const RCP<Array<int> >& cellLIDs,
			RCP<Array<double> >& A
			) const {;}

	/** */
	virtual void postApply(
			const int multiplicity ,
			const int vectorPerCell,
			const CellJacobianBatch& JTrans,
		    const CellJacobianBatch& JVol,
			const Array<int>& facetIndex,
			const RCP<Array<int> >& cellLIDs,
			RCP<Array<double> >& A
			) const {;}

	/** */
	virtual void preApplyTranspose(
			const int multiplicity ,
			const int vectorPerCell,
			const int cellDim,
			const Array<int>& cellLIDs,
			const Array<int>& facetIndex,
			RCP<Array<double> >& A
			) const {;}
};

}

#endif /* SUNDANCETRANSFORMATIONHERMITE_HPP_ */
