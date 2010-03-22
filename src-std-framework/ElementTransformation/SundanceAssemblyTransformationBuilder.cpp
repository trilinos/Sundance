/*
 * SundanceAssemblyTransformationBuilder.cpp
 *
 *  Created on: Mar 16, 2010
 *      Author: benk
 */

#include "SundanceAssemblyTransformationBuilder.hpp"
#include "SundanceHNDoFMapBase.hpp"

namespace Sundance {

/* for Assembly transformations*/
AssemblyTransformationBuilder::AssemblyTransformationBuilder( const RCP<IntegralGroup>& group ,
	int groupIndex ,
	int entryPerCell ,
    CellType cellType ,
    CellType maxCellType ,
	const Array<RCP<DOFMapBase> >& rowMaps ,
	const Array<RCP<DOFMapBase> >& colMaps ,
	const Mesh& mesh):
  verb_(0),
  nrCol_(group->nUnkNodes()) ,
  nrRow_(group->nTestNodes()) ,
  preTransformation_(0),
  postTransformation_(0),
  _myRowDOFMap(0),
  _myColDOFMap(0),
  hasTransformation_(false),
  onlyVectorTransformation_(true)
{
	SUNDANCE_MSG2(verb(), "AssemblyTransformationBuilder::AssemblyTransformationBuilder initialized fields:" << mesh.allowsHangingHodes());
	// make different transformations
    if (mesh.allowsHangingHodes()){

    	// if the integration is not on a MaxCell then no HangingNode Trafo!
    	// caz in some cases this would incease the size of the elem matrix
    	if (cellType != maxCellType){
        	hasTransformation_ = false;
    		return;
    	}

    	// we have HANGING NODES, create a corresponding transformation
    	hasTransformation_ = true;
    	if (group->nUnkNodes() > 0){
    		SUNDANCE_MSG2(verb(), "Matrix trafo colMaps.size:" << colMaps.size());
    		SUNDANCE_MSG2(verb(), "Matrix trafo rowMaps.size:" << rowMaps.size());
    		SUNDANCE_MSG2(verb(), "Matrix trafo group->unkBlock():" << group->unkBlock());
    		SUNDANCE_MSG2(verb(), "Matrix trafo group->testBlock():"<< group->testBlock());
    	    // this means we have to multiply Matrix
    		onlyVectorTransformation_ = false;
    		_myColDOFMap = colMaps[group->unkBlock()[0]].get(); // group->unkBlock(), has alway one element ;-)
    		_myRowDOFMap = rowMaps[group->testBlock()[0]].get(); // group->testBlock(), has alway one element ;-)
    		SUNDANCE_MSG2(verb(), "Trafo the Maps to HN Maps" );
    		// create the NH transformation, for pre and post apply
    		const HNDoFMapBase* colMap
    		    = dynamic_cast<const HNDoFMapBase*>(_myColDOFMap);
    		const HNDoFMapBase* rowMap
    		    = dynamic_cast<const HNDoFMapBase*>(_myRowDOFMap);

    		if (colMap != 0)
    		{
    			postTransformation_ = rcp((TransformationBase*)(new TransformationHN( colMap , nrCol_ , nrRow_ )));
    	    }else{
    	    	SUNDANCE_ERROR("AssemblyTransformationBuilder::AssemblyTransformationBuilder, wrong configuration , no Col DoFMap");
    		}

    		if (rowMap != 0)
    		{
    			preTransformation_ = rcp(new TransformationHN( rowMap , nrCol_ , nrRow_ ));
    		}else{
    	    	SUNDANCE_ERROR("AssemblyTransformationBuilder::AssemblyTransformationBuilder, wrong configuration , no Row DoFMap");
    		}
    	} else{
    	    // this means we have to multiply Vector
            if (group->nTestNodes() > 0){
            	SUNDANCE_MSG2(verb(), "Vector trafo colMaps.size:" << colMaps.size());
            	SUNDANCE_MSG2(verb(), "Vector trafo rowMaps.size:" << rowMaps.size());
            	SUNDANCE_MSG2(verb(), "Vector trafo group->unkBlock():" << group->unkBlock());
            	SUNDANCE_MSG2(verb(), "Vector trafo group->testBlock():"<< group->testBlock());
        		_myRowDOFMap = rowMaps[group->testBlock()[0]].get(); // group->testBlock(), has alway one element ;-)

        		// create the NH transformation, only for pre apply
            	const HNDoFMapBase* rowMap
            		    = dynamic_cast<const HNDoFMapBase*>(_myRowDOFMap);

            	if (rowMap != 0)
            	{
            		preTransformation_ = rcp(new TransformationHN( rowMap , nrCol_ , nrRow_ ));
            	}else{
        	    	SUNDANCE_ERROR(" ");
            	}
            }
            else{
            	hasTransformation_ = false;
            }
    	}
    }
    else{
     // Todo: check for Hermite transformation
     // ask the "group" for their basis, and see if it is Hermite, then create fransformaion
    	hasTransformation_ = false;
    }
}

AssemblyTransformationBuilder::~AssemblyTransformationBuilder() {

}

void AssemblyTransformationBuilder::applyTransformsToAssembly(
		const CellJacobianBatch& JTrans,
	    const CellJacobianBatch& JVol,
	    const Array<int>& facetNum,
	    const RCP<Array<int> >& cellLIDs,
	    RCP<Array<double> >& A){

	if (hasTransformation_){
		// do the transformation with the pre and post transformation object
		if (onlyVectorTransformation_){
			preTransformation_->preApply(
					1,       //multiplicity , not used, not needed now
					1,      // vectorPerCell, not used, not needed now
					JTrans, JVol, facetNum, cellLIDs, A);
		}else{
			preTransformation_->preApply(
					1,       //multiplicity , not used, not needed now
					1,      // vectorPerCell, not used, not needed now
					JTrans, JVol, facetNum, cellLIDs, A);
			postTransformation_->postApply(
					1,       //multiplicity , not used, not needed now
					1,      // vectorPerCell, not used, not needed now
					JTrans, JVol, facetNum, cellLIDs, A);

		}
	}
	SUNDANCE_MSG2( verb() , " AssemblyTransformationBuilder::applyTransformsToAssembly ");
}

}
