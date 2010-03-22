/*
 * SundanceDiscreteSpaceTransfBuilder.cpp
 *
 *  Created on: Mar 21, 2010
 *      Author: benk
 */


#include "SundanceDiscreteSpaceTransfBuilder.hpp"


using namespace Sundance;


DiscreteSpaceTransfBuilder::DiscreteSpaceTransfBuilder():
verb_(0),
basisSize_(0),
hasTransformation_(false) , oneTransformation_(false) ,
singleTransformation_()
//, transformations_()
{
}

DiscreteSpaceTransfBuilder::DiscreteSpaceTransfBuilder(
		const Mesh& mesh, const BasisArray& basis,
		const RCP<DOFMapBase>& map):
        verb_(0),
        basisSize_(0),
        hasTransformation_(false) , oneTransformation_(false) ,
        singleTransformation_()
        //, transformations_()
        {

	if (mesh.allowsHangingHodes()){
	    // in this case we have to create a transformation
   		const HNDoFMapBase* myMap
    		    = dynamic_cast<const HNDoFMapBase*>(map.get());
   		if (myMap != 0){
   		   // create one hanging node transformation
   		   SUNDANCE_MSG2( verb() , "DiscreteSpaceTransfBuilder::DiscreteSpaceTransfBuilder basis.size():" << basis.size());
   		   basisSize_ = basis.size();
   		   singleTransformation_ = rcp((TransformationBase*)(new TransformationHN( myMap , 1 , basis.size() ))); //todo: basis.size is wrong
		   hasTransformation_ = true;
		   oneTransformation_ = true;
   		}
	}
	else
	{ // e.g. create Hermite transformation
	}
}

void DiscreteSpaceTransfBuilder::getDoFsWithTransformation(
		                             const Array<int>& dofs,
		                             const int funcID ,
		                             const int nrDoFs ,
		                 			 const int cellDim ,
		                 			 const Array<int>& cellLIDs,
		                             const RCP<GhostView<double> >& ghostView ,
		                             Array<double>& localValues) const {

	if (oneTransformation_){
		   // todo: We assume here that "funcID" does not matter, but in case of MixedDoFMap it matters

		   // get all the DoFs for this functionID
           int nrFunctions = basisSize_;
           int DoFPerElement = dofs.size() / (basisSize_*cellLIDs.size());
           Array<double> tmpArray(DoFPerElement*cellLIDs.size());
           Array<int> facetIndex(1); //just needed for the interface
           int cellI , elemDof;

           SUNDANCE_MSG2( verb() , "DiscreteSpaceTransfBuilder::getDoFsWithTransformation()  nrFunctions:" << nrFunctions
        		   << " DoFPerElement:" << DoFPerElement <<  "  localValues.size():" << localValues.size());
           SUNDANCE_MSG2( verb() , "DiscreteSpaceTransfBuilder::getDoFsWithTransformation() localValues:" << localValues);

           for (int nf = 0 ; nf < nrFunctions ; nf++)
           {
        	   // copy the elements into the temporary array
        	   for(cellI = 0 ; cellI < cellLIDs.size() ; cellI++) {
        		   for (elemDof = 0 ; elemDof < DoFPerElement ; elemDof++)
        		   tmpArray[cellI*DoFPerElement + elemDof] = localValues[(nrFunctions*cellI+nf)*DoFPerElement + elemDof];
        	   }
        	   SUNDANCE_MSG2( verb() ,"getDoFsWithTransformation() before Transformation:" << tmpArray);
		       // make the transformation for all elements once
			   singleTransformation_->preapplyTranspose(
			   	     1 , //const int multiplicity ,
					 1 , //const int vectorPerCell,
					 cellDim, cellLIDs,
					 facetIndex,
					 tmpArray );
			   // copy the elements back
			   SUNDANCE_MSG2( verb() , "getDoFsWithTransformation() after Transformation:" << tmpArray );
        	   for(cellI = 0 ; cellI < cellLIDs.size() ; cellI++) {
        		   for (elemDof = 0 ; elemDof < DoFPerElement ; elemDof++)
        		   localValues[(nrFunctions*cellI+nf)*DoFPerElement + elemDof] = tmpArray[cellI*DoFPerElement + elemDof];
        	   }
           }
	}
	else{
		// todo: we have might have different transformation for each space
		// index is funcID
	}
}
