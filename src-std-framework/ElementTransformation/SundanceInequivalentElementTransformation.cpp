#include "SundanceInequivalentElementTransformation.hpp"

namespace Sundance
{
  InequivalentElementTransformation::InequivalentElementTransformation( const Mesh& mesh ,
									const MixedDOFMap *map ):
    mesh_( mesh ), map_( map ), chunkBases_( map->nBasisChunks() )
  {
    // extract all the bases from the dof map
    for (int i=0;i<map->nBasisChunks();i++) 
      {
	chunkBases_[i] = dynamic_cast<const BasisFamilyBase *>(map->basis(i).get());
      }

  }

  void InequivalentElementTransformation::preApply( const int funcID,
						    const CellJacobianBatch& JTrans,
						    const CellJacobianBatch& JVol,
						    const Array<int>& facetIndex,
						    const RCP<Array<int> >& cellLIDs,
						    RCP<Array<double> >& A
						    ) const
  
  {
    chunkBases_[map_->chunkForFuncID( funcID )]->preApplyTransformation( mesh_.cellType( mesh_.spatialDim() ) , JVol , A );
  }

  void InequivalentElementTransformation::postApply( const int funcID,
						     const CellJacobianBatch& JTrans,
						     const CellJacobianBatch& JVol,
						     const Array<int>& facetIndex,
						     const RCP<Array<int> >& cellLIDs,
						     RCP<Array<double> >& A
						     ) const
  {
    chunkBases_[map_->chunkForFuncID( funcID )]->postApplyTransformation( mesh_.cellType( mesh_.spatialDim() ) , JVol , A );
  }
  
  void InequivalentElementTransformation::preapplyTranspose( const int cellDim,
							     const int funcID,
							     const Array<int>& cellLIDs,
							     const Array<int>& facetIndex,
							     Array<double>& A ) const
  {
    CellJacobianBatch JVol;
    mesh_.getJacobians( mesh_.spatialDim() , cellLIDs , JVol );
    chunkBases_[map_->chunkForFuncID( funcID )]->preApplyTransformationTranspose( mesh_.cellType( mesh_.spatialDim() ) , JVol , A );
  }
  

}
