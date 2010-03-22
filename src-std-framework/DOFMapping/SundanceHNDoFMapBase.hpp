/*
 * SundanceHNMapBase.hpp
 *
 *  Created on: Mar 18, 2010
 *      Author: benk
 */

#ifndef SUNDANCEHNMAPBASE_HPP_
#define SUNDANCEHNMAPBASE_HPP_

#include "SundanceDefs.hpp"
#include "SundanceSpatiallyHomogeneousDOFMapBase.hpp"
namespace Sundance
{
using namespace Teuchos;

/**
 * The abstract class which extends the functionalities of the DOF map <br>
 * The only additional functionality is that we have a restriction on the DOFs
 * , with the pre-fill transformations these constraints can be build in into the matrix
 *
 */
class HNDoFMapBase : public SpatiallyHomogeneousDOFMapBase
{
public:

	/** Empty Ctor */
	HNDoFMapBase(const Mesh& mesh, int nFuncs,
		    const ParameterList& verbParams ):
		    		SpatiallyHomogeneousDOFMapBase(mesh, nFuncs, verbParams){;}

	/**
	 * @param cellLID [in] the maxCell LID input
	 * @param trafoMatrixSize [in/out]
	 * @param doTransform [out]
	 * @param transfMatrix [out] (we assume that the array is already pre-sized )*/
	  virtual void getTrafoMatrixForCell(
		    int cellLID,
		    int& trafoMatrixSize,
		    bool& doTransform,
		    Array<double>& transfMatrix ) const=0;

	  /** Function used for plotting for hanging node DOFMaps
	   *  Returns for one hanging node (element) the global DoFs which contribute to that hanging local DoF
	   * @param cellDim [in] the dimension
	   * @param cellLID [in] the LID of the cell
	   * @param funcID  [in] the function ID, (to wchich the DOFs belong)
	   * @param dofs    [out] the global DoF s
	   * @param coefs   [out] the coefficient of each global DoF */
	  virtual void getDOFsForCell(
			int cellDim,
			int cellLID,
	        int funcID,
	        Array<int>& dofs ,
	        Array<double>& coefs ) const=0;

	  /** Returns the dimension where the DoF map is defined
	   *  For HN we do transformation only for the maxCell type */
	  int getSpacialMeshDim() const { return mesh().spatialDim();}

protected:

	  /** Function used in the transformation matrix creation, where
	   * we need to insert one node LID into the element list,
	   * If the LID is already in the list then retrurn the position, otherwise add to the end <br>
	   * (functionality similar to the Set)
	   * */
	  static inline void returnInsertionPosition(Array<int>& fLIDs , int fLID , int& insertPos){
		  int pos;
		  for(pos = 0 ; pos < fLIDs.size() ; pos++)
		  {
			if (fLIDs[pos] == fLID){
				insertPos = pos; // we found, nothing to do now
				break;
			}
		  }
		  if (pos >= fLIDs.size()){ // add to the end
			  insertPos = pos;
			  fLIDs.resize(fLIDs.size()+1);
			  fLIDs[insertPos] = fLID;
		  }
	  }

};

}


#endif /* SUNDANCEHNMAPBASE_HPP_ */
