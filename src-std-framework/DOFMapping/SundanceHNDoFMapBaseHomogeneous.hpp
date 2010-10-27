/*
 * SundanceHNMapBase.hpp
 *
 *  Created on: Mar 18, 2010
 *      Author: benk
 */

#ifndef SUNDANCEHNMAPBASEHOMOGENEOUS_HPP_
#define SUNDANCEHNMAPBASEHOMOGENEOUS_HPP_

#include "SundanceDefs.hpp"
#include "SundanceSpatiallyHomogeneousDOFMapBase.hpp"
#include "SundanceHNDoFMapBase.hpp"
namespace Sundance
{
using namespace Teuchos;

/**
 * The abstract class which extends the functionalities of the DOF map <br>
 * The only additional functionality is that we have a restriction on the DOFs
 * , with the pre-fill transformations these constraints can be build in into the matrix
 *
 */
class HNDoFMapBaseHomogeneous : public SpatiallyHomogeneousDOFMapBase , public HNDoFMapBase
{
public:

	/** Empty Ctor */
	HNDoFMapBaseHomogeneous(const Mesh& mesh, int nFuncs, int setupVerb) :
		    		SpatiallyHomogeneousDOFMapBase(mesh, nFuncs, setupVerb),
		    		HNDoFMapBase(mesh, nFuncs, setupVerb) {;}

	virtual ~HNDoFMapBaseHomogeneous() {;}

protected:

};

}


#endif /* SUNDANCEHNMAPBASEHOMOGENEOUS_HPP_ */
