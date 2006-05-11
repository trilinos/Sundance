#ifndef PYSUNDANCE_FIATQUADRATUREADAPTER_H
#define PYSUNDANCE_FIATQUADRATUREADAPTER_H

#include "SundanceQuadratureFamilyBase.hpp"
#include "Python.h"

namespace SundanceStdFwk
{
	using namespace SundanceUtils;
	using namespace SundanceStdMesh;
	using namespace SundanceStdMesh::Internal;
	using namespace Internal;
	using namespace SundanceCore;
	using namespace SundanceCore::Internal;
	using namespace Teuchos;

	class FIATQuadratureAdapter : public QuadratureFamilyBase
	{
	public:
		FIATQuadratureAdapter( PyObject *py_quad_factory , int order );
		virtual ~FIATQuadratureAdapter();
		virtual void getPoints( const CellType& cellType ,
								Array<Point>& quadPoints ,
								Array<double>& quadWeights ) const;
	private:
		// constructor will pretabulate for line, triangle, and tet and
		// store in these arrays
		Array<Array<Point> > pts_;
		Array<Array<double> > wts_;
	};

}

#endif