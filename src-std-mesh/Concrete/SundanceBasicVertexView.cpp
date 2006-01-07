#include "SundanceBasicVertexView.hpp"
#include "Teuchos_Utils.hpp"

using namespace Teuchos;
using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;


string VertexView::toString() const
{
  int* ptr = *base_ +  offset_*length_;
	string rtn="{";
	for (int i=0; i<length_; i++) 
		{
			rtn += Teuchos::toString(ptr[i]);
			if (i < length_-1) rtn += ", ";
		}
	rtn += "}";
	return rtn;
}
