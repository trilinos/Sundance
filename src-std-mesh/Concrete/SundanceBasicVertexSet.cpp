#include "SundanceBasicVertexSet.hpp"
#include "Teuchos_Utils.hpp"

using namespace Teuchos;
using namespace Sundance;
using namespace Sundance::StdMesh::Internal;


string VertexSet::toString() const
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
