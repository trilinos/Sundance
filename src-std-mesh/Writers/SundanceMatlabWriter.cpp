/* @HEADER@ */
/* @HEADER@ */

#include "SundanceMatlabWriter.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "Teuchos_XMLObject.hpp"



using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace Teuchos;
using namespace TSFExtended;



void MatlabWriter::write() const 
{
  ofstream os(filename().c_str());
  int numNodes = mesh().numCells(0);
          
  for (int i=0; i<numNodes; i++)
    {
      os << mesh().nodePosition(i);
      for (int j=0; j<pointScalarFields().size(); j++)
        {
          const RefCountPtr<FieldBase>& expr = pointScalarFields()[j];
          os << " " << expr->getData(0, i, 0);
        }
      os << endl;
    }
}
