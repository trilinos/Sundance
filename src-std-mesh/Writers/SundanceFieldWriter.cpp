/* @HEADER@ */
/* @HEADER@ */

#include "SundanceFieldWriter.hpp"
#include "SundanceExceptions.hpp"

using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace Teuchos;
using namespace TSFExtended;

void FieldWriter::addMesh(const Mesh& mesh) const
{
  ptr()->addMesh(mesh);
}

void FieldWriter::write() const
{
  ptr()->write();
}

void FieldWriter::addField(const string& name, 
                           const Handle<FieldBase>& field)
{
  ptr()->addField(name, field.ptr());
}
