/* @HEADER@ */
/* @HEADER@ */

#include "SundanceFieldWriterBase.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"


using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace Teuchos;
using namespace TSFExtended;

FieldWriterBase::FieldWriterBase(const string& filename) 
  : filename_(filename),
    nProc_(0), 
    myRank_(-1),
    meshID_(-1),
    mesh_(),
    comments_(),
    pointScalarFields_(),
    cellScalarFields_(),
    pointScalarNames_(),
    cellScalarNames_(),
    pointVectorFields_(),
    cellVectorFields_(),
    pointVectorNames_(),
    cellVectorNames_()
{;}


void FieldWriterBase::impersonateParallelProc(int nProc, int rank)
{
  nProc_ = nProc;
  myRank_ = rank;
}

int FieldWriterBase::nProc() const
{
  if (nProc_ < 1) return mesh().comm().getNProc();
  return nProc_;
}

int FieldWriterBase::myRank() const
{
  if (myRank_ < 0) return mesh().comm().getRank();
  return myRank_;
}




void FieldWriterBase::addMesh(const Mesh& mesh) 
{
  if (meshID_ < 0)
    {
      mesh_ = mesh;
      meshID_ = mesh.id();
    }
                     
  TEST_FOR_EXCEPTION(meshID_ != mesh.id(), RuntimeError,
                     "FieldWriterBase::setMesh(): inconsistent meshes: "
                     "existing mesh has meshID=" << meshID_ << ", newly "
                     "added mesh has meshID=" << mesh.id());
}

  void FieldWriterBase::addField(const string& /* name */, 
                               const RefCountPtr<FieldBase>& /* expr */) 
{

//   string fieldName = name;
//   if (name.length() == 0) fieldName = expr.name();

//   if (expr.length() > 1)
//     {
//       /* expr is a list of fields */
//       for (int i=0; i<expr.length(); i++)
//         {
//           string n = name + "_" + TSF::toString(i);
//           addField(n, expr[i]);
//         }
//     }
//   else 
//     {
//       /* expr is a single scalar field */

//       /* get the mesh. Check for consistency with previously obtained meshes */
//       Mesh m;
//       expr.getMesh(m);
//       setMesh(m);

//       /* Now determine whether it is a nodal field or a cell field */
//       BasisFamily basis = expr.getBasis();
//       if (basis.order() > 0)
//         {
//           pointScalarFields_.append(expr);
//           pointScalarNames_.append(fieldName);
//         }
//       else
//         {
//           cellScalarFields_.append(expr);
//           cellScalarNames_.append(fieldName); 
//         }
//     }
  
}

void FieldWriterBase::addCommentLine(const string& line) 
{
  comments_.append(line);
}






