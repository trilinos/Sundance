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

void FieldWriterBase::addField(const string& name, 
                               const RefCountPtr<FieldBase>& expr) 
{

  string fieldName = name;

  if (expr->numElems() > 1)
    {
      TEST_FOR_EXCEPTION(expr->numElems() > 1, RuntimeError,
                         "FieldWriterBase::addField not ready for vector fields");
    }
  else 
    {
      /* expr is a single scalar field */
      
      pointScalarFields_.append(expr);
      pointScalarNames_.append(fieldName);
    }
  
}

void FieldWriterBase::addCommentLine(const string& line) 
{
  comments_.append(line);
}






