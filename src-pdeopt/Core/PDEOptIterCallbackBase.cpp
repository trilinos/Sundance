#include "PDEOptIterCallbackBase.hpp"
#include "PDEOptPDEConstrainedObjBase.hpp"
#include "Sundance.hpp"

namespace Sundance
{

DefaultIterCallback::DefaultIterCallback(
  const std::string& filename, 
  const std::string& type,
  int frequency)
  : type_(type), filename_(filename), frequency_(frequency)
{}

void DefaultIterCallback::call(const PDEConstrainedObjBase* obj, 
  int iter) const
{
  if (iter % frequency_ != 0) return;
 
  string name = filename_ + "-iter-" + Teuchos::toString(iter);
  
  FieldWriter writer;

  if (type_=="VTK")
  {
    writer = new VTKWriter(name);
  }
  else if (type_=="Exodus")
  {
    writer = new ExodusWriter(name);
  }
  else if (type_=="Matlab")
  {
    writer = new MatlabWriter(name);
  }
  else 
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, RuntimeError, 
      "writer type [" << type_ << "] not defined");
  }

  Array<Expr> state = obj->stateVars();
  Array<Expr> adjoint = obj->adjointVars();
  Expr design = obj->designVar();

  writer.addMesh(obj->mesh());

  for (int b=0; b<state.size(); b++)
  {
    for (int i=0; i<state[b].size(); i++)
    {
      string tag = "[" + Teuchos::toString(b)
        + "][" + Teuchos::toString(i) + "]";
      writer.addField("state" + tag, new ExprFieldWrapper(state[b][i]));
      writer.addField("adjoint" + tag, new ExprFieldWrapper(adjoint[b][i]));
    }
  }

  for (int i=0; i<design.size(); i++)
  {
    string tag = "[" + Teuchos::toString(i) + "]";
    writer.addField("design" + tag, new ExprFieldWrapper(design[i]));
  }
  
  writer.write();
}


}
