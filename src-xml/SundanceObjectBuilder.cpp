#include "Sundance.hpp"
#include "SundanceObjectBuilder.hpp"
#include "SundanceExprParser.hpp"

using namespace SundanceXML;
using SundanceCore::List;


void ObjectBuilder::interpret(const XMLObject& xml)
{
  const string& tag = xml.getTag();

  if (tag=="Include")
    {
      const string& filename = xml.getRequired("filename");
      interpret(filename);
    }
  else if (ExprParser::basisTypes().contains(tag))
    {
      createBasis(xml);
    }
  else if (ExprParser::meshTypes().contains(tag))
    {
      createMesh(xml);
    }
  else if (ExprParser::quadTypes().contains(tag))
    {
      createQuadrature(xml);
    }
  else if (ExprParser::cellFilterTypes().contains(tag))
    {
      createMesh(xml);
    }
  else if (ExprParser::exprTypes().contains(tag))
    {
      createExpr(xml);
    }
  else if (tag=="DiscreteSpace")
    {
      createDiscreteSpace(xml);
    }
  else if (tag=="LinearProblem")
    {
      createLinearProblem(xml);
    }
}

Expr ObjectBuilder::createExpr(const XMLObject& xml)
{
  Expr rtn;

  const string& tag = xml.getTag();
  string name;
  if (xml.hasAttribute("name")) name = xml.getRequired("name");
  
  if (tag=="Expr")
    {
      
    }
  else if (tag=="Derivative")
    {
      int direction = xml.getRequiredInt("order");
      rtn = new Derivative(direction);
    }
  else if (tag=="CoordExpr")
    {
      int direction = xml.getRequiredInt("order");
      rtn = new CoordExpr(direction, name);
    }
  else if (tag=="UnknownFunction")
    {
      TEST_FOR_EXCEPTION(xml.hasAttribute("basis") && xml.numChildren() != 0,
                         RuntimeError,
                         "redundant basis specification in " << xml);
      TEST_FOR_EXCEPTION(!xml.hasAttribute("basis") && xml.numChildren() == 0,
                         RuntimeError,
                         "no basis specification in " << xml);
      BasisFamily myBasis;
      if (xml.hasAttribute("basis"))
        {
          const string& basisName = xml.getRequired("basis");
          myBasis = basis(basisName);
        }
      else
        {
          myBasis = createBasis(xml.getChild(0));
        }
      return new UnknownFunction(myBasis, name);
    }
  else if (tag=="TestFunction")
    {
      TEST_FOR_EXCEPTION(xml.hasAttribute("basis") && xml.numChildren() != 0,
                         RuntimeError,
                         "redundant basis specification in " << xml);
      TEST_FOR_EXCEPTION(!xml.hasAttribute("basis") && xml.numChildren() == 0,
                         RuntimeError,
                         "no basis specification in " << xml);
      BasisFamily myBasis;
      if (xml.hasAttribute("basis"))
        {
          const string& basisName = xml.getRequired("basis");
          myBasis = basis(basisName);
        }
      else
        {
          myBasis = createBasis(xml.getChild(0));
        }
      return new TestFunction(myBasis, name);
    }
  else if (tag=="Var")
    {
      rtn = expr(xml.getRequired("name"));
    }
  else if (tag=="Integral")
    {
    }
  else if (tag=="EssentialBC")
    {
    }
  else if (tag=="Sum")
    {
      Expr left = createExpr(xml.getChild(0));
      Expr right = createExpr(xml.getChild(1));
      rtn =  left + right;
    }
  else if (tag=="Product")
    {
      Expr left = createExpr(xml.getChild(0));
      Expr right = createExpr(xml.getChild(1));
      rtn = left * right;
    }
  else if (tag=="Reciprocal")
    {
      Expr arg = createExpr(xml.getChild(0));
      rtn = 1.0/arg;
    }
  else if (tag=="Function")
    {
      //    rtn = createFunction(xml);
    }
  else if (tag=="UnaryMinus")
    {
      Expr arg = createExpr(xml.getChild(0));
      rtn = -arg;
    }
}

LinearProblem ObjectBuilder::createLinearProblem(const XMLObject& xml)
{
  cerr << "dummy LP builder for " << xml << endl;
}

DiscreteSpace ObjectBuilder::createDiscreteSpace(const XMLObject& xml)
{
  cerr << "dummy discrete space builder for " << xml << endl;
}

Mesh ObjectBuilder::createMesh(const XMLObject& xml)
{
  const string& tag = xml.getTag();
  const string& name = xml.getRequired("name");
  
  MeshSource src;

  if (tag=="Exodus")
    {
      const string& filename = xml.getRequired("file");
      src = new ExodusNetCDFMeshReader(filename, 
                                       new BasicSimplicialMeshType());
    }
  else if (tag=="Triangle")
    {
      const string& filename = xml.getRequired("file");
      src = new TriangleMeshReader(filename, 
                                   new BasicSimplicialMeshType());
    }
  else if (tag=="Rectangle")
    {
      double ax = xml.getRequiredDouble("ax");
      double bx = xml.getRequiredDouble("bx");
      double ay = xml.getRequiredDouble("ay");
      double by = xml.getRequiredDouble("by");
      int nx = xml.getRequiredInt("nx");
      int ny = xml.getRequiredInt("ny");
      int npx = xml.getRequiredInt("npx");
      int npy = xml.getRequiredInt("npy");
      
      src = new PartitionedRectangleMesher(ax, bx, nx, npx,
                                           ay, by, ny, npy, 
                                           new BasicSimplicialMeshType());
    }
  else if (tag=="Line")
    {
      double ax = xml.getRequiredDouble("ax");
      double bx = xml.getRequiredDouble("bx");
      int nx = xml.getRequiredInt("nx");
      src = new PartitionedLineMesher(ax, bx, nx, 
                                      new BasicSimplicialMeshType());
    }
  else
    {
      TEST_FOR_EXCEPTION(true, RuntimeError,
                         "unregognized mesh source type " << tag
                         << " in ObjectBuilder::createMesh()");
                         
    }
  
  Mesh mesh = src.getMesh();
  
  mesh_.put(name, mesh);

  return mesh;
}


BasisFamily ObjectBuilder::createBasis(const XMLObject& xml) 
{
  const string& tag = xml.getTag();
  const string& name = xml.getRequired("name");
  
  BasisFamily basis;
  if (tag=="Lagrange")
    {
      int order = xml.getRequiredInt("order");
      basis = new Lagrange(order);
    }
  else
    {
      TEST_FOR_EXCEPTION(true, RuntimeError,
                         "unregognized basis type " << tag
                         << " in ObjectBuilder::createBasis()");
    }
  
  basis_.put(name, basis);
  return basis;
}

QuadratureFamily ObjectBuilder::createQuadrature(const XMLObject& xml) 
{
  const string& tag = xml.getTag();
  const string& name = xml.getRequired("name");
  
  QuadratureFamily quadrature;
  if (tag=="GaussianQuadrature")
    {
      int order = xml.getRequiredInt("order");
      quadrature = new GaussianQuadrature(order);
    }
  else
    {
      TEST_FOR_EXCEPTION(true, RuntimeError,
                         "unregognized quadrature type " << tag
                         << " in ObjectBuilder::createQuadrature()");
    }
  
  quadrature_.put(name, quadrature);
  return quadrature;
}


CellFilter ObjectBuilder::createCellFilter(const XMLObject& xml) 
{
  const string& tag = xml.getTag();
  const string& name = xml.getRequired("name");

  CellFilter filter;

  if (tag=="MaximalCellFilter")
    {
      filter = new MaximalCellFilter();
    }
  else if (tag=="BoundaryCellFilter")
    {
      filter = new BoundaryCellFilter();
    }
  else if (tag=="LabeledCellFilter")
    {
      int label = xml.getRequiredInt("label");
      const string& superName = xml.getRequired("super");
      const CellFilter& super = cellFilter(superName);
      filter = super.labeledSubset(label);
    }
  
  filter_.put(name, filter);
  return filter;
}






