#include "Sundance.hpp"
#include "SundanceExprParser.hpp"

using namespace SundanceXML;
using SundanceCore::List;

int main(int argc, void** argv)
{
  
  try
		{
      Sundance::init(&argc, &argv);

      // Interpreter interpreter;

//       string fqFile = Sundance::searchForFile("test.xml");
  
//       FileInputSource fis(fqFile);
      
//       XMLObject xml = fis.getObject();

//       interpreter.interpret(xml);

      Array<string> lines;
      
      lines.append("d=a+b*c;");

      lines.append("a+b*c;");

      lines.append("u = UnknownFunction(Lagrange(1));");

      lines.append("u = CoordExpr(0);");

      lines.append("mesh = Exodus('fred');");

      lines.append("q4 = GaussianQuadrature(4);");

      lines.append("mesh = Line(0.0, 1.0, 10);");


      for (int i=0; i<lines.size(); i++)
        {
          cerr << endl << endl ;
          cerr << "----------- parsing " << lines[i] << endl;
          cerr << endl << endl 
               << ExprParser::evaluateAssignment(lines[i]).toString() << endl;
        }
      
    }
	catch(std::exception& e)
		{
      cerr << e.what() << endl;
		}
  Sundance::finalize();
}


// void Interpreter::interpret(const XMLObject& xml)
// {
//   const string& tag = xml.getTag();

//   if (tag=="Include")
//     {
//       const string& filename = xml.getRequired("filename");
//       interpret(filename);
//     }
//   else if (tag=="Solver")
//     {
      
//     }
//   else if (tag=="CellFilter")
//     {
      
//     }
//   else if (tag=="Define")
//     {
      
//     }

//   for (int i=0; i<xml.numChildren(); i++)
//     {
//       const XMLObject& child = xml.getChild(i);
      
//     }
// }


// void Interpreter::createMesh(const XMLObject& xml)
// {
//   checkTag(xml, "Mesh");
//   const string& name = xml.getRequired("name");
  
//   TEST_FOR_EXCEPTION(xml.numChildren() != 1, RuntimeError,
//                      "Interpreter::createMesh() expected exactly one child "
//                      "in " << xml.toString());

//   const XMLObject& child = xml.getChild(0);
//   const string& tag = child.getTag();

//   MeshSource src;

//   if (tag=="Exodus")
//     {
//       const string& filename = child.getRequired("file");
//       src = new ExodusNetCDFMeshReader(filename, 
//                                        new BasicSimplicialMeshType());
//     }
//   else if (tag=="Triangle")
//     {
//       const string& filename = child.getRequired("file");
//       src = new TriangleMeshReader(filename, 
//                                    new BasicSimplicialMeshType());
//     }
//   else if (tag=="Rectangle")
//     {
//       double ax = child.getRequiredDouble("ax");
//       double bx = child.getRequiredDouble("bx");
//       double ay = child.getRequiredDouble("ay");
//       double by = child.getRequiredDouble("by");
//       int nx = child.getRequiredInt("nx");
//       int ny = child.getRequiredInt("ny");
//       int npx = child.getRequiredInt("npx");
//       int npy = child.getRequiredInt("npy");
      
//       src = new PartitionedRectangleMesher(ax, bx, nx, npx,
//                                            ay, by, ny, npy, 
//                                            new BasicSimplicialMeshType());
//     }
//   else if (tag=="Line")
//     {
//       double ax = child.getRequiredDouble("ax");
//       double bx = child.getRequiredDouble("bx");
//       int nx = child.getRequiredInt("nx");
//       src = new PartitionedLineMesher(ax, bx, nx, 
//                                       new BasicSimplicialMeshType());
//     }
//   else
//     {
//       TEST_FOR_EXCEPTION(true, RuntimeError,
//                          "unregognized mesh source type " << tag
//                          << " in Interpreter::createMesh()");
                         
//     }
  
//   Mesh mesh = src.getMesh();
  
//   mesh_.put(name, mesh);
// }


// void Interpreter::createBasis(const XMLObject& xml)
// {
//   checkTag(xml, "BasisFamily");
//   const string& name = xml.getRequired("name");
  
//   TEST_FOR_EXCEPTION(xml.numChildren() != 0, RuntimeError,
//                      "Interpreter::createBasis() expected exactly zero "
//                      "children in " << xml.toString());

//   const string& type = xml.getRequired("type");

//   BasisFamily basis;
//   if (type=="Lagrange")
//     {
//       int order = xml.getRequiredInt("order");
//       basis = new Lagrange(order);
//     }
//   else
//     {
//       TEST_FOR_EXCEPTION(true, RuntimeError,
//                          "unregognized basis type " << type
//                          << " in Interpreter::createBasis()");
//     }
  
//   basis_.put(name, basis);
// }


// void Interpreter::createQuad(const XMLObject& xml)
// {
//   checkTag(xml, "QuadatureFamily");
//   const string& name = xml.getRequired("name");
  
//   TEST_FOR_EXCEPTION(xml.numChildren() != 0, RuntimeError,
//                      "Interpreter::createQuad() expected exactly zero "
//                      "children in " << xml.toString());

//   const string& type = xml.getRequired("type");

//   QuadratureFamily quad;
//   if (type=="Gaussian")
//     {
//       int order = xml.getRequiredInt("order");
//       quad = new GaussianQuadrature(order);
//     }
//   else
//     {
//       TEST_FOR_EXCEPTION(true, RuntimeError,
//                          "unregognized quadrature type " << type
//                          << " in Interpreter::createQuad()");
//     }
  
//   quad_.put(name, quad);
// }


// void Interpreter::createMaximalCellFilter(const XMLObject& xml)
// {
//   checkTag(xml, "MaximalCells");
//   const string& name = xml.getRequired("name");
  
//   TEST_FOR_EXCEPTION(xml.numChildren() != 0, RuntimeError,
//                      "Interpreter::createMaximalCellFilter() expected exactly zero "
//                      "children in " << xml.toString());

//   CellFilter filter = new MaximalCellFilter();

//   cellFilter_.put(name, filter);
// }


// void Interpreter::createBoundaryCellFilter(const XMLObject& xml)
// {
//   checkTag(xml, "BoundaryCells");
//   const string& name = xml.getRequired("name");
  
//   TEST_FOR_EXCEPTION(xml.numChildren() != 0, RuntimeError,
//                      "Interpreter::createBoundaryCellFilter() "
//                      "expected exactly zero "
//                      "children in " << xml.toString());

//   CellFilter filter = new BoundaryCellFilter();

//   cellFilter_.put(name, filter);
// }


// void Interpreter::createLabeledCellFilter(const XMLObject& xml)
// {
//   checkTag(xml, "CellFilter");
//   const string& name = xml.getRequired("name");
  
//   TEST_FOR_EXCEPTION(xml.numChildren() != 0, RuntimeError,
//                      "Interpreter::createCellFilter() "
//                      "expected exactly zero "
//                      "children in " << xml.toString());

//   const string& superName = xml.getRequired("superset");

//   TEST_FOR_EXCEPTION(!cellFilter_.containsKey(superName), RuntimeError,
//                      "cell filter [" << superName << "] not found in "
//                      "map " << cellFilter_.toString());

//   CellFilter super = cellFilter_.get(superName);
//   CellFilter filter;

//   if (xml.hasAttribute("label"))
//     {
//       int label = xml.getRequiredInt("label");
//       filter = super.labeledSubset(label);
//     }

//   cellFilter_.put(name, filter);
// }






