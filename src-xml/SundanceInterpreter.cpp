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

      lines.append("bas = List(Lagrange(1), Lagrange(1));");

      lines.append("u = UnknownFunction(bas);");

      lines.append("x = CoordExpr(0);");

      lines.append("dx = Derivative(0);");

      lines.append("dy = Derivative(1);");

      lines.append("grad = List(dx, dy);");

      lines.append("grad = List(Derivative(0), Derivative(1));");

      lines.append("grad = List(Derivative(0), dx, List(dx, dy));");

      lines.append("y = sin(x);");

      lines.append("mesh = Exodus('fred');");

      lines.append("q4 = GaussianQuadrature(4);");

      lines.append("mesh = Line(0.0, 1.0, 10);");

      lines.append("bdry = BoundaryCellFilter();");

      lines.append("interior = MaximalCellFilter();");

      lines.append("left = LabeledCellFilter(bdry, 'left');");

      lines.append("eqn = Integral(interior, (grad*v)*(grad*u) + f*v) "
                   "+ Integral(left, v*sin(x), q4)");

      lines.append("bc = EssentialBC(top, v*(u-1.0))");

      lines.append("prob1 = LinearProblem(mesh, eqn1, bc, v, u)");

      lines.append("prob2 = LinearProblem(mesh, eqn2, bc2, List(vx,vy), List(ux,uy))");


                  


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








