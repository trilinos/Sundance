#include "Sundance.hpp"
#include "SundanceExprParser.hpp"

using namespace SundanceXML;
using SundanceCore::List;

int main(int argc, void** argv)
{
  
  try
		{
      Sundance::init(&argc, &argv);

      Interpreter interpreter;

      string fqFile = Sundance::searchForFile("test.xml");
  
      FileInputSource fis(fqFile);
      
      XMLObject xml = fis.getObject();

      interpreter.interpret(xml);
      
    }
	catch(std::exception& e)
		{
      cerr << e.what() << endl;
		}
  Sundance::finalize();
}


void Interpreter::interpret(const XMLObject& xml)
{
  const string& tag = xml.getTag();

  if (tag=="Include")
    {
      const string& filename = xml.getRequired("filename");
      interpret(filename);
    }
  else if (tag=="Solver")
    {
      
    }
  else if (tag=="CellFilter")
    {
      
    }
  else if (tag=="Define")
    {
      
    }

  for (int i=0; i<xml.numChildren(); i++)
    {
      const XMLObject& child = xml.getChild(i);
      
    }
}
