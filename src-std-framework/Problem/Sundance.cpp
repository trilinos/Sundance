/* @HEADER@ */
/* @HEADER@ */

#include "Sundance.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "SundanceAssembler.hpp"
#include "SundanceQuadratureIntegral.hpp"
#include "SundanceRefIntegral.hpp"
#include "SundanceEvaluator.hpp"
#include "SundanceUnaryFunctor.hpp"
#include "SundanceBruteForceEvaluator.hpp"
#include "SundanceInstructionCachingEvaluator.hpp"
#include "SundanceGrouperBase.hpp"
#include "SundanceSparsityPattern.hpp"


static Time& totalTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("total Sundance time"); 
  return *rtn;
}

void Sundance::init(int* argc, void*** argv)
{
  /* start up MPI. In a serial run, this will be a no-op */
  MPISession::init(argc, argv);

  /* read standard command line flags */
  CommandLineProcessor clp;
  string configFilename = "SundanceConfig.xml";
  bool fpCheck = false;
  int workSetSize = 100;
  clp.setOption("config", &configFilename, "Configuration file");
  clp.setOption("fpcheck", "nofpcheck", &fpCheck, 
                "Check results of math lib calls in expr evals");
  clp.setOption("workset", &workSetSize, 
                "Work set size");


  clp.throwExceptions(false);

  CommandLineProcessor::EParseCommandLineReturn rtn = clp.parse(*argc, (char**) *argv);

  TEST_FOR_EXCEPTION(rtn != CommandLineProcessor::PARSE_SUCCESSFUL,
                     RuntimeError,
                     "Command-line parsing failed");

  /* process the settings file */
  setSettings(configFilename);


  UnaryFunctor::checkResults() = fpCheck;
  Assembler::workSetSize() = workSetSize;
  FunctionalEvaluator::workSetSize() = workSetSize;
  
} 

void Sundance::handleException(std::exception& e)
{
  cerr << "Sundance detected exception: " << endl;
  cerr << e.what() << endl;
}


void Sundance::finalize()
{
  TimeMonitor::summarize();
  MPISession::finalize();
}

string Sundance::searchForFile(const string& name)
{
  string pathSep = "/";
  Array<string> path;
  path.append(".");
  path.append("../../etc");
  path.append("../../../etc");

  for (int i=0; i<path.size(); i++)
    {
      ifstream fileToTry((path[i] + pathSep + name).c_str());
      if (!fileToTry) continue;
      return path[i] + pathSep + name;
    }
}

void Sundance::setSettings(const string& settingsFile)
{
  string fqFile = searchForFile(settingsFile);
  
  FileInputSource fis(fqFile);

  XMLObject xml = fis.getObject();

  for (int i=0; i<xml.numChildren(); i++)
    {
      const XMLObject& child = xml.getChild(i);
      if (child.getTag()=="Parameter")
        {
          const string& name = child.getRequired("name");
          if (name=="Work Set Size")
            {
              int workSetSize = child.getRequiredInt("value");
              Assembler::workSetSize() = workSetSize;
              FunctionalEvaluator::workSetSize() = workSetSize;
            }
          else if (name=="Check for Floating Point Errors")
            {
              UnaryFunctor::checkResults() = child.getRequiredBool("value");
            }
        }
      else if (child.getTag()=="Verbosity")
        {
          const string& context = child.getRequired("context");
          const string& value = child.getRequired("value");
          cerr << context << " verbosity=" << value << endl;
          if (context=="Evaluation")
            {
              SundanceCore::Internal::Evaluator::classVerbosity() = verbosity(value);
            }
          else if (context=="Assembly")
            {
              Internal::Assembler::classVerbosity() = verbosity(value);
            }
          else if (context=="Quadrature")
            {
              Internal::QuadratureIntegral::classVerbosity() = verbosity(value);
            }
          else if (context=="Reference Integration")
            {
              Internal::RefIntegral::classVerbosity() = verbosity(value);
            }
          else if (context=="DOF Mapping")
            {
              Internal::DOFMapBase::classVerbosity() = verbosity(value);
            }
          else if (context=="Symbolic Sparsity Determination")
            {
              SundanceCore::Internal::SparsityPattern::classVerbosity() = verbosity(value);
            }
          else if (context=="Integral Grouping")
            {
              Internal::GrouperBase::classVerbosity() = verbosity(value);
            }
        }
      else if (child.getTag()=="DefaultMesh")
        {
          string type = child.getRequired("type");
          if (type=="BasicSimplicial")
            {
              MeshSource::defaultMeshType() 
                = new BasicSimplicialMeshType();
            }
        }
      else if (child.getTag()=="DefaultQuadrature")
        {
          string type = child.getRequired("type");
          int order = child.getRequiredInt("order");
          if (type=="Gaussian")
            {
              QuadratureFamilyStub::defaultQuadrature() 
                = rcp(new GaussianQuadrature(order));
            }
        }
      else if (child.getTag()=="Evaluation")
        {
          string type = child.getRequired("type");
          if (type=="Instruction Caching")
            {
              EvaluatorFactory::defaultEvaluator() 
                = rcp(new InstructionCachingEvaluatorFactory());
              cerr << "Evaluation method: instruction caching"<< endl;
            }
          else
            {
              EvaluatorFactory::defaultEvaluator() 
                = rcp(new BruteForceEvaluatorFactory());
              cerr << "Evaluation method: brute force"<< endl;
            }
        }
      
      
    }
}

VerbositySetting Sundance::verbosity(const string& str)
{
  if (str=="Low")
    {
      return VerbLow;
    }
  else if (str=="Medium")
    {
      return VerbMedium;
    }
  else if (str=="High")
    {
      return VerbHigh;
    }
  else if (str=="Extreme")
    {
      return VerbExtreme;
    }
  return VerbSilent;
}


bool Sundance::checkTest(double error, double tol)
{
  return error < tol;
}

void Sundance:: passFailTest(double error, double tol)
{
  cerr << "error norm = " << error << endl;
  cerr << "tolerance = " << tol << endl;
  if (checkTest(error, tol))
    {
      cerr << "test PASSED" << endl;
    }
  else
    {
      cerr << "test FAILED" << endl;
    }
}
