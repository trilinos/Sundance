/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#include "Sundance.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "SundanceAssembler.hpp"
#include "SundanceQuadratureIntegral.hpp"
#include "SundanceQuadratureEvalMediator.hpp"
#include "SundanceRefIntegral.hpp"
#include "SundanceEvaluator.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceCellJacobianBatch.hpp"
#include "SundanceUnaryFunctor.hpp"
#include "SundanceGrouperBase.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "SundanceDefaultPath.hpp"
#include "SundanceVersionString.hpp"


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

  Tabs tab;

  cerr << "Simulation built using Sundance version " << VersionString::number() 
       << " (" << VersionString::date() << ")" << endl;

  cerr << "Sundance is copyright (C) 2005 Sandia National Laboratories and is"
       << endl;
  cerr << "licensed under the GNU Lesser General Public License, version 2.1" << endl;
  cerr << tab << endl;
  

  /* read standard command line flags */
  string configFilename = "SundanceConfig.xml";

  bool defaultFpCheck = false;
  bool cmdFpCheck = defaultFpCheck;
  int defaultWorkSetSize = 100;
  int cmdWorkSetSize = defaultWorkSetSize;

  Assembler::workSetSize() = defaultWorkSetSize;

  clp().setOption("config", &configFilename, "Configuration file");
  clp().setOption("fpcheck", "nofpcheck", &cmdFpCheck, 
                  "Check results of math lib calls in expr evals");
  clp().setOption("workset", &cmdWorkSetSize, 
                  "Work set size");


  clp().throwExceptions(false);

  CommandLineProcessor::EParseCommandLineReturn rtn 
    = clp().parse(*argc, (char**) *argv);

  TEST_FOR_EXCEPTION(rtn != CommandLineProcessor::PARSE_SUCCESSFUL,
                     RuntimeError,
                     "Command-line parsing failed");

  /* process the settings file */
  setSettings(configFilename);


  if (cmdWorkSetSize != defaultWorkSetSize)
    {
      Assembler::workSetSize() = cmdWorkSetSize;
    }
} 

void Sundance::setOption(const string& optionName, 
                         int& value, 
                         const string& helpMsg)
{
  clp().setOption(optionName.c_str(), &value, helpMsg.c_str());
}

void Sundance::setOption(const string& optionName, 
                         double& value, 
                         const string& helpMsg)
{
  clp().setOption(optionName.c_str(), &value, helpMsg.c_str());
}

void Sundance::setOption(const string& optionName, 
                         string& value, 
                         const string& helpMsg)
{
  clp().setOption(optionName.c_str(), &value, helpMsg.c_str());
}

void Sundance::setOption(const string& optionTrueName, 
                         const string& optionFalseName,
                         bool& value, 
                         const string& helpMsg)
{
  clp().setOption(optionTrueName.c_str(), 
                  optionFalseName.c_str(),
                  &value, 
                  helpMsg.c_str());
}



void Sundance::handleException(std::exception& e)
{
  cerr << "Sundance detected exception: " << endl;
  cerr << e.what() << endl;
}


void Sundance::finalize()
{
  Tabs tab;
  cerr << tab << "eval vector flops: " << EvalVector::totalFlops() << endl;
  cerr << tab << "quadrature flops: " << QuadratureIntegral::totalFlops() << endl;
  cerr << tab << "ref integration flops: " 
       << RefIntegral::totalFlops() << endl;
  cerr << tab << "cell jacobian batch flops: " << CellJacobianBatch::totalFlops() << endl;
  cerr << tab << "quadrature eval mediator: " << QuadratureEvalMediator::totalFlops() << endl;
  TimeMonitor::summarize();
  MPISession::finalize();
}

string Sundance::searchForFile(const string& name)
{
  string pathSep = "/";
  Array<string> path = parsePathStr();

  for (int i=0; i<path.size(); i++)
    {
      ifstream fileToTry((path[i] + pathSep + name).c_str());
      if (!fileToTry) continue;
      return path[i] + pathSep + name;
    }

  TEST_FOR_EXCEPTION(true, RuntimeError, "could not find file "
                     << name << " in path " << path);
}

string Sundance::getPathStr() 
{
  char* pathEnvStr = getenv("SUNDANCE_PATH");

  if (pathEnvStr == NULL) 
    {
      return defaultSundancePath();
    }
  else
    {
      return pathEnvStr;
    }
}

Array<string> Sundance::parsePathStr() 
{
  string pathStr = getPathStr();
  
  Array<string> rtn;

  int begin;
  int end;
  
  begin = pathStr.find_first_not_of(":");
  
  while (begin < pathStr.length())
    {
      end = pathStr.find_first_of(":", begin);

      rtn.append(pathStr.substr(begin, end-begin));
      begin = pathStr.find_first_not_of(":", end);
    }

  return rtn;
}


void Sundance::setSettings(const XMLObject& xml)
{
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
            }
          else if (name=="Check for Floating Point Errors")
            {
              UnaryFunctor::checkResults() = child.getRequiredBool("value");
            }
          else if (name=="Shadow Calculations with String Values")
            {
              EvalVector::shadowOps() = child.getRequiredBool("value");
            }
        }
      else if (child.getTag()=="Verbosity")
        {
          const string& context = child.getRequired("context");
          const string& value = child.getRequired("value");
          if (context=="Evaluation")
            {
              SundanceCore::Internal::Evaluator::classVerbosity() = verbosity(value);
            }
          else if (context=="Assembly")
            {
              Internal::Assembler::classVerbosity() = verbosity(value);
            }
          else if (context=="Linear Problem")
            {
              LinearProblem::classVerbosity() = verbosity(value);
            }
          else if (context=="Integration Management")
            {
              Internal::IntegralGroup::classVerbosity() = verbosity(value);
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
              SundanceCore::Internal::SparsitySuperset::classVerbosity() = verbosity(value);
              SundanceCore::Internal::EvaluatableExpr::classVerbosity() = verbosity(value);
            }
          else if (context=="Integral Grouping")
            {
              Internal::GrouperBase::classVerbosity() = verbosity(value);
            }
          else if (context=="Mesh Creation")
            {
              SundanceStdMesh::Internal::MeshSourceBase::classVerbosity() 
                = verbosity(value);
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
    }
}

void Sundance::setSettings(const string& settingsFile)
{
  string fqFile = searchForFile(settingsFile);
  
  FileInputSource fis(fqFile);

  XMLObject xml = fis.getObject();

  cerr << "settings are: " << xml.toString() << endl;

  setSettings(xml);

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
