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

int main(int argc, char** argv)
{
  try
  {
    /* Read the XML filename as a command-line option */
    string xmlFilename = "paramExample.xml";
    Sundance::setOption("xml-file", xmlFilename, "XML filename");
      
    /* Initialize */
    Sundance::init(&argc, &argv);

    /* Read a parameter list from the XML file */
    ParameterXMLFileReader reader(xmlFilename);
    ParameterList params = reader.getParameters();

    /* Get the parameters for the "Widget" sublist */
    const ParameterList& widget = params.sublist("Widget");
    Out::root() << "widget region label: " << widget.get<int>("Region") << endl;
    Out::root() << "widget material: " << widget.get<string>("Material") << endl;
    Out::root() << "widget density: " << widget.get<double>("Density") << endl;

    /* Get the parameters for the "Gizmo" sublist */
    const ParameterList& gizmo = params.sublist("Gizmo");
    Out::root() << "gizmo region label: " << gizmo.get<int>("Region") << endl;
    Out::root() << "gizmo material: " << gizmo.get<string>("Material") << endl;
    Out::root() << "gizmo density: " << gizmo.get<double>("Density") << endl;

    Sundance::passFailTest(true);
  }
	catch(std::exception& e)
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); 
  return Sundance::testStatus();
}


