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
    /* Declare variables for the options to be set on command line, and
     * initialize with default values.
     */
    int someInt = 137;
    double someDouble = 3.14159;
    string someString = "blue";
    bool someBool = false;

    Sundance::setOption("integer", someInt, "An integer");
    Sundance::setOption("alpha", someDouble, "A double");
    Sundance::setOption("color", someString, "What is your favorite color?");
    Sundance::setOption("lie", "truth", someBool, "I am lying.");

    /* Now that command-line parsing has been set up, call init */ 
    Sundance::init(&argc, &argv);

    /* Just for the heck of it, do something with the options */
    Out::root() << "User input:" << endl;
    Out::root() << "An integer: " << someInt << endl;
    Out::root() << "A double-precision number: " << someDouble << endl;
    Out::root() << "Favorite color: " << someString << endl;
    Out::root() << "I am lying: " << someBool << endl;

    Sundance::passFailTest(true);
  }
	catch(std::exception& e)
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); 
  return Sundance::testStatus();
}
