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


#include "SundanceHNodeMesher2D.hpp"

using namespace Sundance;
using namespace Teuchos;


HNodeMesher2D::HNodeMesher2D(const ParameterList& params)
: MeshSourceBase(params),
  _position_x(params.get<double>("position_x")),
  _position_y(params.get<double>("position_y")),
  _offset_x(params.get<int>("offset_x")),
  _offset_y(params.get<double>("offset_y")),
  _resolution_x(params.get<double>("resolution_x")),
  _resolution_y(params.get<double>("resolution_y"))
{
	// nothing to do here
}

Mesh HNodeMesher2D::fillMesh() const
{
	// here we create the mesh and return to the Sundance
	HNodeMesh2D *hnodegrid;
	Mesh mesh = createMesh(2);
	// get the pointer
	hnodegrid = (HNodeMesh2D*) mesh.ptr().get();
	hnodegrid->createMesh( _position_x , _position_y , _offset_x , _offset_y , _resolution_x , _resolution_y );
	return (mesh);
}

