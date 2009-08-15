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

#ifndef SUNDANCE_IDENTITYREORDERER_H
#define SUNDANCE_IDENTITYREORDERER_H



#include "SundanceDefs.hpp"
#include "SundanceCellReordererImplemBase.hpp"
#include "SundanceCellReordererBase.hpp"
#include "SundanceHandleable.hpp"

namespace SundanceStdMesh
{
  
namespace Internal
{
using namespace Teuchos;
using namespace SundanceUtils;
using Teuchos::RCP;

/**
 * The identity reorderer walks through cells in whatever
 * order they are numbered the the mesh. 
 */
class IdentityReordererImplem : public CellReordererImplemBase
{
public:
  /** */
  IdentityReordererImplem(const MeshBase* mesh); 
      
  /** */
  virtual ~IdentityReordererImplem(){;}
    
  /** */
  virtual int advance(int currentLID) const {return currentLID+1;}
};
}

/** */
class IdentityReorderer 
  : public Internal::GenericCellReordererFactory<Internal::IdentityReordererImplem>
{
public:
  IdentityReorderer(){;}

  virtual ~IdentityReorderer(){;}

  GET_RCP(Internal::CellReordererFactoryBase);
};
}


#endif
