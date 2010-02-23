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

#ifndef SUNDANCE_FEKETEQUADRATURE_H
#define SUNDANCE_FEKETEQUADRATURE_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceQuadratureFamilyBase.hpp"

namespace Sundance
{

using namespace Teuchos;

/**
 * Family of Fekete approach integration rules, suitable for integration
 * over arbitrary regions
 * (For tensor elements the Gauss-Lobatto points are Fekete points)
 */
class FeketeQuadrature: public QuadratureFamilyBase
{
public:
	/** */
	FeketeQuadrature(int order);

	/** */
	virtual ~FeketeQuadrature()
	{
		;
	}

	/** */
	virtual XMLObject toXML() const;

	/** Describable interface */
	virtual string description() const
	{
		return "FeketeQuadrature[order=" + Teuchos::toString(order())
				+ "]";
	}

	/* handleable boilerplate */
	GET_RCP(QuadratureFamilyStub)
	;

protected:
	/** compute a rule for the reference line cell */
	virtual void getLineRule(Array<Point>& quadPoints,
			Array<double>& quadWeights) const;

	/** compute a rule for the reference triangle cell */
	virtual void getTriangleRule(Array<Point>& quadPoints,
			Array<double>& quadWeights) const;

	/** compute a rule for the reference quad cell */
	virtual void getQuadRule(Array<Point>& quadPoints,
			Array<double>& quadWeights) const;

	/** compute a rule for the reference brick cell */
	virtual void getBrickRule(Array<Point>& quadPoints,
			Array<double>& quadWeights) const;
};
}

#endif
