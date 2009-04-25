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

#ifndef SUNDANCE_BASISFAMILYBASE_H
#define SUNDANCE_BASISFAMILYBASE_H

#include "SundanceDefs.hpp"
#include "SundanceBasisDOFTopologyBase.hpp"
#include "SundanceBasisReferenceEvaluationBase.hpp"
#include "TSFHandleable.hpp"
#include "TSFPrintable.hpp"
#include "TSFObjectWithVerbosity.hpp"
#include "SundanceTypeUtils.hpp"
#include "Teuchos_XMLObject.hpp"

namespace SundanceStdFwk {


using Teuchos::Array;
using SundanceUtils::Point;
using SundanceStdMesh::CellType;
using SundanceCore::Internal::MultiIndex;
using Teuchos::RefCountPtr;
using Teuchos::rcp;


/** 
 *
 */
class BasisFamilyBase
  : public TSFExtended::Handleable<BasisFamilyBase>,
    public TSFExtended::Printable,
    public TSFExtended::ObjectWithVerbosity<BasisFamilyBase>,
    public BasisDOFTopologyBase,
    public BasisReferenceEvaluationBase
{
public:

  /** 
   * \brief Return the polynomial order of the basis functions, for use in
   * determining the quadrature rule required to integrate a product of
   * basis functions exactly. The polynomial order will be the smallest
   * integer for which all mixed partial derivatives vanish exactly.
   *  
   * Note: in H(div) and H(curl) spaces the order of accuracy is not
   * always an integer, and the relationship between the order of accuracy
   * and the return value of the order() method is not necessarily simple
   * (for instance, it can depend on things such as the convexity of the
   * boundary). Thus it is best to think of this method as specifying
   * the required order of quadrature, and neither the order of accuracy of
   * approximation nor the order to which the space is complete.
   */
  virtual int order() const = 0 ;

  /** 
   * \brief Return the dimension of the members of 
   * a vector-valued basis. Return 1 if the basis
   * is scalar-valued. 
   */
  virtual int dim() const = 0 ;

  /** \brief Inform caller as to whether I am a scalar basis. Default
   * implementation returns false. Overridden by ScalarBasis. */
  virtual bool isScalarBasis() const {return false;}

  /** \brief Inform caller as to whether I am a covariant basis. Default
   * implementation returns false. Overridden by CovariantBasis. */
  virtual bool isCovariantBasis() const {return false;}

  /** \brief Inform caller as to whether I am a contravariant basis. Default
   * implementation returns false. Overridden by CovariantBasis.  */
  virtual bool isContravariantBasis() const {return false;}

  /** */
  virtual bool lessThan(const BasisDOFTopologyBase* other) const ;
};





/* ----- Subtypes of BasisFamilyBase that specify transformation type --- */

/**     
 * Base class for scalar-valued basis families. Bases for scalar
 * fields living in, e.g., H1 or L2, should derive from this class. 
 */
class ScalarBasis 
  : public BasisFamilyBase
{
public:
  /** Inform caller that I am a scalar basis */
  bool isScalarBasis() const {return true;}
  
   /** 
    * \brief Return the dimension of the members of a scalar-valued basis  
    */
  int dim() const {return 1;}
};

/** 
 * Base class for bases for covariant vector fields (fields in H(div) 
 * spaces).
 */
class CovariantVectorBasis
  : public BasisFamilyBase
{
public:
  /** Inform caller that I am a covariant basis */
  bool isCovariantBasis() const {return true;}

};


/** 
 * Base class for bases for contravariant vector fields (fields in H(curl) 
 * spaces).
 */
class ContravariantVectorBasis
  : public BasisFamilyBase
{
public:
  /** Inform caller that I am a contravariant basis */
  bool isContravariantBasis() const {return true;}
};



} // namespace SundanceStdFwk


#endif
