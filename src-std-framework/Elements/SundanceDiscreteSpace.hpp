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

#ifndef SUNDANCE_DISCRETESPACE_H
#define SUNDANCE_DISCRETESPACE_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceDOFMapBase.hpp"
#include "TSFVectorType.hpp"
#include "TSFVectorSpace.hpp"
#include "TSFObjectWithVerbosity.hpp"

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace TSFExtended;
  using namespace SundanceStdMesh;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;
  using namespace Internal;

  /** 
   * DiscreteSpace represents a discrete finite-element space (i.e., 
   * a mesh and a basis).
   */
  class DiscreteSpace : public ObjectWithVerbosity<DiscreteSpace>
  {
  public:
    /** */
    DiscreteSpace(){;}
    /** */
    DiscreteSpace(const Mesh& mesh, const BasisFamily& basis,
                  const VectorType<double>& vecType);
    /** */
    DiscreteSpace(const Mesh& mesh, const Array<BasisFamily>& basis,
                  const VectorType<double>& vecType);

    /** */
    DiscreteSpace(const Mesh& mesh, const Array<BasisFamily>& basis,
                  const RefCountPtr<DOFMapBase>& map,
                  const VectorType<double>& vecType);

    /** */
    const RefCountPtr<DOFMapBase>& map() const {return map_;}

    /** return the number of functions */
    int nFunc() const {return basis_.size();}

    /** */
    const Array<BasisFamily>& basis() const {return basis_;}

    /** */
    Vector<double> createVector() const {return vecSpace_.createMember();}

    /** */
    VectorSpace<double> vecSpace() const {return vecSpace_;}

    /** */
    const Mesh& mesh() const {return mesh_;}

    /** */
    const VectorType<double>& vecType() const {return vecType_;}

    /** */
    void importGhosts(const Vector<double>& x,
                      RefCountPtr<GhostView<double> >& ghostView) const ;
  private:
    /** */
    RefCountPtr<DOFMapBase> map_;

    /** */
    Mesh mesh_;

    /** */
    Array<BasisFamily> basis_;

    /** */
    VectorSpace<double> vecSpace_;

    /** */
    VectorType<double> vecType_;

    /** */
    RefCountPtr<GhostImporter<double> > ghostImporter_;
  };

}



#endif
