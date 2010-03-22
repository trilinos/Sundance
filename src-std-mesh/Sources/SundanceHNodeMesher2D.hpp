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

#ifndef SUNDANCE_HNODEMESHER2D_H_
#define SUNDANCE_HNODEMESHER2D_H_

#include "SundanceDefs.hpp"
#include "SundanceMeshSourceBase.hpp"
#include "SundanceHNodeMesh2D.hpp"


namespace Sundance
{
class HNodeMesher2D : public MeshSourceBase
  {
  public:
    /**     */
	HNodeMesher2D(double position_x, double position_y,
                        double offset_x , double offset_y,
                        double resolution_x , double resolution_y ,
                        const MeshType& meshType,
                        const MPIComm& comm = MPIComm::world())
      :
      MeshSourceBase(meshType, comm),
      _position_x(position_x), _position_y(position_y),
      _offset_x(offset_x), _offset_y(offset_y),
      _resolution_x(resolution_x) , _resolution_y(resolution_y) {;}


    /** Create a rectangle mesher from a ParameterList */
	HNodeMesher2D(const ParameterList& params);

    /** */
    virtual ~HNodeMesher2D() {;}

    /** Print a short descriptive string */
    virtual string description() const
    {return "HNodeMesher2D[pos x =" + Teuchos::toString(_position_x)
       + ", pos y=" + Teuchos::toString(_position_y)
       + ", offset x=" + Teuchos::toString(_offset_x) +
       + ", offset y=" + Teuchos::toString(_offset_y)
       + ", resolution_x=" + Teuchos::toString(_resolution_x)
       + ", resolution_y=" + Teuchos::toString(_resolution_y)+"]";}


    /** Return a ref count pointer to self */
    virtual RCP<MeshSourceBase> getRcp() {return rcp(this);}


  protected:

    /** The method which all Mesher should have */
    virtual Mesh fillMesh() const ;

  private:

    /** The X coordinate of the origin point (lower left)*/
    double _position_x;
    /** The Y coordinate of the origin point (lower left)*/
    double _position_y;
    /** The offset (length) of the grid in the X direction*/
    double _offset_x;
    /** The offset (length) of the grid in the Y direction*/
    double _offset_y;
    /** */
    double _resolution_x;
    /** */
    double _resolution_y;

  };
}
#endif /* SUNDANCE_HNODEMESHER2D_H_ */
