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

#ifndef SUNDANCE_REGIONQUADCOMBO_H
#define SUNDANCE_REGIONQUADCOMBO_H


#include "SundanceDefs.hpp"
#include "SundanceMap.hpp"
#include "Teuchos_Utils.hpp"
#include "SundanceWatchFlag.hpp"
#include "SundanceCellFilterStub.hpp"
#include "SundanceQuadratureFamilyStub.hpp"
#include "SundanceOrderedTuple.hpp"
#include "SundanceOrderedHandle.hpp"


namespace SundanceCore
{
using namespace Teuchos;
using namespace SundanceUtils;
using std::string;
using SundanceUtils::Map;

/** */
typedef OrderedTriple<OrderedHandle<CellFilterStub>,
                      OrderedHandle<QuadratureFamilyStub>,
                      WatchFlag> RegTriple;
/** 
 * Expressions may appear in more than one subregions of a problem,
 * for instance in an internal domain and also on a boundary. On
 * those different subregions, a given expression might be subject
 * to different sets of functional derivatives; thus, different
 * evaluation regions might have different sparsity patterns. 
 * It is therefore necessary to build and store
 * sparsity information on a region-by-region basis. 
 *
 * Class RegionQuadCombo is used as an identifier for regions. The
 * only thing it needs to do is to be useable as a key in a STL map.
 */
class RegionQuadCombo 
{
public:
  /** */
  RegionQuadCombo();
  /** */
  RegionQuadCombo(const RefCountPtr<CellFilterStub>& domain,
    const RefCountPtr<QuadratureFamilyStub>& quad,
    const WatchFlag& watch = WatchFlag());

  /** */
  inline bool operator==(const RegionQuadCombo& other) const
    {return id_==other.id_;}

  /** */
  string toString() const ;

  /** */
  bool operator<(const RegionQuadCombo& other) const
    {return id_ < other.id_;}

  /** */
  const RefCountPtr<CellFilterStub>& domain() const {return domain_;}

  /** */
  const RefCountPtr<QuadratureFamilyStub>& quad() const 
    {return quad_;}

  /** */
  const WatchFlag& watch() const 
    {return watch_;}

private:

  /** */
  int id_;

  /** */
  RefCountPtr<CellFilterStub> domain_;

  /** */
  RefCountPtr<QuadratureFamilyStub> quad_;

  /** */
  WatchFlag watch_;
          
  /** */
  static int getID(const RefCountPtr<CellFilterStub>& domain,
    const RefCountPtr<QuadratureFamilyStub>& quad,
    const WatchFlag& watch);

  /** */
  static int topID() {static int rtn=0; return rtn++;}

  /** */
  static Map<RegTriple, int>& domainAndQuadToIDMap() ;
};

}

namespace std
{
/** \relates SundanceCore::RegionQuadCombo*/
inline ostream& operator<<(ostream& os, 
  const SundanceCore::RegionQuadCombo& c)
{
  os << c.toString();
  return os;
}
}

namespace Teuchos
{
using std::string;

/** \relates SundanceCore::RegionQuadCombo */
inline string toString(const SundanceCore::RegionQuadCombo& h)
{return h.toString();}

}

#endif
