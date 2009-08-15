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

#ifndef SUNDANCE_SPECTRALBASIS_H
#define SUNDANCE_SPECTRALBASIS_H

#include "SundanceDefs.hpp"
#include "SundanceHandle.hpp"
#include "SundanceHandleable.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceSpectralBasisBase.hpp"



using namespace std;
using namespace SundanceUtils;

namespace SundanceCore
{
  /** Doxygen doc for SpectralBasis */
  class SpectralBasis : public SundanceUtils::Handle<SpectralBasisBase>
  {
  public:
    /* boilerplate handle ctors */
    HANDLE_CTORS(SpectralBasis, SpectralBasisBase);

    /** Return the dim of the Spectral Basis */
    int getDim() const {return ptr()->getDim();}

    /** Return the order of the Spectral Basis */
    int getOrder() const {return ptr()->getOrder();}

    /** Return the maximum number of terms */
    int nterms() const {return ptr()->nterms();}
    
    /** Return the basis element stored in the basis array index */
    int getElement(int i) const {return ptr()->getElement(i);}
    
    /** expectation operator */
    double expectation(int i, int j, int k) const 
    {return ptr()->expectation(i,j,k);}

    /** Write to a string */
    string toString() const {return ptr()->toString();}
  };
}

namespace std
{
  /** \relates  SundanceCore::SpectralBasis */
  inline ostream& operator<<(ostream& os, const SundanceCore::SpectralBasis& s)
  {
    os << s.toString();
    return os;
  }
}

#endif
