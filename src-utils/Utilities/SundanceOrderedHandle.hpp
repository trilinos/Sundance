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

#ifndef SUNDANCEORDEREDHANDLE_HPP
#define SUNDANCEORDEREDHANDLE_HPP

#include "TSFConfigDefs.hpp"
#include "TSFHandle.hpp"
#include <typeinfo>


#define ORDERED_HANDLE_CTORS(handle, contents) \
/** Empty ctor */ \
handle() : OrderedHandle<contents >() {;} \
/** Construct a #handle with a raw pointer to a #contents */ \
handle(TSFExtended::Handleable<contents >* rawPtr) : OrderedHandle<contents >(rawPtr) {;} \
/** Construct a #handle with a smart pointer to a #contents */ \
handle(const RefCountPtr<contents >& smartPtr) : OrderedHandle<contents >(smartPtr){;}




namespace SundanceUtils
{
  using namespace Teuchos;

  /**
   * Class OrderedHandle is an extension to TSFExtended::Handle that 
   * includes a comparison operator ("<" operator) so that
   * the handle can be used in ordered containers such as STL maps and sets.
   */
  template <class PointerType>
  class OrderedHandle : public TSFExtended::Handle<PointerType>
  {
  public:
    /** empty ctor */
    OrderedHandle() : TSFExtended::Handle<PointerType>() {;}

    /** Construct from a raw ptr */
    OrderedHandle(TSFExtended::Handleable<PointerType>* rawPtr) : TSFExtended::Handle<PointerType>(rawPtr) {;}

    /** Construct from a smart ptr*/
    OrderedHandle(const RefCountPtr<PointerType>& smartPtr) 
      : TSFExtended::Handle<PointerType>(smartPtr) {;}

    /** comparison operator */
    bool operator<(const OrderedHandle<PointerType>& other) const 
    {
      /* first compare types */
      const PointerType* me = this->ptr().get();
      const PointerType* you = other.ptr().get();
      if (me==0 && you==0) 
        {
          return false;
        }
      if (me==0) 
        {
          return true;
        }
      if (you==0) 
        {
          return false;
        }

      if (typeid(*me).before(typeid(*you))) 
        {
          return true;
        }

      if (typeid(*you).before(typeid(*me)))
        {
          return false;
        }
      
      /* if the types are equal, compare values of the contents using
       * the lessThan() method. */
      bool rtn = this->ptr()->lessThan(other.ptr().get());
      return rtn;
    }

    /** */
    bool operator!=(const OrderedHandle<PointerType>& other) const
      {
        return (*this < other) || (other < *this);
      }


    /** */
    bool operator==(const OrderedHandle<PointerType>& other) const
      {
        return !(*this != other);
      }

    
  };
}

namespace std
{
  /** \relates OrderedHandle */
  template <class P> inline 
  std::ostream& operator<<(std::ostream& os, 
                           const SundanceUtils::OrderedHandle<P>& h)
  {
    h.print(os);
    return os;
  }
}




#endif

