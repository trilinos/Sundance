/* @HEADER@ */
/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// **********************************************************************/
/* @HEADER@ */

#ifndef SundanceHANDLE_HPP
#define SundanceHANDLE_HPP

#include "SundanceDefs.hpp"
#include "SundancePrintable.hpp"
#include "Teuchos_Describable.hpp"
#include "SundanceHandleable.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_MPIComm.hpp"


/** \def This helper macro defines boilerplate constructors for classes deriving
 * from Handle.
 *
 * If class MyHandle is a handle to a type MyType, simply 
 * put
 * \code
 * HANDLE_CTORS(MyHandle, MyType);
 * \endcode
 * in the class declaration of MyHandle and the macro will create 
 * an empty ctor, a ctor from a smart ptr, and a ctor from a raw pointer. 
 * The macro will also create appropriate doxygen for the handle ctors */
#define HANDLE_CTORS(handle, contents) \
/** Empty ctor */ \
handle() : SundanceUtils::Handle<contents >() {;} \
/** Construct a #handle with a raw pointer to a #contents */ \
handle(SundanceUtils::Handleable<contents >* rawPtr) : SundanceUtils::Handle<contents >(rawPtr) {;} \
/** Construct a #handle with a smart pointer to a #contents */ \
handle(const Teuchos::RefCountPtr<contents >& smartPtr) : SundanceUtils::Handle<contents >(smartPtr){;}



#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceUtils
{
  using namespace Teuchos;

  /**
   * Class SundanceUtils::Handle provides a general implementation
   * of the common features of reference-counted handles.
   */
  template <class PointerType>
  class Handle : public Describable
  {
  public:
    /** Empty ctor  */
    Handle() : ptr_() {;}

    /** Construct from a smart pointer */
    Handle(const RefCountPtr<PointerType>& _ptr) : ptr_(_ptr) {;}

    /** Construct from a raw pointer to a Handleable.  */
    Handle(Handleable<PointerType>* rawPtr) : ptr_(rawPtr->getRcp()) {;}

    /** Read-only access to the underlying smart pointer. */
    const RefCountPtr<PointerType>& ptr() const {return ptr_;}

    /** Read-write access to the underlying smart pointer. */
    RefCountPtr<PointerType>& ptr() {return ptr_;}

    /** 
     * Print to a stream using the Printable interface. 
     * If the contents of the handle cannot be 
     * downcasted or crosscasted to a Printable*, an exception
     * will be thrown 
     */
    void print(std::ostream& os) const ;

    /** 
     * Return a short descriptive string using the Describable interface.
     * If the contents of the handle cannot be 
     * downcasted or crosscasted to a Describable*, an exception
     * will be thrown. 
     */
    std::string description() const ;


    /** */
    void describe(
      Teuchos::FancyOStream                &out_arg,
      const Teuchos::EVerbosityLevel      verbLevel
      ) const 
      {
        const Describable* p = dynamic_cast<const Describable*>(ptr().get());
        if (p!=0) p->describe(out_arg, verbLevel);
        else out_arg << description();
      }

    /** 
     * Return the verbosity setting using the ObjectWithVerbosity
     * interface. If the contents of the handle cannot be downcasted
     * or crosscasted into an ObjectWithVerbosity, an exception will
     * be thrown. 
     */
    VerbositySetting verbosity() const 
    {
      const ObjectWithVerbosity<PointerType>* v 
        = dynamic_cast<const ObjectWithVerbosity<PointerType>*>(ptr_.get());
      
      TEST_FOR_EXCEPTION(v==0, std::runtime_error,
                         "Attempted to cast non-verbose "
                         "pointer to an ObjectWithVerbosity");
      return v->verbosity();
    }

    /** 
     * Return a writeable reference to 
     * the verbosity setting using the ObjectWithVerbosity
     * interface. If the contents of the handle cannot be downcasted
     * or crosscasted into an ObjectWithVerbosity, an exception will
     * be thrown. 
     */
    VerbositySetting& verbosity() 
    {
      ObjectWithVerbosity<PointerType>* v 
        = dynamic_cast<ObjectWithVerbosity<PointerType>*>(ptr_.get());
      
      TEST_FOR_EXCEPTION(v==0, std::runtime_error,
                         "Attempted to cast non-verbose "
                         "pointer to an ObjectWithVerbosity");
      return v->verbosity();
    }

    /** Writeable access to the class-wide verbosity setting for the handled type */
    static VerbositySetting& classVerbosity() 
    {
      return PointerType::classVerbosity();
    }

  private:
    RefCountPtr<PointerType> ptr_;
  };

  /* implementation of print() */
  template <class PointerType> inline 
  void Handle<PointerType>::print(std::ostream& os) const 
  {
    const Printable* p = dynamic_cast<const Printable*>(ptr_.get());
    const Describable* d = dynamic_cast<const Describable*>(ptr_.get());
      
    TEST_FOR_EXCEPTION(p==0 && d==0, std::runtime_error,
      "Attempt to print a Handle to an object that is neither "
      "Printable nor Describable");

    if (p!=0) p->print(os);
    else if (d!=0) os << description();
  }

  /* implementation of description() */
  template <class PointerType> inline
  std::string Handle<PointerType>::description() const 
  {
    const Describable* p = dynamic_cast<const Describable*>(ptr_.get());
    
    TEST_FOR_EXCEPTION(p==0, std::runtime_error,
                       "Attempted to cast non-describable "
                       "pointer to a Describable");
    return p->description();
  }
}


template <class PointerType> inline
std::ostream& operator<<(std::ostream& os, const SundanceUtils::Handle<PointerType>& h)
{
  h.print(os);
  return os;
}

#define STREAM_OUT(handleType) \
inline ostream& operator<<(ostream& os, const handleType& h) \
{h.print(os); return os;}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif

