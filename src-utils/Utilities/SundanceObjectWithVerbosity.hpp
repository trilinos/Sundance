/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_OBJECTWITHVERBOSITY_H
#define SUNDANCE_OBJECTWITHVERBOSITY_H

#include "SundanceDefs.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceParamUtils.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_FancyOStream.hpp"


namespace SundanceUtils
{
using Teuchos::RefCountPtr;
using Teuchos::rcp;
using Teuchos::FancyOStream;
using Teuchos::ParameterList;
using Teuchos::ParameterEntry;


/** 
 * Defines traits for getting default verbosity parameters from a class
 */
template <class X>
class VerbosityTraits
{
public:
  static RefCountPtr<ParameterList> defaultVerbParams()
    {return X::defaultVerbParams();}
};

/**
 * ObjectWithVerbosity and the related verbosity() method
 * provide an interface for getting/setting
 * verbosity flags for classes or instances. 
 *
 * All objects start out with a verbosity setting of 0.
 * 
 * You can set verbosity for a single instance of a class, or for
 * the whole class. To set for an instance, use the verbosity()
 * member function, for example,
 * \code
 * Mesh mesh1 = reader1.getMesh();
 * Mesh mesh2 = reader2.getMesh();
 * Mesh mesh3 = reader3.getMesh();
 * mesh1.verbosity() = 3;
 * \endcode
 * which sets the verbosity of <tt>mesh1</tt> to 3 and leaves
 * those of <tt>mesh2</tt> and <tt>mesh3</tt> unchanged.
 *
 * Alternatively, you can set a default verbosity for an entire
 * class, for example,
 * \code
 * Mesh mesh1 = reader1.getMesh();
 * Mesh mesh2 = reader2.getMesh();
 * Mesh mesh3 = reader3.getMesh();
 * mesh1.verbosity() = 3;
 * verbosity<Mesh>() = 2;
 * \endcode
 * which sets the default verbosity to 2. Since <tt>mesh1</tt>
 * has its own verbosity setting of 3, 
 * it will use it rather than the
 * default, but <tt>mesh2</tt> and <tt>mesh3</tt> will use 2.
 * 
 */
template <class X>
class ObjectWithVerbosity
{
public:
  /** \deprecated Construct, starting silent */
  ObjectWithVerbosity(int verb=classVerbosity()) : verb_() {;}

  /** */
  int verb() const {return verb_;}

  /** */
  int& verb() {return verb_;}

  /** \deprecated Writeable access to the default verbosity for the class */
  static int& classVerbosity() 
    {
      static int rtn = 0;
      return rtn;
    }

  /** */
  static FancyOStream& os()
    {
      static RefCountPtr<std::ostream> os = rcp(&std::cout, false);
      static RefCountPtr<FancyOStream> rtn = fancyOStream(os);
      static bool first = true;
      if (first)
      {
        rtn->setShowProcRank(true);
        first = false;
      }
      return *rtn;
    }
  

private:
  /** */
  int verb_;
};

/** 
 * \relates ObjectWithVerbosity
 * Global method for setting verbosity of a class
 */
template <class X> int& verbosity() 
{
  return X::classVerbosity();
}

template <class X> 
class ParameterControlledObjectWithVerbosity 
  : public ObjectWithVerbosity<X> 
{
public:
  /** \deprecated Construct, starting silent */
  ParameterControlledObjectWithVerbosity() : ObjectWithVerbosity<X>() {;}

  /** Construct with a parameter list controlling the verbosity settings */
  ParameterControlledObjectWithVerbosity(const std::string& objName, const ParameterList& p)
    : ObjectWithVerbosity<X>(),
    verbControlParams_() 
    {
      RefCountPtr<ParameterList> defaults = VerbosityTraits<X>::defaultVerbParams();
      TEST_FOR_EXCEPTION(defaults->name() != objName, std::runtime_error,
        "mismatched ParameterList names for verbosity control: expected "
        << defaults->name() << ", got " << objName);
      TEST_FOR_EXCEPTION(defaults->name() != p.name(), std::runtime_error,
        "mismatched ParameterList names for verbosity control: expected "
        << defaults->name() << ", got " << p.name());
      verbControlParams_ = rcp(new ParameterList(mergeParams(*defaults, p)));
    }

  /** */
  int verbLevel(const std::string& context) const 
    {
      const ParameterEntry* pe = verbControlParams_->getEntryPtr(context);
      TEST_FOR_EXCEPTION(pe==0, std::runtime_error,
        "parameter with name \"" << context << "\" not found in verbosity "
        "control parameter list " << verbControlParams_);
      TEST_FOR_EXCEPTION(pe->getAny().type() != typeid(int),
        std::runtime_error,
        "context parameter name \"" 
        << context << "\" does not have an integer value in verbosity "
        "control parameter list " << verbControlParams_);

      return Teuchos::any_cast<int>(pe->getAny());
    }

  /** */
  const ParameterList& verbSublist(const std::string& name) const 
    {
      TEST_FOR_EXCEPTION(!verbControlParams_->isSublist(name),
        std::runtime_error,
        "context parameter name \"" 
        << name << "\" is not a sublist in verbosity "
        "control parameter list " << *verbControlParams_);

      return verbControlParams_->sublist(name);
    }

  /** */
  ParameterList mergeParams(const ParameterList& pDef, const ParameterList& pIn) const
    {
      return mergeParamLists(pDef, pIn);
    }
       
  /** */
  const ParameterList params() const {return *verbControlParams_;}

  /** */
  RefCountPtr<ParameterList> modifiableParams() const {return verbControlParams_;}
private:
  RefCountPtr<ParameterList> verbControlParams_;
};



}





#endif
