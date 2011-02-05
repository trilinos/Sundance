/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_GHOST_VIEW_HPP
#define PLAYA_GHOST_VIEW_HPP

#include "PlayaVectorDecl.hpp"

namespace Playa
{
  using namespace Teuchos;

  /**
   * GhostView is an interface for read-only views
   * of vector elements including selected
   * off-processor elements. GhostView has no standard constructor; subclasses
   * should be constructed using the importView() method of GhostImporter.
   */
  template <class Scalar>
  class GhostView : public Printable
  {
  public:
    /** Virtual dtor */
    virtual ~GhostView(){;}
    
    /** Indicate whether the value at the given global index is accessible
     * in this view. */
    virtual bool isAccessible(int globalIndex) const = 0 ;

    
    /** */
    virtual const double& getElement(int globalIndex) const = 0 ;
    
    /** */
    virtual void getElements(const int* globalIndices, int numElems,
      Teuchos::Array<double>& elems) const = 0 ;
    
    /**  */
    virtual void print(std::ostream& os) const = 0 ;

  private:
  };

}

#endif
