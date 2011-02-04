/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_GHOSTIMPORTER_HPP
#define PLAYA_GHOSTIMPORTER_HPP

#include "PlayaVectorDecl.hpp"
#include "PlayaVectorSpaceDecl.hpp"
#include "PlayaGhostView.hpp"

namespace Playa
{
using namespace Teuchos;

/**
 * In many applications it is necessary to view some subset of 
 * off-processor, or "ghost", elements of a vector. In 
 * matrix-vector multiplications,
 * access to off-processor elements is assumed to be handled internally
 * by the apply() method of LinearOp subtypes, so the Playa Vector type
 * does not need explicit accessors for ghost elements. However, in 
 * application interfaces such as finite element assembly engines, 
 * read-only access to ghost elements is sometimes required. The abstract
 * classes GhostImporter and GhostView define flexible interfaces
 * through which a set of required ghosts can be defined, ghost values
 * can be imported, and element values can be accessed through
 * global indices.
 *
 * Class GhostImporter is used to specify the set of ghost elements
 * that must be imported to this processor, and then to carry out the import.
 * It will often be the case that we do many imports with the same
 * set of ghost indices; for example, in a nonlinear problem the
 * import of the same set of ghost indices
 * will be repeated at each function evaluation. Therefore, it makes sense
 * to do the definition of the ghost index set and the import as
 * distinct methods. The definition of the ghost index set should be done
 * in the constructors of GhostImporter subclasses.
 */
template <class Scalar>
class GhostImporter
{
public:
  /** virtual dtor */
  virtual ~GhostImporter(){;}

  /** 
   * Import the ghost elements of the given vector
   * as specified during construction of this object. 
   */
  virtual void importView(const Vector<Scalar>& x,
    RCP<GhostView<Scalar> >& ghostView) const = 0 ;

};

}

#endif
