/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_BASICVERTEXSET_H
#define SUNDANCE_BASICVERTEXSET_H


#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"

namespace Sundance
{
  namespace Internal
  {
    using namespace Teuchos;

    /**
     * VertexSet is a helper class whose use reduces the number
     * of temporary arrays created during the process of face 
     * identification. 
     */
    class VertexSet
    {
    public:
      /** empty ctor */
      VertexSet() : base_(0), offset_(0), length_(0), rotation_(0) {;}
      /** Construct with an offset into an array of vertex triplets */
      VertexSet(int** base, int offset, int length)
        : base_(base), offset_(offset), length_(length), rotation_(0) {;}

      /**
       * Return a hash code for the vertex set. Because the vertices 
       * might be presented in
       * permuted order, we need a hash function that is invariant 
       * under permutations.
       */
      int hashCode() const ;

      bool operator==(const VertexSet& other) const ;

      string toString() const ;

      /** 
       * When a face is already present, a new face with the same 
       * vertices will see the vertices in reversed order since both
       * cells will have ordered their face nodes such that they
       * run counterclockwise with the outward normals. The vertices
       * may also be "rotated", since the two cells might list
       * this face's vertices starting at a different point. For example,
       * if cell A has a face <tt>{1, 2, 3}</tt>, an adjoining cell
       * B sharing that face might list them in any one of the following
       * orders:
       *
       * \code
       * {1,3,2}  // rotation = 0
       * {2,1,3}  // rotation = 1
       * {3,2,1}  // rotation = 2
       * \endcode
       */
      int rotation() const {return rotation_;}

    private:
      int** base_;
      int offset_;
      int length_;
      mutable int rotation_;
    };

    /** 
     * Two vertex sets are equal when their vertices are identical
     * to within a reversal of ordering and possibly a rotation.
     */
    inline bool VertexSet::operator==(const VertexSet& other) const
    {
      int rotation = -1;

      int* p = *base_ + offset_*length_;
      int* op = *(other.base_) + other.offset_*length_;

      for (int i=0; i<length_; i++)
        {
          if (p[i] == op[0])
            {
              rotation = i;
              break;
            }
        }
      if (rotation == -1) return false;

      int n = length_;
      for (int i=0; i<length_; i++)
        {
          if (p[n-((rotation+i+1)%n) - 1] != op[i])
            return false;
        }
      rotation_ = rotation;
      return true;
    }

    /**
     * Return a hash code for the vertex set. Because the vertices 
     * might be presented in
     * permuted order, we need a hash function that is invariant 
     * under permutations.
     */
    inline int VertexSet::hashCode() const
    {
      int rtn = 0;
      int* p = *base_ + offset_*length_;
      int prod = 1;

      for (int i=0; i<length_; i++)
        {
          rtn += p[i];
        }

      return rtn;
    }
  }
}

namespace Teuchos
{
  /** \relates VertexSet */
  inline int hashCode(const Sundance::StdMesh::Internal::VertexSet& v) 
  {return v.hashCode();}

  /** \relates VertexSet */
  inline string toString(const Sundance::StdMesh::Internal::VertexSet& v) 
  {return v.toString();}
}


#endif /* DOXYGEN_DEVELOPER_ONLY */

#endif
