/* @HEADER@ */

/* @HEADER@ */

#ifndef SUNDANCE_VERTEX_SORT_H
#define SUNDANCE_VERTEX_SORT_H


#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"


namespace Sundance
{
using Teuchos::Array;

/** */
template <class T> inline
void insertionSort(Teuchos::Array<T>& A)
{
  int N = A.size();
  for (int i=1; i<N; i++)
  {
    T val = A[i];
    int j = i-1;
    bool done = false;
    while ( !done )
    {
      if (A[j] > val)
      {
        A[j+1]=A[j];
        j=j-1;
        if ( j<0 ) done = true;
      }
      else done = true;
      A[j+1]=val;
    }
  }
}
 
/* Sort, returning the key associated with the permutation */
void vertexSort(Array<int>& verts, int* key);

/** Return the permutation that produces the specified key */
void getKeyedPerm(int key, Array<int>& digits);

/** Compute base^N */
int iPow(int base, int n);

/** */
int exFacetIndexToUFCFacetIndex(int meshDim, int permKey,
  int exFacetID);

/** */
int ufcFacetIndexToExFacetIndex(int meshDim, int ufcFacetID);

/** */
int exVertPosToUFCVertPos(int meshDim, int permKey, int exVertPos);

/** */
int mapExSideToMissingVertex(int dim, int exFaceID);

/** */
Array<int> exSideVertPos(int dim, int f);
/** */
Array<int> ufcSideVertPos(int dim, int f);

/** */
Array<int> ufcSide(int f, const Array<int>& verts);

/** */
Array<int> exSide(int f, const Array<int>& verts);


}


#endif
