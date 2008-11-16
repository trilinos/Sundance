#ifndef SUNDANCERIVARAMESH_HPP
#define SUNDANCERIVARAMESH_HPP

#include "SundanceDefs.hpp"
#include "SundanceRivaraElement.hpp"
#include <stack>
#include "SundanceMap.hpp"
#include "SundanceIncrementallyCreatableMesh.hpp"
#include "SundanceRivaraElementIterator.hpp"


namespace SundanceStdMesh
{
namespace Rivara
{
using SundanceUtils::Map;
using std::stack;


class RivaraMesh 
{
public:
  RivaraMesh(int dim, const MPIComm& comm);

  int addNode(const RefCountPtr<Node>& node);
  int addVertex(int globalIndex, const Point& x, int ownerProcID, int label);

  void addElement(const RefCountPtr<Element>& tri);
  int addElement(int globalIndex, const Array<int>& vertexGIDs, int ownerProc,
    int label);

  RefCountPtr<Edge> tryEdge(const RefCountPtr<Node>& a,
    const RefCountPtr<Node>& b,
    int& edgeSign);

  RefCountPtr<Face> tryFace(const RefCountPtr<Node>& a,
    const RefCountPtr<Node>& b,
    const RefCountPtr<Node>& c);

  const RefCountPtr<Face>& getFace(const RefCountPtr<Node>& a,
    const RefCountPtr<Node>& b,
    const RefCountPtr<Node>& c) const ;

  const RefCountPtr<Node>& node(int i) const {return nodes_[i];}

  int numNodes() const {return nodes_.length();}

  stack<Element*>& refinementSet()
    {return refinementSet_;}

  stack<double>& refinementAreas()
    {return refinementAreas_;}

  void refine();

  ElementIterator iterator() const ;

  friend class ElementIterator;

  RefCountPtr<Element> element(int i) const {return elements_[i];}

  int numRootElements() const {return elements_.length();}

  int numElements() const ;

  int spatialDim() const ;

  int& nextGID() {return nextGID_;}

  int nextGID() const {return nextGID_;}
private:

  int spatialDim_;
  
  int nextGID_;

  Array<RefCountPtr<Node> > nodes_;

  Array<RefCountPtr<Edge> > edges_;

  Array<RefCountPtr<Face> > faces_;

  Array<RefCountPtr<Element> > elements_;

  Array<Map<int, int> > nodeToEdgeMap_;

  Map<FaceNodes, int> faceToLIDMap_;

  stack<Element*> refinementSet_;

  stack<double> refinementAreas_;
};
}
}

#endif
