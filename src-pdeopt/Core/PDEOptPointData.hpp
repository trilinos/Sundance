#ifndef PDEOPT_POINTDATA_H
#define PDEOPT_POINTDATA_H

#include "Sundance.hpp"
#include <fstream>


namespace Sundance
{
/**
 * PointData is a utility class that helps create the Sundance objects required
 * to include point measurements in an inversion problem. 
 */
class PointData 
{
public:
  /** Create a point data object given a set of sensor locations and
   * corresponding sensor readings */
  PointData(const Array<Point>& locations, const Array<double>& values,
    const double& pointComparisonTolerance);

  /** Create a point data object given a set of sensor locations and
   * a field to probe at those points */
  PointData(const Array<Point>& locations, const Expr& field,
    const double& pointComparisonTolerance);

  /** Read point data from an XML file */
  PointData(const XMLObject& xml, const Mesh& mesh);

    

  /** Return an expr that, when called at one of the measurement locations,
   * returns the value of the measurement there. */
  Expr sensorValues() const {return sensorVals_;}

  /** Return a cell filter that identifies the sensor points */
  CellFilter sensorLocations() const {return sensorLocations_;}


  /** Adjust a set of points to the nearest vertices on a mesh */
  static Array<Point> snapToMesh(const Mesh& mesh, const Array<Point>& locations) ;

  /** Find the vertex nearest to a specified point */
  static Point nearestMeshPoint(const Mesh& mesh, const Point& x);
    
    
private:
  void init(const Array<Point>& locations, const Array<double>& values,
    const double& pointComparisonTolerance);

  Expr sensorVals_;
  CellFilter sensorLocations_;
};



/** 
 * Do lexigraphic comparison on points, including a bit of sloppiness 
 */
class SloppyPointComparitor : public std::binary_function<Point, Point, bool>
{
public: 
  SloppyPointComparitor(const double& tol) : tol_(tol) {;}

  bool operator()(const Point& p1, const Point& p2) const ;
private:
  double tol_;
};


/** 
 * This is a functor that allows the sensor readings to be obtained through
 * a Sundance user-defined expression.
 */
class PointDataExprFunctor : public PointwiseUserDefFunctor0
{
public:
  /** */
  PointDataExprFunctor(const Array<Point>& locations, 
    const Array<double>& values,
    const double& pointComparisonTolerance);

  /** */
  virtual void eval0(const double* in, double* out) const ;

  /** */
  int numArgs() const {return dim_;}

private:
  std::map<Point, double, SloppyPointComparitor> pointToValueMap_;
  int dim_;
};
  
/** 
 * This is a functor that allows a positional cell predicate to test
 * against the locations of the sensors.
 */
class PointDataCellPredicateFunctor : public CellPredicateFunctorBase,
                                      public Playa::Handleable<CellPredicateFunctorBase>
{
public:
  /** */
  PointDataCellPredicateFunctor(const Array<Point>& locations, 
    const double& pointComparisonTolerance);
    
  /** */
  virtual bool operator()(const Point& x) const ;

  GET_RCP(CellPredicateFunctorBase);

private:
  std::set<Point, SloppyPointComparitor> pointSet_;
};
  

}

#endif
