/*
 * SundanceGaussLobattoQuadrature.hpp
 *
 *  Created on: Jan 20, 2011
 *      Author: benk
 */

#ifndef SUNDANCEGAUSSLOBATTOQUADRATURE_HPP_
#define SUNDANCEGAUSSLOBATTOQUADRATURE_HPP_

#include "SundanceDefs.hpp"
#include "PlayaTabs.hpp"
#include "SundanceQuadratureFamilyBase.hpp"

namespace Sundance {

class GaussLobattoQuadrature : public QuadratureFamilyBase {

public:

	/** In this case the */
	GaussLobattoQuadrature(int order);

	/** */
	virtual ~GaussLobattoQuadrature()
	{;}

	/** */
	virtual XMLObject toXML() const;

	/** Describable interface */
	virtual std::string description() const
	{
		return "GaussLobattoQuadrature[order=" + Teuchos::toString(order()) + "]";
	}

	/* handleable boilerplate */
	GET_RCP(QuadratureFamilyStub);

protected:
	/** compute a rule for the reference line cell */
	virtual void getLineRule(Array<Point>& quadPointsL,
			Array<double>& quadWeights) const;

	/** compute a rule for the reference triangle cell */
	virtual void getTriangleRule(Array<Point>& quadPoints,
			Array<double>& quadWeights) const;

	/** compute a rule for the reference quad cell */
	virtual void getQuadRule(Array<Point>& quadPoints,
			Array<double>& quadWeights) const;

	/** compute a rule for the reference tet cell */
	virtual void getTetRule(Array<Point>& quadPoints,
			Array<double>& quadWeights) const ;

	/** compute a rule for the reference brick cell */
	virtual void getBrickRule(Array<Point>& quadPoints,
			Array<double>& quadWeights) const;

	/**
	 * Compute adapted weights according to curve
	 *
	 * @param cellType
	 * @param cellDim
	 * @param cellLID
	 * @param facetIndex
	 * @param mesh
	 * @param globalCurve
	 * @param quadPoints
	 * @param quadWeights
	 * @param changedWeights
	 */
	virtual void getAdaptedWeights(const CellType& cellType, int cellDim,
			int cellLID, int facetIndex, const Mesh& mesh,
			const ParametrizedCurve& globalCurve,
			Array<Point>& quadPoints, Array<double>& quadWeights,
			bool &weightsChanged) const;

	/** Get the weights for one quad */
	virtual void getAdaptedQuadWeights(int cellLID, const Mesh& mesh,
			const ParametrizedCurve& globalCurve, Array<Point>& quadPoints,
			Array<double>& quadWeights, bool& weightsChanged) const;

private:

	/** get the triangle quadrature points for the adaptive integration*/
	void getTriangleQuadPoints(Array<Point>& pnt , Array<double>& weight ) const;

	/**
	 * @param i , the i-th point
	 * @param pnt 1D points where we interpolate
	 * @param x the position where we evaluate */
	inline double evalLagrangePol(int i , Array<Point>& pnt , double x) const{
		int m = pnt.size();
		double p = 1.0;
		for (int j = 0 ; (j <= i-1)&&(j<m) ; j++){
		 p = p *( x - pnt[j][0] )/(pnt[i][0] - pnt[j][0]);
		}
		for (int j = i+1 ; j < m ; j++){
		 p = p *( x - pnt[j][0] )/(pnt[i][0] - pnt[j][0]);
		}
	    return p;
	}

	/** this function applyes quadrature to the lagrange functions <br>
	 * all the coordinates are in the reference coordinates */
	inline void makeInterpolantQuad( double px, double py, double ofx, double ofy,
			                         int nr1DPoints , int nr2DPoints ,
			                         Array<Point>& linePoints ,
			                         Array<Point>& quadPoints ,
			                         Array<double>& pointWeights ,
			                         Array<double>& weightsOut ,
			                         double areFakt ) const {
		int xi,yi;
		double xc,yc, intR;

        SUNDANCE_MSG3(verb_, "px="<<px<<",py="<<py<<",ofx="<<ofx<<",ofy="<<ofy <<
        		" , nr2DPoints:" << nr2DPoints << " , pointWeights.size():" << pointWeights.size());

		for (int i = 0 ; i < nr2DPoints ; i++){
			xi = i % nr1DPoints; yi = i / nr1DPoints; intR = 0.0;
			// for each quadrature point
			for (int q = 0 ; q < pointWeights.size() ; q++){
				xc = px + quadPoints[q][0]*ofx;
				yc = py + quadPoints[q][1]*ofy;
				//SUNDANCE_MSG3(verb_, "q="<<q<<",xc="<<xc<<",yc="<<yc<<",pointWeights[q]=" << pointWeights[q] );
				intR = intR + pointWeights[q] * ( evalLagrangePol(yi , linePoints , yc) * evalLagrangePol(xi , linePoints , xc) );
			}
			intR = intR * ::fabs(ofx*ofy/areFakt);
			weightsOut[i] = weightsOut[i] + intR;
		}
	}

	int nrPointin1D_;

	/** the verbosity of the object*/
	int verb_;

	static int quadsEdgesPoints[4][2];
};

}

#endif /* SUNDANCEGAUSSLOBATTOQUADRATURE_HPP_ */
