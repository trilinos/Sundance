/*
 * SundanceTrapesoidQuadrature.hpp
 *
 *  Created on: Jan 20, 2011
 *      Author: benk
 */

#ifndef SUNDANCETRAPESOIDQUADRATURE_HPP_
#define SUNDANCETRAPESOIDQUADRATURE_HPP_

#include "SundanceDefs.hpp"
#include "SundanceQuadratureFamilyBase.hpp"

namespace Sundance {

class TrapesoidQuadrature : public QuadratureFamilyBase {

public:

	/** In this case the */
	TrapesoidQuadrature(int resolution);

	/** */
	virtual ~TrapesoidQuadrature()
	{;}

	/** */
	virtual XMLObject toXML() const;

	/** Describable interface */
	virtual std::string description() const
	{
		return "TrapesoidQuadrature[order=" + Teuchos::toString(order()) + "]";
	}

	/* handleable boilerplate */
	GET_RCP(QuadratureFamilyStub);

protected:
	/** compute a rule for the reference line cell */
	virtual void getLineRule(Array<Point>& quadPoints,
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

private:

	/** the homogeneous resolution */
	int resolution_;

};

}

#endif /* SUNDANCETRAPESOIDQUADRATURE_HPP_ */
