/*
 * SundanceGaussLobattoQuadrature.cpp
 *
 *  Created on: Jan 20, 2011
 *      Author: benk
 */

#include "SundanceGaussLobattoQuadrature.hpp"
#include "PlayaTabs.hpp"

using namespace Sundance;


int GaussLobattoQuadrature::quadsEdgesPoints[4][2] = { {0,1} , {0,2} , {1,3} , {2,3} };


GaussLobattoQuadrature::GaussLobattoQuadrature(int order) :
   QuadratureFamilyBase(order) {
	nrPointin1D_ = order+1;
	verb_ = 0;
}

XMLObject GaussLobattoQuadrature::toXML() const
{
	XMLObject rtn("GaussLobattoQuadrature");
	rtn.addAttribute("order", Teuchos::toString(order()));
	return rtn;
}

void GaussLobattoQuadrature::getLineRule(
		Array<Point>& quadPointsL,
		Array<double>& quadWeights) const {

	int n = order()+1; //todo:

	Array<double> quadPoints;
	quadPoints.resize(n);
	quadPointsL.resize(n);
	quadWeights.resize(n);

	// ================== QUAD POINTS ========================
	switch (n) {
	case 2: { quadPoints[0]=0.0; quadPoints[1] = 1; break; }
	case 3: { quadPoints[0]=0.0; quadPoints[1] = 0.5; quadPoints[2] = 1.0; break; }
	case 4:	{ quadPoints[0]=0.0; quadPoints[1] = 0.5-0.5/::sqrt(5.0); quadPoints[2] = 0.5+0.5/::sqrt(5.0); quadPoints[3] = 1.0; break; }
	case 5:	{ quadPoints[0]=0.0; quadPoints[1] = 0.5-0.5*::sqrt(3.0/7.0); quadPoints[2] = 0.5;
	        quadPoints[3] = 0.5+0.5*::sqrt(3.0/7.0); quadPoints[4] = 1.0; break; }
	case 6: { double t0=::sqrt(7.0);
	        double t1=(7.0+2.0*t0)/21.0;
	        double t2=(7.0-2.0*t0)/21.0;
	        quadPoints[0] = 0; quadPoints[1] = 0.5-0.5*::sqrt(t1); quadPoints[2] = 0.5-0.5*::sqrt(t2);
	        quadPoints[3] = 0.5+0.5*::sqrt(t2); quadPoints[4] = 0.5+0.5*::sqrt(t1); quadPoints[5] = 1.0; break; }
	case 7: {
		quadPoints[0] = 0.00000000000000000000e+00;
		quadPoints[1] = 8.48880518607165179823e-02;
		quadPoints[2] = 2.65575603264642912116e-01;
		quadPoints[3] = 5.00000000000000000000e-01;
		quadPoints[4] = 7.34424396735357087884e-01;
		quadPoints[5] = 9.15111948139283537529e-01;
		quadPoints[6] = 1.00000000000000000000e+00;
	    break; }
	case 8: {
		quadPoints[0] = 0.00000000000000000000e+00;
		quadPoints[1] = 6.41299257451967141819e-02;
		quadPoints[2] = 2.04149909283428854234e-01;
		quadPoints[3] = 3.95350391048760574364e-01;
		quadPoints[4] = 6.04649608951239425636e-01;
		quadPoints[5] = 7.95850090716571090255e-01;
		quadPoints[6] = 9.35870074254803285818e-01;
		quadPoints[7] = 1.00000000000000000000e+00;
		break; }
	case 9: {
		quadPoints[0] = 0.00000000000000000000e+00;
		quadPoints[1] = 5.01210022942699118254e-02;
		quadPoints[2] = 1.61406860244631134016e-01;
		quadPoints[3] = 3.18441268086910922452e-01;
		quadPoints[4] = 5.00000000000000000000e-01;
		quadPoints[5] = 6.81558731913089133059e-01;
		quadPoints[6] = 8.38593139755368865984e-01;
		quadPoints[7] = 9.49878997705730032663e-01;
		quadPoints[8] = 1.00000000000000000000e+00;
	    break; }
	case 10: {
		quadPoints[0] = 0.00000000000000000000e+00;
		quadPoints[1] = 4.02330459167705711820e-02;
		quadPoints[2] = 1.30613067447247432895e-01;
		quadPoints[3] = 2.61037525094777733692e-01;
		quadPoints[4] = 4.17360521166806497373e-01;
		quadPoints[5] = 5.82639478833193447116e-01;
		quadPoints[6] = 7.38962474905222266308e-01;
		quadPoints[7] = 8.69386932552752567105e-01;
		quadPoints[8] = 9.59766954083229428818e-01;
		quadPoints[9] = 1.00000000000000000000e+00;
		break; }
	 case 11: {
		quadPoints[0] = 0.00000000000000000000e+00;
		quadPoints[1] = 3.29992847959704183047e-02;
		quadPoints[2] = 1.07758263168427792511e-01;
		quadPoints[3] = 2.17382336501897477365e-01;
		quadPoints[4] = 3.52120932206530290465e-01;
		quadPoints[5] = 5.00000000000000000000e-01;
		quadPoints[6] = 6.47879067793469709535e-01;
		quadPoints[7] = 7.82617663498102578146e-01;
		quadPoints[8] = 8.92241736831572263000e-01;
		quadPoints[9] = 9.67000715204029637206e-01;
		quadPoints[10] = 1.00000000000000000000e+00;
	    break; }
	 case 12: {
		quadPoints[0] = 0.00000000000000000000e+00;
		quadPoints[1] = 2.75503638885589152707e-02;
		quadPoints[2] = 9.03603391779966846897e-02;
		quadPoints[3] = 1.83561923484069688950e-01;
		quadPoints[4] = 3.00234529517325543502e-01;
		quadPoints[5] = 4.31723533572536233294e-01;
		quadPoints[6] = 5.68276466427463766706e-01;
		quadPoints[7] = 6.99765470482674456498e-01;
		quadPoints[8] = 8.16438076515930255539e-01;
		quadPoints[9] = 9.09639660822003315310e-01;
		quadPoints[10] = 9.72449636111441084729e-01;
		quadPoints[11] = 1.00000000000000000000e+00;
	    break; }
	}

	// transform the array of doubles into an array of Points
	for (int i = 0 ; i < n ; i++){
		quadPointsL[i] = Point(quadPoints[i]);
	}

	// ================== WEIGHTS ========================
	switch (n) {
	case 2:{
		quadWeights[0] = 0.5;
		quadWeights[1] = 0.5;
	    break; }
	case 3:{
		quadWeights[0] = 1.0/6.0; quadWeights[1] = 2.0/3.0; quadWeights[2] = 1.0/6.0;
		break;}
	case 4:{
		quadWeights[0] = 1.0/12.0; quadWeights[1] = 5.0/12.0; quadWeights[2] = 5.0/12.0; quadWeights[3] = 1.0/12.0;
		break;}
	case 5:{
		quadWeights[0] = 0.05; quadWeights[1] = 49.0/180.0; quadWeights[2] = 32.0/90.0; quadWeights[3] = 49.0/180.0; quadWeights[4] = 0.05;
		break;}
	case 6:{
	    double t0=::sqrt(7.0);
	    double t1=(7.0+2.0*t0)/21.0;
	    double t2=(7.0-2.0*t0)/21.0;
	    double k1=(1.0-t0)*(1.0-t0);
	    double k2=(1.0+t0)*(1.0+t0);
	    quadWeights[0] = 1.0/30.0; quadWeights[1] = 0.3/(t1*k1); quadWeights[2] = 0.3/(t2*k2); quadWeights[3] = 0.3/(t2*k2);
	    quadWeights[4] = 0.3/(t1*k1); quadWeights[5] = 1.0/30.0;
	    break;}
	case 7:{
		quadWeights[0] = 2.38095238095238082021e-02;
		quadWeights[1] = 1.38413023680782953928e-01;
		quadWeights[2] = 2.15872690604931305458e-01;
		quadWeights[3] = 2.43809523809523809312e-01;
		quadWeights[4] = 2.15872690604931305458e-01;
		quadWeights[5] = 1.38413023680782953928e-01;
		quadWeights[6] = 2.38095238095238082021e-02;
	    break;}
	case 8:{
		quadWeights[0] = 1.78571428571428561516e-02;
		quadWeights[1] = 1.05352113571753072674e-01;
		quadWeights[2] = 1.70561346241752204156e-01;
		quadWeights[3] = 2.06229397329351860080e-01;
		quadWeights[4] = 2.06229397329351860080e-01;
		quadWeights[5] = 1.70561346241752204156e-01;
		quadWeights[6] = 1.05352113571753072674e-01;
		quadWeights[7] = 1.78571428571428561516e-02;
	    break;}
	case 9:{
		quadWeights[0] = 1.38888888888888881179e-02;
		quadWeights[1] = 8.27476807804027880699e-02;
		quadWeights[2] = 1.37269356250080826198e-01;
		quadWeights[3] = 1.73214255486523083238e-01;
		quadWeights[4] = 1.85759637188208620584e-01;
		quadWeights[5] = 1.73214255486523083238e-01;
		quadWeights[6] = 1.37269356250080826198e-01;
		quadWeights[7] = 8.27476807804027880699e-02;
		quadWeights[8] = 1.38888888888888881179e-02;
	    break;}
	case 10:{
		quadWeights[0] = 1.11111111111111115352e-02;
		quadWeights[1] = 6.66529954255350304271e-02;
		quadWeights[2] = 1.12444671031563220298e-01;
		quadWeights[3] = 1.46021341839841889421e-01;
		quadWeights[4] = 1.63769880591948718829e-01;
		quadWeights[5] = 1.63769880591948718829e-01;
		quadWeights[6] = 1.46021341839841889421e-01;
		quadWeights[7] = 1.12444671031563220298e-01;
		quadWeights[8] = 6.66529954255350304271e-02;
		quadWeights[9] = 1.11111111111111115352e-02;
	    break;}
	case 11:{
		quadWeights[0] = 9.09090909090909046752e-03;
		quadWeights[1] = 5.48061366334974126024e-02;
		quadWeights[2] = 9.35849408901525958715e-02;
		quadWeights[3] = 1.24024052132014145355e-01;
		quadWeights[4] = 1.43439562389503921791e-01;
		quadWeights[5] = 1.50108797727845355574e-01;
		quadWeights[6] = 1.43439562389503921791e-01;
		quadWeights[7] = 1.24024052132014145355e-01;
		quadWeights[8] = 9.35849408901525958715e-02;
		quadWeights[9] = 5.48061366334974126024e-02;
		quadWeights[10] = 9.09090909090909046752e-03;
	    break;}
	case 12:{
		quadWeights[0] = 7.57575757575757596785e-03;
		quadWeights[1] = 4.58422587065981240739e-02;
		quadWeights[2] = 7.89873527821850218711e-02;
		quadWeights[3] = 1.06254208880510653268e-01;
		quadWeights[4] = 1.25637801599600640312e-01;
		quadWeights[5] = 1.35702620455348088591e-01;
		quadWeights[6] = 1.35702620455348088591e-01;
		quadWeights[7] = 1.25637801599600640312e-01;
		quadWeights[8] = 1.06254208880510653268e-01;
		quadWeights[9] = 7.89873527821850218711e-02;
		quadWeights[10] = 4.58422587065981240739e-02;
		quadWeights[11] = 7.57575757575757596785e-03;
	    break;}
   }
}

void GaussLobattoQuadrature::getTriangleRule(
		Array<Point>& quadPoints,
		Array<double>& quadWeights) const {
    // todo: implement this
	SUNDANCE_ERROR("Triangle rule not available for " << toXML());
}

void GaussLobattoQuadrature::getQuadRule(
		Array<Point>& quadPoints,
		Array<double>& quadWeights) const {
	Array<Point> quadPointsLine;
	Array<double> quadWeightsLine;
	// get the line rule
    this->getLineRule( quadPointsLine, quadWeightsLine );

    int nrPointPerAxis = quadPointsLine.size();
    // we simply take the tensor product
    quadPoints.resize(nrPointPerAxis*nrPointPerAxis);
    quadWeights.resize(nrPointPerAxis*nrPointPerAxis);

    int pcount = 0;
    for (int ix = 0 ; ix < nrPointPerAxis ; ix ++){
    	for (int iy = 0 ; iy < nrPointPerAxis ; iy ++){
    		// here we take the tensor product of the
    		quadPoints[pcount] = Point( quadPointsLine[ix][0] , quadPointsLine[iy][0] );
    		quadWeights[pcount] = quadWeightsLine[ix] * quadWeightsLine[iy];
    		pcount++;
    	}
    }
}


void GaussLobattoQuadrature::getTetRule(
		Array<Point>& quadPoints,
		Array<double>& quadWeights) const {
   // todo: implement this
   SUNDANCE_ERROR("Tet rule not available for " << toXML());
}


void GaussLobattoQuadrature::getBrickRule(
		Array<Point>& quadPoints,
		Array<double>& quadWeights) const {
	// todo: implement this
	SUNDANCE_ERROR("Brick rule not available for " << toXML());
}


void GaussLobattoQuadrature::getAdaptedWeights(
		const CellType& cellType, int cellDim,
		int cellLID, int facetIndex, const Mesh& mesh,
		const ParametrizedCurve& globalCurve,
		Array<Point>& quadPoints, Array<double>& quadWeights,
		bool &weightsChanged) const {

	Tabs tabs;

	weightsChanged = true; // todo: this not always might be true to save computations we might test if this is the case
	if (mesh.IsSpecialWeightValid() && mesh.hasSpecialWeight(cellDim, cellLID))
	{
		mesh.getSpecialWeight(cellDim, cellLID, quadWeights);
		//SUNDANCE_MSG3(verb, tabs << "GaussianQuadrature::getAdaptedWeights Found cached weights for cell LID " << cellLID)
		return;
	}
	// if we have no quad points then get them
	if (quadPoints.size() <= 1) getPoints(cellType,  quadPoints, quadWeights);

	// Maximal cell type
	CellType maxCellType = mesh.cellType(mesh.spatialDim());

	switch (maxCellType)
	{
	// We have a triangle based mesh
	case TriangleCell:
	{
		switch (cellType)
		{
		case TriangleCell:
		{
			// todo: implement this
			break;
		}
		default:
#ifndef TRILINOS_7
			SUNDANCE_ERROR ("getAdaptedWeights rule for ACI(Adaptive Cell Integration) not available for submaximal cell " << maxCellType << " in a triangle mesh")
			;
#else
			SUNDANCE_ERROR7("getAdaptedWeights rule for ACI(Adaptive Cell Integration) not available for submaximal cell " << maxCellType << " in a triangle mesh");
#endif
		}
		break;
	}
	case QuadCell:
	{
		switch(cellType)
		{
		case QuadCell:
		{
			// call the method to integrate the Quad
			getAdaptedQuadWeights(cellLID, mesh, globalCurve, quadPoints, quadWeights, weightsChanged);
			break;
		}
		default:
#ifndef TRILINOS_7
		SUNDANCE_ERROR ("getAdaptedWeights rule for ACI(Adaptive Cell Integration) not available for submaximal cell " << maxCellType << " in a triangle mesh")
		;
#else
		SUNDANCE_ERROR7("getAdaptedWeights rule for ACI(Adaptive Cell Integration) not available for submaximal cell " << maxCellType << " in a triangle mesh");
#endif
		}
		break;
	}
	default:
#ifndef TRILINOS_7
		SUNDANCE_ERROR("getAdaptedWeights rule for ACI(Adaptive Cell Integration) not available for cell type " << cellType)
		;
#else
		SUNDANCE_ERROR7("getAdaptedWeights rule for ACI(Adaptive Cell Integration) not available for cell type " << cellType);
#endif
	}

	//store the weights in the mesh
	mesh.setSpecialWeight(cellDim, cellLID, quadWeights);
}


void GaussLobattoQuadrature::getAdaptedQuadWeights(int cellLID, const Mesh& mesh,
		const ParametrizedCurve& globalCurve, Array<Point>& quadPoints,
		Array<double>& quadWeights, bool& weightsChanged) const{

	int maxStack = 1000;
	int maxLevel = 5;
	int stackIndex = 0;
	Tabs tabs;

	Array< Array<Point> > pointStack(maxStack);     // -> the intersection points of each quad cell in the stack (in ref coords)
	Array< Array<Point> > quadPointStack(maxStack); // -> the points which form the quad cell (in ref coords)
	Array< Array<int> > intersectionPointStack(maxStack); // -> which edges are intersected for the actual quad cell
	Array<int> levelStack(maxStack);                      // -> the refinement level of the actual quad cell
	Array<int> refinedStack(maxStack,-1);                 // -> if a quad cell has been replaced then this will be > 0 , -1 otherwise
	Array<int> intersectioncaseStack(maxStack,-1);        // -> the intersection case can be from 0,..,6 , if it is something else then refine
	Array<double> alpha1(maxStack);                       // -> the integration coefficient for the first part of the integral
	Array<double> alpha2(maxStack);                       // -> the integration coefficient for the second part of the integral
	// the actual points of the quad cell and the actual intersection points in real coords
	Array<Point> quadOfPoints(4);
	Array<Point> intesectPoints(8);
	int nrIPoints , tmp;

	// first add the initial quad to the quad stack
	for (int i = 0 ; i < 4 ; i++){
		quadOfPoints[i] = mesh.nodePosition( mesh.facetLID(2,cellLID , 0 , i , nrIPoints ) );
	}

	// first divide the parent cells in smaller quads ,
	pointStack[0].resize(0);
	quadPointStack[0].resize(4);
	quadPointStack[0][0] = Point(0.0,0.0); quadPointStack[0][1] = Point(1.0,0.0);
	quadPointStack[0][2] = Point(0.0,1.0); quadPointStack[0][3] = Point(1.0,1.0);
	intersectionPointStack[0].resize(0);
	levelStack[0] = 0;
	tmp = 0;
	// for each edge test for intersection
	for (int i = 0 ; i < 4 ; i++)
	{
		// for each edge test intersection point
		globalCurve.returnIntersectPoints( quadOfPoints[quadsEdgesPoints[i][0]] , quadOfPoints[quadsEdgesPoints[i][1]],
				 nrIPoints , intesectPoints );
		SUNDANCE_MSG3(verb_, " Test edge: "<< quadOfPoints[quadsEdgesPoints[i][0]] << " -- " << quadOfPoints[quadsEdgesPoints[i][1]] );
		pointStack[0].resize( tmp + nrIPoints );
		intersectionPointStack[0].resize( tmp + nrIPoints );
		for (int j = 0 ; j < nrIPoints ; j++){
			// transform to reference coordinates back
			pointStack[0][tmp+j] = Point( ( intesectPoints[j][0] - quadOfPoints[0][0])/( quadOfPoints[3][0]-quadOfPoints[0][0] ) ,
					                      ( intesectPoints[j][1] - quadOfPoints[0][1])/( quadOfPoints[3][1]-quadOfPoints[0][1] )  );
			intersectionPointStack[0][tmp+j] = i;
			SUNDANCE_MSG3(verb_, " adding Ints Points:"<< tmp+j << " edge:" << intersectionPointStack[0][tmp+j] << " P:"
					<< pointStack[0][tmp+j] << " Real P:" << intesectPoints[j] );
		}
		tmp = tmp + nrIPoints;
	}

	// select the correct case
	intersectioncaseStack[0] = -1;
	if (intersectionPointStack[0].size() == 0){
		intersectioncaseStack[0] = 6; // no intersection point
		alpha1[0] = alpha2[0] = globalCurve.integrationParameter(quadOfPoints[0]);
	}
	if (intersectionPointStack[0].size() == 2){
       if ((intersectionPointStack[0][0] == 0 ) && (intersectionPointStack[0][1] == 1 )) {
    	   intersectioncaseStack[0] = 0;
    	   alpha1[0] = globalCurve.integrationParameter(quadOfPoints[0]);
    	   alpha2[0] = globalCurve.integrationParameter(quadOfPoints[3]);
       }
       if ((intersectionPointStack[0][0] == 0 ) && (intersectionPointStack[0][1] == 2 )) {
    	   intersectioncaseStack[0] = 2;
    	   alpha1[0] = globalCurve.integrationParameter(quadOfPoints[0]);
    	   alpha2[0] = globalCurve.integrationParameter(quadOfPoints[1]);
       }
       if ((intersectionPointStack[0][0] == 0 ) && (intersectionPointStack[0][1] == 3 )) {
    	   alpha1[0] = globalCurve.integrationParameter(quadOfPoints[0]);
    	   alpha2[0] = globalCurve.integrationParameter(quadOfPoints[3]);
    	   if (pointStack[0][0][0] > pointStack[0][1][0]) intersectioncaseStack[0] = 41;
    	   else intersectioncaseStack[0] = 42;
       }
       if ((intersectionPointStack[0][0] == 1 ) && (intersectionPointStack[0][1] == 2 )) {
    	   alpha1[0] = globalCurve.integrationParameter(quadOfPoints[0]);
    	   alpha2[0] = globalCurve.integrationParameter(quadOfPoints[3]);
    	   if (pointStack[0][0][1] > pointStack[0][1][1]) intersectioncaseStack[0] = 51;
    	   else intersectioncaseStack[0] = 52;
       }
       if ((intersectionPointStack[0][0] == 1 ) && (intersectionPointStack[0][1] == 3 )) {
    	   alpha1[0] = globalCurve.integrationParameter(quadOfPoints[0]);
    	   alpha2[0] = globalCurve.integrationParameter(quadOfPoints[2]);
    	   intersectioncaseStack[0] = 1;
       }
       if ((intersectionPointStack[0][0] == 2 ) && (intersectionPointStack[0][1] == 3 )) {
    	   alpha1[0] = globalCurve.integrationParameter(quadOfPoints[0]);
    	   alpha2[0] = globalCurve.integrationParameter(quadOfPoints[3]);
    	   intersectioncaseStack[0] = 3;
       }
	}
	SUNDANCE_MSG3(verb_, " intersectioncaseStack[0]:"<< intersectioncaseStack[0] << " , alpha1[0]:" << alpha1[0]
	                       << " , alpha2[0]" << alpha2[0]);
	stackIndex = 1;
	//quadsEdgesPoints


	// the criterion for division is if the are after the refinement is different significantly
	// then before the division
	// put these quads with the intersection points in one list
	// todo: make here a while cycle, till we arrive


	// get all the possible weights and quadrature points
	Array<Point> linePoints;
	Array<double> lineWeights;
	getLineRule( linePoints , lineWeights);
	Array<Point> quadQuadPoints;
	Array<double> quadQuadWeights;
	getQuadRule( quadQuadPoints , quadQuadWeights);
	Array<Point> trianglePoints;
	Array<double> triangleWeights;
	getTriangleQuadPoints( trianglePoints , triangleWeights);
	double summWeights = 0.0;

	int nr1DPoints = linePoints.size();
	int nr2DPoints = quadQuadPoints.size();

	//those must be equal
	TEST_FOR_EXCEPTION( quadPoints.size() != quadQuadPoints.size() , std::runtime_error,
			"quadPoints.size() != quadQuadPoints.size() , size1:" << quadPoints.size() << " , size2:" << quadQuadPoints.size());
	for (int q = 0; q < nr2DPoints; q++) {
		SUNDANCE_MSG3(verb_, " Quad point quadWeights["<<q<<"]="<<quadWeights[q]);
		summWeights = summWeights + quadWeights[q];
		quadWeights[q] = 0.0;
	}
	SUNDANCE_MSG3(verb_, " Summ old weights = " << summWeights );

    // for all elements in the list make the quadrature
	for (int listI = 0 ; listI < stackIndex ; listI++){
		// look at this quad only when it is not refined
		if (refinedStack[listI] < 0 ){
			SUNDANCE_MSG3(verb_, tabs << "getAdaptedQuadWeights Integrate quad listI:" << listI <<
					",intersectioncaseStack[listI]:" << intersectioncaseStack[listI]);
			// ========== calculate first the integral over the whole quad
			Array<double> tmpWeightsQuad( quadWeights.size() , 0.0 );
			double  ofx , ofy , px , py;
			// integrate the whole quad, which result will be later used
            ofx = (quadPointStack[listI][1][0] - quadPointStack[listI][0][0]);
            ofy = (quadPointStack[listI][2][1] - quadPointStack[listI][0][1]);
            px = quadPointStack[listI][0][0]; py = quadPointStack[listI][0][1];
            // call the function which makes the quadrature
            makeInterpolantQuad( px, py, ofx, ofy, nr1DPoints , nr2DPoints ,linePoints ,
            		             quadQuadPoints , quadQuadWeights ,tmpWeightsQuad , 1.0 );
            SUNDANCE_MSG3(verb_, tabs << " end quadring the whole quad, now quadring the remaining parts " );
			switch (intersectioncaseStack[listI]){
			// =================================================================================================================
			case 0:{
				Array<double> tmpWeightsTriangle( quadWeights.size() , 0.0 );
				// integrate the triangle
                ofx = (pointStack[listI][0][0] - quadPointStack[listI][0][0]);
                ofy = (pointStack[listI][1][1] - quadPointStack[listI][0][1]);
                px = quadPointStack[listI][0][0]; py = quadPointStack[listI][0][1];
                // call the function which makes the quadrature
	            makeInterpolantQuad( px, py, ofx, ofy, nr1DPoints , nr2DPoints ,linePoints ,
	            		trianglePoints , triangleWeights ,tmpWeightsTriangle , 2.0 );
				// now add up the two diferent arreas
				for (int i = 0 ; i < nr2DPoints ; i++){
					SUNDANCE_MSG3(verb_, tabs << "i=" << i << ",Wquad[i]:" << tmpWeightsQuad[i] << " , Wtriag[i]:"<< tmpWeightsTriangle[i] );
					quadWeights[i] = alpha1[listI]*tmpWeightsTriangle[i] + alpha2[listI]*( tmpWeightsQuad[i] - tmpWeightsTriangle[i] );
				}
			    break;}
			// =================================================================================================================
			case 1:{
				Array<double> tmpWeightsTriangle( quadWeights.size() , 0.0 );
				// integrate the triangle
                ofx = (pointStack[listI][1][0] - quadPointStack[listI][0][0]);
                ofy = (pointStack[listI][0][1] - quadPointStack[listI][2][1]);
                px = quadPointStack[listI][2][0]; py = quadPointStack[listI][2][1];
                // call the function which makes the quadrature
	            makeInterpolantQuad( px, py, ofx, ofy, nr1DPoints , nr2DPoints ,linePoints ,
	            		trianglePoints , triangleWeights ,tmpWeightsTriangle , 2.0 );
				// now add up the two different areas
				for (int i = 0 ; i < nr2DPoints ; i++){
					SUNDANCE_MSG3(verb_, tabs << "i=" << i << ",Wquad[i]:" << tmpWeightsQuad[i] << " , Wtriag[i]:"<< tmpWeightsTriangle[i] );
					quadWeights[i] = alpha2[listI]*tmpWeightsTriangle[i] + alpha1[listI]*( tmpWeightsQuad[i] - tmpWeightsTriangle[i] );
				}
			    break;}
			// =================================================================================================================
			case 2:{
				Array<double> tmpWeightsTriangle( quadWeights.size() , 0.0 );
				// integrate the triangle
                ofx = (pointStack[listI][0][0] - quadPointStack[listI][1][0]);
                ofy = (pointStack[listI][1][1] - quadPointStack[listI][1][1]);
                px = quadPointStack[listI][1][0]; py = quadPointStack[listI][1][1];
                // call the function which makes the quadrature
	            makeInterpolantQuad( px, py, ofx, ofy, nr1DPoints , nr2DPoints ,linePoints ,
	            		trianglePoints , triangleWeights ,tmpWeightsTriangle , 2.0 );
				// now add up the two diferent arreas
				for (int i = 0 ; i < nr2DPoints ; i++){
					SUNDANCE_MSG3(verb_, tabs << "i=" << i << ",Wquad[i]:" << tmpWeightsQuad[i] << " , Wtriag[i]:"<< tmpWeightsTriangle[i] );
					quadWeights[i] = alpha2[listI]*tmpWeightsTriangle[i] + alpha1[listI]*( tmpWeightsQuad[i] - tmpWeightsTriangle[i] );
				}
			    break;}
			// =================================================================================================================
			case 3:{
				Array<double> tmpWeightsTriangle( quadWeights.size() , 0.0 );
				// integrate the triangle
                ofx = (pointStack[listI][1][0] - quadPointStack[listI][3][0]);
                ofy = (pointStack[listI][0][1] - quadPointStack[listI][3][1]);
                px = quadPointStack[listI][3][0]; py = quadPointStack[listI][3][1];
                // call the function which makes the quadrature
	            makeInterpolantQuad( px, py, ofx, ofy, nr1DPoints , nr2DPoints ,linePoints ,
	            		trianglePoints , triangleWeights ,tmpWeightsTriangle , 2.0 );
				// now add up the two diferent arreas
				for (int i = 0 ; i < nr2DPoints ; i++){
					SUNDANCE_MSG3(verb_, tabs << "i=" << i << ",Wquad[i]:" << tmpWeightsQuad[i] << " , Wtriag[i]:"<< tmpWeightsTriangle[i] );
					quadWeights[i] = alpha2[listI]*tmpWeightsTriangle[i] + alpha1[listI]*( tmpWeightsQuad[i] - tmpWeightsTriangle[i] );
				}
			    break;}
			// =================================================================================================================
			case 41: case 42:{
				// integrate the quad
				Array<double> tmpWeightsQuad2( quadWeights.size() , 0.0 );
				if (intersectioncaseStack[listI] == 41){
					ofx = ( pointStack[listI][1][0] - quadPointStack[listI][0][0]);
				} else {
					ofx = ( pointStack[listI][0][0] - quadPointStack[listI][0][0]);
				}
	            ofy = (quadPointStack[listI][2][1] - quadPointStack[listI][0][1]);
	            px = quadPointStack[listI][0][0]; py = quadPointStack[listI][0][1];
	            // call the function which makes the quadrature
	            makeInterpolantQuad( px, py, ofx, ofy, nr1DPoints , nr2DPoints ,linePoints ,
	            		quadQuadPoints , quadQuadWeights ,tmpWeightsQuad2 , 1.0 );
				// integrate the triangle
				Array<double> tmpWeightsTriangle( quadWeights.size() , 0.0 );
				if (intersectioncaseStack[listI] == 41){
					ofx = ( pointStack[listI][0][0] - pointStack[listI][1][0]);
	                ofy = ( quadPointStack[listI][3][1] - quadPointStack[listI][0][1] );
					px = pointStack[listI][1][0]; py = pointStack[listI][0][1];
				} else {
					ofx = ( pointStack[listI][1][0] - pointStack[listI][0][0]);
					ofy = -( quadPointStack[listI][3][1] - quadPointStack[listI][0][1] );
					px = pointStack[listI][0][0]; py = pointStack[listI][1][1];
				}
				// call the function which makes the quadrature
	            makeInterpolantQuad( px, py, ofx, ofy, nr1DPoints , nr2DPoints ,linePoints ,
	            		trianglePoints , triangleWeights ,tmpWeightsTriangle , 2.0 );
				// now add up the two different areas
				for (int i = 0 ; i < nr2DPoints ; i++){
					SUNDANCE_MSG3(verb_, tabs << "i=" << i << " , Wquad[i]:" << tmpWeightsQuad[i] <<
							" , Wquad2[i]:" << tmpWeightsQuad2[i] << " , Wtriag[i]:"<< tmpWeightsTriangle[i] );
					quadWeights[i] = alpha1[listI]*(tmpWeightsQuad2[i]+tmpWeightsTriangle[i]) +
							         alpha2[listI]*( tmpWeightsQuad[i] - tmpWeightsQuad2[i] - tmpWeightsTriangle[i]);
				}
			    break;}
			// =================================================================================================================
			case 51: case 52:{
				// integrate the quad
				Array<double> tmpWeightsQuad2( quadWeights.size() , 0.0 );
				if (intersectioncaseStack[listI] == 52){
					ofy = ( pointStack[listI][0][1] - quadPointStack[listI][0][1]);
				} else {
					ofy = ( pointStack[listI][1][1] - quadPointStack[listI][0][1]);
				}
	            ofx = (quadPointStack[listI][1][0] - quadPointStack[listI][0][0]);
	            px = quadPointStack[listI][0][0]; py = quadPointStack[listI][0][1];
	            // call the function which makes the quadrature
	            makeInterpolantQuad( px, py, ofx, ofy, nr1DPoints , nr2DPoints ,linePoints ,
	            		quadQuadPoints , quadQuadWeights ,tmpWeightsQuad2 , 1.0 );
				// integrate the triangle
				Array<double> tmpWeightsTriangle( quadWeights.size() , 0.0 );
				if (intersectioncaseStack[listI] == 51){
	                ofx = ( quadPointStack[listI][3][0] - quadPointStack[listI][0][0] );
					ofy = ( pointStack[listI][0][1] - pointStack[listI][1][1]);
					px = pointStack[listI][0][0]; py = pointStack[listI][1][1];
				} else {
					ofx = -( quadPointStack[listI][3][0] - quadPointStack[listI][0][0] );
					ofy = ( pointStack[listI][1][1] - pointStack[listI][0][1]);
					px = pointStack[listI][1][0]; py = pointStack[listI][0][1];
				}
				// call the function which makes the quadrature
	            makeInterpolantQuad( px, py, ofx, ofy, nr1DPoints , nr2DPoints ,linePoints ,
	            		trianglePoints , triangleWeights ,tmpWeightsTriangle , 2.0 );
				// now add up the two different areas
				for (int i = 0 ; i < nr2DPoints ; i++){
					SUNDANCE_MSG3(verb_, tabs << "i=" << i << " , Wquad[i]:" << tmpWeightsQuad[i] <<
							" , Wquad2[i]:" << tmpWeightsQuad2[i] << " , Wtriag[i]:"<< tmpWeightsTriangle[i] );
					quadWeights[i] = alpha1[listI]*(tmpWeightsQuad2[i]+tmpWeightsTriangle[i]) +
							         alpha2[listI]*(tmpWeightsQuad[i] - tmpWeightsQuad2[i] - tmpWeightsTriangle[i]);
				}
			    break;}
			// =================================================================================================================
			case 6:{
				// no intersection point just integrate the quad elem
				for (int i = 0 ; i < nr2DPoints ; i++){
				   SUNDANCE_MSG3(verb_, tabs << "i=" << i << " , Wquad[i]:" << tmpWeightsQuad[i] );
                   quadWeights[i] = quadWeights[i] + alpha1[listI]*tmpWeightsQuad[i];
				}
			    break;}
			default:{
				// throw error
				TEST_FOR_EXCEPTION( true , std::runtime_error , "Quad cell not integrable:" << intersectioncaseStack[listI]);
			    break; }
			}
		}
	}
	// just print the weights
	summWeights = 0.0;
	for (int q = 0; q < nr2DPoints; q++) {
		summWeights = summWeights + quadWeights[q];
		SUNDANCE_MSG3(verb_, " New weights quadWeights["<<q<<"]="<<quadWeights[q]);
	}
	SUNDANCE_MSG3(verb_, " Summ new weights = " << summWeights );
}

void GaussLobattoQuadrature::getTriangleQuadPoints(Array<Point>& pnt  ,Array<double>& weight ) const{
	// we took this directly from the Matlab code of Prof.Ulbrich
	int order = nrPointin1D_;
	if (order==2){
	  pnt.resize(3); weight.resize(3);
	  pnt[0] = Point(0.50000000000000000000 , 0.00000000000000000000);
	  pnt[1] = Point(0.50000000000000000000 , 0.50000000000000000000);
	  pnt[2] = Point(0.00000000000000000000 , 0.50000000000000000000);
	  weight[0] = 0.33333333333333333333;
	  weight[1] = 0.33333333333333333333;
	  weight[2] = 0.33333333333333333333;
	} else if (order==3){
		pnt.resize(4); weight.resize(4);
		pnt[0] = Point(0.33333333333333333333 , 0.33333333333333333333);
		pnt[1] = Point(0.60000000000000000000 , 0.20000000000000000000);
		pnt[2] = Point(0.20000000000000000000 , 0.60000000000000000000);
		pnt[3] = Point(0.20000000000000000000 , 0.20000000000000000000);
		weight[0] = -0.56250000000000000000;
		weight[1] = 0.52083333333333333333;
		weight[2] = 0.52083333333333333333;
		weight[3] = 0.52083333333333333333;
	} else if (order==4) {
		pnt.resize(6); weight.resize(6);
		pnt[0] = Point(0.816847572980459 , 0.091576213509771);
		pnt[1] = Point(0.091576213509771 , 0.816847572980459);
		pnt[2] = Point(0.091576213509771 , 0.091576213509771);
		pnt[3] = Point(0.108103018168070 , 0.445948490915965);
		pnt[4] = Point(0.445948490915965 , 0.108103018168070);
		pnt[5] = Point(0.445948490915965 , 0.445948490915965);
		weight[0] = 0.109951743655322;
		weight[1] = 0.109951743655322;
		weight[2] = 0.109951743655322;
		weight[3] = 0.223381589678011;
		weight[4] = 0.223381589678011;
		weight[5] = 0.223381589678011;
	} else if (order==5) {
		pnt.resize(7); weight.resize(7);
		pnt[0] = Point(0.33333333333333333 , 0.33333333333333333);
		pnt[1] = Point(0.79742698535308720 , 0.10128650732345633);
		pnt[2] = Point(0.10128650732345633 , 0.79742698535308720);
		pnt[3] = Point(0.10128650732345633 , 0.10128650732345633);
		pnt[4] = Point(0.05971587178976981 , 0.47014206410511505);
		pnt[5] = Point(0.47014206410511505 , 0.05971587178976981);
		pnt[6] = Point(0.47014206410511505 , 0.47014206410511505);
		weight[0] = 0.22500000000000000;
		weight[1] = 0.12593918054482717;
		weight[2] = 0.12593918054482717;
		weight[3] = 0.12593918054482717;
		weight[4] = 0.13239415278850616;
		weight[5] = 0.13239415278850616;
		weight[6] = 0.13239415278850616;
	} else if (order==6) {
		pnt.resize(9); weight.resize(9);
		pnt[0] = Point(0.124949503233232 , 0.437525248383384);
		pnt[1] = Point(0.437525248383384 , 0.124949503233232);
		pnt[2] = Point(0.437525248383384 , 0.437525248383384);
		pnt[3] = Point(0.797112651860071 , 0.165409927389841);
		pnt[4] = Point(0.797112651860071 , 0.037477420750088);
		pnt[5] = Point(0.165409927389841 , 0.797112651860071);
		pnt[6] = Point(0.165409927389841 , 0.037477420750088);
		pnt[7] = Point(0.037477420750088 , 0.797112651860071);
		pnt[8] = Point(0.037477420750088 , 0.165409927389841);
		weight[0] = 0.205950504760887;
		weight[1] = 0.205950504760887;
		weight[2] = 0.205950504760887;
		weight[3] = 0.063691414286223;
		weight[4] = 0.063691414286223;
		weight[5] = 0.063691414286223;
		weight[6] = 0.063691414286223;
		weight[7] = 0.063691414286223;
		weight[8] = 0.063691414286223;
	} else if (order==7) {
		pnt.resize(13); weight.resize(13);
		pnt[0] = Point(0.333333333333333 , 0.333333333333333);
		pnt[1] = Point(0.479308067841923 , 0.260345966079038);
		pnt[2] = Point(0.260345966079038 , 0.479308067841923);
		pnt[3] = Point(0.260345966079038 , 0.260345966079038);
		pnt[4] = Point(0.869739794195568 , 0.065130102902216);
		pnt[5] = Point(0.065130102902216 , 0.869739794195568);
		pnt[6] = Point(0.065130102902216  ,0.065130102902216);
		pnt[7] = Point(0.638444188569809 , 0.312865496004875);
		pnt[8] = Point(0.638444188569809 , 0.048690315425316);
		pnt[9] = Point(0.312865496004875 , 0.638444188569809);
		pnt[10] = Point(0.312865496004875 , 0.048690315425316);
		pnt[11] = Point(0.048690315425316 , 0.638444188569809);
		pnt[12] = Point(0.048690315425316 , 0.312865496004875);
		weight[0] = -0.149570044467670;
		weight[1] = 0.175615257433204;
		weight[2] = 0.175615257433204;
		weight[3] = 0.175615257433204;
		weight[4] = 0.053347235608839;
		weight[5] = 0.053347235608839;
		weight[6] = 0.053347235608839;
		weight[7] = 0.077113760890257;
		weight[8] = 0.077113760890257;
		weight[9] = 0.077113760890257;
		weight[10] = 0.077113760890257;
		weight[11] = 0.077113760890257;
		weight[12] = 0.077113760890257;
    } else if (order==8) {
		pnt.resize(19); weight.resize(19);
		pnt[0] = Point(0.3333333333333333 , 0.3333333333333333);
		pnt[1] = Point(0.7974269853530872 , 0.1012865073234563);
		pnt[2] = Point(0.1012865073234563 , 0.7974269853530872);
		pnt[3] = Point(0.1012865073234563 , 0.1012865073234563);
		pnt[4] = Point(0.0597158717897698 , 0.4701420641051151);
		pnt[5] = Point(0.4701420641051151 , 0.0597158717897698);
		pnt[6] = Point(0.4701420641051151 , 0.4701420641051151);
		pnt[7] = Point(0.5357953464498992 , 0.2321023267750504);
		pnt[8] = Point(0.2321023267750504 , 0.5357953464498992);
		pnt[9] = Point(0.2321023267750504 , 0.2321023267750504);
		pnt[10] = Point(0.9410382782311209 , 0.0294808608844396);
		pnt[11] = Point(0.0294808608844396 , 0.9410382782311209);
		pnt[12] = Point(0.0294808608844396 , 0.0294808608844396);
		pnt[13] = Point(0.7384168123405100 , 0.2321023267750504);
		pnt[14] = Point(0.7384168123405100 , 0.0294808608844396);
		pnt[15] = Point(0.2321023267750504 , 0.7384168123405100);
		pnt[16] = Point(0.2321023267750504 , 0.0294808608844396);
		pnt[17] = Point(0.0294808608844396 , 0.7384168123405100);
		pnt[18] = Point(0.0294808608844396 , 0.2321023267750504);
		weight[0] = 0.0378610912003147;
		weight[1] = 0.0376204254131829;
		weight[2] = 0.0376204254131829;
		weight[3] = 0.0376204254131829;
		weight[4] = 0.0783573522441174;
		weight[5] = 0.0783573522441174;
		weight[6] = 0.0783573522441174;
		weight[7] = 0.1162714796569659;
		weight[8] = 0.1162714796569659;
		weight[9] = 0.1162714796569659;
		weight[10] = 0.0134442673751655;
		weight[11] = 0.0134442673751655;
		weight[12] = 0.0134442673751655;
		weight[13] = 0.0375097224552317;
		weight[14] = 0.0375097224552317;
		weight[15] = 0.0375097224552317;
		weight[16] = 0.0375097224552317;
		weight[17] = 0.0375097224552317;
		weight[18] = 0.0375097224552317;
    } else if (order == 9) {
		pnt.resize(19); weight.resize(19);
		pnt[0] = Point(0.33333333333333331     ,  0.33333333333333331);
		pnt[1] = Point(2.06349616025259287E-002,  0.48968251919873701);
		pnt[2] = Point(0.48968251919873701     ,  2.06349616025259287E-002);
		pnt[3] = Point(0.48968251919873701      , 0.48968251919873701);
		pnt[4] = Point(0.12582081701412900     ,  0.43708959149293553);
		pnt[5] = Point(0.43708959149293553     ,  0.12582081701412900);
		pnt[6] = Point(0.43708959149293553     ,  0.43708959149293553);
		pnt[7] = Point(0.62359292876193562     ,  0.18820353561903219);
		pnt[8] = Point(0.18820353561903219     ,  0.62359292876193562);
		pnt[9] = Point(0.18820353561903219     ,  0.18820353561903219);
		pnt[10] = Point(0.91054097321109406     ,  4.47295133944529688E-002);
		pnt[11] = Point(4.47295133944529688E-002,  0.91054097321109406);
		pnt[12] = Point(4.47295133944529688E-002,  4.47295133944529688E-002);
		pnt[13] = Point(0.74119859878449801     ,  3.68384120547362581E-002);
		pnt[14] = Point(0.74119859878449801     ,  0.22196298916076573);
		pnt[15] = Point(3.68384120547362581E-002,  0.74119859878449801);
		pnt[16] = Point(3.68384120547362581E-002,  0.22196298916076573);
		pnt[17] = Point(0.22196298916076573     ,  0.74119859878449801);
		pnt[18] = Point(0.22196298916076573     ,  3.68384120547362581E-002 );
		weight[0] = 9.71357962827961025E-002;
		weight[1] = 3.13347002271398278E-002;
		weight[2] = 3.13347002271398278E-002;
		weight[3] = 3.13347002271398278E-002;
		weight[4] = 7.78275410047754301E-002;
		weight[5] = 7.78275410047754301E-002;
		weight[6] = 7.78275410047754301E-002;
		weight[7] = 7.96477389272090969E-002;
		weight[8] = 7.96477389272090969E-002;
		weight[9] = 7.96477389272090969E-002;
		weight[10] = 2.55776756586981006E-002;
		weight[11] = 2.55776756586981006E-002;
		weight[12] = 2.55776756586981006E-002;
		weight[13] = 4.32835393772893970E-002;
		weight[14] = 4.32835393772893970E-002;
		weight[15] = 4.32835393772893970E-002;
		weight[16] = 4.32835393772893970E-002;
		weight[17] = 4.32835393772893970E-002;
		weight[18] = 4.32835393772893970E-002;
    } else if (order<=11) {
		pnt.resize(28); weight.resize(28);
		pnt[0] = Point(0.33333333333333333 , 0.333333333333333333);
		pnt[1] = Point(0.9480217181434233  , 0.02598914092828833);
		pnt[2] = Point(0.02598914092828833 , 0.9480217181434233);
		pnt[3] = Point(0.02598914092828833 , 0.02598914092828833);
		pnt[4] = Point(0.8114249947041546  , 0.09428750264792270);
		pnt[5] = Point(0.09428750264792270 , 0.8114249947041546);
		pnt[6] = Point(0.09428750264792270 , 0.09428750264792270);
		pnt[7] = Point(0.01072644996557060 , 0.4946367750172147);
		pnt[8] = Point(0.4946367750172147  , 0.01072644996557060);
		pnt[9] = Point(0.4946367750172147  , 0.4946367750172147);
		pnt[10] = Point(0.5853132347709715  , 0.2073433826145142);
		pnt[11] = Point(0.2073433826145142  , 0.5853132347709715);
		pnt[12] = Point(0.2073433826145142  , 0.2073433826145142);
		pnt[13] = Point(0.1221843885990187  , 0.4389078057004907);
		pnt[14] = Point(0.4389078057004907  , 0.1221843885990187);
		pnt[15] = Point(0.4389078057004907  , 0.4389078057004907);
		pnt[16] = Point(0.6779376548825902  , 0.04484167758913055);
		pnt[17] = Point(0.6779376548825902  , 0.27722066752827925);
		pnt[18] = Point(0.04484167758913055 , 0.6779376548825902);
		pnt[19] = Point(0.04484167758913055 , 0.27722066752827925);
		pnt[20] = Point(0.27722066752827925 , 0.6779376548825902);
		pnt[21] = Point(0.27722066752827925 , 0.04484167758913055);
		pnt[22] = Point(0.8588702812826364  , 0.00000000000000000);
		pnt[23] = Point(0.8588702812826364  , 0.1411297187173636);
		pnt[24] = Point(0.0000000000000000  , 0.8588702812826364);
		pnt[25] = Point(0.0000000000000000  , 0.1411297187173636);
		pnt[26] = Point(0.1411297187173636  , 0.8588702812826364);
		pnt[27] = Point(0.1411297187173636  , 0.0000000000000000);
		weight[0] = 0.08797730116222190;
		weight[1] = 0.008744311553736190;
		weight[2] = 0.008744311553736190;
		weight[3] = 0.008744311553736190;
		weight[4] = 0.03808157199393533;
		weight[5] = 0.03808157199393533;
		weight[6] = 0.03808157199393533;
		weight[7] = 0.01885544805613125;
		weight[8] = 0.01885544805613125;
		weight[9] = 0.01885544805613125;
		weight[10] = 0.07215969754474100;
		weight[11] = 0.07215969754474100;
		weight[12] = 0.07215969754474100;
		weight[13] = 0.06932913870553720;
		weight[14] = 0.06932913870553720;
		weight[15] = 0.06932913870553720;
		weight[16] = 0.04105631542928860;
		weight[17] = 0.04105631542928860;
		weight[18] = 0.04105631542928860;
		weight[19] = 0.04105631542928860;
		weight[20] = 0.04105631542928860;
		weight[21] = 0.04105631542928860;
		weight[22] = 0.007362383783300573;
		weight[23] = 0.007362383783300573;
		weight[24] = 0.007362383783300573;
		weight[25] = 0.007362383783300573;
		weight[26] = 0.007362383783300573;
		weight[27] = 0.007362383783300573;
    } else if (order<=13) {
    	//
    }
}
