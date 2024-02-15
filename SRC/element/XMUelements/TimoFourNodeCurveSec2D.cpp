/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.38 $
// $Date: 2010-06-01 23:41:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/dispBeamColumn/TimodispBeamColumn.cpp,v $

// Written: MHS
// Created: Feb 2001
//
// Description: This file contains the class definition for TimoFourNodeCurveSec2D.

#include <TimoFourNodeCurveSec2D.h>
#include <Node.h>
#include <SectionForceDeformation.h>
#include <CrdTransf.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <CompositeResponse.h>
#include <math.h>
#include <ElementalLoad.h>
#include <ElementAPI.h>
#include <LegendreBeamIntegration.h>

Matrix TimoFourNodeCurveSec2D::B(4, 12);
Matrix TimoFourNodeCurveSec2D::Btran(12, 4);
Matrix TimoFourNodeCurveSec2D::K(12, 12);
Vector TimoFourNodeCurveSec2D::P(12);
//double TimoFourNodeCurveSec2D::workArea[100];

void *
OPS_TimoFourNodeCurveSec2D(void)
{
	if (OPS_GetNumRemainingInputArgs() < 7) {
		opserr << "insufficient arguments: eleTag, iNode, jNode, kNode, lNode, numSec, secTag \n";
		return 0;
	}

	Element* theEle = 0;

	// get the ele tag
	int iData[7];
	int numData = 7;
	
	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "TimoFourNodeCurveSec2D::WARNING: invalid integer inputs \n";
		return 0;
	}

	int eleTag, iNode, jNode, kNode, lNode, numSec, secTag, numInt;

	eleTag = iData[0];
	iNode = iData[1];
	jNode = iData[2];
	kNode = iData[3];
	lNode = iData[4];
	numSec = iData[5];
	secTag = iData[6];

	SectionForceDeformation* sec = OPS_getSectionForceDeformation(secTag);
	if (sec == 0) {
		opserr << "section " << "not found\n";
		delete sec;
		return 0;
	}
	SectionForceDeformation** secs = new SectionForceDeformation* [numSec];
	for (int i = 0; i < numSec; i++)
		secs[i] = sec;

	BeamIntegration* beamIntegr = new LegendreBeamIntegration();

	theEle = new TimoFourNodeCurveSec2D(eleTag, iNode, jNode, kNode, lNode, numSec, secs, *beamIntegr);

	delete[] secs;

	return theEle;
}

TimoFourNodeCurveSec2D::TimoFourNodeCurveSec2D(int tag, int nd1, int nd2, int nd3, int nd4,
	int numSec, SectionForceDeformation** s,
	BeamIntegration& bi)
	:Element(tag, ELE_TAG_TimoFourNodeCurveSec2D),
	numSections(numSec), theSections(0), beamInt(0), L(0), costheta(0), sintheta(0),
	connectedExternalNodes(4), T(12, 12), Ttran(12, 12)
{
	// Allocate arrays of pointers to SectionForceDeformations
	theSections = new SectionForceDeformation * [numSections];

	if (theSections == 0) {
		opserr << "TimoFourNodeCurveSec2D::TimoFourNodeCurveSec2D - failed to allocate section model pointer\n";
		exit(-1);
	}

	for (int i = 0; i < numSections; i++) {

		// Get copies of the material model for each integration point
		theSections[i] = s[i]->getCopy();

		// Check allocation
		if (theSections[i] == 0) {
			opserr << "TimoFourNodeCurveSec2D::TimoFourNodeCurveSec2D -- failed to get a copy of section model\n";
			exit(-1);
		}
	}

	beamInt = bi.getCopy();

	if (beamInt == 0) {
		opserr << "TimoFourNodeCurveSec2D::TimoFourNodeCurveSec2D - failed to copy beam integration\n";
		exit(-1);
	}

	// Set connected external node IDs
	connectedExternalNodes(0) = nd1;
	connectedExternalNodes(1) = nd2;
	connectedExternalNodes(2) = nd3;
	connectedExternalNodes(3) = nd4;

	theNodes[0] = 0;
	theNodes[1] = 0;
	theNodes[2] = 0;
	theNodes[3] = 0;
}

TimoFourNodeCurveSec2D::TimoFourNodeCurveSec2D()
	:Element(0, ELE_TAG_TimoFourNodeCurveSec2D),
	numSections(0), theSections(0), beamInt(0), L(0), costheta(0), sintheta(0), 
	connectedExternalNodes(4), T(12, 12), Ttran(12, 12) 
{
	theNodes[0] = 0;
	theNodes[1] = 0;
	theNodes[2] = 0;
	theNodes[3] = 0;
}

TimoFourNodeCurveSec2D::~TimoFourNodeCurveSec2D()
{
	for (int i = 0; i < numSections; i++) {
		if (theSections[i])
			delete theSections[i];
	}

	// Delete the array of pointers to SectionForceDeformation pointer arrays
	if (theSections)
		delete[] theSections;

	if (beamInt != 0)
		delete beamInt;
}

int
TimoFourNodeCurveSec2D::getNumExternalNodes() const
{
	return 4;
}

const ID&
TimoFourNodeCurveSec2D::getExternalNodes()
{
	return connectedExternalNodes;
}

Node**
TimoFourNodeCurveSec2D::getNodePtrs()
{
	return theNodes;
}

int
TimoFourNodeCurveSec2D::getNumDOF()
{
	return 12;
}

void
TimoFourNodeCurveSec2D::setDomain(Domain* theDomain)
{
	// Check Domain is not null - invoked when object removed from a domain
	if (theDomain == 0) {
		theNodes[0] = 0;
		theNodes[1] = 0;
		theNodes[2] = 0;
		theNodes[3] = 0;
		return;
	}

	int Nd1 = connectedExternalNodes(0);
	int Nd2 = connectedExternalNodes(1);
	int Nd3 = connectedExternalNodes(2);
	int Nd4 = connectedExternalNodes(3);

	theNodes[0] = theDomain->getNode(Nd1);
	theNodes[1] = theDomain->getNode(Nd2);
	theNodes[2] = theDomain->getNode(Nd3);
	theNodes[3] = theDomain->getNode(Nd4);

	if (theNodes[0] == 0 || theNodes[1] == 0 || theNodes[2] == 0 || theNodes[3] == 0) {
		opserr << "WARNING TimoFourNodeCurveSec2D (tag: %d), node not found in domain" << this->getTag() << endln;;
		return;
	}

	int dofNd1 = theNodes[0]->getNumberDOF();
	int dofNd2 = theNodes[1]->getNumberDOF();
	int dofNd3 = theNodes[2]->getNumberDOF();
	int dofNd4 = theNodes[3]->getNumberDOF();

	if (dofNd1 != 3 || dofNd2 != 3 || dofNd3 != 3 || dofNd4 != 3) {
		opserr << "the node dof should be ndf1 = 3, ndf2 = 2, ndf3 = 3" << endln;
		return;
	}

	Vector direct = theNodes[3]->getCrds() - theNodes[0]->getCrds();
	//	double L = crdTransf->getInitialLength();

	L = direct.Norm();
	costheta = direct(0) / L;
	sintheta = direct(1) / L;

	T(0, 0) = costheta;
	T(1, 0) = -sintheta;
	T(0, 1) = sintheta;
	T(1, 1) = costheta;
	T(2, 2) = 1;
	T(3, 3) = costheta;
	T(4, 3) = -sintheta;
	T(3, 4) = sintheta;
	T(4, 4) = costheta;
	T(5, 5) = 1;
	T(6, 6) = costheta;
	T(7, 6) = -sintheta;
	T(6, 7) = sintheta;
	T(7, 7) = costheta;
	T(8, 8) = 1;
	T(9, 9) = costheta;
	T(10, 9) = -sintheta;
	T(9, 10) = sintheta;
	T(10, 10) = costheta;
	T(11, 11) = 1;


	Ttran.addMatrixTranspose(0, T, 1.0);

	if (L == 0.0) {
		opserr << "WARNING TimoFourNodeCurveSec2D (tag: %d), the length should not be 0" << this->getTag() << endln;;
		return;
	}

	this->DomainComponent::setDomain(theDomain);

	this->update();
}

int
TimoFourNodeCurveSec2D::commitState()
{
	int retVal = 0;

	// call element commitState to do any base class stuff
	if ((retVal = this->Element::commitState()) != 0) {
		opserr << "TimoFourNodeCurveSec2D::commitState () - failed in base class";
	}

	// Loop over the integration points and commit the material states
	for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->commitState();

	return retVal;
}

int
TimoFourNodeCurveSec2D::revertToLastCommit()
{
	int retVal = 0;

	// Loop over the integration points and revert to last committed state
	for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToLastCommit();

	return retVal;
}

int
TimoFourNodeCurveSec2D::revertToStart()
{
	int retVal = 0;

	// Loop over the integration points and revert states to start
	for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToStart();

	return retVal;
}

int
TimoFourNodeCurveSec2D::update(void)
{
	int err = 0;

	Vector displ(12);
	for (int i = 0; i < 4; i++)
	{
		const Vector& dispi = theNodes[i]->getTrialDisp();
		for (int j = 0; j < 3; j++)
			displ(3 * i + j) = dispi(j);
	}

	//const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
	double xi[maxNumSections];
	beamInt->getSectionLocations(numSections, L, xi);

	// Loop over the integration points
	for (int i = 0; i < numSections; i++) {

		calculateB(xi[i]);
	
		Vector e = B * T * displ;

		err += theSections[i]->setTrialSectionDeformation(e);
	}

	if (err != 0) {
		opserr << "TimoFourNodeCurveSec2D::update() - failed setTrialSectionDeformations()\n";
		return err;
	}

	return 0;
}

const Matrix&
TimoFourNodeCurveSec2D::getTangentStiff()
{
	double xi[maxNumSections];
	beamInt->getSectionLocations(numSections, L, xi);
	double wt[maxNumSections];
	beamInt->getSectionWeights(numSections, L, wt);

	// Loop over the integration points
	Matrix Kl(12, 12);
	for (int i = 0; i < numSections; i++) 
	{
		const Matrix& ks = theSections[i]->getSectionTangent();

		calculateB(xi[i]);

		Kl += Btran * ks * B * wt[i] * L;
	}

	K = Ttran * Kl * T;

	return K;
}

const Vector&
TimoFourNodeCurveSec2D::getResistingForce()
{
	double xi[maxNumSections];
	beamInt->getSectionLocations(numSections, L, xi);
	double wt[maxNumSections];
	beamInt->getSectionWeights(numSections, L, wt);

	// Loop over the integration points
	Vector Pl(12);
	for (int i = 0; i < numSections; i++) 
	{
		const Vector& s = theSections[i]->getStressResultant();

		calculateB(xi[i]);

		Pl += Btran * s  * wt[i] * L;
	}

	P = Ttran * Pl;

	return P;
}

void
TimoFourNodeCurveSec2D::calculateB(double xi)
{
	B.Zero();
	B(0, 0) = -(11 - 36 * xi + 27 * xi * xi) / (2 * L);
	B(0, 3) = 9 * (2 - 10 * xi + 9 * xi * xi) / (2 * L);
	B(0, 6) = -9 * (1 - 8 * xi + 9 * xi * xi) / (2 * L);
	B(0, 9) = (2 - 18 * xi + 27 * xi * xi) / (2 * L);

	B(1, 2) = -(11 - 36 * xi + 27 * xi * xi) / (2 * L);
	B(1, 5) = 9 * (2 - 10 * xi + 9 * xi * xi) / (2 * L);
	B(1, 8) = -9 * (1 - 8 * xi + 9 * xi * xi) / (2 * L);
	B(1, 11) = (2 - 18 * xi + 27 * xi * xi) / (2 * L);

	B(2, 2) = (11 - 36 * xi + 27 * xi * xi) / (2 * L);
	B(2, 5) = -9 * (2 - 10 * xi + 9 * xi * xi) / (2 * L);
	B(2, 8) = 9 * (1 - 8 * xi + 9 * xi * xi) / (2 * L);
	B(2, 11) = -(2 - 18 * xi + 27 * xi * xi) / (2 * L);
	B(2, 1) = 9 * (2 - 3 * xi) / (L * L);
	B(2, 4) = -9 * (5 - 9 * xi) / (L * L);
	B(2, 7) = 9 * (4 - 9 * xi) / (L * L);
	B(2, 10) = -9 * (1 - 3 * xi) / (L * L);

	B(3, 1) = -(11 - 36 * xi + 27 * xi * xi) / (2 * L);
	B(3, 4) = 9 * (2 - 10 * xi + 9 * xi * xi) / (2 * L);
	B(3, 7) = -9 * (1 - 8 * xi + 9 * xi * xi) / (2 * L);
	B(3, 10) = (2 - 18 * xi + 27 * xi * xi) / (2 * L);
	B(3, 2) = -(1 - 3 * xi) * (2 - 3 * xi) * (1 - xi) /2;
	B(3, 5) = -9 * (2 - 3 * xi) * (1 - xi) * xi / 2;
	B(3, 8) = 9 * (1 - 3 * xi) * (1 - xi) * xi / 2;
	B(3, 11) = -(1 - 3 * xi) * (2 - 3 * xi) * xi / 2;

	Btran.Zero();
	Btran.addMatrixTranspose(0.0, B, 1.0);
}

const Matrix&
TimoFourNodeCurveSec2D::getInitialStiff()
{
	return K;
}

const Matrix&
TimoFourNodeCurveSec2D::getMass()
{
	return K;
}

void
TimoFourNodeCurveSec2D::zeroLoad(void)
{

}

int
TimoFourNodeCurveSec2D::addLoad(ElementalLoad* theLoad, double loadFactor)
{
	return 0;
}

int
TimoFourNodeCurveSec2D::addInertiaLoadToUnbalance(const Vector& accel)
{
	return 0;
}

const Vector&
TimoFourNodeCurveSec2D::getResistingForceIncInertia()
{

	this->getResistingForce();

	return P;
}

int
TimoFourNodeCurveSec2D::sendSelf(int commitTag, Channel& theChannel)
{
	return 0;
}

int
TimoFourNodeCurveSec2D::recvSelf(int commitTag, Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
	return 0;
}

void
TimoFourNodeCurveSec2D::Print(OPS_Stream& s, int flag)
{
	s << "\nTimoFourNodeCurveSec2D, element id:  " << this->getTag() << endln;
	s << "\tConnected external nodes:  " << connectedExternalNodes;

	beamInt->Print(s, flag);

	for (int i = 0; i < numSections; i++)
		theSections[i]->Print(s, flag);
}


int
TimoFourNodeCurveSec2D::displaySelf(Renderer& theViewer, int displayMode, float fact)
{
	// first determine the end points of the quad based on
	// the display factor (a measure of the distorted image)
	const Vector& end1Crd = theNodes[0]->getCrds();
	const Vector& end2Crd = theNodes[1]->getCrds();

	static Vector v1(3);
	static Vector v2(3);

	if (displayMode >= 0) {
		const Vector& end1Disp = theNodes[0]->getDisp();
		const Vector& end2Disp = theNodes[1]->getDisp();

		for (int i = 0; i < 2; i++) {
			v1(i) = end1Crd(i) + end1Disp(i) * fact;
			v2(i) = end2Crd(i) + end2Disp(i) * fact;
		}
	}
	else {
		int mode = displayMode * -1;
		const Matrix& eigen1 = theNodes[0]->getEigenvectors();
		const Matrix& eigen2 = theNodes[1]->getEigenvectors();
		if (eigen1.noCols() >= mode) {
			for (int i = 0; i < 2; i++) {
				v1(i) = end1Crd(i) + eigen1(i, mode - 1) * fact;
				v2(i) = end2Crd(i) + eigen2(i, mode - 1) * fact;
			}
		}
		else {
			for (int i = 0; i < 2; i++) {
				v1(i) = end1Crd(i);
				v2(i) = end2Crd(i);
			}
		}
	}

	return theViewer.drawLine(v1, v2, 1.0, 1.0);
}

Response*
TimoFourNodeCurveSec2D::setResponse(const char** argv, int argc,
	OPS_Stream& output)
{
	Response* theResponse = 0;

	output.tag("ElementOutput");
	output.attr("eleType", "TimoFourNodeCurveSec2D");
	output.attr("eleTag", this->getTag());
	output.attr("node1", connectedExternalNodes[0]);
	output.attr("node2", connectedExternalNodes[1]);

	// global force - 
	if (strcmp(argv[0], "forces") == 0 || strcmp(argv[0], "force") == 0
		|| strcmp(argv[0], "globalForce") == 0 || strcmp(argv[0], "globalForces") == 0) {

		output.tag("ResponseType", "Px_1");
		output.tag("ResponseType", "Py_1");
		output.tag("ResponseType", "Mz_1");
		output.tag("ResponseType", "Px_2");
		output.tag("ResponseType", "Py_2");
		output.tag("ResponseType", "Mz_2");

		theResponse = new ElementResponse(this, 1, P);


		// local force -
	}
	else if (strcmp(argv[0], "localForce") == 0 || strcmp(argv[0], "localForces") == 0) {

		output.tag("ResponseType", "N1");
		output.tag("ResponseType", "V1");
		output.tag("ResponseType", "M1");
		output.tag("ResponseType", "N2");
		output.tag("ResponseType", "V2");
		output.tag("ResponseType", "M2");

		theResponse = new ElementResponse(this, 2, P);


		// basic force -
	}
	else if (strcmp(argv[0], "basicForce") == 0 || strcmp(argv[0], "basicForces") == 0) {

		output.tag("ResponseType", "N");
		output.tag("ResponseType", "M1");
		output.tag("ResponseType", "M2");

		theResponse = new ElementResponse(this, 9, Vector(3));

		// chord rotation -
	}
	else if (strcmp(argv[0], "chordRotation") == 0 || strcmp(argv[0], "chordDeformation") == 0
		|| strcmp(argv[0], "basicDeformation") == 0) {

		output.tag("ResponseType", "eps");
		output.tag("ResponseType", "theta1");
		output.tag("ResponseType", "theta2");

		theResponse = new ElementResponse(this, 3, Vector(3));

		// plastic rotation -
	}
	else if (strcmp(argv[0], "plasticRotation") == 0 || strcmp(argv[0], "plasticDeformation") == 0) {

		output.tag("ResponseType", "epsP");
		output.tag("ResponseType", "theta1P");
		output.tag("ResponseType", "theta2P");

		theResponse = new ElementResponse(this, 4, Vector(3));

	}
	else if (strcmp(argv[0], "RayleighForces") == 0 || strcmp(argv[0], "rayleighForces") == 0) {

		theResponse = new ElementResponse(this, 12, P);
	}

	// section response -
	else if (strstr(argv[0], "sectionX") != 0) {
		if (argc > 2) {
			float sectionLoc = atof(argv[1]);

			double xi[maxNumSections];
			beamInt->getSectionLocations(numSections, L, xi);

			sectionLoc /= L;

			float minDistance = fabs(xi[0] - sectionLoc);
			int sectionNum = 0;
			for (int i = 1; i < numSections; i++) {
				if (fabs(xi[i] - sectionLoc) < minDistance) {
					minDistance = fabs(xi[i] - sectionLoc);
					sectionNum = i;
				}
			}

			output.tag("GaussPointOutput");
			output.attr("number", sectionNum + 1);
			output.attr("eta", xi[sectionNum] * L);

			theResponse = theSections[sectionNum]->setResponse(&argv[2], argc - 2, output);
		}
	}
	else if (strstr(argv[0], "section") != 0) {

		if (argc > 1) {

			int sectionNum = atoi(argv[1]);

			if (sectionNum > 0 && sectionNum <= numSections && argc > 2) {

				output.tag("GaussPointOutput");
				output.attr("number", sectionNum);
				double xi[maxNumSections];
				beamInt->getSectionLocations(numSections, L, xi);
				output.attr("eta", xi[sectionNum - 1] * L);

				theResponse = theSections[sectionNum - 1]->setResponse(&argv[2], argc - 2, output);

				output.endTag();

			}
			else if (sectionNum == 0) { // argv[1] was not an int, we want all sections, 

				CompositeResponse* theCResponse = new CompositeResponse();
				int numResponse = 0;
				double xi[maxNumSections];
				beamInt->getSectionLocations(numSections, L, xi);

				for (int i = 0; i < numSections; i++) {

					output.tag("GaussPointOutput");
					output.attr("number", i + 1);
					output.attr("eta", xi[i] * L);

					Response* theSectionResponse = theSections[i]->setResponse(&argv[1], argc - 1, output);

					output.endTag();

					if (theSectionResponse != 0) {
						numResponse = theCResponse->addResponse(theSectionResponse);
					}
				}

				if (numResponse == 0) // no valid responses found
					delete theCResponse;
				else
					theResponse = theCResponse;
			}
		}
	}

	// curvature sensitivity along element length
	else if (strcmp(argv[0], "dcurvdh") == 0)
		return new ElementResponse(this, 5, Vector(numSections));

	// basic deformation sensitivity
	else if (strcmp(argv[0], "dvdh") == 0)
		return new ElementResponse(this, 6, Vector(3));

	else if (strcmp(argv[0], "integrationPoints") == 0)
		return new ElementResponse(this, 7, Vector(numSections));

	else if (strcmp(argv[0], "integrationWeights") == 0)
		return new ElementResponse(this, 8, Vector(numSections));

	output.endTag();
	return theResponse;
}

int
TimoFourNodeCurveSec2D::getResponse(int responseID, Information& eleInfo)
{
	double V;

	if (responseID == 1)
		return eleInfo.setVector(this->getResistingForce());

	else if (responseID == 12)
		return eleInfo.setVector(this->getRayleighDampingForces());

	else if (responseID == 2) {
		return eleInfo.setVector(P);
	}

	else if (responseID == 7) {
		//const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
		double xi[maxNumSections];
		beamInt->getSectionLocations(numSections, L, xi);
		Vector locs(numSections);
		for (int i = 0; i < numSections; i++)
			locs(i) = xi[i] * L;
		return eleInfo.setVector(locs);
	}

	else if (responseID == 8) {
		//const Vector &wts = quadRule.getIntegrPointWeights(numSections);
		double wt[maxNumSections];
		beamInt->getSectionWeights(numSections, L, wt);
		Vector weights(numSections);
		for (int i = 0; i < numSections; i++)
			weights(i) = wt[i] * L;
		return eleInfo.setVector(weights);
	}

	else
		return -1;
}
