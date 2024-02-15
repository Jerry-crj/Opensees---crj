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
// Description: This file contains the class definition for TimoshenkoModified.

#include "TimoshenkoModified.h"
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
#include <elementAPI.h>
#include <string>
#include <direct.h>

#define ELE_TAG_TimoshenkoModified 12314153

Matrix TimoshenkoModified::K(6, 6);
Vector TimoshenkoModified::P(6);

int theObjectiveTimeStep = 1500;    // step 1.0e-3;
int numEles = 2;
int timeStep = 0;

void *
OPS_TimoshenkoModified()
{
	if (OPS_GetNumRemainingInputArgs() < 4) {
		opserr << "insufficient arguments: eleTag, iNode, jNode, sectionTag, <-mass mass> <-cmass>\n";
		return 0;
	}

	// inputs: 
	int iData[4];
	int numData = 4;
	if (OPS_GetIntInput(&numData, &iData[0]) < 0) {
		opserr << "WARNING: invalid integer inputs\n";
		return 0;
	}

	// options
	double mass = 0.0;
	numData = 1;
	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* type = OPS_GetString();
		if (strcmp(type, "-mass") == 0) {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetDoubleInput(&numData, &mass) < 0) {
					opserr << "WARNING: invalid mass\n";
					return 0;
				}
			}
		}
		else
			break;
	}

	int secTag = iData[3];
	// check sections
	SectionForceDeformation *section = OPS_getSectionForceDeformation(secTag);
	if (section == 0) {
		opserr << "section " << "not found\n";
		delete section;
		return 0;
	}

	Element *theEle = new TimoshenkoModified(iData[0], iData[1], iData[2], section, mass);
	//delete[] section;
	return theEle;
}

TimoshenkoModified::TimoshenkoModified(int tag, int nd1, int nd2,
	SectionForceDeformation *s, double r)
	:Element(tag, ELE_TAG_TimoshenkoModified),
	theSection(0), connectedExternalNodes(2),
	Q(6), q(3), rho(r), e(3)
{
	// Get copies of the material model for each integration point
	theSection = s->getCopy();

	// Check allocation
	if (theSection == 0) {
		opserr << "TimoshenkoModified::TimoshenkoModified -- failed to get a copy of section model\n";
		exit(-1);
	}

	// Set connected external node IDs
	connectedExternalNodes(0) = nd1;
	connectedExternalNodes(1) = nd2;

	theNodes[0] = 0;
	theNodes[1] = 0;

	q0[0] = 0.0;
	q0[1] = 0.0;
	q0[2] = 0.0;

	p0[0] = 0.0;
	p0[1] = 0.0;
	p0[2] = 0.0;
}

TimoshenkoModified::TimoshenkoModified()
	:Element(0, ELE_TAG_TimoshenkoModified),
	numSections(0), theSection(0), connectedExternalNodes(2),
	Q(6), q(3), rho(0.0), e(3)
{
	q0[0] = 0.0;
	q0[1] = 0.0;
	q0[2] = 0.0;

	p0[0] = 0.0;
	p0[1] = 0.0;
	p0[2] = 0.0;

	theNodes[0] = 0;
	theNodes[1] = 0;
}

TimoshenkoModified::~TimoshenkoModified()
{
	if (theSection)
		delete[] theSection;

}

int
TimoshenkoModified::getNumExternalNodes() const
{
	return 2;
}

const ID&
TimoshenkoModified::getExternalNodes()
{
	return connectedExternalNodes;
}

Node **
TimoshenkoModified::getNodePtrs()
{
	return theNodes;
}

int
TimoshenkoModified::getNumDOF()
{
	return 6;
}

void
TimoshenkoModified::setDomain(Domain *theDomain)
{
	// Check Domain is not null - invoked when object removed from a domain
	if (theDomain == 0) {
		theNodes[0] = 0;
		theNodes[1] = 0;
		return;
	}

	int Nd1 = connectedExternalNodes(0);
	int Nd2 = connectedExternalNodes(1);

	theNodes[0] = theDomain->getNode(Nd1);
	theNodes[1] = theDomain->getNode(Nd2);

	if (theNodes[0] == 0 || theNodes[1] == 0) {
		opserr << "WARNING TimoshenkoModified (tag: %d), node not found in domain" << this->getTag() << endln;;
		return;
	}

	L = (theNodes[0]->getCrds() - theNodes[1]->getCrds()).Norm();

	if (L <1.e-12)
	{
		opserr << "the Length of TimosikoModifid Beam Column element should not be 0" << endln;
	}
	cosX = (theNodes[1]->getCrds() - theNodes[0]->getCrds())(0) / L;
	sinX = (theNodes[1]->getCrds() - theNodes[0]->getCrds())(1) / L;
	
	int dofNd1 = theNodes[0]->getNumberDOF();
	int dofNd2 = theNodes[1]->getNumberDOF();

	if (dofNd1 != 3 || dofNd2 != 3) {
		//opserr << "FATAL ERROR TimoshenkoModified (tag: %d), has differing number of DOFs at its nodes",
		//	this->getTag());

		return;
	}

	this->DomainComponent::setDomain(theDomain);

	this->update();
}

int
TimoshenkoModified::commitState()
{
	int retVal = 0;
	if ((retVal = this->Element::commitState()) != 0) {
		opserr << "TimoshenkoModified::commitState () - failed in base class";
	}
	retVal += theSection->commitState();

	return retVal;
}

int
TimoshenkoModified::revertToLastCommit()
{
	int retVal = 0;

	// Loop over the integration points and revert to last committed state
	retVal += theSection->revertToLastCommit();

	return retVal;
}

int
TimoshenkoModified::revertToStart()
{
	int retVal = 0;

	// Loop over the integration points and revert states to start
	retVal += theSection->revertToStart();

	return retVal;
}

int
TimoshenkoModified::update(void)
{
	int err = 0;

	Vector disp(6);
	Matrix B(3, 6), transf(6, 6);

	for (int i = 0; i < 3; i++)
	{
		disp(i) = theNodes[0]->getTrialDisp()[i];
		disp(i + 3) = theNodes[1]->getTrialDisp()[i];
	}
	// compute coeficient phi
	B = getB();
	transf = getTransf();
	e = B * transf * disp;

	//if (timeStep == theObjectiveTimeStep)
	//{
	//	static int trialNum1 = 0;
	//	static bool outputToTxt = true;
	//	if (trialNum1 >= 1)
	//	{
	//		int currentEle = (trialNum1 - 1) % callTimesEachStep + 1;
	//		if (currentEle == 2 && outputToTxt)
	//		{
	//			int error = mkdir(".\\nodeInfSurf");
	//			ofstream file0;
	//			file0.precision(12);

	//			double pertDisp = 5e-3;
	//			double pertAngle = 5e-5;
	//			int halfPertSteps = 20;
	//			for (int ii = 0; ii < 6; ii++)
	//			{
	//				Vector dispPert = disp;
	//				double pertValue;
	//				if (ii != 2 && ii != 5)
	//					pertValue = pertDisp;
	//				else
	//					pertValue = pertAngle;

	//				for (int k = 0; k < 2 * halfPertSteps; k++)
	//				{
	//					dispPert(ii) = disp(ii) - halfPertSteps * pertValue + k * pertValue;
	//					Vector epsPert = B * transf * dispPert;
	//					err += theSection->setTrialSectionDeformation(epsPert);
	//					Vector forcePert = this->getResistingForceIncInertia();
	//					for (int jj = 0; jj < 6; jj++)
	//					{
	//						std::string sss;
	//						sss = ".\\nodeInfSurf\\pertDisp" + std::to_string(ii + 1) + "force" + std::to_string(jj + 1) + ".txt";
	//						file0.open(sss, std::ios::app);
	//						file0 << dispPert(ii) << "    " << forcePert(jj) << endln;
	//						file0.close();
	//					}
	//				}
	//			}
	//			outputToTxt = false;
	//		}
	//	}
	//	trialNum1++;
	//}


	// Set the section deformations
	err += theSection->setTrialSectionDeformation(e);

	if (err != 0) {
		opserr << "TimoshenkoModified::update() - failed setTrialSectionDeformations()\n";
		return err;
	}

#ifdef DEBUGELEINFO
	if (timeStep == theObjectiveTimeStep)
	{
		static int trialNum = 0;

		if (trialNum >= 1)
		{
			int currentEle = (trialNum - 1) % callTimesEachStep + 1;

			ofstream file;
			file.precision(12);
			std::string s;
			s = ".\\failedStepInfo\\nodeDispEle" + std::to_string(currentEle) + ".txt";
			file.open(s, std::ios::app);
			file << disp(0) << "    " << disp(1) << "    " << disp(2) << "    " << disp(3) << "    " << disp(4) << "    " << disp(5) << endln;
			file.close();
		}
		trialNum++;
	}
#endif

	return 0;
}

const Matrix&
TimoshenkoModified::getTangentStiff()
{
	Matrix  Kl(6, 6), B(3, 6), BTran(6, 3), transf(6, 6), transfTran(6, 6);

	B = getB();

	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 3; j++)
			BTran(i, j) = B(j, i);

	transf = getTransf();

	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
			transfTran(i, j) = transf(j, i);

	// Get the section tangent stiffness and stress resultant
	const Matrix &ks = theSection->getSectionTangent();

	double wt = L;

	Kl = BTran * ks * B * wt;

	K = transfTran * Kl * transf;

#ifdef DEBUGELEINFO
	if (timeStep == theObjectiveTimeStep)
	{
		static int trialNum = 0;
		int currentEle = trialNum % callTimesEachStep + 1;

		ofstream file;
		file.precision(12);
		std::string s;
		s = ".\\failedStepInfo\\eleStiffEle" + std::to_string(currentEle) + ".txt";
		file.open(s, std::ios::app);
		file << K(0, 0) << "    " << K(0, 1) << "    " << K(0, 2) << "    " << K(0, 3) << "    " << K(0, 4) << "    " << K(0, 5) << endln;
		file << K(1, 0) << "    " << K(1, 1) << "    " << K(1, 2) << "    " << K(1, 3) << "    " << K(1, 4) << "    " << K(1, 5) << endln;
		file << K(2, 0) << "    " << K(2, 1) << "    " << K(2, 2) << "    " << K(2, 3) << "    " << K(2, 4) << "    " << K(2, 5) << endln;
		file << K(3, 0) << "    " << K(3, 1) << "    " << K(3, 2) << "    " << K(3, 3) << "    " << K(3, 4) << "    " << K(3, 5) << endln;
		file << K(4, 0) << "    " << K(4, 1) << "    " << K(4, 2) << "    " << K(4, 3) << "    " << K(4, 4) << "    " << K(4, 5) << endln;
		file << K(5, 0) << "    " << K(5, 1) << "    " << K(5, 2) << "    " << K(5, 3) << "    " << K(5, 4) << "    " << K(5, 5) << endln;
		file << endln;
		file.close();

		trialNum++;
	}
#endif

	return K;
}

const Matrix&
TimoshenkoModified::getInitialStiff()
{
	Matrix  Kl(6, 6), B(3, 6), BTran(6, 3), transf(6, 6), transfTran(6, 6);

	B = getB();

	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 3; j++)
			BTran(i, j) = B(j, i);

	transf = getTransf();

	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
			transfTran(i, j) = transf(j, i);

	// Get the section tangent stiffness and stress resultant
	const Matrix &ks = theSection->getInitialTangent();

	double wt = L;

	Kl = BTran * ks * B * wt;

	K = transfTran * Kl * transf;

	return K;
}

const Matrix&
TimoshenkoModified::getMass()
{
	K.Zero();

	if (rho == 0.0)
		return K;

	double m = 0.5 * rho * L;

	K(0, 0) = K(1, 1) = K(3, 3) = K(4, 4) = m;

	return K;
}

void
TimoshenkoModified::zeroLoad(void)
{
	Q.Zero();

	q0[0] = 0.0;
	q0[1] = 0.0;
	q0[2] = 0.0;

	p0[0] = 0.0;
	p0[1] = 0.0;
	p0[2] = 0.0;

	return;
}

int
TimoshenkoModified::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	int type;
	const Vector &data = theLoad->getData(type, loadFactor);

	if (type == LOAD_TAG_Beam2dUniformLoad) {
		double wt = data(0)*loadFactor;  // Transverse (+ve upward)
		double wa = data(1)*loadFactor;  // Axial (+ve from node I to J)

		double V = 0.5*wt*L;
		double M = V * L / 6.0; // wt*L*L/12
		double P = wa * L;

		// Reactions in basic system
		p0[0] -= P;
		p0[1] -= V;
		p0[2] -= V;

		// Fixed end forces in basic system
		q0[0] -= 0.5*P;
		q0[1] -= M;
		q0[2] += M;
	}
	else if (type == LOAD_TAG_Beam2dPointLoad) {
		double P = data(0)*loadFactor;
		double N = data(1)*loadFactor;
		double aOverL = data(2);

		if (aOverL < 0.0 || aOverL > 1.0)
			return 0;

		double a = aOverL * L;
		double b = L - a;

		// Reactions in basic system
		p0[0] -= N;
		double V1 = P * (1.0 - aOverL);
		double V2 = P * aOverL;
		p0[1] -= V1;
		p0[2] -= V2;

		double L2 = 1.0 / (L*L);
		double a2 = a * a;
		double b2 = b * b;

		// Fixed end forces in basic system
		q0[0] -= N * aOverL;
		double M1 = -a * b2 * P * L2;
		double M2 = a2 * b * P * L2;
		q0[1] += M1;
		q0[2] += M2;
	}
	else {
		opserr << "TimoshenkoModified::TimoshenkoModified -- load type unknown for element with tag: "
			<< this->getTag() << "TimoshenkoModified::addLoad()\n";

		return -1;
	}

	return 0;
}

int
TimoshenkoModified::addInertiaLoadToUnbalance(const Vector &accel)
{
	// Check for a quick return
	if (rho == 0.0)
		return 0;

	// Get R * accel from the nodes
	const Vector &Raccel1 = theNodes[0]->getRV(accel);
	const Vector &Raccel2 = theNodes[1]->getRV(accel);

	if (3 != Raccel1.Size() || 3 != Raccel2.Size()) {
		opserr << "TimoshenkoModified::addInertiaLoadToUnbalance matrix and vector sizes are incompatable\n";
		return -1;
	}

	double m = 0.5*rho*L;

	// Want to add ( - fact * M R * accel ) to unbalance
	// Take advantage of lumped mass matrix
	Q(0) -= m * Raccel1(0);
	Q(1) -= m * Raccel1(1);
	Q(3) -= m * Raccel2(0);
	Q(4) -= m * Raccel2(1);

	return 0;
}

const Vector&
TimoshenkoModified::getResistingForce()
{
	Matrix B(3, 6), BTran(6, 3), transf(6, 6), transfTran(6, 6);

	B = getB();

	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 3; j++)
			BTran(i, j) = B(j, i);

	transf = getTransf();

	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
			transfTran(i, j) = transf(j, i);

	const Vector &s = theSection->getStressResultant();
	
	double wt = L;

	Vector tmp(6);
	tmp = BTran * s * wt;
	P = transfTran * BTran * s * wt;

#ifdef DEBUGELEINFO
	if (timeStep == theObjectiveTimeStep)
	{
		static int trialNum = -5;
		if (trialNum > 0)
		{
			int currentEle = (trialNum-1) % callTimesEachStep + 1;

			ofstream file;
			file.precision(12);
			std::string s;
			s = ".\\failedStepInfo\\nodeForceEle" + std::to_string(currentEle) + ".txt";
			file.open(s, std::ios::app);
			file << P(0) << "    " << P(1) << "    " << P(2) << "    " << P(3) << "    " << P(4) << "    " << P(5) << endln;
			file.close();
		}
		trialNum++;
	}
#endif

	return P;
}

const Vector&
TimoshenkoModified::getResistingForceIncInertia()
{

	this->getResistingForce();

	if (rho != 0.0) {
		const Vector &accel1 = theNodes[0]->getTrialAccel();
		const Vector &accel2 = theNodes[1]->getTrialAccel();

		// Compute the current resisting force
		this->getResistingForce();

		double m = 0.5*rho*L;

		P(0) += m * accel1(0);
		P(1) += m * accel1(1);
		P(3) += m * accel2(0);
		P(4) += m * accel2(1);

		// add the damping forces if rayleigh damping
		if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
			P += this->getRayleighDampingForces();

	}
	else {

		// add the damping forces if rayleigh damping
		if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
			P += this->getRayleighDampingForces();
	}

	return P;
}

int
TimoshenkoModified::sendSelf(int commitTag, Channel &theChannel)
{
	return 0;
}

int
TimoshenkoModified::recvSelf(int commitTag, Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{
	return 0;
}

void
TimoshenkoModified::Print(OPS_Stream &s, int flag)
{
	s << "\nTimoshenkoModified, element id:  " << this->getTag() << endln;
	s << "\tConnected external nodes:  " << connectedExternalNodes;
	s << "\tmass density:  " << rho << endln;

	s << "\tEnd 1 Forces (N, V, M): (" << P(0) << ",   " << P(2) << ",   " << P(1) << ")" << endln;
	s << "\tEnd 2 Forces (N, V, M): (" << P(3) << ",   " << P(5) << ",   " << P(4) << ")" << endln;

	theSection->Print(s, flag);
}


int
TimoshenkoModified::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
	return 0;
}

Response*
TimoshenkoModified::setResponse(const char **argv, int argc,
	OPS_Stream &output)
{
	Response *theResponse = 0;
	if (strcmp(argv[0], "section") == 0) {
		if (argc > 1) {
			theResponse = theSection->setResponse(&argv[1], argc - 1, output);
		}
	}
	if (strcmp(argv[0], "strain") == 0) {
		theResponse = new ElementResponse(this, 1, e);
	}
	else if (strcmp(argv[0], "force") == 0)
	{
		Vector res(6);
		theResponse = new ElementResponse(this, 2, res);
	}

	return theResponse;
}

int
TimoshenkoModified::getResponse(int responseID, Information &eleInfo)
{
	switch (responseID) {
	case 1:
		return eleInfo.setVector(this->e);
	case 2:
		return eleInfo.setVector(this->getResistingForce());
	default:
		return -1;
	}
}


Matrix 
TimoshenkoModified::getB(void)
{
	Matrix temp(3, 6);
	temp(0, 0) = -1 / L; 
	temp(0, 3) = 1 / L;
	temp(1, 2) = 1 / L;
	temp(1, 5) = -1 / L;
	temp(2, 1) = -1 / L;
	temp(2, 2) = -0.5;
	temp(2, 4) = 1 / L;
	temp(2, 5) = -0.5;
	return temp;
}


Matrix
TimoshenkoModified::getTransf(void)
{
	Matrix temp(6, 6);
	temp(0, 0) = cosX;
	temp(0, 1) = sinX;
	temp(1, 0) = -sinX;
	temp(1, 1) = cosX;
	temp(2, 2) = 1;
	temp(3, 3) = cosX;
	temp(3, 4) = sinX;
	temp(4, 3) = -sinX;
	temp(4, 4) = cosX;
	temp(5, 5) = 1;
	return temp;
}

