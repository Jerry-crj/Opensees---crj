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

// $Revision: 1.10 $
// $Date: 2011/03/10 22:51:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/shell/PDTrussWithShear2DFinitDef.cpp,v $

// Written: Leopoldo Tesser, Diego A. Talledo, VÈronique Le Corvec
//
// Bathe MITC 4 four node shell element with membrane and drill
// Ref: Dvorkin,Bathe, A continuum mechanics based four node shell
//      element for general nonlinear analysis,
//      Eng.Comput.,1,77-88,1984

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <SectionForceDeformation.h>
#include <Domain.h>
#include <ErrorHandler.h>
#include "PDTrussWithShear2DFinitDef.h"
#include <R3vectors.h>
#include <Renderer.h>
#include <ElementResponse.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <map>

#include <string.h>
#include <fstream>
#include <TclModelBuilder.h>


#define DEBUG_PDTrussWithShear2DFinitDef
#define min(a,b) ( (a)<(b) ? (a):(b) )

//static int numPDTrussWithShear2DFinitDef = 0;

//static data
Matrix  PDTrussWithShear2DFinitDef::stiff(4, 4);
Vector  PDTrussWithShear2DFinitDef::resid(4);
Matrix  PDTrussWithShear2DFinitDef::mass(4, 4);
Matrix  PDTrussWithShear2DFinitDef::B(2, 4);
Matrix  PDTrussWithShear2DFinitDef::Btran(4, 2);

void *
OPS_PDTrussWithShear2DFinitDef(void)
{
	if (OPS_GetNumRemainingInputArgs() < 5) {
		opserr << "WARNING insufficient arguments\n";
		opserr << "want: eleTag? Node1? Node2? A? matTag? \n";
		return 0;
	}

	// get the id and end nodes 
	int iData[3];
	int num = 3;
	if (OPS_GetIntInput(&num, iData) < 0) {
		opserr << "WARNING: invalid integer data\n";
		return 0;
	}

	double area;
	num = 1;

	if (OPS_GetDouble(&num, &area) < 0) {
		opserr << "WARNING: invalid integer data\n";
		return 0;
	}

	int matTag;
	num = 1;
	if (OPS_GetIntInput(&num, &matTag) < 0) {
		opserr << "WARNING: invalid integer data\n";
		return 0;
	}

	NDMaterial* theMat = OPS_GetNDMaterial(matTag);
	if (theMat == 0) {
		opserr << "WARNING section not found\n";
		opserr << "Section: " << theMat;
		return 0;
	}

	Element *theElement = new PDTrussWithShear2DFinitDef(iData[0], iData[1], iData[2], area, theMat);

//	delete[] theSec;
	return theElement;
}


//null constructor
PDTrussWithShear2DFinitDef::PDTrussWithShear2DFinitDef() :
	Element(0, ELE_TAG_PDTrussWithShear2DFinitDef), theMaterial(0),
	connectedExternalNodes(2), A(0), L0(0), 
	Ki(0), load(0), n(2)
{
	nodePointers[0] = 0;
	nodePointers[1] = 0;

	connectedExternalNodes(0) = 0;
	connectedExternalNodes(1) = 0;
}

//*********************************************************************
//full constructor

PDTrussWithShear2DFinitDef::PDTrussWithShear2DFinitDef(int tag, int node1, int node2, double a, NDMaterial* theMat) :
	Element(tag, ELE_TAG_PDTrussWithShear2DFinitDef), A(a), connectedExternalNodes(2), load(0), Ki(0),
	L0(0), n(2)
{
	connectedExternalNodes(0) = node1;
	connectedExternalNodes(1) = node2;

	for (int i = 0; i < 2; i++)
	{
		nodePointers[i] = 0;
	}

	theMaterial = theMat->getCopy();
	if (theMaterial == 0) {
		opserr << "PDTrussWithShear2DFinitDef::PDTrussWithShear2DFinitDef - failed to allocate sectionN model pointer\n";
		exit(-1);
	}
}
//******************************************************************

//destructor 
PDTrussWithShear2DFinitDef::~PDTrussWithShear2DFinitDef()
{
	delete theMaterial;

	if (Ki != 0)
		delete Ki;
}
//**************************************************************************

//set domain
void  PDTrussWithShear2DFinitDef::setDomain(Domain *theDomain)
{
	for (int i = 0; i < 2; i++)
		nodePointers[i] = theDomain->getNode(connectedExternalNodes(i));

	const Vector &crdNode0 = nodePointers[0]->getCrds();
	const Vector &crdNode1 = nodePointers[1]->getCrds();

	int dofNd1 = nodePointers[0]->getNumberDOF();
	int dofNd2 = nodePointers[1]->getNumberDOF();

	if (dofNd1 != 2 || dofNd2 != 2)
	{
		opserr << "WARNING PDTrussWithShear2DFinitDef (tag: %d), the node dof should be ndf1 = 2, ndf2 = 2" << this->getTag() << endln;
		return;
	}

	n = crdNode1 - crdNode0;

	L0 = n.Norm();

	n /= L0;

	if (L0 == 0) {
		opserr << "WARNING PDTrussWithShear2DFinitDef (tag: %d), the length of element should not be 0" << this->getTag() << endln;
		return;
	}

	B(0, 0) = -1 / L0;
	B(1, 1) = -1 / L0;
	B(0, 2) = 1 / L0;
	B(1, 3) = 1 / L0;

	Btran.addMatrixTranspose(0.0, B, 1.0);

	this->DomainComponent::setDomain(theDomain);
}

//get the number of external nodes
int  PDTrussWithShear2DFinitDef::getNumExternalNodes() const
{
	return 2;
}

//return connected external nodes
const ID&  PDTrussWithShear2DFinitDef::getExternalNodes()
{
	return connectedExternalNodes;
}

Node **
PDTrussWithShear2DFinitDef::getNodePtrs(void)
{
	return nodePointers;
}

//return number of dofs
int PDTrussWithShear2DFinitDef::getNumDOF()
{
	return 4;
}

//commit state
int PDTrussWithShear2DFinitDef::commitState()
{
	int success = 0;

	// call element commitState to do any base class stuff
	if ((success = this->Element::commitState()) != 0) {
		opserr << "PDTrussWithShear2DFinitDef::commitState () - failed in base class";
	}

	success += theMaterial->commitState();

	n.Zero();
	n = nodePointers[1]->getCrds() - nodePointers[0]->getCrds();
	n += nodePointers[1]->getDisp() - nodePointers[0]->getDisp();
		
	L0 = n.Norm();

	n /= L0;
	
	calculateB();

	return success;
}

//revert to last commit 
int  PDTrussWithShear2DFinitDef::revertToLastCommit()
{
	int success = 0;

	success += theMaterial->revertToLastCommit();

	n.Zero();
	n = nodePointers[1]->getCrds() - nodePointers[0]->getCrds();
	n += nodePointers[1]->getDisp() - nodePointers[0]->getDisp();

	L0 = n.Norm();

	n /= L0;

	calculateB();

	return success;
}


//revert to start 
int PDTrussWithShear2DFinitDef::revertToStart()
{
	int success = 0;

	success += theMaterial->revertToStart();

	return success;
}

//print out element data
void PDTrussWithShear2DFinitDef::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_CURRENTSTATE) {
		s << endln;
		s << "Element type: Shear Wall Shell, tag: " << this->getTag() << endln;
		s << "\tConnected external nodes: " << connectedExternalNodes;
		s << "\tNumber of sections : 1" << endln;
		s << "\nSection Information: " << endln;

		theMaterial->Print(s, OPS_PRINT_PRINTMODEL_MATERIAL);

		s << endln;
	}
}

Response*
PDTrussWithShear2DFinitDef::setResponse(const char **argv, int argc, OPS_Stream &output)
{
	Response *theResponse = 0;

	if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "strain") == 0)
	{
		theResponse = theMaterial->setResponse(argv, argc, output);
	}
	else
		opserr<<"this recorder type is not supported in PDTrussWithShear2DFinitDef "<<endln;

	return theResponse;
}

int
PDTrussWithShear2DFinitDef::getResponse(int responseID, Information &eleInfo)
{
	switch (responseID) {
	case 1: // global forces
		return theMaterial->getResponse(responseID, eleInfo);
	default:
		return -1;
	}
}


//get residual
const Vector&  PDTrussWithShear2DFinitDef::getResistingForce()
{
	Vector stress = theMaterial->getStress();

	resid.Zero();

	resid = Ttran * Btran * stress * L0;

	// subtract external loads 
	if (load != 0)
		resid -= *load;

	return resid;
}


//return stiffness matrix 
const Matrix& PDTrussWithShear2DFinitDef::getTangentStiff()
{

	Matrix tangent = theMaterial->getTangent();

	stiff.Zero();
	
	stiff = Ttran * Btran * tangent * B * T * L0;

	return stiff;
}


//return mass matrix
const Matrix&  PDTrussWithShear2DFinitDef::getMass()
{
	return mass;
}


//return secant matrix 
const Matrix& PDTrussWithShear2DFinitDef::getInitialStiff()
{
	if (Ki != 0)
		return *Ki;

	Matrix tangent = theMaterial->getInitialTangent();

	Matrix Ttran(2, 2);
	Ttran.addMatrixTranspose(0.0, T, 1.0);
	stiff = (Ttran * tangent * T) * A;

	Ki = new Matrix(stiff);

	return stiff;
}


void  PDTrussWithShear2DFinitDef::zeroLoad()
{
	if (load != 0)
		load->Zero();

	return;
}


int
PDTrussWithShear2DFinitDef::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	opserr << "PDTrussWithShear2DFinitDef::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
	return -1;
}



int
PDTrussWithShear2DFinitDef::addInertiaLoadToUnbalance(const Vector &accel)
{

	static Vector r(24);

	int i;

	int allRhoZero = 0;
	if (theMaterial->getRho() != 0.0)
		allRhoZero = 1;

	if (allRhoZero == 0)
		return 0;

	//formInertiaLoad();

	int count = 0;
	for (i = 0; i < 2; i++) {
		const Vector &Raccel = nodePointers[i]->getRV(accel);
		for (int j = 0; j < 6; j++)
			r(count++) = Raccel(j);
	}

	if (load == 0)
		load = new Vector(24);

	load->addMatrixVector(1.0, mass, r, -1.0);

	return 0;
}


//get residual with inertia terms
const Vector&  PDTrussWithShear2DFinitDef::getResistingForceIncInertia()
{
	static Vector res(2);

	//formResistingForce();
	//formInertiaLoad();

	res = resid;
	// add the damping forces if rayleigh damping
	if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
		res += this->getRayleighDampingForces();

	// subtract external loads 
	if (load != 0)
		res -= *load;

	return res;
}

// º∆À„”¶±‰
int
PDTrussWithShear2DFinitDef::update(void)
{
	int success = 0;

	Vector detug(2);
	for (int i = 0; i < 2; i++)
	{
		Vector nug = nodePointers[i]->getTrialDisp();
		detug(i) = nug(1) - nug(0);
	}

	double Ln;

	for (int i = 0; i < 2; i++)
		Ln += detug(i) * detug(i);

	Ln = pow(Ln, 0.5);

	Vector e(2);
	e(0) = (Ln - L0) / L0;
	e(1) = detug - (Ln - L0) * n;

	success += theMaterial->setTrialStrain(e);

	return 0;
}


//**********************************************************************

int  PDTrussWithShear2DFinitDef::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	return res;
}

int  PDTrussWithShear2DFinitDef::recvSelf(int commitTag,
	Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{
	int res = 0;

	return res;
}
//**************************************************************************

int
PDTrussWithShear2DFinitDef::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{

	int error = 0;

	return error;
}


void 
PDTrussWithShear2DFinitDef::calculateB()
{
	B(0, 0) = -1 / L0;
	B(1, 1) = -1 / L0;
	B(0, 2) = 1 / L0;
	B(1, 3) = 1 / L0;

	Btran.addMatrixTranspose(0.0, B, 1.0);
}
