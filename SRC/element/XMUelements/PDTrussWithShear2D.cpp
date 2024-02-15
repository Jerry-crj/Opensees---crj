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
// $Source: /usr/local/cvs/OpenSees/SRC/element/shell/PDTrussWithShear2D.cpp,v $

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
#include "PDTrussWithShear2D.h"
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


#define DEBUG_PDTrussWithShear2D
#define min(a,b) ( (a)<(b) ? (a):(b) )

//static int numPDTrussWithShear2D = 0;

//static data
Matrix  PDTrussWithShear2D::stiff(4, 4);
Vector  PDTrussWithShear2D::resid(4);
Matrix  PDTrussWithShear2D::mass(4, 4);

void *
OPS_PDTrussWithShear2D(void)
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

	Element *theElement = new PDTrussWithShear2D(iData[0], iData[1], iData[2], area, theMat);

//	delete[] theSec;
	return theElement;
}


//null constructor
PDTrussWithShear2D::PDTrussWithShear2D() :
	Element(0, ELE_TAG_PDTrussWithShear2D), theMaterial(0),
	connectedExternalNodes(2), parallel(0), perpendicular(2), A(0), L(0),
	Ki(0), load(0), T(4, 4), Ttran(4, 4), B(2, 4), Btran(4, 2)
{
	nodePointers[0] = 0;
	nodePointers[1] = 0;

	connectedExternalNodes(0) = 0;
	connectedExternalNodes(1) = 0;
}

//*********************************************************************
//full constructor

PDTrussWithShear2D::PDTrussWithShear2D(int tag, int node1, int node2, double a, NDMaterial* theMat) :
	Element(tag, ELE_TAG_PDTrussWithShear2D), A(a), connectedExternalNodes(2), load(0), Ki(0),
	parallel(2), perpendicular(2), L(0), T(4, 4), Ttran(4, 4), B(2, 4), Btran(4, 2)
{
	connectedExternalNodes(0) = node1;
	connectedExternalNodes(1) = node2;

	for (int i = 0; i < 2; i++)
	{
		nodePointers[i] = 0;
	}

	theMaterial = theMat->getCopy();
	if (theMaterial == 0) {
		opserr << "PDTrussWithShear2D::PDTrussWithShear2D - failed to allocate sectionN model pointer\n";
		exit(-1);
	}
}
//******************************************************************

//destructor 
PDTrussWithShear2D::~PDTrussWithShear2D()
{
	delete theMaterial;

	if (Ki != 0)
		delete Ki;
}
//**************************************************************************

//set domain
void  PDTrussWithShear2D::setDomain(Domain *theDomain)
{
	for (int i = 0; i < 2; i++)
		nodePointers[i] = theDomain->getNode(connectedExternalNodes(i));

	const Vector &crdNode0 = nodePointers[0]->getCrds();
	const Vector &crdNode1 = nodePointers[1]->getCrds();

	int dofNd1 = nodePointers[0]->getNumberDOF();
	int dofNd2 = nodePointers[1]->getNumberDOF();

	if (dofNd1 != 2 || dofNd2 != 2)
	{
		opserr << "WARNING PDTrussWithShear2D (tag: %d), the node dof should be ndf1 = 2, ndf2 = 2" << this->getTag() << endln;
		return;
	}

	parallel = crdNode1 - crdNode0;
	L = parallel.Norm();
	
	if (L == 0) {
		opserr << "WARNING PDTrussWithShear2D (tag: %d), the length of element should not be 0" << this->getTag() << endln;
		return;
	}

	parallel = parallel / L;
	perpendicular(0) = -parallel(1);
	perpendicular(1) = parallel(0);

	T(0, 0) = parallel(0);
	T(1, 0) = -parallel(1);
	T(0, 1) = parallel(1);
	T(1, 1) = parallel(0);

	T(2, 2) = parallel(0);
	T(3, 2) = -parallel(1);
	T(2, 3) = parallel(1);
	T(3, 3) = parallel(0);

	Ttran.addMatrixTranspose(0.0, T, 1.0);

	B(0, 0) = -1 / L;
	B(1, 1) = -1 / L;
	B(0, 2) = 1 / L;
	B(1, 3) = 1 / L;

	Btran.addMatrixTranspose(0.0, B, 1.0);

	this->DomainComponent::setDomain(theDomain);
}

//get the number of external nodes
int  PDTrussWithShear2D::getNumExternalNodes() const
{
	return 2;
}

//return connected external nodes
const ID&  PDTrussWithShear2D::getExternalNodes()
{
	return connectedExternalNodes;
}

Node **
PDTrussWithShear2D::getNodePtrs(void)
{
	return nodePointers;
}

//return number of dofs
int PDTrussWithShear2D::getNumDOF()
{
	return 4;
}

//commit state
int PDTrussWithShear2D::commitState()
{
	int success = 0;

	// call element commitState to do any base class stuff
	if ((success = this->Element::commitState()) != 0) {
		opserr << "PDTrussWithShear2D::commitState () - failed in base class";
	}

	success += theMaterial->commitState();

	return success;
}

//revert to last commit 
int  PDTrussWithShear2D::revertToLastCommit()
{
	int success = 0;

	success += theMaterial->revertToLastCommit();

	return success;
}


//revert to start 
int PDTrussWithShear2D::revertToStart()
{
	int success = 0;

	success += theMaterial->revertToStart();

	return success;
}

//print out element data
void PDTrussWithShear2D::Print(OPS_Stream &s, int flag)
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
PDTrussWithShear2D::setResponse(const char **argv, int argc, OPS_Stream &output)
{
	Response *theResponse = 0;

	if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "strain") == 0)
	{
		theResponse = theMaterial->setResponse(argv, argc, output);
	}
	else
		opserr<<"this recorder type is not supported in PDTrussWithShear2D "<<endln;

	return theResponse;
}

int
PDTrussWithShear2D::getResponse(int responseID, Information &eleInfo)
{
	switch (responseID) {
	case 1: // global forces
		return theMaterial->getResponse(responseID, eleInfo);
	default:
		return -1;
	}
}


//get residual
const Vector&  PDTrussWithShear2D::getResistingForce()
{
	Vector stress = theMaterial->getStress();

	resid.Zero();

	resid = Ttran * Btran * stress * L;

	// subtract external loads 
	if (load != 0)
		resid -= *load;

	return resid;
}


//return stiffness matrix 
const Matrix& PDTrussWithShear2D::getTangentStiff()
{

	Matrix tangent = theMaterial->getTangent();

	stiff.Zero();
	
	stiff = Ttran * Btran * tangent * B * T * L;

	return stiff;
}


//return mass matrix
const Matrix&  PDTrussWithShear2D::getMass()
{
	return mass;
}


//return secant matrix 
const Matrix& PDTrussWithShear2D::getInitialStiff()
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


void  PDTrussWithShear2D::zeroLoad()
{
	if (load != 0)
		load->Zero();

	return;
}


int
PDTrussWithShear2D::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	opserr << "PDTrussWithShear2D::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
	return -1;
}



int
PDTrussWithShear2D::addInertiaLoadToUnbalance(const Vector &accel)
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
const Vector&  PDTrussWithShear2D::getResistingForceIncInertia()
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
PDTrussWithShear2D::update(void)
{
	int success = 0;

	Vector displ(4);
	for (int i = 0; i < 2; i++)
	{
		Vector nodeDisp = nodePointers[i]->getTrialDisp();
		for (int j = 0; j < 2; j++)
			displ(2 * i + j) = nodeDisp(j);
	}

	Vector e = B * T * displ;

	success += theMaterial->setTrialStrain(e);

	return 0;
}


//**********************************************************************

int  PDTrussWithShear2D::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	return res;
}

int  PDTrussWithShear2D::recvSelf(int commitTag,
	Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{
	int res = 0;

	return res;
}
//**************************************************************************

int
PDTrussWithShear2D::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{

	int error = 0;

	return error;
}

