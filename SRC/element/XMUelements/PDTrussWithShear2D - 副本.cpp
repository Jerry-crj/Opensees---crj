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

// Written: Leopoldo Tesser, Diego A. Talledo, Véronique Le Corvec
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
		opserr << "want: eleTag? Node1? Node2? Node3? Node4? secTag?\n";
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
	connectedExternalNodes(2), initdVec(0), initnVec(2), A(0), L(0),
	initL(0), Ki(0), load(0), T(2, 2), Ttran(2, 2)
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
	initdVec(2), initnVec(2), L(0), initL(0), T(2, 2), Ttran(2, 2)
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

	initdVec = crdNode1 - crdNode0;
	L = initdVec.Norm();
	initL = L;

	if (L < 1e-13) {
		opserr << "PDTrussWithShear2D::PDTrussWithShear2D - failed to get height or length\n";
		exit(-1);
	}

	initdVec = initdVec / L;
	initnVec(0) = -initdVec(1);
	initnVec(1) = initdVec(0);

	Ttran(0, 0) = initdVec(0);
	Ttran(1, 0) = initdVec(1);
	Ttran(0, 1) = -initdVec(1);
	Ttran(1, 1) = initdVec(0);

	T.addMatrixTranspose(0.0, Ttran, 1.0);


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

	if (strcmp(argv[0], "section") == 0)
	{
		if (argc > 2)
		{
			theResponse = theMaterial->setResponse(&argv[1], argc - 1, output);
		}
		else
		{
			opserr << "this recorder type is not supported in SW-element " << endln;
		}
	}
	else
	{
		opserr<<"this recorder type is not supported in SW-element "<<endln;
	}

	return theResponse;
}

int
PDTrussWithShear2D::getResponse(int responseID, Information &eleInfo)
{
	int cnt = 0;
	static Vector stresses(32);
	static Vector strains(32);

	switch (responseID) {
	case 1: // global forces
		return eleInfo.setVector(this->getResistingForce());
		break;

	case 2: // stresses
		for (int i = 0; i < 4; i++) {

			// Get material stress response
			const Vector &sigma = theMaterial->getStress();
			stresses(cnt) = sigma(0);
			stresses(cnt + 1) = sigma(1);
			stresses(cnt + 2) = sigma(2);
			stresses(cnt + 3) = sigma(3);
			stresses(cnt + 4) = sigma(4);
			stresses(cnt + 5) = sigma(5);
			stresses(cnt + 6) = sigma(6);
			stresses(cnt + 7) = sigma(7);
			cnt += 8;
		}
		return eleInfo.setVector(stresses);
		break;
	case 3: //strain
		for (int i = 0; i < 4; i++) {

			// Get section deformation
			const Vector &deformation = theMaterial->getStrain();
			strains(cnt) = deformation(0);
			strains(cnt + 1) = deformation(1);
			strains(cnt + 2) = deformation(2);
			strains(cnt + 3) = deformation(3);
			strains(cnt + 4) = deformation(4);
			strains(cnt + 5) = deformation(5);
			strains(cnt + 6) = deformation(6);
			strains(cnt + 7) = deformation(7);
			cnt += 8;
		}
		return eleInfo.setVector(strains);
		break;
	default:
		return -1;
	}
}


//get residual
const Vector&  PDTrussWithShear2D::getResistingForce()
{
	Vector crd1, crd2;
	crd1 = nodePointers[0]->getCrds();
	crd2 = nodePointers[1]->getCrds();

	Vector stress = theMaterial->getStress();

	Vector fb(2);
	fb = (Ttran * stress) * A;

	resid.Zero();
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
		{
			if (i == 0)
				resid(2 * i + j) = -fb(j);
			else
				resid(2 * i + j) = fb(j);
		}
	//{
	//	if (i == 0)
	//	{
	//		resid(2 * i) = -fb(0);
	//		resid(2 * i + 1) = fb(1);
	//	}
	//	else
	//	{
	//		resid(2 * i) = fb(0);
	//		resid(2 * i + 1) = -fb(1);
	//	}
	//}

	// subtract external loads 
	if (load != 0)
		resid -= *load;

	return resid;
}


//return stiffness matrix 
const Matrix& PDTrussWithShear2D::getTangentStiff()
{

	Matrix tangent = theMaterial->getTangent();

	Matrix kb(2, 2);
	kb = (Ttran * tangent * T) * A / initL;

	stiff.Zero();
	for (int i=0;i<2;i++)
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 2; k++)
				for (int l = 0; l < 2; l++)
				{
					if ((i == 0 && j == 0) || (i == 1 && j == 1))
						stiff(2 * i + k, 2 * j + l) = kb(k, l);
					else
						stiff(2 * i + k, 2 * j + l) = -kb(k, l);
				}
		}
	
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

// ¼ÆËãÓ¦±ä
int
PDTrussWithShear2D::update(void)
{
	int success = 0;

	Vector node1Disp, node2Disp, deltaDisp;
	node1Disp = nodePointers[0]->getTrialDisp();
	node2Disp = nodePointers[1]->getTrialDisp();

	deltaDisp = node2Disp - node1Disp;

	Vector strain(2);

	strain(0) = (deltaDisp ^ initdVec) / L;
	strain(1) = (deltaDisp / L - strain(0) * initdVec) ^ initnVec;

	Vector crd1, crd2;
	crd1 = nodePointers[0]->getCrds();
	crd2 = nodePointers[1]->getCrds();


	success += theMaterial->setTrialStrain(strain);

	return 0;
}


//**********************************************************************

int  PDTrussWithShear2D::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	//// note: we don't check for dataTag == 0 for Element
	//// objects as that is taken care of in a commit by the Domain
	//// object - don't want to have to do the check if sending data
	//int dataTag = this->getDbTag();


	//// Now quad sends the ids of its materials
	//int matDbTag;

	//static ID idData(14);

	//int i;
	//for (i = 0; i < 4; i++) {
	//	idData(i) = theSections->getClassTag();
	//	matDbTag = theSections->getDbTag();
	//	// NOTE: we do have to ensure that the material has a database
	//	// tag if we are sending to a database channel.
	//	if (matDbTag == 0) {
	//		matDbTag = theChannel.getDbTag();
	//		if (matDbTag != 0)
	//			theSections->setDbTag(matDbTag);
	//	}
	//	idData(i + 4) = matDbTag;
	//}

	//idData(8) = this->getTag();
	//idData(9) = connectedExternalNodes(0);
	//idData(10) = connectedExternalNodes(1);
	//idData(11) = connectedExternalNodes(2);
	//idData(12) = connectedExternalNodes(3);
	//if (doUpdateBasis == true)
	//	idData(13) = 0;
	//else
	//	idData(13) = 1;

	//res += theChannel.sendID(dataTag, commitTag, idData);
	//if (res < 0) {
	//	opserr << "WARNING PDTrussWithShear2D::sendSelf() - " << this->getTag() << " failed to send ID\n";
	//	return res;
	//}

	//static Vector vectData(4);
	//vectData(0) = alphaM;
	//vectData(1) = betaK;
	//vectData(2) = betaK0;
	//vectData(3) = betaKc;

	//res += theChannel.sendVector(dataTag, commitTag, vectData);
	//if (res < 0) {
	//	opserr << "WARNING PDTrussWithShear2D::sendSelf() - " << this->getTag() << " failed to send ID\n";
	//	return res;
	//}

	//// Finally, quad asks its material objects to send themselves
	//for (i = 0; i < 4; i++) {
	//	res += theSectionsN->sendSelf(commitTag, theChannel);
	//	if (res < 0) {
	//		opserr << "WARNING PDTrussWithShear2D::sendSelf() - " << this->getTag() << " failed to send its Material\n";
	//		return res;
	//	}
	//}

	return res;
}

int  PDTrussWithShear2D::recvSelf(int commitTag,
	Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{
	int res = 0;

	//int dataTag = this->getDbTag();

	//static ID idData(14);
	//// Quad now receives the tags of its four external nodes
	//res += theChannel.recvID(dataTag, commitTag, idData);
	//if (res < 0) {
	//	opserr << "WARNING PDTrussWithShear2D::recvSelf() - " << this->getTag() << " failed to receive ID\n";
	//	return res;
	//}

	//this->setTag(idData(8));
	//connectedExternalNodes(0) = idData(9);
	//connectedExternalNodes(1) = idData(10);
	//connectedExternalNodes(2) = idData(11);
	//connectedExternalNodes(3) = idData(12);
	//if (idData(13) == 0)
	//	doUpdateBasis = true;
	//else
	//	doUpdateBasis = false;

	//static Vector vectData(4);
	//res += theChannel.recvVector(dataTag, commitTag, vectData);
	//if (res < 0) {
	//	opserr << "WARNING PDTrussWithShear2D::sendSelf() - " << this->getTag() << " failed to send ID\n";
	//	return res;
	//}

	//alphaM = vectData(0);
	//betaK = vectData(1);
	//betaK0 = vectData(2);
	//betaKc = vectData(3);

	//int i;

	//if (theSectionsN == 0) {
	//	for (i = 0; i < 4; i++) {
	//		int matClassTag = idData(i);
	//		int matDbTag = idData(i + 4);
	//		// Allocate new material with the sent class tag
	//		theSectionsN = theBroker.getNewSection(matClassTag);
	//		if (theSectionsN == 0) {
	//			opserr << "PDTrussWithShear2D::recvSelf() - Broker could not create NDMaterial of class type" << matClassTag << endln;;
	//			return -1;
	//		}
	//		// Now receive materials into the newly allocated space
	//		theSectionsN->setDbTag(matDbTag);
	//		res += theSectionsN->recvSelf(commitTag, theChannel, theBroker);
	//		if (res < 0) {
	//			opserr << "PDTrussWithShear2D::recvSelf() - material " << i << "failed to recv itself\n";
	//			return res;
	//		}
	//	}
	//}
	//// Number of materials is the same, receive materials into current space
	//else {
	//	for (i = 0; i < 4; i++) {
	//		int matClassTag = idData(i);
	//		int matDbTag = idData(i + 4);
	//		// Check that material is of the right type; if not,
	//		// delete it and create a new one of the right type
	//		if (theSectionsN->getClassTag() != matClassTag) {
	//			delete theSectionsN;
	//			theSectionsN = theBroker.getNewSection(matClassTag);
	//			if (theSectionsN == 0) {
	//				opserr << "PDTrussWithShear2D::recvSelf() - Broker could not create NDMaterial of class type" << matClassTag << endln;
	//				exit(-1);
	//			}
	//		}
	//		// Receive the material
	//		theSectionsN->setDbTag(matDbTag);
	//		res += theSectionsN->recvSelf(commitTag, theChannel, theBroker);
	//		if (res < 0) {
	//			opserr << "PDTrussWithShear2D::recvSelf() - material " << i << "failed to recv itself\n";
	//			return res;
	//		}
	//	}
	//}

	return res;
}
//**************************************************************************

int
PDTrussWithShear2D::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
	//// first determine the end points of the quad based on
	//// the display factor (a measure of the distorted image)
	//// store this information in 4 3d vectors v1 through v4
	//const Vector &end1Crd = nodePointers[0]->getCrds();
	//const Vector &end2Crd = nodePointers[1]->getCrds();
	//const Vector &end3Crd = nodePointers[2]->getCrds();
	//const Vector &end4Crd = nodePointers[3]->getCrds();

	//static Matrix coords(4, 3);
	//static Vector values(4);
	//static Vector P(24);

	//for (int j = 0; j < 4; j++)
	//	values(j) = 0.0;

	//if (displayMode >= 0) {
	//	// Display mode is positive:
	//	// display mode = 0 -> plot no contour
	//	// display mode = 1-8 -> plot 1-8 stress resultant

	//	// Get nodal displacements
	//	const Vector &end1Disp = nodePointers[0]->getDisp();
	//	const Vector &end2Disp = nodePointers[1]->getDisp();
	//	const Vector &end3Disp = nodePointers[2]->getDisp();
	//	const Vector &end4Disp = nodePointers[3]->getDisp();

	//	// Get stress resultants
	//	if (displayMode <= 8 && displayMode > 0) {
	//		for (int i = 0; i < 4; i++) {
	//			const Vector &stress = theSections->getStressResultant();
	//			values(i) = stress(displayMode - 1);
	//		}
	//	}

	//	// Get nodal absolute position = OriginalPosition + (Displacement*factor)
	//	for (int i = 0; i < 3; i++) {
	//		coords(0, i) = end1Crd(i) + end1Disp(i)*fact;
	//		coords(1, i) = end2Crd(i) + end2Disp(i)*fact;
	//		coords(2, i) = end3Crd(i) + end3Disp(i)*fact;
	//		coords(3, i) = end4Crd(i) + end4Disp(i)*fact;
	//	}
	//}
	//else {
	//	// Display mode is negative.
	//	// Plot eigenvectors
	//	int mode = displayMode * -1;
	//	const Matrix &eigen1 = nodePointers[0]->getEigenvectors();
	//	const Matrix &eigen2 = nodePointers[1]->getEigenvectors();
	//	const Matrix &eigen3 = nodePointers[2]->getEigenvectors();
	//	const Matrix &eigen4 = nodePointers[3]->getEigenvectors();
	//	if (eigen1.noCols() >= mode) {
	//		for (int i = 0; i < 3; i++) {
	//			coords(0, i) = end1Crd(i) + eigen1(i, mode - 1)*fact;
	//			coords(1, i) = end2Crd(i) + eigen2(i, mode - 1)*fact;
	//			coords(2, i) = end3Crd(i) + eigen3(i, mode - 1)*fact;
	//			coords(3, i) = end4Crd(i) + eigen4(i, mode - 1)*fact;
	//		}
	//	}
	//	else {
	//		for (int i = 0; i < 3; i++) {
	//			coords(0, i) = end1Crd(i);
	//			coords(1, i) = end2Crd(i);
	//			coords(2, i) = end3Crd(i);
	//			coords(3, i) = end4Crd(i);
	//		}
	//	}
	//}

	int error = 0;

	//// Draw a poligon with coordinates coords and values (colors) corresponding to values vector
	//error += theViewer.drawPolygon(coords, values);

	return error;
}

