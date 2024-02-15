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
// $Source: /usr/local/cvs/OpenSees/SRC/element/shell/PDTrussWithShear3D.cpp,v $

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
#include "PDTrussWithShear3D.h"
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


#define DEBUG_PDTrussWithShear3D
#define min(a,b) ( (a)<(b) ? (a):(b) )

//static int numPDTrussWithShear3D = 0;

//static data
Matrix  PDTrussWithShear3D::stiff(6, 6);
Vector  PDTrussWithShear3D::resid(6);
Matrix  PDTrussWithShear3D::mass(6, 6);

void*
OPS_PDTrussWithShear3D(void)
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

	Element* theElement = new PDTrussWithShear3D(iData[0], iData[1], iData[2], area, theMat);

	//	delete[] theSec;
	return theElement;
}


//null constructor
PDTrussWithShear3D::PDTrussWithShear3D() :
	Element(0, ELE_TAG_PDTrussWithShear3D), theMaterial(0),
	connectedExternalNodes(2), parallel(3), perpendicular1(3), perpendicular2(3), A(0), L(0),
	Ki(0), load(0), T(3, 3), Ttran(3, 3), B(3, 6), dirct(2)
{
	nodePointers[0] = 0;
	nodePointers[1] = 0;

	connectedExternalNodes(0) = 0;
	connectedExternalNodes(1) = 0;
}

//*********************************************************************
//full constructor

PDTrussWithShear3D::PDTrussWithShear3D(int tag, int node1, int node2, double a, NDMaterial* theMat) :
	Element(tag, ELE_TAG_PDTrussWithShear3D), A(a), connectedExternalNodes(2), load(0), Ki(0),
	parallel(3), perpendicular1(3), perpendicular2(3), L(0), T(3, 3), Ttran(3, 3), B(3, 6), dirct(2)
{
	connectedExternalNodes(0) = node1;
	connectedExternalNodes(1) = node2;

	for (int i = 0; i < 2; i++)
	{
		nodePointers[i] = 0;
	}

	theMaterial = theMat->getCopy();
	if (theMaterial == 0) {
		opserr << "PDTrussWithShear3D::PDTrussWithShear3D - failed to allocate sectionN model pointer\n";
		exit(-1);
	}
}
//******************************************************************

//destructor 
PDTrussWithShear3D::~PDTrussWithShear3D()
{
	delete theMaterial;

	if (Ki != 0)
		delete Ki;
}
//**************************************************************************

//set domain
void  PDTrussWithShear3D::setDomain(Domain* theDomain)
{
	for (int i = 0; i < 2; i++)
		nodePointers[i] = theDomain->getNode(connectedExternalNodes(i));

	const Vector& crdNode0 = nodePointers[0]->getCrds();
	const Vector& crdNode1 = nodePointers[1]->getCrds();

	parallel = crdNode1 - crdNode0;
	L = parallel.Norm();

	if (L < 1e-13) {
		opserr << "PDTrussWithShear3D::PDTrussWithShear3D - failed to get height or length\n";
		exit(-1);
	}

	parallel = parallel / L;

	if (abs(parallel(0)) > 1.e-14)
	{
		perpendicular1(0) = -parallel(1);
		perpendicular1(1) = parallel(0);
	}
	else if (abs(parallel(1)) > 1.e-14)
	{
		perpendicular1(1) = -parallel(2);
		perpendicular1(2) = parallel(1);
	}
	else
	{
		perpendicular1(0) = -parallel(2);
		perpendicular1(2) = parallel(0);
	}
	perpendicular2(0) = parallel(1) * perpendicular1(2) - parallel(2) * perpendicular1(1);
	perpendicular2(1) = parallel(2) * perpendicular1(0) - parallel(0) * perpendicular1(2);
	perpendicular2(2) = parallel(0) * perpendicular1(1) - parallel(1) * perpendicular1(0);

	for (int i = 0; i < 3; i++)
	{
		T(0, i) = parallel(i);
		T(1, i) = perpendicular1(i);
		T(2, i) = perpendicular2(i);
	}

	Ttran.addMatrixTranspose(0.0, T, 1.0);

	B(0, 0) = -1 / L;
	B(1, 1) = -1 / L;
	B(2, 2) = -1 / L;
	B(0, 3) = 1 / L;
	B(1, 4) = 1 / L;
	B(2, 5) = 1 / L;

	this->DomainComponent::setDomain(theDomain);
}

//get the number of external nodes
int  PDTrussWithShear3D::getNumExternalNodes() const
{
	return 2;
}

//return connected external nodes
const ID& PDTrussWithShear3D::getExternalNodes()
{
	return connectedExternalNodes;
}

Node**
PDTrussWithShear3D::getNodePtrs(void)
{
	return nodePointers;
}

//return number of dofs
int PDTrussWithShear3D::getNumDOF()
{
	return 6;
}

//commit state
int PDTrussWithShear3D::commitState()
{
	int success = 0;

	// call element commitState to do any base class stuff
	if ((success = this->Element::commitState()) != 0) {
		opserr << "PDTrussWithShear3D::commitState () - failed in base class";
	}

	success += theMaterial->commitState();

	return success;
}

//revert to last commit 
int  PDTrussWithShear3D::revertToLastCommit()
{
	int success = 0;

	success += theMaterial->revertToLastCommit();

	return success;
}


//revert to start 
int PDTrussWithShear3D::revertToStart()
{
	int success = 0;

	success += theMaterial->revertToStart();

	return success;
}

//print out element data
void PDTrussWithShear3D::Print(OPS_Stream& s, int flag)
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
PDTrussWithShear3D::setResponse(const char** argv, int argc, OPS_Stream& output)
{
	Response* theResponse = 0;

	if (strcmp(argv[0], "-material") == 0 || strcmp(argv[0], "material") == 0)
	{
		theResponse = theMaterial->setResponse(&argv[1], argc - 1, output);
	}
	else
	{
		opserr << "this recorder type is not supported in PDTrussWithShear3D " << endln;
	}

	return theResponse;
}

int
PDTrussWithShear3D::getResponse(int responseID, Information& eleInfo)
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
			const Vector& sigma = theMaterial->getStress();
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
			const Vector& deformation = theMaterial->getStrain();
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
const Vector& PDTrussWithShear3D::getResistingForce()
{
	Vector stress = theMaterial->getStress();

	Vector stressEx(3);
	stressEx(0) = stress(0);
	stressEx(1) = stress(1) * dirct(0);
	stressEx(2) = stress(1) * dirct(1);

	Vector fb(3);
	fb = (Ttran * stressEx) * A;

	resid.Zero();
	for (int i = 0; i < 2; i++)
	{
		Matrix BItran(3, 3);
		for (int ii = 0; ii < 3; ii++)
			for (int jj = 0; jj < 3; jj++)
				BItran(jj, ii) = B(ii, i * 3 + jj);

		Vector Fb = BItran * fb * L;
		for (int j = 0; j < 3; j++)
			resid(3 * i + j) = Fb(j);
	}

	if (load != 0)
		resid -= *load;

	return resid;
}


//return stiffness matrix 
const Matrix& PDTrussWithShear3D::getTangentStiff()
{
	Matrix tangent = theMaterial->getTangent();

	Matrix tangentEx(3, 3);
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			tangentEx(i, j) = tangent(i, j);

	tangentEx(0, 2) = tangentEx(0, 1);
	tangentEx(2, 0) = tangentEx(1, 0);
	tangentEx(2, 2) = tangentEx(1, 1);

	Matrix kb(3, 3);
	kb = (Ttran * tangentEx * T) * A;

	stiff.Zero();
	for (int i = 0; i < 2; i++)
	{
		Matrix BI(3, 3);
		for (int ii = 0; ii < 3; ii++)
			for (int jj = 0; jj < 3; jj++)
				BI(ii, jj) = B(ii, i * 3 + jj);

		for (int j = 0; j < 2; j++)
		{
			Matrix BJtran(3, 3);
			for (int ii = 0; ii < 3; ii++)
				for (int jj = 0; jj < 3; jj++)
					BJtran(jj, ii) = B(ii, j * 3 + jj);

			Matrix KIJ = BJtran * kb * BI * L;

			for (int k = 0; k < 3; k++)
				for (int l = 0; l < 3; l++)
					stiff(3 * i + k, 3 * j + l) = KIJ(k, l);
		}
	}

	return stiff;
}


//return mass matrix
const Matrix& PDTrussWithShear3D::getMass()
{
	return mass;
}


//return secant matrix 
const Matrix& PDTrussWithShear3D::getInitialStiff()
{
	if (Ki != 0)
		return *Ki;

	Matrix tangent = theMaterial->getInitialTangent();

	Matrix Ttran(3, 3);

	Matrix tangentEx(3, 3);
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			tangentEx(i, j) = tangent(i, j);

	Ttran.addMatrixTranspose(0.0, T, 1.0);
	stiff = (Ttran * tangentEx * T) * A;

	Ki = new Matrix(stiff);

	return stiff;
}


void  PDTrussWithShear3D::zeroLoad()
{
	if (load != 0)
		load->Zero();

	return;
}


int
PDTrussWithShear3D::addLoad(ElementalLoad* theLoad, double loadFactor)
{
	opserr << "PDTrussWithShear3D::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
	return -1;
}



int
PDTrussWithShear3D::addInertiaLoadToUnbalance(const Vector& accel)
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
		const Vector& Raccel = nodePointers[i]->getRV(accel);
		for (int j = 0; j < 6; j++)
			r(count++) = Raccel(j);
	}

	if (load == 0)
		load = new Vector(24);

	load->addMatrixVector(1.0, mass, r, -1.0);

	return 0;
}


//get residual with inertia terms
const Vector& PDTrussWithShear3D::getResistingForceIncInertia()
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
PDTrussWithShear3D::update(void)
{
	int success = 0;

	Vector node1Disp, node2Disp, deltaDisp, detDispParallel(3), detDispPerpendicular(3);
	node1Disp = nodePointers[0]->getTrialDisp();
	node2Disp = nodePointers[1]->getTrialDisp();

	deltaDisp = node2Disp - node1Disp;
	double detDispParallelNorm = deltaDisp ^ parallel;
	detDispPerpendicular = deltaDisp - detDispParallelNorm * parallel;
	double detDispPerpendicularNorm = detDispPerpendicular.Norm();

	Matrix B2(3, 3);
	for (int ii = 0; ii < 3; ii++)
		for (int jj = 0; jj < 3; jj++)
			B2(ii, jj) = B(ii, 3 + jj);

	Vector strain(2), strainEx;

	Vector tmp = T * deltaDisp;

	strainEx = B2 * T * deltaDisp;
	strain(0) = strainEx(0);
	strain(1) = sqrt(strainEx(1) * strainEx(1) + strainEx(2) * strainEx(2));

	Vector crd1, crd2;
	crd1 = nodePointers[0]->getCrds();
	crd2 = nodePointers[1]->getCrds();

	success += theMaterial->setTrialStrain(strain);

	if (strain(1) > 1.e-14)
	{
		dirct(0) = strainEx(1) / strain(1);
		dirct(1) = strainEx(2) / strain(1);
	}

	Vector stress = theMaterial->getStress();
	Vector stressEx(3);
	stressEx(0) = stress(0);
	stressEx(1) = stress(1) * dirct(0);
	stressEx(2) = stress(1) * dirct(1);

	Vector fb(3);
	fb = (Ttran * stressEx) * A;

	return 0;
}


int  PDTrussWithShear3D::sendSelf(int commitTag, Channel& theChannel)
{
	int res = 0;

	return res;
}

int  PDTrussWithShear3D::recvSelf(int commitTag,
	Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
	int res = 0;

	return res;
}

int
PDTrussWithShear3D::displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes, int numMode)
{

	int error = 0;

	return error;
}

