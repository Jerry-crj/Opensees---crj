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
// $Source: /usr/local/cvs/OpenSees/SRC/element/shell/ShearWall.cpp,v $

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
#include "ShearWall.h"
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


#define DEBUG_ShearWall
#define min(a,b) ( (a)<(b) ? (a):(b) )

//static int numShearWall = 0;

//static data
Matrix  ShearWall::stiff(24, 24);
Vector  ShearWall::resid(24);
Matrix  ShearWall::mass(24, 24);
static Matrix B(8, 4);
//static double workArea[100];

//extern void printCommand(int argc, TCL_Char **argv);

void *
OPS_ShearWall(void)
{
	if (OPS_GetNumRemainingInputArgs() < 6) {
		opserr << "WARNING insufficient arguments\n";
		opserr << "want: eleTag? Node1? Node2? Node3? Node4? secTag?\n";
		return 0;
	}

	// get the id and end nodes 
	int idata[6];
	int num = 6;
	if (OPS_GetIntInput(&num, idata) < 0) {
		opserr << "WARNING: invalid integer data\n";
		return 0;
	}

//	int numSec = 1;
//	int numSecN = 1;

	int sectionID = idata[5];

	SectionForceDeformation *theSec = OPS_getSectionForceDeformation(sectionID);
	if (theSec == 0) {
		opserr << "WARNING section not found\n";
		opserr << "Section: " << sectionID;
		return 0;
	}

	BeamIntegration *bi = new LegendreBeamIntegration();

	Element *theElement = new ShearWall(idata[0], idata[1], idata[2], idata[3], idata[4],
		theSec, bi);

//	delete[] theSec;
	return theElement;
}


//null constructor
ShearWall::ShearWall() :
	Element(0, ELE_TAG_ShearWall), theSections(0),
	connectedExternalNodes(4), applyLoad(0), load(0), H(0), Ki(0)
{
	B.Zero();
	//for (int i = 0; i < 4; i++)
	//	theSections[i] = 0;
	connectedExternalNodes(0) = 0;
	connectedExternalNodes(1) = 0;
	connectedExternalNodes(2) = 0;
	connectedExternalNodes(3) = 0;

}

//*********************************************************************
//full constructor

ShearWall::ShearWall(int tag, int node1, int node2, int node3, int node4,
	SectionForceDeformation *s,BeamIntegration *bi) :
	Element(tag, ELE_TAG_ShearWall), theSections(s), 
	connectedExternalNodes(4), applyLoad(0), load(0), Ki(0)
{

	B.Zero();

	connectedExternalNodes(0) = node1;
	connectedExternalNodes(1) = node2;
	connectedExternalNodes(2) = node3;
	connectedExternalNodes(3) = node4;

	for (int i = 0; i < 4; i++)
	{
		nodePointers[i] = 0;
	}

	theSections = s->getCopy();
	if (theSections == 0) {
		opserr << "ShearWall::ShearWall - failed to allocate sectionN model pointer\n";
		exit(-1);
	}

	beamInt = bi->getCopy();

	if (beamInt == 0) {
		opserr << "ShearWall::ShearWall - failed to copy beam integration\n";
		exit(-1);
	}

}
//******************************************************************

//destructor 
ShearWall::~ShearWall()
{
	delete theSections;

	if (crdTransf)
		delete crdTransf;

	if (crdTransfTran)
		delete crdTransfTran;

	if (Ki != 0)
		delete Ki;
}
//**************************************************************************

//set domain
void  ShearWall::setDomain(Domain *theDomain)
{
	//node pointers
	for (int i = 0; i < 4; i++)
		nodePointers[i] = theDomain->getNode(connectedExternalNodes(i));

	const Vector &crdNode0 = nodePointers[0]->getCrds();
	const Vector &crdNode1 = nodePointers[1]->getCrds();
	const Vector &crdNode2 = nodePointers[2]->getCrds();
	Vector x = crdNode1 - crdNode0;
	Vector y = crdNode2 - crdNode1;
	L = x.Norm();
	H = y.Norm();

	if ((H < 1e-13 && H>-1e-13) || (L < 1e-13 && L>-1e-13)) {
		opserr << "ShearWall::ShearWall - failed to get height or length\n";
		exit(-1);
	}

	this->DomainComponent::setDomain(theDomain);
}

//get the number of external nodes
int  ShearWall::getNumExternalNodes() const
{
	return 4;
}

//return connected external nodes
const ID&  ShearWall::getExternalNodes()
{
	return connectedExternalNodes;
}

Node **
ShearWall::getNodePtrs(void)
{
	return nodePointers;
}

//return number of dofs
int ShearWall::getNumDOF()
{
	return 24;
}

//commit state
int ShearWall::commitState()
{
	int success = 0;

	// call element commitState to do any base class stuff
	if ((success = this->Element::commitState()) != 0) {
		opserr << "ShearWall::commitState () - failed in base class";
	}

	success += theSections->commitState();

	return success;
}

//revert to last commit 
int  ShearWall::revertToLastCommit()
{
	int success = 0;

	success += theSections->revertToLastCommit();

	return success;
}


//revert to start 
int ShearWall::revertToStart()
{
	int success = 0;

	success += theSections->revertToStart();

	return success;
}

//print out element data
void ShearWall::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_CURRENTSTATE) {
		s << endln;
		s << "Element type: Shear Wall Shell, tag: " << this->getTag() << endln;
		s << "\tConnected external nodes: " << connectedExternalNodes;
		s << "\tNumber of sections : 1" << endln;
		s << "\nSection Information: " << endln;

		theSections->Print(s, OPS_PRINT_PRINTMODEL_MATERIAL);

		s << endln;
	}
}

Response*
ShearWall::setResponse(const char **argv, int argc, OPS_Stream &output)
{
	Response *theResponse = 0;

	if (strcmp(argv[0], "section") == 0)
	{
		if (argc > 2)
		{
			theResponse = theSections->setResponse(&argv[1], argc - 1, output);
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
ShearWall::getResponse(int responseID, Information &eleInfo)
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
			const Vector &sigma = theSections->getStressResultant();
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
			const Vector &deformation = theSections->getSectionDeformation();
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
const Vector&  ShearWall::getResistingForce()
{
	formResistingForce();

	// subtract external loads 
	if (load != 0)
		resid -= *load;

	return resid;
}


//return stiffness matrix 
const Matrix& ShearWall::getTangentStiff()
{	//do tangent and residual here

	formTangentStiff();

	return stiff;
}


//return mass matrix
const Matrix&  ShearWall::getMass()
{

	formMass();

	return mass;
}


//return secant matrix 
const Matrix& ShearWall::getInitialStiff()
{
	if (Ki != 0)
		return *Ki;

	//strains ordered : du/dx  dv/dx  du/dy+dv/dy  dty/dx 
	//				 -dtx/dy   dty/dy-dtx/dx  dw/dy-tx  dw/dx+ty
	static const int ndm = 3;
	static const int ndf = 6;
	static const int nstress = 8;
	static const int numberNodes = 4;
	static const int numberVirNode = 2; // maybe bug here
	static const int nShape = 2; // maybe bug here

	double xi = 0.5;
	double dvol = 1.0;

	static double detJacobi = H;  // determinant of jacaobian matrix 
	static double shp[nShape][numberVirNode];  //shape functions at a gauss point
	//static double Shape[nShape][numberVirNode][numSectionsN + numSectionsS]; //all the shape functions

	static Vector residJ(ndf - 2); //nodeJ residual 
	static Matrix stiffJK(ndf - 2, ndf - 2); //nodeJK stiffness 
	static Vector secResut(nstress);  //stress
//	static Matrix secStiff(nstress, nstress);  //material tangent

	//---------B-matrices------------------------------------
	static Matrix BJ(nstress, ndf - 2);      // B matrix node J
	static Matrix BJtran(ndf - 2, nstress);
	static Matrix BK(nstress, ndf - 2);      // B matrix node k
	static Matrix BJtranD(ndf - 2, nstress);

	int p, q;

	stiff.Zero();

	computeBasis();

	//calculate shape functions
	double gaussCoor = xi;
	getShapeFnc(gaussCoor, shp);

	const Matrix& secStiff = theSections->getInitialTangent();

	int jj = 0;
	for (int j = 0; j < numberNodes; j++)
	{
		BJ = computeB(j, shp);

		//multiply by volume element
		//secStiff *= dvol[i];

		for (int p = 0; p < nstress; p++)
		{
			for (int q = 0; q < ndf; q++)
			{
				BJtran(q, p) = BJ(p, q);
			}
		}

		BJtranD = BJtran * secStiff;

		int kk = 0;
		for (int k = 0; k < numberNodes; k++)
		{
			BK = computeB(k, shp);
			stiffJK = BJtranD * BK;

			for (int p = 0; p < nstress; p++)
			{
				for (int q = 0; q < ndf; q++)
				{
					stiff(jj + p, kk + q) = stiffJK(p, q)*dvol * H;
				}
			}
			kk += ndf; // end j
		}
		jj += ndf; //end i
	}

	int ii = 0;
	for (int i = 0; i < numberNodes * 2; i++) // disp and rotation per node
	{
		int jj = 0;
		for (int j = 0; j < numberNodes * 2; j++)
		{
			Matrix stiffSub(3, 3), stiffTransfSub(3, 3);

			for (p = 0; p < 3; p++)
			{
				for (q = 0; q < 3; q++)
				{
					stiffSub(p, q) = stiff(ii + p, jj + q);
				}
			}

			stiffTransfSub = (*crdTransfTran) * stiffSub * (*crdTransf);

			for (p = 0; p < 3; p++)
			{
				for (q = 0; q < 3; q++)
				{
					stiff(ii + p, jj + q) = stiffTransfSub(p, q);
				}
			}

			jj += ndf / 2;
		} // end j
		ii += ndf / 2;
	} //end i

	Ki = new Matrix(stiff);

	return stiff;
}


void  ShearWall::zeroLoad()
{
	if (load != 0)
		load->Zero();

	applyLoad = 0;

	appliedB[0] = 0.0;
	appliedB[1] = 0.0;
	appliedB[2] = 0.0;
	appliedB[3] = 0.0;
	appliedB[4] = 0.0;
	appliedB[5] = 0.0;

	return;
}


int
ShearWall::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	opserr << "ShearWall::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
	return -1;
}



int
ShearWall::addInertiaLoadToUnbalance(const Vector &accel)
{

	static Vector r(24);

	int i;

	int allRhoZero = 0;
	if (theSections->getRho() != 0.0)
		allRhoZero = 1;

	if (allRhoZero == 0)
		return 0;

	formInertiaLoad();

	int count = 0;
	for (i = 0; i < 4; i++) {
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
const Vector&  ShearWall::getResistingForceIncInertia()
{
	static Vector res(24);

	formResistingForce();
	formInertiaLoad();

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
ShearWall::update(void)
{

	//strains ordered : du/dx  dv/dx  du/dy+dv/dy  dty/dx 
	//				 -dtx/dy   dty/dy-dtx/dx  dw/dy-tx  dw/dx+ty
	static const int ndm = 3;
	static const int ndf = 6;
	static const int nstress = 8;
	static const int numberNodes = 4;
	static const int numberVirNodes = 2;
	static const int nShape = 2; // maybe bug here

	//--------------------- normal ----------------

	double xi = 0.5;
	double dvol = 1.0;

	static double shp[nShape][numberVirNodes];  //shape functions at a gauss point


	//---------B-matrices------------------------------------
	static Matrix BJ(nstress, ndf);      // B matrix node J

	int p, q;

	computeBasis();

	int success;

	static Vector strain(nstress);

	success = 0;
	//calculate shape functions
	double gaussCoor = xi;
	getShapeFnc(gaussCoor, shp);

	strain.Zero();

	for (int j = 0; j < numberNodes; j++)
	{
		const Vector &ul = nodePointers[j]->getTrialDisp();
		Vector ulDisp(3), ulRot(3), ulDispBasic(3), ulRotBasic(3), ulBasic(4);
		ulDisp(0) = ul(0); ulDisp(1) = ul(1); ulDisp(2) = ul(2);
		ulRot(0) = ul(3); ulRot(1) = ul(4); ulRot(2) = ul(5);

		ulDispBasic = (*crdTransf) * ulDisp;
		ulRotBasic = (*crdTransf) * ulRot;

		ulBasic(0) = ulDispBasic(0);
		ulBasic(1) = ulDispBasic(1);
		ulBasic(2) = ulDispBasic(2);
		ulBasic(3) = ulRotBasic(0);

		BJ = computeB(j, shp);

		Vector BJu(nstress);
		BJu = BJ * ulBasic;
		strain += BJu;
	}
	success += theSections->setTrialSectionDeformation(strain);

	return 0;
}


//**********************************************************************

int  ShearWall::sendSelf(int commitTag, Channel &theChannel)
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
	//	opserr << "WARNING ShearWall::sendSelf() - " << this->getTag() << " failed to send ID\n";
	//	return res;
	//}

	//static Vector vectData(4);
	//vectData(0) = alphaM;
	//vectData(1) = betaK;
	//vectData(2) = betaK0;
	//vectData(3) = betaKc;

	//res += theChannel.sendVector(dataTag, commitTag, vectData);
	//if (res < 0) {
	//	opserr << "WARNING ShearWall::sendSelf() - " << this->getTag() << " failed to send ID\n";
	//	return res;
	//}

	//// Finally, quad asks its material objects to send themselves
	//for (i = 0; i < 4; i++) {
	//	res += theSectionsN->sendSelf(commitTag, theChannel);
	//	if (res < 0) {
	//		opserr << "WARNING ShearWall::sendSelf() - " << this->getTag() << " failed to send its Material\n";
	//		return res;
	//	}
	//}

	return res;
}

int  ShearWall::recvSelf(int commitTag,
	Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{
	int res = 0;

	//int dataTag = this->getDbTag();

	//static ID idData(14);
	//// Quad now receives the tags of its four external nodes
	//res += theChannel.recvID(dataTag, commitTag, idData);
	//if (res < 0) {
	//	opserr << "WARNING ShearWall::recvSelf() - " << this->getTag() << " failed to receive ID\n";
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
	//	opserr << "WARNING ShearWall::sendSelf() - " << this->getTag() << " failed to send ID\n";
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
	//			opserr << "ShearWall::recvSelf() - Broker could not create NDMaterial of class type" << matClassTag << endln;;
	//			return -1;
	//		}
	//		// Now receive materials into the newly allocated space
	//		theSectionsN->setDbTag(matDbTag);
	//		res += theSectionsN->recvSelf(commitTag, theChannel, theBroker);
	//		if (res < 0) {
	//			opserr << "ShearWall::recvSelf() - material " << i << "failed to recv itself\n";
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
	//				opserr << "ShearWall::recvSelf() - Broker could not create NDMaterial of class type" << matClassTag << endln;
	//				exit(-1);
	//			}
	//		}
	//		// Receive the material
	//		theSectionsN->setDbTag(matDbTag);
	//		res += theSectionsN->recvSelf(commitTag, theChannel, theBroker);
	//		if (res < 0) {
	//			opserr << "ShearWall::recvSelf() - material " << i << "failed to recv itself\n";
	//			return res;
	//		}
	//	}
	//}

	return res;
}
//**************************************************************************

int
ShearWall::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
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

void
ShearWall::formResistingForce()
{

	//strains ordered : du/dx  dv/dx  du/dy+dv/dy  dty/dx 
	//				 -dtx/dy   dty/dy-dtx/dx  dw/dy-tx  dw/dx+ty
	static const int ndm = 3;
	static const int ndf = 6;
	static const int nstress = 8;
	static const int numberNodes = 4;
	static const int numberVirNodes = 2; // maybe bug here
	static const int nShape = 2; // maybe bug here

	// ------------------- normal -------------------------
	double xi = 0.5;
	double dvol = 1;

	static double detJacobi = H;  // determinant of jacaobian matrix 
	static double shp[nShape][numberVirNodes];  //shape functions at a gauss point

	static Vector residJ(ndf - 2);  //nodeJ residual 
	static Matrix BJ(nstress, ndf - 2);  // B matrix node J
	static Matrix BJtran(ndf - 2, nstress);

	int p, q;

	resid.Zero();

	computeBasis();

	int count = 0;

	//calculate shape functions
	double gaussCoor = xi;
	getShapeFnc(gaussCoor, shp);

	//compute the stress

	const Vector &secResult = theSections->getStressResultant();

	int jj = 0;
	for (int j = 0; j < numberNodes; j++)
	{

		BJ = computeB(j, shp);

		for (int p = 0; p < nstress; p++)
			for (int q = 0; q < ndf - 2; q++)
				BJtran(q, p) = BJ(p, q);

		residJ.Zero();
		residJ.addMatrixVector(1.0, BJtran, secResult, 1.0);
		for (int p = 0; p < ndf - 2; p++)
			resid(jj + p) += residJ(p) * dvol * H;

		jj += ndf;
	}

	// coordinate transformation (from local to global)
	int ii = 0;
	for (int i = 0; i < numberNodes; i++)
	{
		Vector forceLocal(3), forceGlobal(3), momentLocal(3), momentGlobal(3);

		forceLocal(0) = resid(ii);
		forceLocal(1) = resid(ii + 1);
		forceLocal(2) = resid(ii + 2);
		momentLocal(0) = resid(ii + 3);
		momentLocal(1) = 0;
		momentLocal(2) = 0;

		forceGlobal = (*crdTransfTran) * forceLocal;
		momentGlobal = (*crdTransfTran) * momentLocal;

		resid(ii) = forceGlobal(0);
		resid(ii + 1) = forceGlobal(1);
		resid(ii + 2) = forceGlobal(2);
		resid(ii + 3) = momentGlobal(0);
		resid(ii + 4) = momentGlobal(1);
		resid(ii + 5) = momentGlobal(2);

		ii += ndf;
	}
	return;
}


void
ShearWall::formTangentStiff()
{

	//strains ordered : du/dx  dv/dx  du/dy+dv/dy  dty/dx 
	//				 -dtx/dy   dty/dy-dtx/dx  dw/dy-tx  dw/dx+ty
	static const int ndm = 3;
	static const int ndf = 6;
	static const int nstress = 8;
	static const int numberNodes = 4;
	static const int numberVirNodes = 2; // maybe bug here
	static const int nShape = 2; // maybe bug h
	// --------------- normal ------------------

	double xi = 0.5;
	double dvol = 1.0;

	static double detJacobi = H;  // determinant of jacaobian matrix 
	static double shp[nShape][numberVirNodes];  //shape functions at a gauss point

	static Matrix stiffJK(ndf - 2, ndf - 2); //nodeJK stiffness 
	static Matrix BJ(nstress, ndf - 2);      // B matrix node J
	static Matrix BJtran(ndf - 2, nstress);
	static Matrix BK(nstress, ndf - 2);      // B matrix node k
	static Matrix BJtranD(ndf - 2, nstress);

	int p, q;

	stiff.Zero();

	computeBasis();

	//calculate shape functions
	double gaussCoor = xi;
	getShapeFnc(gaussCoor, shp);

	const Matrix& secStiff = theSections->getSectionTangent();

	int	jj = 0;
	for (int j = 0; j < numberNodes; j++)
	{

		BJ = computeB(j, shp);

		for (int p = 0; p < nstress; p++)
			for (int q = 0; q < ndf; q++)
				BJtran(q, p) = BJ(p, q);

		BJtranD = BJtran * secStiff;

		int kk = 0;
		for (int k = 0; k < numberNodes; k++)
		{
			BK = computeB(k, shp);

			stiffJK = BJtranD * BK;

			for (int p = 0; p < ndf - 2; p++)
				for (int q = 0; q < ndf - 2; q++)
					stiff(jj + p, kk + q) += stiffJK(p, q)*dvol * H;

			kk += ndf; // end j
		}   // end k
		jj += ndf; //end i
	}  // end j

	double err = 1e-10;
	stiff(0, 0) *= 1 + err; stiff(0, 6) *= 1 - err; stiff(6, 0) *= 1 - err; stiff(6, 6) *= 1 + err;
	stiff(12, 12) *= 1 + err; stiff(12, 18) *= 1 - err; stiff(18, 12) *= 1 - err; stiff(18, 18) *= 1 + err;

	stiff(4, 4) = err; stiff(4, 10) = -err; stiff(10, 4) = -err; stiff(10, 10) = err;
	stiff(16, 16) = err; stiff(16, 22) = -err; stiff(22, 16) = -err; stiff(22, 22) = err;

	stiff(5, 5) = err; stiff(5, 11) = -err; stiff(11, 5) = -err; stiff(11, 11) = err;
	stiff(17, 17) = err; stiff(17, 23) = -err; stiff(23, 17) = -err; stiff(23, 23) = err;

	int ii = 0;
	for (int i = 0; i < numberNodes * 2; i++) // disp and rotation per node
	{
		int jj = 0;
		for (int j = 0; j < numberNodes * 2; j++)
		{
			Matrix stiffSub(3, 3), stiffTransfSub(3, 3);

			for (p = 0; p < 3; p++)
			{
				for (q = 0; q < 3; q++)
				{
					stiffSub(p, q) = stiff(ii + p, jj + q);
				}
			}

			stiffTransfSub = (*crdTransfTran) * stiffSub * (*crdTransf);

			for (p = 0; p < 3; p++)
			{
				for (q = 0; q < 3; q++)
				{
					stiff(ii + p, jj + q) = stiffTransfSub(p, q);
				}
			}

			jj += ndf / 2;
		} // end j
		ii += ndf / 2;
	} //end i
	return;
}


void ShearWall::formMass()
{

}


void ShearWall::formInertiaLoad()
{

}


const Matrix &
ShearWall::computeB(int iNode, double shp[2][2])
{
	switch (iNode)
	{
	case 0:
		B(0, 1) = +shp[1][0] / 2;
		B(1, 1) = -shp[1][0] / L;
		B(2, 3) = +shp[1][0] / 2;
		B(3, 3) = -shp[1][0] / L;
		B(4, 0) = +shp[1][0] / 2;
		B(4, 1) = -shp[0][0] / L;
		B(5, 3) = -shp[0][0] / L;
		B(6, 2) = +shp[1][0] / 2;
		B(6, 3) = -shp[0][0] / 2;
		B(7, 2) = +shp[1][0] / L;
		B(7, 3) = -shp[0][0] / L;
		break;
	case 1:
		B(0, 1) = +shp[1][0] / 2;
		B(1, 1) = +shp[1][0] / L;
		B(2, 3) = +shp[1][0] / 2;
		B(3, 3) = +shp[1][0] / L;
		B(4, 0) = +shp[1][0] / 2;
		B(4, 1) = +shp[0][0] / L;
		B(5, 3) = +shp[0][0] / L;
		B(6, 2) = +shp[1][0] / 2;
		B(6, 3) = -shp[0][0] / 2;
		B(7, 2) = -shp[1][0] / L;
		B(7, 3) = +shp[0][0] / L;
		break;
	case 2:
		B(0, 1) = +shp[1][1] / 2;
		B(1, 1) = +shp[1][1] / L;
		B(2, 3) = +shp[1][1] / 2;
		B(3, 3) = +shp[1][1] / L;
		B(4, 0) = +shp[1][1] / 2;
		B(4, 1) = +shp[0][1] / L;
		B(5, 3) = +shp[0][1] / L;
		B(6, 2) = +shp[1][1] / 2;
		B(6, 3) = -shp[0][1] / 2;
		B(7, 2) = -shp[1][1] / L;
		B(7, 3) = +shp[0][1] / L;
		break;
	case 3:
		B(0, 1) = +shp[1][1] / 2;
		B(1, 1) = -shp[1][1] / L;
		B(2, 3) = +shp[1][1] / 2;
		B(3, 3) = -shp[1][1] / L;
		B(4, 0) = +shp[1][1] / 2;
		B(4, 1) = -shp[0][1] / L;
		B(5, 3) = -shp[0][1] / L;
		B(6, 2) = +shp[1][1] / 2;
		B(6, 3) = -shp[0][1] / 2;
		B(7, 2) = +shp[1][1] / L;
		B(7, 3) = -shp[0][1] / L;
		break;
	}
	return B;
};

void
ShearWall::getShapeFnc(double loc, double shp[2][2])
{
	// the information to calculate B
	shp[0][0] = 1 - loc;
	shp[0][1] = loc;
	shp[1][0] = -1 / H;
	shp[1][1] = 1 / H;
}


void
ShearWall::computeBasis()
{
	// compute the global - local transform matrix 
	const Vector &crdNode0 = nodePointers[0]->getCrds();
	const Vector &crdNode1 = nodePointers[1]->getCrds();
	const Vector &crdNode2 = nodePointers[2]->getCrds();
	Vector x = crdNode1 - crdNode0;
	Vector y = crdNode2 - crdNode1;
	double xLen = x.Norm();
	double yLen = y.Norm();

	for (int i = 0; i < 3; i++)
	{
		x(i) = x(i) / xLen;
		y(i) = y(i) / yLen;
	}

	Vector z(3);
	z(0) = x(1)*y(2) - x(2)*y(1);
	z(1) = x(2)*y(0) - x(0)*y(2);
	z(2) = x(0)*y(1) - x(1)*y(0);

	crdTransf = new Matrix(3, 3);
	for (int i = 0; i < 3; i++)
	{
		(*crdTransf)(0, i) = x(i);
		(*crdTransf)(1, i) = y(i);
		(*crdTransf)(2, i) = z(i);
	}

	crdTransfTran = new Matrix(3, 3);
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			(*crdTransfTran)(i, j) = (*crdTransf)(j, i);
		}
	}

	Vector temp(3), rst(3);
	temp.Zero();
	temp(1) = 1;
	rst = (*crdTransf)*temp;

}


