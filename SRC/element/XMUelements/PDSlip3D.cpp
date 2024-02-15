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

// $Revision: 1.35 $
// $Date: 2009-10-13 21:14:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/PDSlip/PDSlip.cpp,v $

// Written: MHS
// Created: Feb 2000
// Revised: Dec 2000 for efficiency
//
// Description: This file contains the class definition for PDSlip.

#include <PDSlip3D.h>
#include <Node.h>
#include <UniaxialMaterial.h>
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
#include <ElementalLoad.h>
#include <elementAPI.h>

void* OPS_PDSlip3D()
{
	int ndm = OPS_GetNDM();
	int ndf = OPS_GetNDF();

	if (ndm != 3 || ndf != 3) {
		opserr << "WARNING -- model dimensions and/or nodal DOF not compatible with quad element\n";
		return 0;
	}

	if (OPS_GetNumRemainingInputArgs() < 10) {
		opserr << "WARNING insufficient arguments\n";
		opserr << "Want: element PDSlip3D eleTag? iNode? jNode? kNode? lNode? thk? type? matTag? <pressure? rho? b1? b2?>\n";
		return 0;
	}

	// PDSlip3DId, iNode, jNode, kNode, lNode
	int idata[10];
	int num = 10;
	if (OPS_GetIntInput(&num, idata) < 0) {
		opserr << "WARNING: invalid integer inputs\n";
		return 0;
	}

	int matTag;
	num = 1;
	if (OPS_GetIntInput(&num, &matTag) < 0) {
		opserr << "WARNING: invalid matTag\n";
		return 0;
	}

	UniaxialMaterial* mat1 = OPS_getUniaxialMaterial(matTag);
	if (mat1 == 0) {
		opserr << "WARNING material not found\n";
		opserr << "Material: " << matTag;
		opserr << "\nPDSlip3D element: " << idata[0] << endln;
		return 0;
	}

	if (OPS_GetIntInput(&num, &matTag) < 0) {
		opserr << "WARNING: invalid matTag\n";
		return 0;
	}

	UniaxialMaterial* mat2 = OPS_getUniaxialMaterial(matTag);
	if (mat2 == 0) {
		opserr << "WARNING material not found\n";
		opserr << "Material: " << matTag;
		opserr << "\nPDSlip3D element: " << idata[0] << endln;
		return 0;
	}

	double data[3];
	Vector d(data, 3);
	num = 3;

	if (OPS_GetDoubleInput(&num, data) < 0) {
		opserr << "WARNING: invalid double inputs\n";
		return 0;
	}

	return new PDSlip3D(idata[0], idata[1], idata[2], idata[3], idata[4], idata[5],
		idata[6], idata[7], idata[8], idata[9], *mat1, *mat2, data[0], data[1], data[2]);
}


double PDSlip3D::matrixData[729];
Matrix PDSlip3D::K(matrixData, 27, 27);
Vector PDSlip3D::P(27);
double PDSlip3D::shp[8];

PDSlip3D::PDSlip3D(int tag, int nd1, int nd2, int nd3, int nd4, int nd5,
	int nd6, int nd7, int nd8, int nd9, 
	UniaxialMaterial& m1, UniaxialMaterial& m2, double p1, double p2, double p3)
	:Element(tag, ELE_TAG_PDSlip3D),
	theSprMaterial(0), connectedExternalNodes(9),
	Ki(0), dv(3), T(3, 3), Ttran(3, 3), pst(3), deform(0), force(0)
{
	// Allocate arrays of pointers to NDMaterials
	theSprMaterial = new UniaxialMaterial * [3];

	if (theSprMaterial == 0) {
		opserr << "PDSlip3D::PDSlip3D - failed allocate material model pointer\n";
		exit(-1);
	}

	theSprMaterial[0] = m1.getCopy();
	theSprMaterial[1] = m2.getCopy();
	theSprMaterial[2] = m2.getCopy();

	// Check allocation
	if ((theSprMaterial[0] == 0) || (theSprMaterial[1] == 0) || (theSprMaterial[2] == 0)) {
		opserr << "PDSlip3D::PDSlip3D -- failed to get a copy of material model\n";
		exit(-1);
	}

	// Set connected external node IDs
	connectedExternalNodes(0) = nd1;
	connectedExternalNodes(1) = nd2;
	connectedExternalNodes(2) = nd3;
	connectedExternalNodes(3) = nd4;
	connectedExternalNodes(4) = nd5;
	connectedExternalNodes(5) = nd6;
	connectedExternalNodes(6) = nd7;
	connectedExternalNodes(7) = nd8;
	connectedExternalNodes(8) = nd9;

	for (int i = 0; i < 9; i++)
		theNodes[i] = 0;

	dv(0) = p1;
	dv(1) = p2;
	dv(2) = p3;

	if (dv.Norm() < 1e-10)
		opserr << "PDSlip3D::PDSlip3D -- Dirction Vector.Norm is zero \n";

	dv = dv / dv.Norm();

	this->getTransf();

}

PDSlip3D::PDSlip3D()
	:Element(0, ELE_TAG_PDSlip3D),
	theSprMaterial(0), connectedExternalNodes(9),
	Ki(0), dv(3), T(3, 3), Ttran(3, 3), pst(3), deform(0), force(0)
{
	for (int i = 0; i < 9; i++)
		theNodes[i] = 0;
}

PDSlip3D::~PDSlip3D()
{
	for (int i = 0; i < 3; i++) {
		if (theSprMaterial[i])
			delete theSprMaterial[i];
	}

	// Delete the array of pointers to NDMaterial pointer arrays
	if (theSprMaterial)
		delete[] theSprMaterial;

	if (Ki != 0)
		delete Ki;
}

int
PDSlip3D::getNumExternalNodes() const
{
	return 9;
}

const ID&
PDSlip3D::getExternalNodes()
{
	return connectedExternalNodes;
}


Node**
PDSlip3D::getNodePtrs(void)
{
	return theNodes;
}

int
PDSlip3D::getNumDOF()
{
	return 27;
}

void
PDSlip3D::setDomain(Domain* theDomain)
{
	// Check Domain is not null - invoked when object removed from a domain
	if (theDomain == 0) {
		theNodes[0] = 0;
		theNodes[1] = 0;
		theNodes[2] = 0;
		theNodes[3] = 0;
		theNodes[4] = 0;
		theNodes[5] = 0;
		theNodes[6] = 0;
		theNodes[7] = 0;
		theNodes[8] = 0;
		return;
	}

	int Nd1 = connectedExternalNodes(0);
	int Nd2 = connectedExternalNodes(1);
	int Nd3 = connectedExternalNodes(2);
	int Nd4 = connectedExternalNodes(3);
	int Nd5 = connectedExternalNodes(4);
	int Nd6 = connectedExternalNodes(5);
	int Nd7 = connectedExternalNodes(6);
	int Nd8 = connectedExternalNodes(7);
	int Nd9 = connectedExternalNodes(8);

	theNodes[0] = theDomain->getNode(Nd1);
	theNodes[1] = theDomain->getNode(Nd2);
	theNodes[2] = theDomain->getNode(Nd3);
	theNodes[3] = theDomain->getNode(Nd4);
	theNodes[4] = theDomain->getNode(Nd5);
	theNodes[5] = theDomain->getNode(Nd6);
	theNodes[6] = theDomain->getNode(Nd7);
	theNodes[7] = theDomain->getNode(Nd8);
	theNodes[8] = theDomain->getNode(Nd9);

	if (theNodes[0] == 0 || theNodes[1] == 0 || theNodes[2] == 0 || theNodes[3] == 0 || theNodes[4] == 0
		|| theNodes[5] == 0 || theNodes[6] == 0 || theNodes[7] == 0 || theNodes[8] == 0) {
		//opserr << "FATAL ERROR PDSlip3D (tag: %d), node not found in domain",
		//	this->getTag());

		return;
	}

	int dofNd1 = theNodes[0]->getNumberDOF();
	int dofNd2 = theNodes[1]->getNumberDOF();
	int dofNd3 = theNodes[2]->getNumberDOF();
	int dofNd4 = theNodes[3]->getNumberDOF();
	int dofNd5 = theNodes[4]->getNumberDOF();
	int dofNd6 = theNodes[5]->getNumberDOF();
	int dofNd7 = theNodes[6]->getNumberDOF();
	int dofNd8 = theNodes[7]->getNumberDOF();
	int dofNd9 = theNodes[8]->getNumberDOF();

	if (dofNd1 != 3 || dofNd2 != 3 || dofNd3 != 3 || dofNd4 != 3 || dofNd5 != 3
		|| dofNd6 != 3 || dofNd7 != 3 || dofNd8 != 3 || dofNd9 != 3) {
		//opserr << "FATAL ERROR PDSlip3D (tag: %d), has differing number of DOFs at its nodes",
		//	this->getTag());

		return;
	}

	this->DomainComponent::setDomain(theDomain);

	this->determinPst();

	if (pst(0) < -1 || pst(0) > 1 || pst(1) < -1 || pst(1) > 1 || pst(2) < -1 || pst(2) > 1)
	{
		opserr << "PDSlip3D::PDSlip3D -- Wrong range of coordinates in parameter space \n";
		exit(-1);
	}
}

int
PDSlip3D::commitState()
{
	int retVal = 0;

	// call element commitState to do any base class stuff
	if ((retVal = this->Element::commitState()) != 0) {
		opserr << "PDSlip3D::commitState () - failed in base class";
	}

	// Loop over the integration points and commit the material states

	for (int i = 0; i < 3; i++)
		retVal += theSprMaterial[i]->commitState();

	return retVal;
}

int
PDSlip3D::revertToLastCommit()
{
	int retVal = 0;

	// Loop over the integration points and revert to last committed state

	for (int i = 0; i < 3; i++)
		retVal += theSprMaterial[i]->revertToLastCommit();

	return retVal;
}

int
PDSlip3D::revertToStart()
{
	int retVal = 0;

	// Loop over the integration points and revert states to start

	for (int i = 0; i < 3; i++)
		retVal += theSprMaterial[i]->revertToStart();

	return retVal;
}


int
PDSlip3D::update()
{
	static Vector disp;

	static double u[3][8];

	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 3; j++) {
			disp = theNodes[i]->getTrialDisp();
			u[j][i] = disp(j);
		}
	}

	// Determine Jacobian for this integration point
	this->shapeFunction(pst(0), pst(1), pst(2));

	static Vector dispRef(3);

	dispRef.Zero();
	for (int beta = 0; beta < 8; beta++) {
		dispRef(0) += shp[beta] * u[0][beta];
		dispRef(1) += shp[beta] * u[1][beta];
		dispRef(2) += shp[beta] * u[2][beta];
	}

	disp = theNodes[8]->getTrialDisp();

	Vector du(3);
	du(0) = 1;
	du = T * (disp - dispRef);

	int ret = 0;

	// Set the material strain
	ret += theSprMaterial[0]->setTrialStrain(du(0));
	ret += theSprMaterial[1]->setTrialStrain(du(1));
	ret += theSprMaterial[2]->setTrialStrain(du(2));

	deform = du(0);

	return ret;
}


const Matrix&
PDSlip3D::getTangentStiff()
{

	K.Zero();

	Matrix ksg(3, 3);
	ksg(0, 0) = theSprMaterial[0]->getTangent();
	ksg(1, 1) = theSprMaterial[1]->getTangent();
	ksg(2, 2) = theSprMaterial[2]->getTangent();

	ksg = Ttran * ksg * T;

	this->shapeFunction(pst(0), pst(1), pst(2));


	for (int alpha = 0, ia = 0; alpha < 8; alpha++, ia += 3) {
		for (int beta = 0, ib = 0; beta < 8; beta++, ib += 3) {
			for (int ii = 0; ii < 3; ii++) {
				for (int jj = 0; jj < 3; jj++) {
					K(ia + ii, ib + jj) = shp[alpha] * shp[beta] * ksg(ii, jj);
				}
			}
		}
	}

	for (int alpha = 0, ia = 0; alpha < 8; alpha++, ia += 3) {
		for (int ii = 0; ii < 3; ii++) {
			for (int jj = 0; jj < 3; jj++) {
				K(ia + ii, 24 + jj) = -shp[alpha] * ksg(ii, jj);
				K(24 + ii, ia + jj) = -shp[alpha] * ksg(ii, jj);
			}
		}
	}

	for (int ii = 0; ii < 3; ii++) {
		for (int jj = 0; jj < 3; jj++) {
			K(24 + ii, 24 + jj) = ksg(ii, jj);
		}
	}

	return K;
}


const Matrix&
PDSlip3D::getInitialStiff()
{
	if (Ki != 0)
		return *Ki;

	K.Zero();

	Ki = new Matrix(K);

	Matrix ksg(3, 3);
	ksg(0, 0) = theSprMaterial[0]->getInitialTangent();
	ksg(1, 1) = theSprMaterial[1]->getInitialTangent();
	ksg(2, 2) = theSprMaterial[2]->getInitialTangent();

	ksg = Ttran * ksg * T;

	this->shapeFunction(pst(0), pst(1), pst(2));

	for (int alpha = 0, ia = 0; alpha < 8; alpha++, ia += 3) {
		for (int beta = 0, ib = 0; beta < 8; beta++, ib += 3) {
			for (int ii = 0; ii < 3; ii++) {
				for (int jj = 0; jj < 3; jj++) {
					(*Ki)(ia + ii, ib + jj) = shp[alpha] * shp[beta] * ksg(ii, jj);
				}
			}
		}
	}

	for (int alpha = 0, ia = 0; alpha < 8; alpha++, ia += 3) {
		for (int ii = 0; ii < 3; ii++) {
			for (int jj = 0; jj < 3; jj++) {
				(*Ki)(ia + ii, 24 + jj) = -shp[alpha] * ksg(ii, jj);
				(*Ki)(24 + ii, ia + jj) = -shp[alpha] * ksg(ii, jj);
			}
		}
	}

	for (int ii = 0; ii < 3; ii++) {
		for (int jj = 0; jj < 3; jj++) {
			(*Ki)(24 + ii, 24 + jj) = ksg(ii, jj);
		}
	}

	return *Ki;
}

const Matrix&
PDSlip3D::getMass()
{
	K.Zero();

	return K;
}

void
PDSlip3D::zeroLoad(void)
{
	return;
}

int
PDSlip3D::addLoad(ElementalLoad* theLoad, double loadFactor)
{
	// Added option for applying body forces in load pattern: C.McGann, U.Washington
	int type;
	const Vector& data = theLoad->getData(type, loadFactor);

	if (type == LOAD_TAG_SelfWeight) {
		return 0;
	}
	else {
		opserr << "PDSlip3D::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
		return -1;
	}

	return -1;
}

int
PDSlip3D::addInertiaLoadToUnbalance(const Vector& accel)
{
	return 0;
}

const Vector&
PDSlip3D::getResistingForce()
{
	P.Zero();

	Vector sig(3), sigg(3);
	sig(0) = theSprMaterial[0]->getStress();
	sig(1) = theSprMaterial[1]->getStress();
	sig(2) = theSprMaterial[2]->getStress();

	sigg = Ttran * sig;

	this->shapeFunction(pst(0), pst(1), pst(2));

	for (int alpha = 0, ia = 0; alpha < 8; alpha++, ia += 3) {
		for (int ii = 0; ii < 3; ii++) {
			P(ia+ii) = -shp[alpha] * sigg(ii);
		}
	}

	for (int ii = 0; ii < 3; ii++)
		P(24 + ii) = sigg(ii);

	force = sig(0);

	return P;
}

const Vector&
PDSlip3D::getResistingForceIncInertia()
{
	return P;
}

int
PDSlip3D::sendSelf(int commitTag, Channel& theChannel)
{
	int res = 0;

	// note: we don't check for dataTag == 0 for Element
	// objects as that is taken care of in a commit by the Domain
	// object - don't want to have to do the check if sending data
	int dataTag = this->getDbTag();

	// Quad packs its data into a Vector and sends this to theChannel
	// along with its dbTag and the commitTag passed in the arguments
	static Vector data(5);
	data(0) = this->getTag();

	data(1) = alphaM;
	data(2) = betaK;
	data(3) = betaK0;
	data(4) = betaKc;

	res += theChannel.sendVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING PDSlip3D::sendSelf() - " << this->getTag() << " failed to send Vector\n";
		return res;
	}

	// Now quad sends the ids of its materials
	int matDbTag;

	static ID idData(12);

	for (int i = 0; i < 2; i++)
	{
		idData(i) = theSprMaterial[i]->getClassTag();
		matDbTag = theSprMaterial[i]->getDbTag();

		if (matDbTag == 0)
		{
			matDbTag = theChannel.getDbTag();
			if (matDbTag != 0)
				theSprMaterial[i]->setDbTag(matDbTag);
		}

		idData[i + 2] = matDbTag;
	}

	idData(4) = connectedExternalNodes(0);
	idData(5) = connectedExternalNodes(1);
	idData(6) = connectedExternalNodes(2);
	idData(7) = connectedExternalNodes(3);
	idData(8) = connectedExternalNodes(4);

	res += theChannel.sendID(dataTag, commitTag, idData);

	if (res < 0) {
		opserr << "WARNING FourNodeQuad::sendSelf() - " << this->getTag() << " failed to send ID\n";
		return res;
	}

	// Finally, quad asks its material objects to send themselves
	for (int i = 0; i < 2; i++) {
		res += theSprMaterial[i]->sendSelf(commitTag, theChannel);
		if (res < 0) {
			opserr << "WARNING FourNodeQuad::sendSelf() - " << this->getTag() << " failed to send its Material\n";
			return res;
		}
	}

	opserr << "sendSelf() is called !!!!!!!" << endln;

	return res;
}

int
PDSlip3D::recvSelf(int commitTag, Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
	int res = 0;

	int dataTag = this->getDbTag();

	// Quad creates a Vector, receives the Vector and then sets the 
	// internal data with the data in the Vector
	static Vector data(9);
	res += theChannel.recvVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING PDSlip3D::recvSelf() - failed to receive Vector\n";
		return res;
	}

	this->setTag((int)data(0));

	alphaM = data(5);
	betaK = data(6);
	betaK0 = data(7);
	betaKc = data(8);

	static ID idData(12);
	// Quad now receives the tags of its four external nodes
	res += theChannel.recvID(dataTag, commitTag, idData);
	if (res < 0) {
		opserr << "WARNING PDSlip3D::recvSelf() - " << this->getTag() << " failed to receive ID\n";
		return res;
	}

	connectedExternalNodes(0) = idData(8);
	connectedExternalNodes(1) = idData(9);
	connectedExternalNodes(2) = idData(10);
	connectedExternalNodes(3) = idData(11);

	return res;
}

void
PDSlip3D::Print(OPS_Stream& s, int flag)
{
	if (flag == 2) {

		s << "#PDSlip3D\n";

		int i;
		const int numNodes = 4;
		const int nstress = 3;

		for (i = 0; i < numNodes; i++) {
			const Vector& nodeCrd = theNodes[i]->getCrds();
			const Vector& nodeDisp = theNodes[i]->getDisp();
			s << "#NODE " << nodeCrd(0) << " " << nodeCrd(1) << " " << endln;
		}

		// spit out the section location & invoke print on the scetion
		const int numMaterials = 4;

		static Vector avgStress(nstress);
		static Vector avgStrain(nstress);
		avgStress.Zero();
		avgStrain.Zero();

		avgStress /= numMaterials;
		avgStrain /= numMaterials;

		s << "#AVERAGE_STRESS ";
		for (i = 0; i < nstress; i++)
			s << avgStress(i) << " ";
		s << endln;

		s << "#AVERAGE_STRAIN ";
		for (i = 0; i < nstress; i++)
			s << avgStrain(i) << " ";
		s << endln;
	}
}

int
PDSlip3D::displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes, int numMode)
{

	// first set the quantity to be displayed at the nodes;
	// if displayMode is 1 through 3 we will plot material stresses otherwise 0.0

	static Vector values(4);

	for (int j = 0; j < 4; j++)
		values(j) = 0.0;


	// now  determine the end points of the quad based on
	// the display factor (a measure of the distorted image)
	// store this information in 4 3d vectors v1 through v4
	const Vector& end1Crd = theNodes[0]->getCrds();
	const Vector& end2Crd = theNodes[1]->getCrds();
	const Vector& end3Crd = theNodes[2]->getCrds();
	const Vector& end4Crd = theNodes[3]->getCrds();

	static Matrix coords(4, 3);

	if (displayMode >= 0) {

		const Vector& end1Disp = theNodes[0]->getDisp();
		const Vector& end2Disp = theNodes[1]->getDisp();
		const Vector& end3Disp = theNodes[2]->getDisp();
		const Vector& end4Disp = theNodes[3]->getDisp();

		for (int i = 0; i < 2; i++) {
			coords(0, i) = end1Crd(i) + end1Disp(i) * fact;
			coords(1, i) = end2Crd(i) + end2Disp(i) * fact;
			coords(2, i) = end3Crd(i) + end3Disp(i) * fact;
			coords(3, i) = end4Crd(i) + end4Disp(i) * fact;
		}
	}
	else {
		int mode = displayMode * -1;
		const Matrix& eigen1 = theNodes[0]->getEigenvectors();
		const Matrix& eigen2 = theNodes[1]->getEigenvectors();
		const Matrix& eigen3 = theNodes[2]->getEigenvectors();
		const Matrix& eigen4 = theNodes[3]->getEigenvectors();
		if (eigen1.noCols() >= mode) {
			for (int i = 0; i < 2; i++) {
				coords(0, i) = end1Crd(i) + eigen1(i, mode - 1) * fact;
				coords(1, i) = end2Crd(i) + eigen2(i, mode - 1) * fact;
				coords(2, i) = end3Crd(i) + eigen3(i, mode - 1) * fact;
				coords(3, i) = end4Crd(i) + eigen4(i, mode - 1) * fact;
			}
		}
		else {
			for (int i = 0; i < 2; i++) {
				coords(0, i) = end1Crd(i);
				coords(1, i) = end2Crd(i);
				coords(2, i) = end3Crd(i);
				coords(3, i) = end4Crd(i);
			}
		}
	}

	int error = 0;

	// finally we  the element using drawPolygon
	error += theViewer.drawPolygon(coords, values, this->getTag());

	return error;
}

Response*
PDSlip3D::setResponse(const char** argv, int argc,
	OPS_Stream& output)
{
	Response* theResponse = 0;

	if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0) {
		Vector res(1);
		theResponse = new ElementResponse(this, 1, res);
	}
	else if (strcmp(argv[0], "deform") == 0) {
		Vector res(1);
		theResponse = new ElementResponse(this, 2, res);
	}

	return theResponse;
}

int
PDSlip3D::getResponse(int responseID, Information& eleInfo)
{
	if (responseID == 1) {

		Vector res(1);
		res(0) = force;

		return eleInfo.setVector(res);

	}
	else if (responseID == 2) {

		Vector res(1);
		res(0) = deform;

		return eleInfo.setVector(res);

	}
	else
		return -1;
}

void
PDSlip3D::shapeFunction(double xi, double eta, double zeta)
{
	double onePlusxi = 1.0 + xi;
	double oneMinusxi = 1.0 - xi;
	double onePluseta = 1.0 + eta;
	double oneMinuseta = 1.0 - eta;
	double onePluszeta = 1.0 + zeta;
	double oneMinuszeta = 1.0 - zeta;

	shp[0] = 0.125 * oneMinusxi * oneMinuseta * oneMinuszeta;	// N_1
	shp[1] = 0.125 * onePlusxi * oneMinuseta * oneMinuszeta;	// N_2
	shp[2] = 0.125 * onePlusxi * onePluseta * oneMinuszeta;		// N_3
	shp[3] = 0.125 * oneMinusxi * onePluseta * oneMinuszeta;		// N_4
	shp[4] = 0.125 * oneMinusxi * oneMinuseta * onePluszeta;	// N_5
	shp[5] = 0.125 * onePlusxi * oneMinuseta * onePluszeta;	// N_6
	shp[6] = 0.125 * onePlusxi * onePluseta * onePluszeta;		// N_7
	shp[7] = 0.125 * oneMinusxi * onePluseta * onePluszeta;		// N_8

}


void
PDSlip3D::determinPst()
{
	const Vector& nd1Crds = theNodes[0]->getCrds();

	double xMin, yMin, zMin, xMax, yMax, zMax;
	xMin = nd1Crds(0);
	yMin = nd1Crds(1);
	zMin = nd1Crds(2);
	xMax = nd1Crds(0);
	yMax = nd1Crds(1);
	zMax = nd1Crds(2);

	for (int i = 1; i < 8; i++)
	{
		const Vector& ndiCrd = theNodes[i]->getCrds();

		if (xMin > ndiCrd(0))
			xMin = ndiCrd(0);
		if (yMin > ndiCrd(1))
			yMin = ndiCrd(1);
		if (zMin > ndiCrd(2))
			zMin = ndiCrd(2);
		if (xMax < ndiCrd(0))
			xMax = ndiCrd(0);
		if (yMax < ndiCrd(1))
			yMax = ndiCrd(1);
		if (zMax < ndiCrd(2))
			zMax = ndiCrd(2);
	}

	const Vector& nd9Crd = theNodes[8]->getCrds();

	pst(0) = 2 * (nd9Crd(0) - xMin) / (xMax - xMin) - 1;
	pst(1) = 2 * (nd9Crd(1) - yMin) / (yMax - yMin) - 1;
	pst(2) = 2 * (nd9Crd(2) - zMin) / (zMax - zMin) - 1;

}

void 
PDSlip3D::getTransf()
{
	double m = dv(0);
	int loc = 0;

	for (int i = 1; i < 3; i++)
	{
		if (abs(m) < abs(dv(i)))
		{
			m = dv(i);
			loc = i;
		}
	}

	Vector n(3);
	n(loc) = -dv((loc + 1) % 3);
	n((loc + 1) % 3) = m;
	n = n / n.Norm();

	T(0, 0) = dv(0);
	T(0, 1) = dv(1);
	T(0, 2) = dv(2);
	T(1, 0) = n(0);
	T(1, 1) = n(1);
	T(1, 2) = n(2);
	T(2, 0) = T(0, 1) * T(1, 2) - T(1, 1) * T(0, 2);
	T(2, 1) = T(1, 0) * T(0, 2) - T(0, 0) * T(1, 2);
	T(2, 2) = T(0, 0) * T(1, 1) - T(0, 1) * T(1, 0);

	Ttran.addMatrixTranspose(1.0, T, 1.0);
}
