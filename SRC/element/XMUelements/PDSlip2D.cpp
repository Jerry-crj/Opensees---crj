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

#include <PDSlip2D.h>
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

void* OPS_PDSlip2D()
{
	int ndm = OPS_GetNDM();
	int ndf = OPS_GetNDF();

	if (ndm != 2 || ndf != 2) {
		opserr << "WARNING -- model dimensions and/or nodal DOF not compatible with quad element\n";
		return 0;
	}

	if (OPS_GetNumRemainingInputArgs() < 10) {
		opserr << "WARNING insufficient arguments\n";
		opserr << "Want: element PDSlip2D eleTag? iNode? jNode? kNode? lNode? thk? type? matTag? <pressure? rho? b1? b2?>\n";
		return 0;
	}

	// PDSlip2DId, iNode, jNode, kNode, lNode
	int idata[6];
	int num = 6;
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
		opserr << "\nPDSlip2D element: " << idata[0] << endln;
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
		opserr << "\nPDSlip2D element: " << idata[0] << endln;
		return 0;
	}

	double data[2];
	Vector d(data, 2);
	num = 2;

	if (OPS_GetDoubleInput(&num, data) < 0) {
		opserr << "WARNING: invalid double inputs\n";
		return 0;
	}

	return new PDSlip2D(idata[0], idata[1], idata[2], idata[3], idata[4], idata[5],
		*mat1, *mat2, data[0], data[1]);
}


double PDSlip2D::matrixData[100];
Matrix PDSlip2D::K(matrixData, 10, 10);
Vector PDSlip2D::P(10);
double PDSlip2D::shp[4];

PDSlip2D::PDSlip2D(int tag, int nd1, int nd2, int nd3, int nd4, int nd5,
	UniaxialMaterial& m1, UniaxialMaterial& m2,
	double p1, double p2)
	:Element(tag, ELE_TAG_PDSlip2D),
	theSprMaterial(0), connectedExternalNodes(5),
	Ki(0), dv(2), T(2, 2), Ttran(2, 2), pst(2), deform(0), force(0)
{
	// Allocate arrays of pointers to NDMaterials
	theSprMaterial = new UniaxialMaterial * [2];

	if (theSprMaterial == 0) {
		opserr << "PDSlip2D::PDSlip2D - failed allocate material model pointer\n";
		exit(-1);
	}

	theSprMaterial[0] = m1.getCopy();
	theSprMaterial[1] = m2.getCopy();

	// Check allocation
	if ((theSprMaterial[0] == 0) || (theSprMaterial[1] == 0)) {
		opserr << "PDSlip2D::PDSlip2D -- failed to get a copy of material model\n";
		exit(-1);
	}

	// Set connected external node IDs
	connectedExternalNodes(0) = nd1;
	connectedExternalNodes(1) = nd2;
	connectedExternalNodes(2) = nd3;
	connectedExternalNodes(3) = nd4;
	connectedExternalNodes(4) = nd5;

	for (int i = 0; i < 5; i++)
		theNodes[i] = 0;

	dv(0) = p1;
	dv(1) = p2;

	if (dv.Norm() < 1e-10)
		opserr << "PDSlip2D::PDSlip2D -- Dirction Vector.Norm is zero \n";

	dv = dv / dv.Norm();

	this->getTransf();

}

PDSlip2D::PDSlip2D()
	:Element(0, ELE_TAG_PDSlip2D),
	theSprMaterial(0), connectedExternalNodes(5),
	Ki(0), dv(2), T(2, 2), Ttran(2, 2), pst(2), deform(0), force(0)
{
	for (int i = 0; i < 5; i++)
		theNodes[i] = 0;
}

PDSlip2D::~PDSlip2D()
{
	for (int i = 0; i < 2; i++) {
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
PDSlip2D::getNumExternalNodes() const
{
	return 5;
}

const ID&
PDSlip2D::getExternalNodes()
{
	return connectedExternalNodes;
}


Node**
PDSlip2D::getNodePtrs(void)
{
	return theNodes;
}

int
PDSlip2D::getNumDOF()
{
	return 10;
}

void
PDSlip2D::setDomain(Domain* theDomain)
{
	// Check Domain is not null - invoked when object removed from a domain
	if (theDomain == 0) {
		theNodes[0] = 0;
		theNodes[1] = 0;
		theNodes[2] = 0;
		theNodes[3] = 0;
		theNodes[4] = 0;
		return;
	}

	int Nd1 = connectedExternalNodes(0);
	int Nd2 = connectedExternalNodes(1);
	int Nd3 = connectedExternalNodes(2);
	int Nd4 = connectedExternalNodes(3);
	int Nd5 = connectedExternalNodes(4);

	theNodes[0] = theDomain->getNode(Nd1);
	theNodes[1] = theDomain->getNode(Nd2);
	theNodes[2] = theDomain->getNode(Nd3);
	theNodes[3] = theDomain->getNode(Nd4);
	theNodes[4] = theDomain->getNode(Nd5);

	if (theNodes[0] == 0 || theNodes[1] == 0 || theNodes[2] == 0 ||
		theNodes[3] == 0 || theNodes[4] == 0) {
		//opserr << "FATAL ERROR PDSlip2D (tag: %d), node not found in domain",
		//	this->getTag());

		return;
	}

	int dofNd1 = theNodes[0]->getNumberDOF();
	int dofNd2 = theNodes[1]->getNumberDOF();
	int dofNd3 = theNodes[2]->getNumberDOF();
	int dofNd4 = theNodes[3]->getNumberDOF();
	int dofNd5 = theNodes[4]->getNumberDOF();

	if (dofNd1 != 2 || dofNd2 != 2 || dofNd3 != 2 || dofNd4 != 2 || dofNd5 != 2) {
		//opserr << "FATAL ERROR PDSlip2D (tag: %d), has differing number of DOFs at its nodes",
		//	this->getTag());

		return;
	}

	this->DomainComponent::setDomain(theDomain);

	this->determinPst();

	if (pst(0) < -1 || pst(0) > 1 || pst(1) < -1 || pst(1) > 1)
	{
		opserr << "PDSlip2D::PDSlip2D -- Wrong range of coordinates in parameter space \n";
		exit(-1);
	}
}

int
PDSlip2D::commitState()
{
	int retVal = 0;

	// call element commitState to do any base class stuff
	if ((retVal = this->Element::commitState()) != 0) {
		opserr << "PDSlip2D::commitState () - failed in base class";
	}

	// Loop over the integration points and commit the material states

	for (int i = 0; i < 2; i++)
		retVal += theSprMaterial[i]->commitState();

	return retVal;
}

int
PDSlip2D::revertToLastCommit()
{
	int retVal = 0;

	// Loop over the integration points and revert to last committed state

	for (int i = 0; i < 2; i++)
		retVal += theSprMaterial[i]->revertToLastCommit();

	return retVal;
}

int
PDSlip2D::revertToStart()
{
	int retVal = 0;

	// Loop over the integration points and revert states to start

	for (int i = 0; i < 2; i++)
		retVal += theSprMaterial[i]->revertToStart();

	return retVal;
}


int
PDSlip2D::update()
{
	const Vector& disp1 = theNodes[0]->getTrialDisp();
	const Vector& disp2 = theNodes[1]->getTrialDisp();
	const Vector& disp3 = theNodes[2]->getTrialDisp();
	const Vector& disp4 = theNodes[3]->getTrialDisp();
	const Vector& disp5 = theNodes[4]->getTrialDisp();

	static double u[2][4];

	u[0][0] = disp1(0);
	u[1][0] = disp1(1);
	u[0][1] = disp2(0);
	u[1][1] = disp2(1);
	u[0][2] = disp3(0);
	u[1][2] = disp3(1);
	u[0][3] = disp4(0);
	u[1][3] = disp4(1);

	static Vector disp5ref(2);

	int ret = 0;

	// Determine Jacobian for this integration point
	this->shapeFunction(pst(0), pst(1));

	disp5ref.Zero();
	for (int beta = 0; beta < 4; beta++) {
		disp5ref(0) += shp[beta] * u[0][beta];
		disp5ref(1) += shp[beta] * u[1][beta];
	}

	Vector du(2);
	du = T * (disp5 - disp5ref);

	// Set the material strain
	ret += theSprMaterial[0]->setTrialStrain(du(0));
	ret += theSprMaterial[1]->setTrialStrain(du(1));

	deform = du(0);

	return ret;
}


const Matrix&
PDSlip2D::getTangentStiff()
{

	K.Zero();

	Matrix ksg(2, 2);
	ksg(0, 0) = theSprMaterial[0]->getTangent();
	ksg(1, 1) = theSprMaterial[1]->getTangent();

	ksg = Ttran * ksg * T;

	this->shapeFunction(pst(0), pst(1));

	for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia += 2) {
		for (int beta = 0, ib = 0; beta < 4; beta++, ib += 2) {
			K(ia, ib) = shp[alpha] * shp[beta] * ksg(0, 0);
			K(ia, ib + 1) = shp[alpha] * shp[beta] * ksg(0, 1);
			K(ia + 1, ib) = shp[alpha] * shp[beta] * ksg(1, 0);
			K(ia + 1, ib + 1) = shp[alpha] * shp[beta] * ksg(1, 1);
		}
	}

	for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia += 2) {
		K(ia, 8) = -shp[alpha] * ksg(0, 0);
		K(ia, 9) = -shp[alpha] * ksg(0, 1);
		K(ia + 1, 8) = -shp[alpha] * ksg(1, 0);
		K(ia + 1, 9) = -shp[alpha] * ksg(1, 1);
		K(8, ia) = K(ia, 8);
		K(9, ia) = K(ia, 9);
		K(8, ia + 1) = K(ia + 1, 8);
		K(9, ia + 1) = K(ia + 1, 9);
	}

	K(8, 8) = ksg(0, 0);
	K(8, 9) = ksg(0, 1);
	K(9, 8) = ksg(1, 0);
	K(9, 9) = ksg(1, 1);

	return K;
}


const Matrix&
PDSlip2D::getInitialStiff()
{
	if (Ki != 0)
		return *Ki;

	K.Zero();

	Ki = new Matrix(K);

	Matrix ksg(2, 2);
	ksg(0, 0) = theSprMaterial[0]->getInitialTangent();
	ksg(1, 1) = theSprMaterial[1]->getInitialTangent();

	ksg = Ttran * ksg * T;

	this->shapeFunction(pst(0), pst(1));

	for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia += 2) {
		for (int beta = 0, ib = 0; beta < 4; beta++, ib += 2) {
			(*Ki)(ia, ib) = shp[alpha] * shp[beta] * ksg(0, 0);
			(*Ki)(ia, ib + 1) = shp[alpha] * shp[beta] * ksg(0, 1);
			(*Ki)(ia + 1, ib) = shp[alpha] * shp[beta] * ksg(1, 0);
			(*Ki)(ia + 1, ib + 1) = shp[alpha] * shp[beta] * ksg(1, 1);
		}
	}

	for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia += 2) {
		(*Ki)(ia, 8) = -shp[alpha] * ksg(0, 0);
		(*Ki)(ia, 9) = -shp[alpha] * ksg(0, 1);
		(*Ki)(ia + 1, 8) = -shp[alpha] * ksg(1, 0);
		(*Ki)(ia + 1, 9) = -shp[alpha] * ksg(1, 1);
		(*Ki)(8, ia) = (*Ki)(ia, 8);
		(*Ki)(9, ia) = (*Ki)(ia, 9);
		(*Ki)(8, ia + 1) = (*Ki)(ia + 1, 8);
		(*Ki)(9, ia + 1) = (*Ki)(ia + 1, 9);
	}

	(*Ki)(8, 8) = ksg(0, 0);
	(*Ki)(8, 9) = ksg(0, 1);
	(*Ki)(9, 8) = ksg(1, 0);
	(*Ki)(9, 9) = ksg(1, 1);

	return *Ki;
}

const Matrix&
PDSlip2D::getMass()
{
	K.Zero();

	return K;
}

void
PDSlip2D::zeroLoad(void)
{
	return;
}

int
PDSlip2D::addLoad(ElementalLoad* theLoad, double loadFactor)
{
	// Added option for applying body forces in load pattern: C.McGann, U.Washington
	int type;
	const Vector& data = theLoad->getData(type, loadFactor);

	if (type == LOAD_TAG_SelfWeight) {
		return 0;
	}
	else {
		opserr << "PDSlip2D::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
		return -1;
	}

	return -1;
}

int
PDSlip2D::addInertiaLoadToUnbalance(const Vector& accel)
{
	return 0;
}

const Vector&
PDSlip2D::getResistingForce()
{
	P.Zero();

	Vector sig(2), sigg(2);
	sig(0) = theSprMaterial[0]->getStress();
	sig(1) = theSprMaterial[1]->getStress();

	sigg = Ttran * sig;

	this->shapeFunction(pst(0), pst(1));

	for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia += 2) {
		P(ia) =  -shp[alpha] * sigg(0);
		P(ia + 1) = -shp[alpha] * sigg(1);
	}

	P(8) = sigg(0);
	P(9) = sigg(1);

	force = sigg(0);

	return P;
}

const Vector&
PDSlip2D::getResistingForceIncInertia()
{
	return P;
}

int
PDSlip2D::sendSelf(int commitTag, Channel& theChannel)
{
	int res = 0;

	// note: we don't check for dataTag == 0 for Element
	// objects as that is taken care of in a commit by the Domain
	// object - don't want to have to do the check if sending data
	int dataTag = this->getDbTag();

	// Quad packs its data into a Vector and sends this to theChannel
	// along with its dbTag and the commitTag passed in the arguments
	static Vector data(9);
	data(0) = this->getTag();

	data(5) = alphaM;
	data(6) = betaK;
	data(7) = betaK0;
	data(8) = betaKc;

	res += theChannel.sendVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING PDSlip2D::sendSelf() - " << this->getTag() << " failed to send Vector\n";
		return res;
	}


	// Now quad sends the ids of its materials
	int matDbTag;

	static ID idData(12);

	return res;
}

int
PDSlip2D::recvSelf(int commitTag, Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
	int res = 0;

	int dataTag = this->getDbTag();

	// Quad creates a Vector, receives the Vector and then sets the 
	// internal data with the data in the Vector
	static Vector data(9);
	res += theChannel.recvVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING PDSlip2D::recvSelf() - failed to receive Vector\n";
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
		opserr << "WARNING PDSlip2D::recvSelf() - " << this->getTag() << " failed to receive ID\n";
		return res;
	}

	connectedExternalNodes(0) = idData(8);
	connectedExternalNodes(1) = idData(9);
	connectedExternalNodes(2) = idData(10);
	connectedExternalNodes(3) = idData(11);

	return res;
}

void
PDSlip2D::Print(OPS_Stream& s, int flag)
{
	if (flag == 2) {

		s << "#PDSlip2D\n";

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
PDSlip2D::displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes, int numMode)
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
PDSlip2D::setResponse(const char** argv, int argc,
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
PDSlip2D::getResponse(int responseID, Information& eleInfo)
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
PDSlip2D::shapeFunction(double xi, double eta)
{
	double oneMinuseta = 1.0 - eta;
	double onePluseta = 1.0 + eta;
	double oneMinusxi = 1.0 - xi;
	double onePlusxi = 1.0 + xi;

	shp[0] = 0.25 * oneMinusxi * oneMinuseta;	// N_1
	shp[1] = 0.25 * onePlusxi * oneMinuseta;	// N_2
	shp[2] = 0.25 * onePlusxi * onePluseta;		// N_3
	shp[3] = 0.25 * oneMinusxi * onePluseta;	// N_4

}


void
PDSlip2D::determinPst()
{
	const Vector& nd1Crds = theNodes[0]->getCrds();
	const Vector& nd2Crds = theNodes[1]->getCrds();
	const Vector& nd3Crds = theNodes[2]->getCrds();
	const Vector& nd4Crds = theNodes[3]->getCrds();

	double xMin, yMin, xMax, yMax;
	xMin = nd1Crds(0);
	yMin = nd1Crds(1);
	xMax = nd1Crds(0);
	yMax = nd1Crds(1);

	for (int i = 1; i < 4; i++)
	{
		const Vector& ndiCrd = theNodes[i]->getCrds();

		if (xMin > ndiCrd(0))
			xMin = ndiCrd(0);
		if (yMin > ndiCrd(1))
			yMin = ndiCrd(1);
		if (xMax < ndiCrd(0))
			xMax = ndiCrd(0);
		if (yMax < ndiCrd(1))
			yMax = ndiCrd(1);
	}

	const Vector& nd5Crd = theNodes[4]->getCrds();

	pst(0) = 2 * (nd5Crd(0) - xMin) / (xMax - xMin) - 1;
	pst(1) = 2 * (nd5Crd(1) - yMin) / (yMax - yMin) - 1;

}

void 
PDSlip2D::getTransf()
{
	T(0, 0) = dv(0);
	T(0, 1) = dv(1);
	T(1, 0) = -dv(1);
	T(1, 1) = dv(0);

	Ttran.addMatrixTranspose(1.0, T, 1.0);
}
